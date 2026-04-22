
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "pc/numericutils.h"
#include "pc/distance.h"
#include "pc/symdync.h"
#include "pc/projection.h"
#include <RcppThread.h>


inline pc::symdync::PatternCausalityRes patcaus(
    const std::vector<std::vector<double>>& Mx,
    const std::vector<std::vector<double>>& My,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    size_t num_neighbors = 0,
    size_t zero_tolerance = 0,
    size_t h = 0,
    std::string dist_metric = "euclidean",
    bool relative = true,
    bool weighted = true,
    size_t threads = 1
){
  // Configure threads (cap at hardware concurrency)
  threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads);

  const size_t n_obs = Mx.size();

  // Initialize distance matrix (n_obs × n_obs) filled with NaN
  std::vector<std::vector<double>> Dx(
      n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

  // --------------------------------------------------------------------------
  // Step 1: Compute pairwise distances between prediction and library indices
  // --------------------------------------------------------------------------
  auto compute_distance = [&](size_t p) {
    size_t pi = pred_indices[p];
    for (size_t li : lib_indices) {
      double dist = pc::distance::distance(Mx[pi], Mx[li], dist_metric, true);
      if (!std::isnan(dist)) {
        Dx[pi][li] = dist;  // assign distance; no mirroring required
      }
    }
  };

  // Parallel or serial execution depending on thread configuration
  if (threads <= 1) {
    for (size_t p = 0; p < pred_indices.size(); ++p)
      compute_distance(p);
  } else {
    RcppThread::parallelFor(0, pred_indices.size(), compute_distance, threads);
  }

  // --------------------------------------------------------------------------
  // Step 2: Generate signature spaces for Mx and My
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> SMx = pc::symdync::genSignatureSpace(Mx, relative);
  std::vector<std::vector<double>> SMy = pc::symdync::genSignatureSpace(My, relative);

  // --------------------------------------------------------------------------
  // Step 3: Predict target signatures for My using local projections
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> PredSMy = pc::projection::projection(
    SMy, Dx, lib_indices, pred_indices, num_neighbors, zero_tolerance, h, threads);

  // --------------------------------------------------------------------------
  // Step 4: Compute pattern-based causality using symbolic pattern comparison
  // --------------------------------------------------------------------------
  pc::symdync::PatternCausalityRes res = 
    pc::symdync::computePatternCausality(SMx, SMy, PredSMy, weighted);

  return res;
}

inline std::vector<std::vector<std::vector<double>>> patcaus(
    const std::vector<std::vector<double>>& Mx,
    const std::vector<std::vector<double>>& My,
    const std::vector<size_t>& libsizes,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    size_t num_neighbors = 0,
    size_t zero_tolerance = 0,
    size_t h = 0,
    std::string dist_metric = "euclidean",
    size_t boot = 99,
    bool random_sample = true,
    unsigned long long seed = 42,
    bool relative = true,
    bool weighted = true,
    size_t threads = 1,
    size_t parallel_level = 0,
    bool verbose = false
){
  // --------------------------------------------------------------------------
  // Step 1: Configure threads and random generators
  // --------------------------------------------------------------------------
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Enforce boot = 1 for deterministic sampling
  if (!random_sample) boot = 1;

  // Prebuild 64-bit RNG pool for reproducibility
  std::vector<std::mt19937_64> rng_pool(boot);
  for (int i = 0; i < boot; ++i) {
    std::seed_seq seq{static_cast<uint64_t>(seed), static_cast<uint64_t>(i)};
    rng_pool[i] = std::mt19937_64(seq);
  }

  // --------------------------------------------------------------------------
  // Step 2: Compute pairwise distances once
  // --------------------------------------------------------------------------
  const size_t n_obs = Mx.size();
  std::vector<std::vector<double>> Dx(
      n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

  bool L1norm = (dist_metric == 1);

  auto compute_distance = [&](size_t p) {
    size_t pi = pred_indices[p];
    for (size_t li : lib_indices) {
      double dist = CppDistance(Mx[pi], Mx[li], L1norm, true);
      if (!std::isnan(dist)) Dx[pi][li] = dist;
    }
  };

  if (threads_sizet != 1)
    for (size_t p = 0; p < pred_indices.size(); ++p) compute_distance(p);
  else
    RcppThread::parallelFor(0, pred_indices.size(), compute_distance, threads_sizet);

  // --------------------------------------------------------------------------
  // Step 3: Generate signature spaces
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> SMx = GenSignatureSpace(Mx, relative);
  std::vector<std::vector<double>> SMy = GenSignatureSpace(My, relative);

  // --------------------------------------------------------------------------
  // Step 4: Initialize results container [3][libsizes][boot]
  // --------------------------------------------------------------------------
  const size_t n_libsizes = libsizes.size();
  std::vector<std::vector<std::vector<double>>> all_results(
      3, std::vector<std::vector<double>>(n_libsizes, std::vector<double>(boot, std::numeric_limits<double>::quiet_NaN())));

  // Optional progress bar
  std::unique_ptr<RcppThread::ProgressBar> bar;
  if (progressbar)
    bar = std::make_unique<RcppThread::ProgressBar>(n_libsizes, 1);

  // --------------------------------------------------------------------------
  // Step 5: Iterate over library sizes
  // --------------------------------------------------------------------------
  for (size_t li = 0; li < n_libsizes; ++li) {
    size_t L = libsizes[li];

    auto process_boot = [&](int b) {
      std::vector<size_t> sampled_lib, sampled_pred;

      if (random_sample) {
        std::vector<size_t> shuffled_lib = lib_indices;
        std::shuffle(shuffled_lib.begin(), shuffled_lib.end(), rng_pool[b]);
        sampled_lib.assign(shuffled_lib.begin(), shuffled_lib.begin() + L);
        // sampled_pred = sampled_lib;
      } else {
        sampled_lib.assign(lib_indices.begin(), lib_indices.begin() + L);
        // sampled_pred = sampled_lib;
      }

      std::vector<std::vector<double>> PredSMy;
      if (parallel_level == 0)
        PredSMy = SignatureProjection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, threads_sizet);
      else
        PredSMy = SignatureProjection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, 1);

      PatternCausalityRes res = GenPatternCausality(SMx, SMy, PredSMy, weighted);

      all_results[0][li][b] = res.TotalPos;
      all_results[1][li][b] = res.TotalNeg;
      all_results[2][li][b] = res.TotalDark;
    };

    if (parallel_level != 0)
      RcppThread::parallelFor(0, boot, process_boot, threads_sizet);
    else
      for (int b = 0; b < boot; ++b) process_boot(b);

    if (progressbar) (*bar)++;
  }

  return all_results;
}