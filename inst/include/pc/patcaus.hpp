
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


pc::symdync::PatternCausalityRes patcaus(
    const std::vector<std::vector<double>>& Mx,
    const std::vector<std::vector<double>>& My,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    size_t num_neighbors = 0,
    size_t zero_tolerance = 0,
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
    SMy, Dx, lib_indices, pred_indices, num_neighbors, zero_tolerance, threads_sizet);

  // --------------------------------------------------------------------------
  // Step 4: Compute pattern-based causality using symbolic pattern comparison
  // --------------------------------------------------------------------------
  PatternCausalityRes res = pc::symdync::computePatternCausality(
    SMx, SMy, PredSMy, weighted);

  return res;
}
