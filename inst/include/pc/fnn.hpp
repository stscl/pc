#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include "pc/numericutils.hpp"
#include "pc/distance.hpp"
#include <RcppThread.h>

double singlefnn(const std::vector<std::vector<double>>& embedding,
                 const std::vector<size_t>& lib,
                 const std::vector<size_t>& pred,
                 size_t E1,
                 size_t E2,
                 const std::string& dist_metric = "euclidean",
                 double Rtol = 10.0,
                 double Atol = 2.0,
                 size_t threads = 1) 
{
    if (embedding.empty() || embedding[0].size() < E2) 
    {
        return std::numeric_limits<double>::quiet_NaN();  // Invalid dimensions
    }
    size_t N = embedding.size();

    std::vector<int> false_flags(pred.size(), -1); // -1 means skip or invalid, 0 means not a false neighbor, 1 means false neighbor

    // --------------------------------------------------------------------------
    // Step 1: Compute false nearest neighbor flag between E1 and E2 dimension
    // --------------------------------------------------------------------------
    auto compute_flag = [&](size_t p) {
        size_t pidx = pred[p];

        double min_dist = std::numeric_limits<double>::max();
        int nn_idx = -1;

        for (size_t j = 0; j < lib.size(); ++j) 
        {
            size_t lidx = lib[j];
            if (pidx == lidx) continue;

            // Compute distance using only the first E1 dimensions
            std::vector<double> xi(embedding[pidx].begin(), embedding[pidx].begin() + E1);
            std::vector<double> xj(embedding[lidx].begin(), embedding[lidx].begin() + E1);
            double dist = pc::distance::distance(xi, xj, dist_metric, true); 

            if (dist < min_dist) 
            {
                min_dist = dist;
                nn_idx = static_cast<int>(lidx); 
            }
        }

        // Skip if no neighbor found or minimum distance is zero
        if (nn_idx == -1 || pc::numericutils::doubleNearlyEqual(min_dist,0.0)) return;

        // Compare the E2-th dimension to check for false neighbors
        double diff = std::abs(embedding[pidx][E2 - 1] - embedding[nn_idx][E2 - 1]);
        double ratio = diff / min_dist;

        // Determine if this is a false neighbor
        if (ratio > Rtol || diff > Atol) 
        {
            false_flags[i] = 1;
        } 
        else 
        {
            false_flags[i] = 0;
        }
    }

    if (threads <= 1)
    {
        for (size_t p = 0; p < pred.size(); ++p)
            compute_flag(p);
    }
    else 
    {   
        RcppThread::parallelFor(0, pred.size(), compute_flag, threads);
    }

    // Aggregate results
    size_t false_count = 0, total = 0;
    for (int flag : false_flags) 
    {
      if (flag >= 0) 
      {
        total++;
        if (flag == 1) false_count++;
      }
    }

    if (total > 0) 
    {
      return static_cast<double>(false_count) / static_cast<double>(total);
    } 
    else 
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
}

/*
 * Compute False Nearest Neighbor (FNN) ratios across multiple embedding dimensions
 *
 * For a given embedding matrix (with each row representing a spatial unit and
 * each column an embedding dimension), this function evaluates the proportion
 * of false nearest neighbors (FNN) as the embedding dimension increases.
 *
 * It iteratively calls `CppSingleFNN` for each embedding dimension pair (E1, E2),
 * where E1 ranges from 1 to D - 1 (D = number of columns), and E2 = E1 + 1.
 * The FNN ratio measures how often a nearest neighbor in dimension E1 becomes
 * distant in dimension E2, suggesting that E1 is insufficient for reconstructing
 * the system.
 *
 * If `parallel_level == 0`, the function executes in serial;
 * otherwise, it uses multithreading to compute FNN ratios for each (E1, E2) pair
 * in parallel.
 *
 * Parameters:
 * - embedding: A vector of vectors where each row is a unit’s embedding.
 *              Must have at least 2 columns (dimensions).
 * - lib: A vector of indices indicating the library set (0-based).
 * - pred: A vector of indices indicating the prediction set (0-based).
 * - Rtol: A vector of relative distance thresholds (one per E1).
 * - Atol: A vector of absolute distance thresholds (one per E1).
 * - L1norm: If true, use L1 (Manhattan) distance; otherwise, use L2 (Euclidean).
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: 0 for serial loop over E1, >0 for parallel loop over E1.
 *
 * Returns:
 * - A vector of FNN ratios corresponding to each E1 from 1 to D - 1.
 *   If not computable for a given E1, NaN is returned at that position.
 */
std::vector<double> CppFNN(const std::vector<std::vector<double>>& embedding,
                           const std::vector<size_t>& lib,
                           const std::vector<size_t>& pred,
                           const std::vector<double>& Rtol,
                           const std::vector<double>& Atol,
                           const std::string& dist_metric = "euclidean",
                           size_t threads = 1,
                           size_t parallel_level = 0) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  size_t max_E2 = embedding[0].size();
  std::vector<double> results(max_E2 - 1, std::numeric_limits<double>::quiet_NaN());

  if (embedding.empty() || embedding[0].size() < 2) {
    return results;  // Not enough dimensions to compute FNN
  }

  if (parallel_level == 0){
    // Loop through E1 = 1 to max_E2 - 1
    for (size_t E1 = 1; E1 < max_E2; ++E1) {
      size_t E2 = E1 + 1;
      double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2, threads_sizet,
                                      parallel_level, Rtol[E1 - 1], Atol[E1 - 1], L1norm);
      results[E1 - 1] = fnn_ratio;
    }
  } else {
    // Parallel computation
    RcppThread::parallelFor(1, max_E2, [&](size_t E1) {
      size_t E2 = E1 + 1;
      double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2, threads_sizet,
                                      parallel_level, Rtol[E1 - 1], Atol[E1 - 1], L1norm);
      results[E1 - 1] = fnn_ratio;
    }, threads_sizet);
  }

  return results;
}