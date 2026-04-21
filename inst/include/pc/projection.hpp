#ifndef PC_PROJECTION_HPP
#define PC_PROJECTION_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "pc/numericutils.h"
#include <RcppThread.h>

namespace pc
{

namespace projection
{ 

/**
 * Predicts signature vectors for a subset of target points using weighted nearest neighbors.
 *
 * This function performs local weighted prediction in the signature space as follows:
 *   1. For each prediction index `p` in `pred_indices`, find its `num_neighbors` nearest neighbors
 *      among `lib_indices` based on distances in `Dx[p][*]`, ignoring NaN distances.
 *   2. Compute exponential weights scaled by the total distance sum to emphasize close points.
 *      If all distances are zero, uniform weights are used instead.
 *   3. For each dimension of the signature space:
 *        - Count how many neighbor signatures are exactly zero.
 *        - If the zero count exceeds `zero_tolerance`, set the predicted value to 0.
 *        - Otherwise, compute a weighted average of valid (non-NaN) neighbor signatures.
 *   4. Predictions are stored and updated only for indices in `pred_indices`; other entries remain undefined (NaN).
 *
 * Parallelization:
 *   - Controlled by the parameter `threads`.
 *   - If `threads <= 1`, computation is serial (standard for-loop).
 *   - Otherwise, the loop over prediction indices is executed in parallel via RcppThread::parallelFor.
 *
 * @param SMy            Signature space of the target variable Y. Shape: (N_obs, E−1)
 * @param Dx             Distance matrix from prediction points to library points. Shape: (SMy.size(), SMy.size())
 * @param lib_indices    Indices of valid library points used for neighbor search (subset of [0, SMy.size())).
 * @param pred_indices   Indices of points to predict (subset of [0, SMy.size())).
 * @param num_neighbors  Number of nearest neighbors to use. If <= 0, defaults to E+1.
 * @param zero_tolerance Maximum allowed zero values per dimension before forcing prediction to zero.
 *                       If <= 0, defaults to E−1.
 * @param h              Prediction horizon (time shift). Defines how far ahead in time the prediction is performed.
 *                       For each base index p, nearest neighbors are identified at time p, and their future states 
 *                       at time (lib_row + h) are used to predict the target state at time (p + h).
 * @param threads        Number of threads to use. If <= 1, runs serially; otherwise runs parallel.
 *
 * @return A matrix of predicted signature vectors, sized SMy.size() × (E−1).
 */
std::vector<std::vector<double>> projection(
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& Dx,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = 0,   /* = std::numeric_limits<int>::min() */
    int zero_tolerance = 0,  /* = std::numeric_limits<int>::max() */
    size_t h = 0,
    size_t threads = 1
) {
  const size_t n_obs = SMy.size();
  const size_t n_sig_dim = SMy[0].size(); // E−1

  // Infer embedding dimension E
  const int E = static_cast<int>(n_sig_dim + 1);

  // Set defaults if sentinel values are used
  if (num_neighbors <= 0) {
    num_neighbors = E + 1; // E+1 neighbors
  }
  if (zero_tolerance <= 0) {
    zero_tolerance = E - 1; // default: E−1
  }

  // Output containers: same size as SMy
  std::vector<std::vector<double>> pred_signatures(
      n_obs,std::vector<double>(n_sig_dim, std::numeric_limits<double>::quiet_NaN()));

  if (SMy.empty() || Dx.empty() || lib_indices.empty() || pred_indices.empty()) {
    return pred_signatures;
  }

  // // standard for-loop for serial computation
  // for (size_t pi = 0; pi < pred_indices.size(); ++pi) {
  //   size_t p = pred_indices[pi];
  //
  //   // Step 1: Collect valid (non-NaN) distances and corresponding library indices
  //   std::vector<double> distances;
  //   std::vector<size_t> valid_libs;
  //   distances.reserve(lib_indices.size());
  //   valid_libs.reserve(lib_indices.size());
  //
  //   for (size_t lib_idx : lib_indices) {
  //     // // Bounds check for Dx[p][lib_idx]
  //     // if (p >= Dx.size() || lib_idx >= Dx[p].size()) {
  //     //   continue;
  //     // }
  //     double d = Dx[p][lib_idx];
  //     if (!std::isnan(d)) {
  //       distances.push_back(d);
  //       valid_libs.push_back(lib_idx);
  //     }
  //   }
  //
  //   // If no valid distances found, prediction is unchanged (NaNs)
  //   if (distances.empty()) {
  //     continue;
  //   }
  //
  //   // Step 2: Determine actual k (cannot exceed available valid neighbors)
  //   size_t k = std::min(static_cast<size_t>(num_neighbors), distances.size());
  //
  //   // Step 3: Find k nearest neighbors via partial sort
  //   std::vector<size_t> neighbor_indices(distances.size());
  //   std::iota(neighbor_indices.begin(), neighbor_indices.end(), 0);
  //
  //   std::partial_sort(
  //     neighbor_indices.begin(),
  //     neighbor_indices.begin() + k,
  //     neighbor_indices.end(),
  //     [&](size_t a, size_t b) {
  //       return (distances[a] < distances[b]) ||
  //         (pc::numericutils::doubleNearlyEqual(distances[a],distances[b]) && a < b);
  //     }
  //   );
  //
  //   // Step 4: Compute exponential weights
  //   // Method 1: traditional for-loop (slightly faster — minimal overhead, best for performance)
  //   double total_dist = 0.0;
  //   for (size_t i = 0; i < k; ++i) {
  //     total_dist += distances[neighbor_indices[i]];
  //   }
  //   // // Method 2: using std::accumulate with a lambda (cleaner, but slightly slower due to lambda call overhead)
  //   // double total_dist = std::accumulate(
  //   //   neighbor_indices.begin(), neighbor_indices.begin() + k, 0.0,
  //   //   [&](double acc, size_t idx) {
  //   //     return acc + distances[idx];
  //   //   });
  //
  //   std::vector<double> weights(k);
  //   if (pc::numericutils::doubleNearlyEqual(total_dist,0.0)) {
  //     // All distances zero → uniform weights
  //     std::fill(weights.begin(), weights.end(), 1.0);
  //   } else {
  //     for (size_t i = 0; i < k; ++i) {
  //       weights[i] = std::exp(-distances[neighbor_indices[i]] / total_dist);
  //       // // Enforce minimum weight to avoid underflow to zero
  //       // double w = std::exp(-distances[neighbor_indices[i]] / total_dist);
  //       // weights[i] = std::max(w, 1e-6);
  //     }
  //   }
  //
  //   double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
  //
  //   // Step 5: Predict each dimension of signature
  //   for (size_t dim = 0; dim < n_sig_dim; ++dim) {
  //     // Count exact zeros in this dimension among neighbors
  //     int zero_count = 0;
  //     double weighted_sum = 0.0;
  //     bool has_valid = false;
  //
  //     for (size_t i = 0; i < k; ++i) {
  //       size_t lib_row = valid_libs[neighbor_indices[i]];
  //       // // Bounds check
  //       // if (lib_row >= SMy.size() || dim >= SMy[lib_row].size()) {
  //       //   continue;
  //       // }
  //       double val = SMy[lib_row][dim];
  //       if (std::isnan(val)) {
  //         continue;
  //       }
  //       if (pc::numericutils::doubleNearlyEqual(val,0.0)) {
  //         zero_count++;
  //       }
  //       weighted_sum += val * weights[i];
  //       has_valid = true;
  //     }
  //
  //     // Apply zero-tolerance rule
  //     if (zero_count > zero_tolerance) {
  //       pred_signatures[p][dim] = 0.0;
  //     } else if (has_valid && total_weight > 0.0) {
  //       pred_signatures[p][dim] = weighted_sum / total_weight; // normalize by total weight
  //     }
  //     // else: remains NaN (handled by initialization)
  //   }
  // }

  // Define per-prediction computation (single index)
  auto predict_fn = [&](size_t pi) {
    size_t p = pred_indices[pi];

    std::vector<double> distances;
    std::vector<size_t> valid_libs;
    distances.reserve(lib_indices.size());
    valid_libs.reserve(lib_indices.size());

    for (size_t lib_idx : lib_indices) {
      double d = Dx[p][lib_idx];
      if (!std::isnan(d)) {
        distances.push_back(d);
        valid_libs.push_back(lib_idx);
      }
    }
    if (distances.empty()) return;

    size_t k = std::min(static_cast<size_t>(num_neighbors), distances.size());
    std::vector<size_t> neighbor_indices(distances.size());
    std::iota(neighbor_indices.begin(), neighbor_indices.end(), 0);
    std::partial_sort(
      neighbor_indices.begin(),
      neighbor_indices.begin() + k,
      neighbor_indices.end(),
      [&](size_t a, size_t b) {
        if (!pc::numericutils::doubleNearlyEqual(distances[a], distances[b])) {
          return distances[a] < distances[b];
        } else {
          return a < b;
        }
      });

    double total_dist = 0.0;
    for (size_t i = 0; i < k; ++i) total_dist += distances[neighbor_indices[i]];
    std::vector<double> weights(k);
    if (pc::numericutils::doubleNearlyEqual(total_dist,0.0)) {
      std::fill(weights.begin(), weights.end(), 1.0);
    } else {
      for (size_t i = 0; i < k; ++i) {
        weights[i] = std::exp(-distances[neighbor_indices[i]] / total_dist);
        // // Enforce minimum weight to avoid underflow to zero
        // double w = std::exp(-distances[neighbor_indices[i]] / total_dist);
        // weights[i] = std::max(w, 1e-6);
      }
    }
    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    for (size_t dim = 0; dim < n_sig_dim; ++dim) {
      int zero_count = 0;
      double weighted_sum = 0.0;
      bool has_valid = false;

      for (size_t i = 0; i < k; ++i) {
        size_t lib_row = valid_libs[neighbor_indices[i]];
        double val = SMy[lib_row][dim];
        if (std::isnan(val)) continue;
        if (pc::numericutils::doubleNearlyEqual(val,0.0)) zero_count++;
        weighted_sum += val * weights[i];
        has_valid = true;
      }

      if (zero_count > zero_tolerance) {
        pred_signatures[p][dim] = 0.0;
      } else if (has_valid && total_weight > 0.0) {
        pred_signatures[p][dim] = weighted_sum / total_weight;
      }
    }
  };

  // Parallel or serial execution
  if (threads <= 1) {
    for (size_t pi = 0; pi < pred_indices.size(); ++pi) predict_fn(pi);
  } else {
    RcppThread::parallelFor(0, pred_indices.size(), predict_fn, threads);
  }

  return pred_signatures;
}

} // namespace projection

}

#endif // PC_PROJECTION_HPP
