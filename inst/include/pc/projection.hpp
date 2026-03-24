/**************************************************************************
 * File: projection.hpp
 *
 * Nonparametric time-delay projection using weighted nearest neighbors.
 *
 * Provides functionality for:
 *   - Predicting future signature vectors in a reconstructed space.
 *   - Performing local neighbor-based forecasting with distance weighting.
 *   - Supporting time-delay prediction via a configurable horizon (h).
 *
 * The method:
 *   - Selects nearest neighbors based on a distance matrix.
 *   - Applies exponential weighting to emphasize closer neighbors.
 *   - Projects neighbor states forward in time (t + h).
 *   - Computes robust weighted averages with NaN handling and zero filtering.
 *
 * Designed for:
 *   - Nonlinear time series prediction
 *   - State-space reconstruction workflows
 *   - Causality analysis (e.g., cross mapping, symbolic methods)
 *
 * Parallel execution is supported via RcppThread.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 *************************************************************************/

#ifndef PC_PROJECTION_HPP
#define PC_PROJECTION_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "pc/numericutils.hpp"
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
     * @param num_neighbors  Number of nearest neighbors to use. If == 0, defaults to E+1.
     * @param zero_tolerance Maximum allowed zero values per dimension before forcing prediction to zero.
     *                       If == 0, defaults to E−1.
     * @param h              Prediction horizon (time shift). Defines how far ahead in time the prediction is performed.
     *                       For each base index p, nearest neighbors are identified at time p, and their future states 
     *                       at time (lib_row + h) are used to predict the target state at time (p + h).
     * @param threads        Number of threads to use. If <= 1, runs serially; otherwise runs parallel.
     *
     * @return A matrix of predicted signature vectors, sized SMy.size() × (E−1).
     */
    inline std::vector<std::vector<double>> projection(
        const std::vector<std::vector<double>>& SMy,
        const std::vector<std::vector<double>>& Dx,
        const std::vector<size_t>& lib_indices,
        const std::vector<size_t>& pred_indices,
        size_t num_neighbors = 0,  
        size_t zero_tolerance = 0,  
        size_t h = 0,
        size_t threads = 1
    ) {
      const size_t n_obs = SMy.size();
      const size_t n_sig_dim = SMy[0].size(); // E−1

      // Infer embedding dimension E
      const size_t E = n_sig_dim + 1;

      // Set defaults if sentinel values are used
      if (num_neighbors == 0) {
        num_neighbors = E + 1; // E+1 neighbors
      }
      if (zero_tolerance == 0) {
        zero_tolerance = E - 1; // default: E−1
      }

      // Output containers: same size as SMy
      std::vector<std::vector<double>> pred_signatures(
          n_obs, std::vector<double>(n_sig_dim, std::numeric_limits<double>::quiet_NaN()));

      if (SMy.empty() || Dx.empty() || lib_indices.empty() || pred_indices.empty()) {
        return pred_signatures;
      }

      // Define per-prediction computation (single index)
      auto predict_fn = [&](size_t pi) {
        size_t p = pred_indices[pi];
        if (p + h >= n_obs) return;

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

        size_t k = std::min(num_neighbors, distances.size());
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

        if (h == 0) { // no projection horizon
          double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

          for (size_t dim = 0; dim < n_sig_dim; ++dim) {
            size_t zero_count = 0;
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
        } else {
          for (size_t dim = 0; dim < n_sig_dim; ++dim) {
            size_t zero_count = 0;
            double weighted_sum = 0.0;
            double used_weight = 0.0;
            bool has_valid = false;

            for (size_t i = 0; i < k; ++i) {
              size_t lib_row = valid_libs[neighbor_indices[i]];
              size_t target_row = lib_row + h;
              if (target_row >= n_obs) continue;
              double val = SMy[target_row][dim];
              if (std::isnan(val)) continue;
              if (pc::numericutils::doubleNearlyEqual(val,0.0)) zero_count++;
              weighted_sum += val * weights[i];
              used_weight += weights[i];
              has_valid = true;
            }

            if (zero_count > zero_tolerance) {
              pred_signatures[p + h][dim] = 0.0;
            } else if (has_valid && used_weight > 0.0) {
              pred_signatures[p + h][dim] = weighted_sum / used_weight;
            }
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
