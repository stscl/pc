/************************************************************************
 *  File: fnn.hpp
 *
 *  False Nearest Neighbor (FNN) computation
 *  for embedding dimension selection in nonlinear systems.
 *
 *  This module provides efficient routines to evaluate the proportion
 *  of false nearest neighbors as embedding dimension increases, a key
 *  criterion for phase space reconstruction.
 *
 *  Core functionality:
 *      - Evaluate FNN ratio between successive dimensions (E1 → E2)
 *      - Compute FNN curve across all embedding dimensions
 *
 *  Method:
 *      For each point, the nearest neighbor is identified in dimension E1.
 *      The stability of this neighbor is then tested in dimension E2.
 *      A neighbor is considered false if:
 *
 *          - Relative distance increase exceeds Rtol, or
 *          - Absolute distance exceeds Atol
 *
 *  Parallelization:
 *      - Supports parallel computation over prediction points
 *      - Optional parallel evaluation across embedding dimensions
 *
 *  Distance metrics:
 *      "euclidean"  : sqrt(sum((x - y)^2))
 *      "maximum"    : max(|x - y|)
 *      "manhattan"  : sum(|x - y|)
 *
 *  Notes:
 *      - Input embedding must have at least 2 dimensions
 *      - Library and prediction sets must be valid index subsets
 *      - Rtol and Atol must match the number of evaluated dimensions
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 ************************************************************************/

#ifndef PC_FNN_HPP
#define PC_FNN_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <algorithm>
#include "pc/numericutils.hpp"
#include "pc/distance.hpp"
#include <RcppThread.h>

namespace pc
{

namespace fnn
{  
    /***************************************************************************************
     *  Compute False Nearest Neighbor (FNN) ratio between two embedding dimensions
     *
     *  This function evaluates the proportion of false nearest neighbors when
     *  increasing the embedding dimension from E1 to E2 (E2 = E1 + 1).
     *
     *  For each prediction point, the k nearest neighbors are identified using
     *  the first E1 dimensions. The neighbor relationships are then re-evaluated
     *  in the higher dimension E2. A neighbor is classified as "false" if:
     *
     *      1. The relative distance increase exceeds Rtol, or
     *      2. The absolute distance in the added dimension exceeds Atol
     *
     *  For each point, the fraction of false neighbors within its k-nearest
     *  neighborhood is computed, and a majority rule is applied to determine
     *  whether the point exhibits false nearest neighbors.
     *
     *  The final FNN ratio is defined as the proportion of points classified
     *  as having false nearest neighbors, indicating insufficient embedding
     *  in lower dimensions.
     *
     *  Parallelization:
     *      Computation over prediction indices can be parallelized using threads.
     *
     *  Parameters:
     *      embedding   : Matrix-like structure (row = observation, col = dimension)
     *      lib         : Indices of library points (0-based)
     *      pred        : Indices of prediction points (0-based)
     *      E1          : Lower embedding dimension
     *      E2          : Higher embedding dimension (E2 is suggested to be E1 + 1)
     *      dist_metric : Distance metric ("euclidean", "manhattan", "maximum")
     *      k           : Number of nearest neighbors used for evaluation
     *      Rtol        : Relative distance threshold
     *      Atol        : Absolute distance threshold
     *      threads     : Number of threads for parallel computation
     *
     *  Returns:
     *      Proportion of false nearest neighbors (double)
     *      Returns NaN if computation is not feasible
     ***************************************************************************************/
    inline double singlefnn(
        const std::vector<std::vector<double>>& embedding,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t E1,
        size_t E2,
        const std::string& dist_metric = "euclidean",
        size_t k = 3,
        double Rtol = 10.0,
        double Atol = 2.0,
        size_t threads = 1) 
    {
        if (embedding.empty() || embedding[0].size() < E2) 
        {
            return std::numeric_limits<double>::quiet_NaN();  // Invalid dimensions
        }

        std::vector<int> false_flags(pred.size(), -1); // -1 means skip or invalid, 0 means not a false neighbor, 1 means false neighbor

        // --------------------------------------------------------------------------
        // Utility to compute false nearest neighbor flag between E1 and E2 dimension
        // --------------------------------------------------------------------------
        auto compute_flag = [&](size_t p) {
            size_t pidx = pred[p];
            
            // // Legacy 1-NN implementation retained for reference.
            // double min_dist = std::numeric_limits<double>::max();
            // size_t nn_idx = embedding.size(); // invalid index placeholder

            // for (size_t j = 0; j < lib.size(); ++j) 
            // {
            //     size_t lidx = lib[j];
            //     if (pidx == lidx) continue;

            //     // Compute distance using only the first E1 dimensions
            //     std::vector<double> xi(embedding[pidx].begin(), embedding[pidx].begin() + E1);
            //     std::vector<double> xj(embedding[lidx].begin(), embedding[lidx].begin() + E1);
            //     double dist = pc::distance::distance(xi, xj, dist_metric, true); 

            //     if (dist < min_dist) 
            //     {
            //         min_dist = dist;
            //         nn_idx = lidx; 
            //     }
            // }

            // // Skip if no neighbor found or minimum distance is zero
            // if (nn_idx == embedding.size() || pc::numericutils::doubleNearlyEqual(min_dist, 0.0)) return;

            // // Compare the E2-th dimension to check for false neighbors
            // double diff = std::abs(embedding[pidx][E2 - 1] - embedding[nn_idx][E2 - 1]);
            // double ratio = diff / min_dist;

            // // Determine if this is a false neighbor
            // if (ratio > Rtol || diff > Atol) 
            // {
            //     false_flags[p] = 1;
            // } 
            // else 
            // {
            //     false_flags[p] = 0;
            // }

            // --------------------------------------------------
            // Collect distances to all library points (E1 space)
            // --------------------------------------------------
            std::vector<std::pair<double, size_t>> dists;
            dists.reserve(lib.size());

            for (size_t j = 0; j < lib.size(); ++j) 
            {
                size_t lidx = lib[j];
                if (pidx == lidx) continue;

                std::vector<double> xi(embedding[pidx].begin(), embedding[pidx].begin() + E1);
                std::vector<double> xj(embedding[lidx].begin(), embedding[lidx].begin() + E1);

                double dist = pc::distance::distance(xi, xj, dist_metric, true);
                if (!std::isnan(dist)) dists.emplace_back(dist, lidx);
            }

            if (dists.empty()) return;

            // --------------------------------------------------
            // Select k nearest neighbors (partial sort)
            // --------------------------------------------------
            size_t k_use = std::min(k, dists.size());

            std::nth_element(
                dists.begin(),
                dists.begin() + k_use,
                dists.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; }
            );

            // --------------------------------------------------
            // Evaluate false neighbor condition over k neighbors
            // --------------------------------------------------
            size_t false_k = 0;
            size_t valid_k = 0;

            for (size_t i = 0; i < k_use; ++i)
            {
                double dist = dists[i].first;
                size_t idx = dists[i].second;

                double diff = std::abs(embedding[pidx][E2 - 1] - embedding[idx][E2 - 1]);

                if (pc::numericutils::doubleNearlyEqual(dist, 0.0))
                {
                    valid_k++;
                    if (diff > Atol) false_k++;
                    continue;
                }

                double ratio = diff / dist;

                if (ratio > Rtol || diff > Atol)
                    false_k++;

                valid_k++;
            }

            // --------------------------------------------------
            // Majority rule at point level
            // --------------------------------------------------
            if (valid_k > 0)
            {
                double frac = static_cast<double>(false_k) / valid_k;
                false_flags[p] = (frac > 0.5) ? 1 : 0;
            }
        };

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

    /************************************************************************************ 
     * Compute False Nearest Neighbor (FNN) ratios across multiple embedding dimensions
     *
     * For a given embedding matrix (with each row representing a sample unit and
     * each column an embedding dimension), this function evaluates the proportion
     * of false nearest neighbors (FNN) as the embedding dimension increases.
     *
     * For each embedding dimension pair (E1, E2), where E2 = E1 + 1, the FNN ratio
     * is computed using a k-nearest neighbor criterion. For each prediction point,
     * the k nearest neighbors are identified in the E1-dimensional space, and their
     * relationships are re-evaluated in the higher dimension E2.
     *
     * A neighbor is classified as "false" if the relative or absolute distance
     * increase exceeds predefined thresholds. For each point, the fraction of false
     * neighbors within its local neighborhood is computed, and a majority rule is
     * applied to determine whether the point exhibits false nearest neighbors.
     *
     * The FNN ratio measures the proportion of points whose neighborhood structure
     * is not preserved when increasing dimensionality, indicating that the lower
     * dimension E1 is insufficient for reconstructing the system.
     *
     * Setting k = 1 reduces the method to the standard single nearest neighbor variant.
     *
     * Parameters:
     * - embedding: A vector of vectors where each row is a unit’s embedding.
     *              Must have at least 2 columns (dimensions).
     * - lib: A vector of indices indicating the library set (0-based).
     * - pred: A vector of indices indicating the prediction set (0-based).
     * - Rtol: A vector of relative distance thresholds (one per E1).
     * - Atol: A vector of absolute distance thresholds (one per E1).
     * - dist_metric: Distance metric ("euclidean", "manhattan", "maximum")
     * - k: Number of nearest neighbors used for evaluation
     * - threads: Number of threads to use for parallel computation.
     * - parallel_level: Parallelization strategy
     *                        0 = prediction point-level parallelism
     *                        1 = embedding dimension-level parallelism
     *
     * Returns:
     * - A vector of FNN ratios corresponding to each E1 from 1 to D - 1.
     *   If not computable for a given E1, NaN is returned at that position.
     *************************************************************************************/
    inline std::vector<double> fnn(
        const std::vector<std::vector<double>>& embedding,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        const std::vector<double>& Rtol,
        const std::vector<double>& Atol,
        const std::string& dist_metric = "euclidean",
        size_t k = 3,
        size_t threads = 1,
        size_t parallel_level = 0) 
    {
        // Configure threads
        threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads);

        size_t max_E2 = embedding[0].size();
        std::vector<double> results(max_E2 - 1, std::numeric_limits<double>::quiet_NaN());

        if (embedding.empty() || embedding[0].size() < 2) 
        {
            return results;  // Not enough dimensions to compute FNN
        }

        if (parallel_level == 0)
        {
            // Loop through E1 = 1 to max_E2 - 1
            for (size_t E1 = 1; E1 < max_E2; ++E1) 
            {
                size_t E2 = E1 + 1;
                double fnn_ratio = singlefnn(embedding, lib, pred, E1, E2, dist_metric,
                                             k, Rtol[E1 - 1], Atol[E1 - 1], threads);
                results[E1 - 1] = fnn_ratio;
            }
        } 
        else 
        {
            RcppThread::parallelFor(1, max_E2, [&](size_t E1) {
                size_t E2 = E1 + 1;
                double fnn_ratio = singlefnn(embedding, lib, pred, E1, E2, dist_metric,
                                             k, Rtol[E1 - 1], Atol[E1 - 1], 1);
                results[E1 - 1] = fnn_ratio;
            }, threads);
        }

        return results;
    }

} // namespace fnn

}

#endif // PC_FNN_HPP
