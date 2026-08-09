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

#ifndef PC_DMI_HPP
#define PC_DMI_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "pc/ksginfo.hpp"
#include <RcppThread.h>

namespace pc
{

namespace dmi
{ 
    /**
     * Predicts signature vectors for a subset of target points using weighted nearest neighbors.
     *
     * This function performs local weighted prediction in the signature space as follows:
     *   1. For each prediction index `p` in `pred_indices`, find its `num_neighbors` nearest neighbors
     *      among `lib_indices` based on distances in `Dx[p][*]`, ignoring NaN distances and `p` itself.
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
    inline std::vector<double> dmi(
        const std::vector<double>& vec,
        const std::vector<size_t>& pred,
        const std::vector<size_t>& tau,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false,
        size_t threads = 1) 
    {
        // Result vector, initialized with NaN
        std::vector<double> result(tau.size(),
            std::numeric_limits<double>::quiet_NaN());

        const size_t n = pred.size();
        const size_t m = tau.size();

        if (vec.empty() || pred.empty() || tau.empty()) 
        {
            return result;
        }

        // Configure threads
        threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads);

        // Matrix: rows = tau.size() + 1, cols = pred.size()
        // Row 0: original values
        // Row i: lagged values with lag = tau[i-1]
        std::vector<std::vector<double>> mat(
            m + 1, std::vector<double>(n, std::numeric_limits<double>::quiet_NaN())
        );

        // Fill the matrix
        for (size_t col = 0; col < n; ++col) 
        {
            size_t idx = pred[col];

            // Column 0: current value
            if (idx < vec.size()) 
            {
                mat[0][col] = vec[idx];
            }

            // Lagged columns
            for (size_t i = 0; i < m; ++i) 
            {
                size_t lag = tau[i];

                // Check if lag is valid
                if (idx >= lag && (idx - lag) < vec.size()) 
                {
                    mat[i + 1][col] = vec[idx - lag];
                } else 
                {
                    // leave as NaN
                    mat[i + 1][col] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }

        if (threads <= 1)
        {
            // Compute MI for each lag
            for (size_t i = 1; i <= m; ++i) 
            {
                result[i - 1] = pc::ksginfo::mi(
                    mat,
                    0,      // reference column (current)
                    i,      // lagged column
                    k,
                    alg,
                    base,
                    normalize
                );
            }
        } 
        else 
        {
            RcppThread::parallelFor(1, m, [&](size_t i) {
                result[i - 1] = pc::ksginfo::mi(
                    mat,
                    0,      // reference column (current)
                    i,      // lagged column
                    k,
                    alg,
                    base,
                    normalize
                );
            }, threads);
        }

        return result;

    } // namespace dmi

}

#endif // PC_DMI_HPP
