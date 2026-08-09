/*************************************************************************
 *  File: dmi.hpp
 *
 *  Delayed mutual information (DMI) using k-nearest neighbor estimators.
 *
 *  This header computes mutual information between a reference series
 *  and its lagged versions across multiple delay steps.
 *
 *  The dmi function constructs a matrix:
 *
 *      mat[var][obs]
 *
 *      where:
 *          var = 0        -> reference values (indexed by pred)
 *          var = i >= 1   -> lagged values with lag tau[i - 1]
 *
 *      Missing lagged observations are filled with NaN.
 *
 *  Data layout:
 *      Series = std::vector<double>
 *      Index  = std::vector<size_t>
 *      Matrix = std::vector<std::vector<double>> // mat[var][obs]
 *
 *  Computation:
 *
 *      For each lag tau[i], mutual information is computed as:
 *
 *          MI( mat[0], mat[i] )
 *
 *      using pc::ksginfo::mi.
 *
 *  Estimator variants (alg parameter):
 *
 *      alg = 0
 *            Kraskov–Stögbauer–Grassberger estimator I (KSG1)
 *
 *      alg = 1
 *            Kraskov–Stögbauer–Grassberger estimator II (KSG2)
 *
 *  Parallelization:
 *
 *      Independent lag computations can be parallelized
 *      using RcppThread::parallelFor.
 *
 *  Dependencies:
 *
 *      pc::ksginfo::mi
 *      RcppThread
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
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
    /*
     * Delayed Mutual Information (DMI) for Multiple Lags
     *
     * This function computes delayed mutual information between a reference
     * variable and its lagged values using a KSG-based estimator. The input
     * series is indexed by a set of prediction locations, and lagged values
     * are constructed according to the specified lag vector.
     *
     * Internally, a matrix of size (m + 1) x n is built, where:
     *   - m = tau.size() is the number of lag steps
     *   - n = pred.size() is the number of samples
     *   - Row 0 stores the reference values (indexed by pred)
     *   - Row i (i >= 1) stores values lagged by tau[i - 1]
     *
     * Missing lagged values (due to boundary issues) are filled with NaN.
     * Mutual information is then computed between the reference row (row 0)
     * and each lagged row using pc::ksginfo::mi.
     *
     * Optionally, computations across lags can be parallelized.
     *
     * Parameters:
     *   vec        - Input numeric vector representing the ordered series
     *   pred       - Indices defining the sample positions (from past to present)
     *   tau        - Vector of lag steps (non-negative integers)
     *   k          - Number of nearest neighbors for KSG estimator (default: 3)
     *   alg        - Algorithm variant for KSG estimator (default: 0)
     *   base       - Logarithm base for mutual information (default: 2.0)
     *   normalize  - Whether to normalize the MI values (default: false)
     *   threads    - Number of threads for parallel execution (default: 1)
     *
     * Returns:
     *   A vector of length tau.size(), where each element corresponds to the
     *   estimated mutual information at the given lag. Values are initialized
     *   as NaN and remain NaN if computation is not feasible.
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

        // Fill row 0 first
        for (size_t col = 0; col < n; ++col)
        {   
            mat[0][col] = pred[col];
        }

        // Fill lagged rows
        for (size_t i = 0; i < m; ++i) 
        {
            size_t lag = tau[i];
            for (size_t col = 0; col < n; ++col) 
            {
                size_t idx = pred[col];
                if (idx >= lag && (idx - lag) < vec.size()) 
                {
                    mat[i + 1][col] = vec[idx - lag];
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
    }

} // namespace dmi

}

#endif // PC_DMI_HPP
