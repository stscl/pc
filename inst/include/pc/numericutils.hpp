/********************************************************************************
 * File: numericutils.hpp
 *
 * Provides helper functions for:
 *   - Floating-point comparison with combined relative and absolute tolerance.
 *   - Portable numeric constants (epsilon and tolerance).
 *   - Special mathematical functions such as the digamma function ψ(x).
 *   - Basic statistical utilities (mean, quantiles).
 *
 * Statistical features:
 *   mean(vec)
 *
 *   quantile(vec, probs)
 *       Computes sample quantiles using R's Type 7 method:
 *           h = 1 + (n - 1) * p
 *
 *       Supports multiple probability inputs and linear interpolation.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ********************************************************************************/

#ifndef PC_NUMERICUTILS_HPP
#define PC_NUMERICUTILS_HPP

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <initializer_list>

namespace pc
{

namespace numericutils
{

    // ==============================
    // Common numeric constants
    // ==============================

    constexpr double DOUBLE_EPS     = std::numeric_limits<double>::epsilon(); // ≈ 2.22e-16
    constexpr double DOUBLE_TOL_ABS = 1.5e-16;  // Absolute tolerance
    constexpr double DOUBLE_TOL_REL = 1.5e-8;   // Relative tolerance

    constexpr double PI = 3.141592653589793238462643383279502884; // π

    // ==============================
    // Floating-point comparison
    // ==============================

    /**
     * @brief Compare two double values with combined relative and absolute tolerance.
     *
     * Implements a numerically stable test for near equality:
     *
     * |x - y| <= max(rel_tol * max(|x|, |y|, 1.0), abs_tol)
     *
     * @param x First value
     * @param y Second value
     * @param rel_tol Relative tolerance (default DOUBLE_TOL_REL)
     * @param abs_tol Absolute tolerance (default DOUBLE_TOL_ABS)
     * @return true if x and y are considered equal within tolerance
     */
    inline bool doubleNearlyEqual(double x,
                                  double y,
                                  double rel_tol = DOUBLE_TOL_REL,
                                  double abs_tol = DOUBLE_TOL_ABS) noexcept
    {
        double diff  = std::fabs(x - y);
        double scale = std::max({1.0, std::fabs(x), std::fabs(y)});
        return diff <= std::max(rel_tol * scale, abs_tol);
    }

    // ==============================
    // Special mathematical functions
    // ==============================

    /**
     * @brief Computes the digamma function ψ(x) with high numerical accuracy.
     *
     * The digamma function is defined as
     *
     *      ψ(x) = d/dx log Γ(x)
     *
     * Implementation strategy
     *
     * 1. Reflection formula for negative arguments
     *
     *      ψ(1 - x) = ψ(x) + π cot(πx)
     *
     *    This allows evaluation for x < 0 except at poles
     *    (0, -1, -2, ...).
     *
     * 2. Recurrence relation
     *
     *      ψ(x) = ψ(x + 1) − 1/x
     *
     *    Used to shift small x values upward into a stable region.
     *
     * 3. Asymptotic Bernoulli expansion for large x
     *
     *      ψ(x) ≈ log(x) − 1/(2x)
     *             − 1/(12x²)
     *             + 1/(120x⁴)
     *             − 1/(252x⁶)
     *             + ...
     *
     * Accuracy is typically around 1e-15 for double precision inputs.
     *
     * @param x Input value
     * @return digamma value ψ(x)
     */
    inline double digamma(double x) noexcept
    {
        // Handle NaN
        if (std::isnan(x))
            return std::numeric_limits<double>::quiet_NaN();

        // Pole at non-positive integer
        if (std::floor(x) == x && x <= 0)
            return std::numeric_limits<double>::infinity();

        // Reflection formula for negative arguments
        if (x < 0.0)
        {
            return digamma(1.0 - x) - PI / std::tan(PI * x);
        }

        double result = 0.0;

        // Recurrence to shift x into asymptotic region
        while (x < 8.0)
        {
            result -= 1.0 / x;
            x += 1.0;
        }

        double inv_x  = 1.0 / x;
        double inv_x2 = inv_x * inv_x;

        // Bernoulli asymptotic expansion
        double series =
            inv_x2 * (-1.0 / 12.0 +
            inv_x2 * ( 1.0 / 120.0 +
            inv_x2 * (-1.0 / 252.0 +
            inv_x2 * ( 1.0 / 240.0 +
            inv_x2 * (-5.0 / 660.0 +
            inv_x2 * ( 691.0 / 32760.0 +
            inv_x2 * (-1.0 / 12.0 +
                      3617.0 / 8160.0)))))));

        return result + std::log(x) - 0.5 * inv_x + series;
    }

    // ==============================
    // Basic statistical functions
    // ==============================

    /**
     * @brief Compute mean of a numeric vector.
     *
     * @param vec Input vector
     * @return Mean value (NaN if no valid data)
     */
    inline double mean(const std::vector<double>& vec) noexcept
    {
        if (vec.empty())
            return std::numeric_limits<double>::quiet_NaN();

        double sum = 0.0;
        size_t count = 0;

        for (double v : vec)
        {
            if (!std::isnan(v))
            {
                sum += v;
                ++count;
            }
        }

        if (count == 0)
            return std::numeric_limits<double>::quiet_NaN();

        return sum / static_cast<double>(count);
    }

    /**
     * @brief Compute quantiles of a numeric vector (R Type 7).
     *
     * Uses linear interpolation:
     *
     *   h = 1 + (n - 1) * p
     *
     * @param vec Input data vector
     * @param probs Probabilities (default: {0.05, 0.5, 0.95})
     * @return Vector of quantiles
     */
    inline std::vector<double> quantile(
        const std::vector<double>& vec,
        const std::vector<double>& probs = {0.05, 0.5, 0.95})
    {
        std::vector<double> clean_vec;

        clean_vec.reserve(vec.size());
        for (double v : vec)
        {
            if (!std::isnan(v))
                clean_vec.push_back(v);
        }

        if (clean_vec.empty())
        {
            return std::vector<double>(
                probs.size(),
                std::numeric_limits<double>::quiet_NaN());
        }

        std::sort(clean_vec.begin(), clean_vec.end());

        size_t n = clean_vec.size();
        std::vector<double> results;
        results.reserve(probs.size());

        for (double p : probs)
        {
            if (p < 0.0 || p > 1.0)
                throw std::out_of_range("probabilities must be between 0 and 1");

            if (n == 1)
            {
                results.push_back(clean_vec[0]);
                continue;
            }

            double h = 1.0 + (n - 1) * p;
            size_t hf = static_cast<size_t>(std::floor(h));
            double gamma = h - hf;

            size_t idx_lower = (hf <= 1) ? 0 : hf - 1;
            size_t idx_upper = (hf >= n) ? n - 1 : hf;

            double q = (1.0 - gamma) * clean_vec[idx_lower] +
                       gamma * clean_vec[idx_upper];

            results.push_back(q);
        }

        return results;
    }

} // namespace numericutils

}

#endif // PC_NUMERICUTILS_HPP
