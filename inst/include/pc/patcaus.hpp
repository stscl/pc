/************************************************************************
 *  File: patcaus.hpp
 *
 *  Pattern Causality implementation based on
 *  symbolic dynamics and local projection methods.
 *
 *  This module provides high-performance routines for detecting
 *  directional causality between two multivariate time series
 *  using pattern-based symbolic representations.
 *
 *  Core workflow:
 *      1. Compute pairwise distances in reconstructed state space
 *      2. Transform trajectories into symbolic signature space
 *      3. Predict target signatures via nearest-neighbor projection
 *      4. Evaluate causality via pattern matching statistics
 *
 *  Notes:
 *      - Distance matrix is computed only for (pred × lib) pairs
 *      - Parallelism can be controlled at projection or bootstrap level
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 ************************************************************************/

#ifndef PC_PATCAUS_HPP
#define PC_PATCAUS_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cstdint>
#include <iterator>
#include <random>
#include <memory>
#include "pc/distance.hpp"
#include "pc/symdync.hpp"
#include "pc/projection.hpp"
#include <RcppThread.h>

namespace pc
{

namespace patcaus
{   
    /**********************************************************************************
     *  Compute pattern-based causality from Mx to My using a fixed
     *  library and prediction set.
     *
     *  Parameters:
     *      Mx, My          : Reconstructed state-space matrices
     *      lib_indices     : Indices used as library (reference points)
     *      pred_indices    : Indices used for prediction
     *      num_neighbors   : Number of nearest neighbors (0 = auto)
     *      zero_tolerance  : Threshold for zero-distance handling
     *      h               : Prediction horizon
     *      dist_metric     : Distance metric ("euclidean", "manhattan", "maximum".)
     *      relative        : Use relative symbolic encoding
     *      weighted        : Use weighted pattern comparison
     *      threads         : Number of threads for parallel execution
     *      save_detail     : Make per-sample causality output
     **********************************************************************************/
    inline pc::symdync::PatternCausalityRes patcaus(
        const std::vector<std::vector<double>>& Mx,
        const std::vector<std::vector<double>>& My,
        const std::vector<size_t>& lib_indices,
        const std::vector<size_t>& pred_indices,
        const size_t& num_neighbors = 0,
        const size_t& zero_tolerance = 0,
        const size_t& h = 0,
        const std::string& dist_metric = "euclidean",
        bool relative = true,
        bool weighted = true,
        size_t threads = 1,
        bool save_detail = true)
    {
    // Configure threads (cap at hardware concurrency)
    threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads);

    const size_t n_obs = Mx.size();

    // Initialize distance matrix (n_obs × n_obs) filled with NaN
    std::vector<std::vector<double>> Dx(
        n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

    // --------------------------------------------------------------------------
    // Step 1: Compute pairwise distances between prediction and library indices
    // --------------------------------------------------------------------------
    auto compute_distance = [&](size_t p) 
    {
        size_t pi = pred_indices[p];
        for (size_t li : lib_indices) 
        {   
            if (pi == li) continue;

            double dist = pc::distance::distance(Mx[pi], Mx[li], dist_metric, true);
            if (!std::isnan(dist)) 
            {
                Dx[pi][li] = dist;  // assign distance; no mirroring required
            }
        }
    };

    // Parallel or serial execution depending on thread configuration
    if (threads <= 1) 
    {
        for (size_t p = 0; p < pred_indices.size(); ++p)
            compute_distance(p);
    } 
    else 
    {
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
    pc::symdync::PatternCausalityRes res;
    bool use_subset = (pred_indices.size() < Mx.size());

    if (!use_subset)
    {
        // --- Full data: no slicing needed ---
        res = pc::symdync::computePatternCausality(
            SMx, SMy, PredSMy, weighted, save_detail);
    }
    else
    {
        // --- Slice Mx and My ---
        std::vector<std::vector<double>> SMx_sub;
        std::vector<std::vector<double>> SMy_sub;
        std::vector<std::vector<double>> PredSMy_sub;

        SMx_sub.reserve(pred_indices.size());
        SMy_sub.reserve(pred_indices.size());
        PredSMy_sub.reserve(pred_indices.size());

        for (size_t i = 0; i < pred_indices.size(); ++i)
        {
            size_t idx = pred_indices[i];
            SMx_sub.push_back(SMx[idx]);
            SMy_sub.push_back(SMy[idx]);
            PredSMy_sub.push_back(PredSMy[idx]);
        }

        // --- Compute patcaus on subset ---
        res = pc::symdync::computePatternCausality(
            SMx_sub, SMy_sub, PredSMy_sub, weighted, save_detail);
    }

    return res;
    }
    
    /************************************************************************
     *  Perform pattern causality analysis across multiple library sizes
     *  and bootstrap replicates.
     *
     *  Parameters:
     *      Mx, My          : Reconstructed state-space matrices
     *      libsizes        : Vector of library sizes to evaluate
     *      lib_indices     : Full set of candidate library indices
     *      pred_indices    : Prediction indices
     *      num_neighbors   : Number of nearest neighbors
     *      zero_tolerance  : Zero-distance handling threshold
     *      h               : Prediction horizon
     *      dist_metric     : Distance metric
     *      boot            : Number of bootstrap replicates
     *      random_sample   : Whether to randomly sample library
     *      seed            : Random seed for reproducibility
     *      relative        : Relative symbolic encoding
     *      weighted        : Weighted pattern comparison
     *      threads         : Number of threads
     *      parallel_level  : Parallelization strategy
     *                        0 = projection-level parallelism
     *                        1 = bootstrap-level parallelism
     *      verbose         : Show progress bar
     *
     *  Returns:
     *      3D array [3 × libsizes × boot]:
     *          [0] Positive causality
     *          [1] Negative causality
     *          [2] Dark causality
     *
     *  Description:
     *      This function evaluates causality robustness by:
     *          - Varying library sizes
     *          - Applying bootstrap resampling
     *
     *      For each library size:
     *          1. Sample library indices (random or sequential)
     *          2. Perform projection-based prediction
     *          3. Compute symbolic pattern causality
     *
     *      Parallelization can be applied either:
     *          - Inside projection (default)
     *          - Across bootstrap replicates
     ************************************************************************/
    inline std::vector<std::vector<std::vector<double>>> patcaus(
        const std::vector<std::vector<double>>& Mx,
        const std::vector<std::vector<double>>& My,
        const std::vector<size_t>& libsizes,
        const std::vector<size_t>& lib_indices,
        const std::vector<size_t>& pred_indices,
        const size_t& num_neighbors = 0,
        const size_t& zero_tolerance = 0,
        const size_t& h = 0,
        const std::string& dist_metric = "euclidean",
        size_t boot = 99,
        bool replace_sampling = true,
        unsigned long long seed = 42,
        bool relative = true,
        bool weighted = true,
        size_t threads = 1,
        size_t parallel_level = 0,
        bool verbose = false)
    {
    // --------------------------------------------------------------------------
    // Step 1: Configure threads and random generators
    // --------------------------------------------------------------------------
    threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads);

    // Prebuild 64-bit RNG pool for reproducibility
    std::vector<std::mt19937_64> rng_pool(boot);
    for (size_t i = 0; i < boot; ++i) {
        std::seed_seq seq{static_cast<uint64_t>(seed), static_cast<uint64_t>(i)};
        rng_pool[i] = std::mt19937_64(seq);
    }

    // --------------------------------------------------------------------------
    // Step 2: Compute pairwise distances once
    // --------------------------------------------------------------------------
    const size_t n_obs = Mx.size();
    std::vector<std::vector<double>> Dx(
        n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

    auto compute_distance = [&](size_t p) 
    {
        size_t pi = pred_indices[p];
        for (size_t li : lib_indices) 
        {   
            if (pi == li) continue;

            double dist = pc::distance::distance(Mx[pi], Mx[li], dist_metric, true);
            if (!std::isnan(dist)) 
            {
                Dx[pi][li] = dist;  // assign distance; no mirroring required
            }
        }
    };

    if (threads <= 1)
        for (size_t p = 0; p < pred_indices.size(); ++p) compute_distance(p);
    else
        RcppThread::parallelFor(0, pred_indices.size(), compute_distance, threads);

    // --------------------------------------------------------------------------
    // Step 3: Generate signature spaces
    // --------------------------------------------------------------------------
    std::vector<std::vector<double>> SMx = pc::symdync::genSignatureSpace(Mx, relative);
    std::vector<std::vector<double>> SMy = pc::symdync::genSignatureSpace(My, relative);

    // --------------------------------------------------------------------------
    // Step 4: Initialize results container [3][libsizes][boot]
    // --------------------------------------------------------------------------
    const size_t n_libsizes = libsizes.size();
    std::vector<std::vector<std::vector<double>>> all_results(
        3, std::vector<std::vector<double>>(n_libsizes, 
            std::vector<double>(boot, std::numeric_limits<double>::quiet_NaN())));

    // Optional progress bar
    std::unique_ptr<RcppThread::ProgressBar> bar;
    if (verbose)
        bar = std::make_unique<RcppThread::ProgressBar>(n_libsizes, 1);

    // --- Check if full set is used ---
    bool use_subset = (pred_indices.size() < Mx.size());

    // --------------------------------------------------------------------------
    // Step 5: Iterate over library sizes
    // --------------------------------------------------------------------------
    for (size_t li = 0; li < n_libsizes; ++li) 
    {
        size_t L = libsizes[li];

        auto process_boot = [&](size_t b) 
        {
            std::vector<size_t> sampled_lib;
            // std::vector<size_t> sampled_pred;
            
            if (boot == 1)
            {
                sampled_lib.assign(lib_indices.begin(), lib_indices.begin() + L);
                // sampled_pred = sampled_lib;
            }
            else if (replace_sampling) 
            {   
                sampled_lib.resize(L);
                std::uniform_int_distribution<size_t> dist(0, lib_indices.size() - 1);

                for (size_t i = 0; i < L; ++i)
                {
                    sampled_lib[i] = lib_indices[dist(rng_pool[b])];
                }
                // sampled_pred = sampled_lib;
            } 
            else 
            {
                std::vector<size_t> shuffled_lib = lib_indices;
                std::shuffle(shuffled_lib.begin(), shuffled_lib.end(), rng_pool[b]);
                sampled_lib.assign(shuffled_lib.begin(), shuffled_lib.begin() + L);
                // sampled_pred = sampled_lib;
            }

            std::vector<std::vector<double>> PredSMy;
            if (parallel_level == 0)
                PredSMy = pc::projection::projection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, h, threads);
            else
                PredSMy = pc::projection::projection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, h, 1);

            pc::symdync::PatternCausalityRes res;
            if (!use_subset)
            {
                // --- Full data: no slicing needed ---
                res = pc::symdync::computePatternCausality(
                    SMx, SMy, PredSMy, weighted, false);
            }
            else
            {
                // --- Slice Mx and My ---
                std::vector<std::vector<double>> SMx_sub;
                std::vector<std::vector<double>> SMy_sub;
                std::vector<std::vector<double>> PredSMy_sub;

                SMx_sub.reserve(pred_indices.size());
                SMy_sub.reserve(pred_indices.size());
                PredSMy_sub.reserve(pred_indices.size());

                for (size_t i = 0; i < pred_indices.size(); ++i)
                {
                    size_t idx = pred_indices[i];
                    SMx_sub.push_back(SMx[idx]);
                    SMy_sub.push_back(SMy[idx]);
                    PredSMy_sub.push_back(PredSMy[idx]);
                }

                // --- Compute patcaus on subset ---
                res = pc::symdync::computePatternCausality(
                    SMx_sub, SMy_sub, PredSMy_sub, weighted, false);
            }

            all_results[0][li][b] = res.TotalPos;
            all_results[1][li][b] = res.TotalNeg;
            all_results[2][li][b] = res.TotalDark;
        };

        if (parallel_level != 0)
            RcppThread::parallelFor(0, boot, process_boot, threads);
        else
            for (size_t b = 0; b < boot; ++b) process_boot(b);

        if (verbose) (*bar)++;
    }

    return all_results;
    }

} // namespace patcaus

}

#endif // PC_PATCAUS_HPP
