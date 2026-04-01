#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include "pc.h"

// Wrapper function to perform false nearest neighbor analysis
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppFNN(
    const Rcpp::NumericVector& target,
    const Rcpp::NumericVector& rt,
    const Rcpp::NumericVector& eps,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    int tau = 1,
    int style = 0,
    const std::string& dist_metric = "euclidean",
    int threads = 1,
    int parallel_level = 0,
    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
    Rcpp::Nullable<int> nrows = R_NilValue)
{
    // --- Input Conversion and Validation ---
    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    const size_t n_obs = tg.size();

    // Convert library indices (R 1-based → C++ 0-based)
    std::vector<size_t> lib_std = Rcpp::as<std::vector<size_t>>(lib);
    for (auto& idx : lib_std) 
    {
        if (idx < 1 || idx > n_obs) 
        {
            Rcpp::stop("lib index %d out of bounds [1, %d]",
                       static_cast<int>(idx),
                       static_cast<int>(n_obs));
        }
        idx -= 1;
    }

    // Convert prediction indices (R 1-based → C++ 0-based)
    std::vector<size_t> pred_std = Rcpp::as<std::vector<size_t>>(pred);
    for (auto& idx : pred_std) 
    {
        if (idx < 1 || idx > n_obs) 
        {
            Rcpp::stop("pred index %d out of bounds [1, %d]",
                       static_cast<int>(idx),
                       static_cast<int>(n_obs));
        }
        idx -= 1;
    }

    // Construct embedding dimension E
    std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
    for (auto& singleE : E_std) 
    {
        if (singleE < 0) singleE = std::abs(singleE);
    }
    size_t max_E = static_cast<size_t>(*std::max_element(E_std.begin(), E_std.end()));
    if (max_E < 2) max_E = 2;

    // Convert rt and eps
    std::vector<double> rt_vec = Rcpp::as<std::vector<double>>(rt);
    std::vector<double> eps_vec = Rcpp::as<std::vector<double>>(eps);

    // Expand rt and eps
    std::vector<double> rt_std(max_E);
    std::vector<double> eps_std(max_E);

    // ---- rt ----
    if (rt_vec.size() == 1) 
    {
        std::fill(rt_std.begin(), rt_std.end(), rt_vec[0]);
    } 
    else 
    {
        size_t src_len = rt_vec.size();
        for (size_t i = 0; i < rt_std.size(); ++i) 
        {
            rt_std[i] = rt_vec[i % src_len];
        }
    }

    // ---- eps ----
    if (eps_vec.size() == 1) 
    {
        std::fill(eps_std.begin(), eps_std.end(), eps_vec[0]);
    } 
    else 
    {
        size_t src_len = eps_vec.size();
        for (size_t i = 0; i < eps_std.size(); ++i) {
            eps_std[i] = eps_vec[i % src_len];
        }
    }

    // --- Embedding Construction ---
    std::vector<std::vector<double>> Mx;

    if (nb.isNotNull()) 
    {
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        std::vector<std::vector<size_t>> nb_std = pc::convert::nb2std(nb.get());
        Mx = pc::embed::embed(
            tg, nb_std, max_E, 
            static_cast<size_t>(std::abs(tau)), 
            static_cast<size_t>(std::abs(style)));
    } 
    else if (nrows.isNotNull())
    {   
        size_t n_rows = static_cast<size_t>(std::abs(Rcpp::as<int>(nrows)));

        std::vector<std::vector<double>> tm = 
            pc::embed::gridVec2Mat(tg, n_rows);
        Mx = pc::embed::embed(
            tm, max_E, 
            static_cast<size_t>(std::abs(tau)), 
            static_cast<size_t>(std::abs(style)));
    }
    else  
    {
        Mx = pc::embed::embed(
            tg, max_E, 
            static_cast<size_t>(std::abs(tau)), 
            static_cast<size_t>(std::abs(style)));

        size_t max_lag = (tau == 0) 
            ? (max_E - 1)
            : ((max_E - 1) * static_cast<size_t>(std::abs(tau)));

        lib_std.erase(
            std::remove_if(lib_std.begin(), lib_std.end(), 
                [&](size_t idx){ return idx + 1 < max_lag; }),
            lib_std.end()
        );

        pred_std.erase(
            std::remove_if(pred_std.begin(), pred_std.end(), 
                [&](size_t idx){ return idx + 1 < max_lag; }),
            pred_std.end()
        );
    }

    // ---- sort + unique lib/pred ----
    std::sort(lib_std.begin(), lib_std.end());
    lib_std.erase(
        std::unique(lib_std.begin(), lib_std.end()),
        lib_std.end()
    );

    std::sort(pred_std.begin(), pred_std.end());
    pred_std.erase(
        std::unique(pred_std.begin(), pred_std.end()),
        pred_std.end()
    );

    // ---- filter lib/pred (remove NaN in target/source) ----
    size_t write = 0;
    for (size_t i = 0; i < lib_std.size(); ++i)
    {
        size_t idx = lib_std[i];
        if (!std::isnan(tg[idx]))
        {
            lib_std[write++] = idx;
        }
    }
    lib_std.resize(write);

    write = 0;
    for (size_t i = 0; i < pred_std.size(); ++i)
    {
        size_t idx = pred_std[i];
        if (!std::isnan(tg[idx]))
        {
            pred_std[write++] = idx;
        }
    }
    pred_std.resize(write);
    
    // --- Prepare for data slicing ---
    std::vector<size_t> selected_indices;
    selected_indices.reserve(lib_std.size() + pred_std.size());
    for (size_t i = 0; i < lib_std.size(); ++i)
        selected_indices.push_back(lib_std[i]);
    for (size_t i = 0; i < pred_std.size(); ++i)
        selected_indices.push_back(pred_std[i]);
    std::sort(selected_indices.begin(), selected_indices.end());
    selected_indices.erase(
        std::unique(selected_indices.begin(), selected_indices.end()),
        selected_indices.end()
    );

    // --- Check if full set is used ---
    bool use_subset = (selected_indices.size() < Mx.size());

    // --- Perform Pattern Causality Analysis ---
    std::vector<double> res;

    if (!use_subset)
    {
        // --- Full data: no slicing needed ---
        res = pc::fnn::fnn(
            Mx, lib_std, pred_std, rt_std, eps_std, dist_metric,
            static_cast<size_t>(std::abs(threads)), 
            static_cast<size_t>(std::abs(parallel_level)));
    }
    else
    {   
        // --- Slice Mx  ---
        std::vector<std::vector<double>> Mx_sub;
        Mx_sub.reserve(selected_indices.size());

        for (size_t i = 0; i < selected_indices.size(); ++i)
        {
            size_t idx = selected_indices[i];
            Mx_sub.push_back(Mx[idx]);
        }

        // --- Subset mode: build index map ---
        std::unordered_map<size_t, size_t> index_map;
        index_map.reserve(selected_indices.size());

        for (size_t i = 0; i < selected_indices.size(); ++i)
        {
            index_map[selected_indices[i]] = i;
        }

        // --- Remap lib indices ---
        for (size_t i = 0; i < lib_std.size(); ++i)
        {
            lib_std[i] = index_map[lib_std[i]];
        }

        // --- Remap pred indices ---
        for (size_t i = 0; i < pred_std.size(); ++i)
        {
            pred_std[i] = index_map[pred_std[i]];
        }

        // --- Run patcaus on subset ---
        res = pc::fnn::fnn(
            Mx_sub, lib_std, pred_std, rt_std, eps_std, dist_metric,
            static_cast<size_t>(std::abs(threads)), 
            static_cast<size_t>(std::abs(parallel_level)));
    }

    // Convert the result back to Rcpp::NumericVector and set names as "E:1", "E:2", ..., "E:n"
    Rcpp::NumericVector result = Rcpp::wrap(res);
    Rcpp::CharacterVector resnames(result.size());
    for (int i = 0; i < result.size(); ++i) {
        resnames[i] = "E:" + std::to_string(i + 1);
    }
    result.names() = resnames;

    // Terminal-friendly hint (one-time, non-intrusive)
    Rcpp::Rcout << "[fnn] Output 'E:i' corresponds to the i-th valid embedding dimension.\n"
                << "[fnn] Input E values exceeding max embeddable dimension were truncated.\n"
                << "[fnn] Please map output indices to original E inputs before interpretation.\n";

    return result;
}
