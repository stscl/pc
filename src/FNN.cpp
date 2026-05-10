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
Rcpp::List RcppFNN(
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
        if (singleE < 0) singleE = std::abs(singleE)
    }
    size_t max_E = static_cast<size_t>(*std::max_element(E_std.begin(), E_std.end()));
    max_E = std::max(2, max_E);

    // Convert rt and eps
    std::vector<double> rt_vec = Rcpp::as<std::vector<double>>(rt);
    std::vector<double> eps_vec = Rcpp::as<std::vector<double>>(eps);

    // Expand rt and eps
    std::vector<double> rt_std(max_E);
    std::vector<double> eps_std(max_E);

    // ---- E ----
    if (E_vec.size() == 1) 
    {
        std::fill(E_std.begin(), E_std.end(), E_vec[0]);
    } 
    else 
    {
        E_std[0] = E_vec[0];
        E_std[1] = E_vec[1];
    }
    
    // Make sure each E is greater than 2
    for (auto& singleE : E_std) {
        if (singleE < 2) singleE = 2;
    }

    // ---- tau ----
    if (tau_vec.size() == 1) 
    {
        std::fill(tau_std.begin(), tau_std.end(), tau_vec[0]);
    } 
    else 
    {
        tau_std[0] = tau_vec[0];
        tau_std[1] = tau_vec[1];
    }

    // --- Embedding Construction ---
    std::vector<std::vector<double>> Mx;
    std::vector<std::vector<double>> My;

    if (nb.isNotNull()) 
    {
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        std::vector<std::vector<size_t>> nb_std = pc::convert::nb2std(nb.get());
        Mx = pc::embed::embed(
            tg, nb_std, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            sg, nb_std, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));
    } 
    else if (nrows.isNotNull())
    {   
        size_t n_rows = static_cast<size_t>(std::abs(Rcpp::as<int>(nrows)));

        std::vector<std::vector<double>> tm = 
            pc::embed::gridVec2Mat(tg, n_rows);
        Mx = pc::embed::embed(
            tm, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));

        std::vector<std::vector<double>> sm = 
            pc::embed::gridVec2Mat(sg, n_rows);
        My = pc::embed::embed(
            sm, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));
    }
    else  
    {
        Mx = pc::embed::embed(
            tg, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            sg, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));

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
        if (!(std::isnan(tg[idx]) || std::isnan(sg[idx])))
        {
            lib_std[write++] = idx;
        }
    }
    lib_std.resize(write);

    write = 0;
    for (size_t i = 0; i < pred_std.size(); ++i)
    {
        size_t idx = pred_std[i];
        if (!(std::isnan(tg[idx]) || std::isnan(sg[idx])))
        {
            pred_std[write++] = idx;
        }
    }
    pred_std.resize(write);
    // Copy prediction indices for mapping back to original dataset
    std::vector<size_t> pred_indices = pred_std;
    
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
    pc::symdync::PatternCausalityRes res;

    if (!use_subset)
    {
        // --- Full data: no slicing needed ---
        res = pc::patcaus::patcaus(
            Mx, My, lib_std, pred_std, 
            static_cast<size_t>(std::abs(num_neighbors)),
            static_cast<size_t>(std::abs(zero_tolerance)),
            static_cast<size_t>(std::abs(h)),
            dist_metric, relative, weighted,
            static_cast<size_t>(std::abs(threads)), true);
    }
    else
    {   
        // --- Slice Mx and My ---
        std::vector<std::vector<double>> Mx_sub;
        std::vector<std::vector<double>> My_sub;

        Mx_sub.reserve(selected_indices.size());
        My_sub.reserve(selected_indices.size());

        for (size_t i = 0; i < selected_indices.size(); ++i)
        {
            size_t idx = selected_indices[i];
            Mx_sub.push_back(Mx[idx]);
            My_sub.push_back(My[idx]);
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
        res = pc::patcaus::patcaus(
            Mx_sub, My_sub, lib_std, pred_std, 
            static_cast<size_t>(std::abs(num_neighbors)),
            static_cast<size_t>(std::abs(zero_tolerance)),
            static_cast<size_t>(std::abs(h)),
            dist_metric, relative, weighted,
            static_cast<size_t>(std::abs(threads)), true);
    }

    // --- Create DataFrame for per-sample causality ---
    const size_t n_samples = res.PatternTypes.size();

    // Allocate vectors
    Rcpp::NumericVector no(n_samples);
    Rcpp::NumericVector positive(n_samples);
    Rcpp::NumericVector negative(n_samples);
    Rcpp::NumericVector dark(n_samples);

    Rcpp::IntegerVector real_index(n_samples);  // original indices (+1 for R)
    Rcpp::CharacterVector pattern_labels(n_samples);

    // Fill values (direct alignment with res)
    for (size_t i = 0; i < n_samples; ++i)
    {
        // Restore original index (+1 for R)
        real_index[i] = static_cast<int>(pred_indices[i] + 1);

        // Directly use i (NOT idx)
        no[i]       = res.NoCausality[i];
        positive[i] = res.PositiveCausality[i];
        negative[i] = res.NegativeCausality[i];
        dark[i]     = res.DarkCausality[i];

        // Pattern type mapping
        switch (res.PatternTypes[i])
        {
            case 0: pattern_labels[i] = "no"; break;
            case 1: pattern_labels[i] = "positive"; break;
            case 2: pattern_labels[i] = "negative"; break;
            case 3: pattern_labels[i] = "dark"; break;
            default: pattern_labels[i] = "unknown"; break;
        }
    }

    // Build DataFrame
    Rcpp::DataFrame causality_df = Rcpp::DataFrame::create(
        Rcpp::Named("index")    = real_index,
        Rcpp::Named("no")       = no,
        Rcpp::Named("positive") = positive,
        Rcpp::Named("negative") = negative,
        Rcpp::Named("dark")     = dark,
        Rcpp::Named("type")     = pattern_labels
    );

    // --- Create summary DataFrame for causal strengths ---

    Rcpp::CharacterVector causal_type = Rcpp::CharacterVector::create("positive", "negative", "dark");
    Rcpp::NumericVector causal_strength = Rcpp::NumericVector::create(res.TotalPos, res.TotalNeg, res.TotalDark);

    Rcpp::DataFrame summary_df = Rcpp::DataFrame::create(
        Rcpp::Named("type") = causal_type,
        Rcpp::Named("strength") = causal_strength
    );

    // --- Return structured results ---

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("causality") = causality_df,
        Rcpp::Named("summary") = summary_df
    );
    out.attr("class") = Rcpp::CharacterVector::create("pc_single");

    return out;
}
