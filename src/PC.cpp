#include <vector>
#include <cmath>
#include <tuple>
#include <limits>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include "pc.h"

// Wrapper function to perform pattern causality analysis
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppPC(
    const Rcpp::NumericVector& target,
    const Rcpp::NumericVector& source,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int style = 0,
    int num_neighbors = 4,
    int zero_tolerance = 0,
    const std::string& dist_metric = "euclidean",
    bool relative = true,
    bool weighted = true,
    int threads = 1,
    int h = 0,
    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
    Rcpp::Nullable<int> nrows = R_NilValue)
{
    // --- Input Conversion and Validation ---
    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    std::vector<double> sg = Rcpp::as<std::vector<double>>(source);
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

    // Convert Rcpp IntegerVector to std::vector<size_t>
    std::vector<size_t> E_vec = Rcpp::as<std::vector<size_t>>(E);
    std::vector<size_t> tau_vec = Rcpp::as<std::vector<size_t>>(tau);

    // Expand E and tau
    std::vector<size_t> E_std(2);
    std::vector<size_t> tau_std(2);

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
            sg, nb_std, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            tg, nb_std, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));
    } 
    else if (nrows.isNotNull())
    {   
        size_t n_rows = static_cast<size_t>(std::abs(Rcpp::as<int>(nrows)));

        std::vector<std::vector<double>> sm = 
            pc::embed::gridVec2Mat(sg, n_rows);
        Mx = pc::embed::embed(
            sm, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));

        std::vector<std::vector<double>> tm = 
            pc::embed::gridVec2Mat(tg, n_rows);
        My = pc::embed::embed(
            tm, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
    }
    else  
    {
        Mx = pc::embed::embed(
            sg, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            tg, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));

        size_t max_E = *std::max_element(E_std.begin(), E_std.end());
        size_t max_tau = *std::max_element(tau_std.begin(), tau_std.end());
        size_t max_lag = (max_tau == 0) 
            ? (max_E - 1)
            : ((max_E - 1) * max_tau);

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

// Wrapper function to perform bootstrapped pattern causality analysis
// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame RcppPCboot(
    const Rcpp::NumericVector& target,
    const Rcpp::NumericVector& source,
    const Rcpp::IntegerVector& libsizes,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int style = 0,
    int num_neighbors = 4,
    int zero_tolerance = 0,
    const std::string& dist_metric = "euclidean",
    int boot = 99,
    bool random_sample = true,
    int seed = 42,
    bool relative = true,
    bool weighted = true,
    int threads = 1,
    int parallel_level = 0,
    bool verbose = false,
    int h = 0,
    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
    Rcpp::Nullable<int> nrows = R_NilValue)
{
    // --- Input Conversion and Validation --------------------------------------
    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    std::vector<double> sg = Rcpp::as<std::vector<double>>(source);
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

    // Convert Rcpp IntegerVector to std::vector<size_t>
    std::vector<size_t> E_vec = Rcpp::as<std::vector<size_t>>(E);
    std::vector<size_t> tau_vec = Rcpp::as<std::vector<size_t>>(tau);

    // Expand E and tau
    std::vector<size_t> E_std(2);
    std::vector<size_t> tau_std(2);

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

    // --- Embedding Construction ------------------------------------------------
    std::vector<std::vector<double>> Mx;
    std::vector<std::vector<double>> My;

    if (nb.isNotNull()) 
    {
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        std::vector<std::vector<size_t>> nb_std = pc::convert::nb2std(nb.get());
        Mx = pc::embed::embed(
            sg, nb_std, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            tg, nb_std, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));
    } 
    else if (nrows.isNotNull())
    {   
        size_t n_rows = static_cast<size_t>(std::abs(Rcpp::as<int>(nrows)));

        std::vector<std::vector<double>> sm = 
            pc::embed::gridVec2Mat(sg, n_rows);
        Mx = pc::embed::embed(
            sm, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));

        std::vector<std::vector<double>> tm = 
            pc::embed::gridVec2Mat(tg, n_rows);
        My = pc::embed::embed(
            tm, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));
    }
    else  
    {
        Mx = pc::embed::embed(
            sg, E_std[0], tau_std[0], static_cast<size_t>(std::abs(style)));
        My = pc::embed::embed(
            tg, E_std[1], tau_std[1], static_cast<size_t>(std::abs(style)));

        size_t max_E = *std::max_element(E_std.begin(), E_std.end());
        size_t max_tau = *std::max_element(tau_std.begin(), tau_std.end());
        size_t max_lag = (max_tau == 0) 
            ? (max_E - 1)
            : ((max_E - 1) * max_tau);

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

    // Validate and preprocess library sizes
    std::vector<size_t> libsizes_std = Rcpp::as<std::vector<size_t>>(libsizes);
    std::vector<size_t> valid_libsizes;
    valid_libsizes.reserve(libsizes_std.size());
    for (size_t s : libsizes_std) {
        if (s > static_cast<size_t>(std::abs(num_neighbors)) && s <= lib_std.size())
        valid_libsizes.push_back(s);
    }

    std::sort(valid_libsizes.begin(), valid_libsizes.end());
    valid_libsizes.erase(std::unique(valid_libsizes.begin(), valid_libsizes.end()), valid_libsizes.end());

    if (valid_libsizes.empty()) {
        Rcpp::warning("[Warning] No valid libsizes after filtering. Using full library size as fallback.");
        valid_libsizes.push_back(lib_std.size());
    }

    // --- Check if full set is used ---
    bool use_subset = (selected_indices.size() < Mx.size());

    // --- Perform Bootstrapped Pattern Causality Analysis -------------------------
    std::vector<std::vector<std::vector<double>>> res;

    if (!use_subset)
    {
        // --- Full data: no slicing needed ---
        res = pc::patcaus::patcaus(
            Mx, My, libsizes_std, lib_std, pred_std, 
            static_cast<size_t>(std::abs(num_neighbors)),
            static_cast<size_t>(std::abs(zero_tolerance)),
            static_cast<size_t>(std::abs(h)), dist_metric, 
            static_cast<size_t>(std::abs(boot)), random_sample, 
            static_cast<unsigned long long>(std::abs(seed)),
            relative, weighted, static_cast<size_t>(std::abs(threads)),
            static_cast<size_t>(std::abs(parallel_level)), verbose);
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
            Mx_sub, My_sub, libsizes_std, lib_std, pred_std, 
            static_cast<size_t>(std::abs(num_neighbors)),
            static_cast<size_t>(std::abs(zero_tolerance)),
            static_cast<size_t>(std::abs(h)), dist_metric, 
            static_cast<size_t>(std::abs(boot)), random_sample, 
            static_cast<unsigned long long>(std::abs(seed)),
            relative, weighted, static_cast<size_t>(std::abs(threads)),
            static_cast<size_t>(std::abs(parallel_level)), verbose);
    }    

    // --- Result Processing -----------------------------------------------------

    // res structure: [3][valid_libsizes][boot]
    // dimension 0: metric type (0=TotalPos,1=TotalNeg,2=TotalDark)
    // dimension 1: libsizes index
    // dimension 2: bootstrap replicates

    int n_types = 3;
    int n_libsizes = static_cast<int>(valid_libsizes.size());
    int n_boot = static_cast<int>(res[0][0].size());

    // Prepare vectors to hold dataframe columns
    std::vector<size_t> df_libsizes;
    std::vector<std::string> df_type;
    std::vector<double> df_causality;
    std::vector<double> df_q05, df_q50, df_q95;  // For quantiles if boot > 1

    bool has_bootstrap = (n_boot > 1);
    const std::string types[3] = {"positive", "negative", "dark"};
    
    Rcpp::DataFrame causality_out;

    if (!has_bootstrap) 
    {
        // boot == 1, simple long format: columns = libsizes, type, causality
        df_libsizes.reserve(n_types * n_libsizes);
        df_type.reserve(n_types * n_libsizes);
        df_causality.reserve(n_types * n_libsizes);

        for (int t = 0; t < n_types; ++t) 
        {
            for (int l = 0; l < n_libsizes; ++l) 
            {
                df_libsizes.push_back(valid_libsizes[l]);
                df_type.push_back(types[t]);
                if(std::isnan(res[t][l][0])) res[t][l][0] = 0; // replace nan causal strength with 0
                df_causality.push_back(res[t][l][0]);
            }
        }

        causality_out = Rcpp::DataFrame::create(
            Rcpp::Named("libsizes") = df_libsizes,
            Rcpp::Named("type") = df_type,
            Rcpp::Named("causality") = df_causality
        );
    } 
    else 
    {
        // boot > 1, summary with mean and quantiles
        df_libsizes.reserve(n_types * n_libsizes);
        df_type.reserve(n_types * n_libsizes);
        df_causality.reserve(n_types * n_libsizes);
        df_q05.reserve(n_types * n_libsizes);
        df_q50.reserve(n_types * n_libsizes);
        df_q95.reserve(n_types * n_libsizes);

        for (int t = 0; t < n_types; ++t)
        {
            for (int l = 0; l < n_libsizes; ++l) 
            {
                const std::vector<double>& boot_vals = res[t][l];
                double mean_val = pc::numericutils::mean(boot_vals);
                std::vector<double> qs = pc::numericutils::quantile(boot_vals, {0.05, 0.5, 0.95});

                // replace nan causal strength with 0
                if(std::isnan(mean_val)) mean_val = 0;
                for(double& q : qs)
                {
                    if(std::isnan(q)) q = 0;
                }

                df_libsizes.push_back(valid_libsizes[l]);
                df_type.push_back(types[t]);
                df_causality.push_back(mean_val);
                df_q05.push_back(qs[0]);
                df_q50.push_back(qs[1]);
                df_q95.push_back(qs[2]);
            }
        }

        causality_out = Rcpp::DataFrame::create(
            Rcpp::Named("libsizes") = df_libsizes,
            Rcpp::Named("type") = df_type,
            Rcpp::Named("mean") = df_causality,
            Rcpp::Named("q05") = df_q05,
            Rcpp::Named("q50") = df_q50,
            Rcpp::Named("q95") = df_q95
        );
    }

    // --- Return structured results --------------------------------------------

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("causality") = causality_out
    );
    out.attr("class") = Rcpp::CharacterVector::create("pc_boot");

    return out;
}

/**
 * Select optimal embedding parameters (E k tau) for pattern causality metrics
 * using a global scan and full tie tracking across three causality measures.
 *
 * The input matrix must contain six columns in the following order:
 *   1. E      embedding dimension
 *   2. k      number of nearest neighbors
 *   3. tau    lag parameter
 *   4. pos    positive causality score
 *   5. neg    negative causality score
 *   6. dark   dark causality score
 *
 * The argument `maximize` specifies which causality metric should be given
 * primary priority. The three metrics are compared in a hierarchical order:
 *
 *   maximize = "positive":  pos then dark then neg
 *   maximize = "negative":  neg then dark then pos
 *   maximize = "dark":      dark then pos then neg
 *
 * The function performs a single pass global scan. During the scan it collects
 * all rows that achieve the joint optimum across the three prioritized metrics
 * within numerical tolerance. When multiple rows tie for the optimum a final
 * deterministic choice is made by selecting the smallest E then the smallest
 * tau then the smallest k. A warning is issued when this final tie breaking
 * procedure is required.
 *
 * @param Emat A numeric matrix with columns: E k tau pos neg dark.
 * @param maximize A string specifying which metric to prioritize.
 * @return IntegerVector containing E k tau in this order.
 */
Rcpp::IntegerVector OptPCparm(const Rcpp::NumericMatrix& Emat,
                              const std::string& maximize = "positive") {

  if (Emat.ncol() != 6) 
  {
    Rcpp::stop("Input matrix must have exactly six columns: E k tau pos neg dark.");
  }
  int n = Emat.nrow();
  if (n == 0) 
  {
    Rcpp::stop("Input matrix must not be empty.");
  }

  if (maximize != "positive" && maximize != "negative" && maximize != "dark") 
  {
    Rcpp::stop("maximize must be one of positive negative or dark.");
  }

  // establish metric priority order
  std::vector<int> priority(3);
  if (maximize == "positive") 
  {
    priority = {3, 5, 4};  // pos dark neg
  } 
  else if (maximize == "negative") 
  {
    priority = {4, 5, 3};  // neg dark pos
  } 
  else 
  {
    priority = {5, 3, 4};  // dark pos neg
  }

  // helper to get metric value
  auto get_metric = [&](int row, int idx) {
    return Emat(row, idx);
  };

  // record of the best metrics for comparison
  std::vector<double> best_vals(3);
  for (int j = 0; j < 3; ++j) 
  {
    best_vals[j] = get_metric(0, priority[j]);
  }

  std::vector<int> best_rows;
  best_rows.push_back(0);

  // global scan
  for (int i = 1; i < n; ++i) 
  {
    bool better = false;
    bool equal_all = true;

    for (int p = 0; p < 3; ++p) 
    {
      double a = get_metric(i, priority[p]);
      double b = best_vals[p];

      if (!pc::numericutils::doubleNearlyEqual(a, b)) 
      {
        if (a > b) 
        {
          better = true;
        }
        equal_all = false;
        break;
      }
    }

    if (better) 
    {
      best_rows.clear();
      best_rows.push_back(i);
      for (int p = 0; p < 3; ++p) 
      {
        best_vals[p] = get_metric(i, priority[p]);
      }
    }
    else if (equal_all) 
    {
      best_rows.push_back(i);
    }
  }

  // if only one globally optimal row return directly
  if (best_rows.size() == 1) 
  {
    int row = best_rows[0];
    return Rcpp::IntegerVector::create(
      static_cast<int>(Emat(row, 0)),
      static_cast<int>(Emat(row, 1)),
      static_cast<int>(Emat(row, 2))
    );
  }

  // issue tie warning
  Rcpp::warning("Multiple parameter sets share the global optimum. The final choice is determined by smallest E then tau then k.");

  // tie break by E then tau then k
  int best_idx = best_rows[0];
  int bestE = Emat(best_idx, 0);
  int bestTau = Emat(best_idx, 2);
  int bestK = Emat(best_idx, 1);

  for (int idx : best_rows) 
  {
    int E = Emat(idx, 0);
    int tau = Emat(idx, 2);
    int k = Emat(idx, 1);

    bool better =
      (E < bestE) ||
      (E == bestE && tau < bestTau) ||
      (E == bestE && tau == bestTau && k < bestK);

    if (better) 
    {
      best_idx = idx;
      bestE = E;
      bestTau = tau;
      bestK = k;
    }
  }

  return Rcpp::IntegerVector::create(bestE, bestK, bestTau);
}

// Wrapper function to perform optimal parameters search for pattern causality analysis
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppPCops(
    const Rcpp::NumericVector& target,
    const Rcpp::NumericVector& source,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    const Rcpp::IntegerVector& k,
    const std::string& maximize = "positive",
    int style = 0,
    int zero_tolerance = 0,
    const std::string& dist_metric = "euclidean",
    bool relative = true,
    bool weighted = true,
    int threads = 1,
    int parallel_level = 0,
    int h = 0,
    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
    Rcpp::Nullable<int> nrows = R_NilValue)
{
    // --- Input Conversion and Validation --------------------------------------

    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    std::vector<double> sg = Rcpp::as<std::vector<double>>(source);
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

    // Unique sorted embedding dimensions, neighbor values, and tau values
    std::vector<size_t> Es = Rcpp::as<std::vector<size_t>>(E);
    // Make sure each E is greater than 2
    for (auto& singleE : Es) 
    {
        if (singleE < 2) singleE = 2;
    }
    std::sort(Es.begin(), Es.end());
    Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

    std::vector<size_t> ks = Rcpp::as<std::vector<size_t>>(k);
    std::sort(ks.begin(), ks.end());
    ks.erase(std::unique(ks.begin(), ks.end()), ks.end());

    std::vector<size_t> taus = Rcpp::as<std::vector<size_t>>(tau);
    std::sort(taus.begin(), taus.end());
    taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

    // Generate unique (E, k, tau) combinations
    std::vector<std::tuple<size_t, size_t, size_t>> unique_EkTau;
    for (size_t ee : Es)
    {
        for (size_t kk : ks)
        {   
            if (kk < ee) continue;

            for (size_t tt : taus)
                unique_EkTau.emplace_back(ee, kk, tt);

        }
    }

    // Process necessay data
    std::vector<std::vector<size_t>> nb_std;
    std::vector<std::vector<double>> tm;
    std::vector<std::vector<double>> sm;
    if (nb.isNotNull()) 
    {   
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        nb_std = pc::convert::nb2std(nb.get());
    }
    else if (nrows.isNotNull())
    {
        size_t n_rows = static_cast<size_t>(std::abs(Rcpp::as<int>(nrows)));
        tm = pc::embed::gridVec2Mat(tg, n_rows);
        sm = pc::embed::gridVec2Mat(sg, n_rows);
    }
    else  
    {
        size_t max_E = *std::max_element(Es.begin(), Es.end());
        size_t max_tau = *std::max_element(taus.begin(), taus.end());
        size_t max_lag = (max_tau == 0) 
            ? (max_E - 1)
            : ((max_E - 1) * max_tau);

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
    bool use_subset = (selected_indices.size() < tg.size());

    if (use_subset)
    {
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
    }
    
    // --- Perform Pattern Causality Analysis -------------------------
    std::vector<std::vector<double>> result(unique_EkTau.size(), std::vector<double>(6));

    if (parallel_level == 0) 
    {
        for (size_t i = 0; i < unique_EkTau.size(); ++i) 
        {
            const size_t Ei   = std::get<0>(unique_EkTau[i]);
            const size_t ki   = std::get<1>(unique_EkTau[i]);
            const size_t taui = std::get<2>(unique_EkTau[i]);
            // auto [Ei, ki, taui] = unique_EkTau[i]; // C++17 structured binding

            // --- Embedding Construction ------------------------------------------------
            std::vector<std::vector<double>> Mx;
            std::vector<std::vector<double>> My;

            if (nb.isNotNull()) 
            {
                Mx = pc::embed::embed(
                    sg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    tg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
            } 
            else if (nrows.isNotNull())
            {   
                Mx = pc::embed::embed(
                    sm, Ei, taui, static_cast<size_t>(std::abs(style)));

                My = pc::embed::embed(
                    tm, Ei, taui, static_cast<size_t>(std::abs(style)));
            }
            else  
            {
                Mx = pc::embed::embed(
                    sg, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    tg, Ei, taui, static_cast<size_t>(std::abs(style)));
            }

            // --- Perform Pattern Causality Analysis ---
            pc::symdync::PatternCausalityRes res;

            if (!use_subset)
            {
                // --- Full data: no slicing needed ---
                res = pc::patcaus::patcaus(
                    Mx, My, lib_std, pred_std, ki,
                    static_cast<size_t>(std::abs(zero_tolerance)),
                    static_cast<size_t>(std::abs(h)),
                    dist_metric, relative, weighted,
                    static_cast<size_t>(std::abs(threads)), false);
            }
            else
            {
                // --- Subset mode: slice Mx and My ---
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

                // --- Run patcaus on subset ---
                res = pc::patcaus::patcaus(
                    Mx_sub, My_sub, lib_std, pred_std, ki,
                    static_cast<size_t>(std::abs(zero_tolerance)),
                    static_cast<size_t>(std::abs(h)),
                    dist_metric, relative, weighted,
                    static_cast<size_t>(std::abs(threads)), false);
            }

            result[i][0] = Ei;
            result[i][1] = ki;
            result[i][2] = taui;
            result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
            result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
            result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
        }
    } 
    else 
    {
        // Configure threads
        size_t threads_sizet = static_cast<size_t>(std::abs(threads));
        threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

        RcppThread::parallelFor(0, unique_EkTau.size(), [&](size_t i) {
            const size_t Ei   = std::get<0>(unique_EkTau[i]);
            const size_t ki   = std::get<1>(unique_EkTau[i]);
            const size_t taui = std::get<2>(unique_EkTau[i]);
            // auto [Ei, ki, taui] = unique_EkTau[i]; // C++17 structured binding

            // --- Embedding Construction ------------------------------------------------
            std::vector<std::vector<double>> Mx;
            std::vector<std::vector<double>> My;

            if (nb.isNotNull()) 
            {
                Mx = pc::embed::embed(
                    sg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    tg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
            } 
            else if (nrows.isNotNull())
            {   
                Mx = pc::embed::embed(
                    sm, Ei, taui, static_cast<size_t>(std::abs(style)));

                My = pc::embed::embed(
                    tm, Ei, taui, static_cast<size_t>(std::abs(style)));
            }
            else  
            {
                Mx = pc::embed::embed(
                    sg, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    tg, Ei, taui, static_cast<size_t>(std::abs(style)));
            }
            
            pc::symdync::PatternCausalityRes res;

            if (!use_subset)
            {
                // --- Full data: no slicing needed ---
                res = pc::patcaus::patcaus(
                    Mx, My, lib_std, pred_std, ki,
                    static_cast<size_t>(std::abs(zero_tolerance)),
                    static_cast<size_t>(std::abs(h)),
                    dist_metric, relative, weighted, 1, false);
            }
            else
            {
                // --- Subset mode: slice Mx and My ---
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

                // --- Run patcaus on subset ---
                res = pc::patcaus::patcaus(
                    Mx_sub, My_sub, lib_std, pred_std, ki,
                    static_cast<size_t>(std::abs(zero_tolerance)),
                    static_cast<size_t>(std::abs(h)),
                    dist_metric, relative, weighted, 1, false);
            }

            result[i][0] = Ei;
            result[i][1] = ki;
            result[i][2] = taui;
            result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
            result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
            result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
        }, threads_sizet);
    }
    
    // --- Convert performance matrix to r--------------------------
    Rcpp::NumericMatrix pmat = pc::convert::mat_std2r(result, true);

    // --- Select optimal parameters ---------------------
    Rcpp::IntegerVector pvec = OptPCparm(pmat, maximize);

    // --- Proper data.frame construction ----------------
    Rcpp::DataFrame pdf = Rcpp::DataFrame::create(
        Rcpp::Named("E")        = pmat(Rcpp::_, 0),
        Rcpp::Named("k")        = pmat(Rcpp::_, 1),
        Rcpp::Named("tau")      = pmat(Rcpp::_, 2),
        Rcpp::Named("Positive") = pmat(Rcpp::_, 3),
        Rcpp::Named("Negative") = pmat(Rcpp::_, 4),
        Rcpp::Named("Dark")     = pmat(Rcpp::_, 5)
    );

    // --- Return structured results --------------------------------------------

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("param") = pvec,
        Rcpp::Named("xmap") = pdf
    );
    out.attr("class") = Rcpp::CharacterVector::create("pc_ops");

    return out;
}
