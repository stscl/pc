#include <vector>
#include <cmath>
#include <tuple>
#include <limits>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
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
    // --- Input Conversion and Validation --------------------------------------

    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    std::vector<double> sg = Rcpp::as<std::vector<double>>(source);
    const size_t n_obs = tg.size();

    // Convert library indices (R 1-based → C++ 0-based)
    std::vector<size_t> lib_std = Rcpp::as<std::vector<size_t>>(lib);
    for (auto& idx : lib_std) 
    {
        if (idx < 1 || idx > n_obs) {
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
        if (idx < 1 || idx > n_obs) {
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

    // --- Perform Pattern Causality Analysis -------------------------
    pc::symdync::PatternCausalityRes res = pc::patcaus::patcaus(
        Mx, My, lib_std, pred_std, 
        static_cast<size_t>(std::abs(num_neighbors)),
        static_cast<size_t>(std::abs(zero_tolerance)),
        static_cast<size_t>(std::abs(h)),
        dist_metric, relative, weighted,
        static_cast<size_t>(std::abs(threads)));

    // --- Create DataFrame for per-sample causality ----------------------------

    const size_t n_samples = res.NoCausality.size();
    Rcpp::LogicalVector real_loop(n_samples, false);
    Rcpp::CharacterVector pattern_labels(n_samples, "no");

    for (size_t rl = 0; rl < res.RealLoop.size(); ++rl) 
    {
        size_t idx = res.RealLoop[rl];
        if (idx < n_samples) 
        {
            // Record validated samples
            real_loop[idx] = true;
            // Map pattern_types (0–3) → descriptive string labels
            switch (res.PatternTypes[rl]) 
            {
                case 0: pattern_labels[idx]  = "no"; break;
                case 1: pattern_labels[idx]  = "positive"; break;
                case 2: pattern_labels[idx]  = "negative"; break;
                case 3: pattern_labels[idx]  = "dark"; break;
                default: pattern_labels[idx] = "unknown"; break;
            }
        }
    }

    Rcpp::DataFrame causality_df = Rcpp::DataFrame::create(
        Rcpp::Named("no") = Rcpp::NumericVector(res.NoCausality.begin(), res.NoCausality.end()),
        Rcpp::Named("positive") = Rcpp::NumericVector(res.PositiveCausality.begin(), res.PositiveCausality.end()),
        Rcpp::Named("negative") = Rcpp::NumericVector(res.NegativeCausality.begin(), res.NegativeCausality.end()),
        Rcpp::Named("dark") = Rcpp::NumericVector(res.DarkCausality.begin(), res.DarkCausality.end()),
        Rcpp::Named("type") = pattern_labels,
        Rcpp::Named("valid") = real_loop
    );

    // --- Create summary DataFrame for causal strengths ------------------------

    Rcpp::CharacterVector causal_type = Rcpp::CharacterVector::create("positive", "negative", "dark");
    Rcpp::NumericVector causal_strength = Rcpp::NumericVector::create(res.TotalPos, res.TotalNeg, res.TotalDark);

    Rcpp::DataFrame summary_df = Rcpp::DataFrame::create(
        Rcpp::Named("type") = causal_type,
        Rcpp::Named("strength") = causal_strength
    );

    // --- Return structured results --------------------------------------------

    return Rcpp::List::create(
        Rcpp::Named("causality") = causality_df,
        Rcpp::Named("summary") = summary_df
    );
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
        if (idx < 1 || idx > n_obs) {
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
        if (idx < 1 || idx > n_obs) {
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

    // --- Perform Bootstrapped Pattern Causality Analysis -------------------------
    std::vector<std::vector<std::vector<double>>> res = pc::patcaus::patcaus(
        Mx, My, libsizes_std, lib_std, pred_std, 
        static_cast<size_t>(std::abs(num_neighbors)),
        static_cast<size_t>(std::abs(zero_tolerance)),
        static_cast<size_t>(std::abs(h)), dist_metric, 
        static_cast<size_t>(std::abs(boot)), random_sample, 
        static_cast<unsigned long long>(std::abs(seed)),
        relative, weighted, static_cast<size_t>(std::abs(threads)),
        static_cast<size_t>(std::abs(parallel_level)), verbose);

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

        return Rcpp::DataFrame::create(
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

        return Rcpp::DataFrame::create(
            Rcpp::Named("libsizes") = df_libsizes,
            Rcpp::Named("type") = df_type,
            Rcpp::Named("mean") = df_causality,
            Rcpp::Named("q05") = df_q05,
            Rcpp::Named("q50") = df_q50,
            Rcpp::Named("q95") = df_q95
        );
    }
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
        if (idx < 1 || idx > n_obs) {
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
        if (idx < 1 || idx > n_obs) {
            Rcpp::stop("pred index %d out of bounds [1, %d]",
                       static_cast<int>(idx),
                       static_cast<int>(n_obs));
        }
        idx -= 1;
    }

    // Unique sorted embedding dimensions, neighbor values, and tau values
    std::vector<size_t> Es = Rcpp::as<std::vector<size_t>>(E);
    std::sort(Es.begin(), Es.end());
    Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

    std::vector<size_t> ks = Rcpp::as<std::vector<size_t>>(k);
    std::sort(ks.begin(), ks.end());
    ks.erase(std::unique(ks.begin(), ks.end()), ks.end());

    std::vector<size_t> taus = Rcpp::as<std::vector<size_t>>(tau);
    std::sort(taus.begin(), taus.end());
    taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

    // Generate unique (E, b, tau) combinations
    std::vector<std::tuple<size_t, size_t, size_t>> unique_EkTau;
    for (size_t ee : Es)
        for (size_t kk : ks)
        for (size_t tt : taus)
            unique_EkTau.emplace_back(ee, kk, tt);

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
    
    // --- Perform Pattern Causality Analysis -------------------------
    std::vector<std::vector<double>> result(unique_EkTau.size(), std::vector<double>(6));

    if (parallel_level == 0) {
        for (size_t i = 0; i < unique_EkTau.size(); ++i) {
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
                    tg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    sg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
            } 
            else if (nrows.isNotNull())
            {   
                Mx = pc::embed::embed(
                    tm, Ei, taui, static_cast<size_t>(std::abs(style)));

                My = pc::embed::embed(
                    sm, Ei, taui, static_cast<size_t>(std::abs(style)));
            }
            else  
            {
                Mx = pc::embed::embed(
                    tg, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    sg, Ei, taui, static_cast<size_t>(std::abs(style)));
            }

            pc::symdync::PatternCausalityRes res = pc::patcaus::patcaus(
                Mx, My, lib_std, pred_std, ki,
                static_cast<size_t>(std::abs(zero_tolerance)),
                static_cast<size_t>(std::abs(h)),
                dist_metric, relative, weighted,
                static_cast<size_t>(std::abs(threads)));

            result[i][0] = Ei;
            result[i][1] = ki;
            result[i][2] = taui;
            result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
            result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
            result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
        }
    } else {
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
                    tg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    sg, nb_std, Ei, taui, static_cast<size_t>(std::abs(style)));
            } 
            else if (nrows.isNotNull())
            {   
                Mx = pc::embed::embed(
                    tm, Ei, taui, static_cast<size_t>(std::abs(style)));

                My = pc::embed::embed(
                    sm, Ei, taui, static_cast<size_t>(std::abs(style)));
            }
            else  
            {
                Mx = pc::embed::embed(
                    tg, Ei, taui, static_cast<size_t>(std::abs(style)));
                My = pc::embed::embed(
                    sg, Ei, taui, static_cast<size_t>(std::abs(style)));
            }

            pc::symdync::PatternCausalityRes res = pc::patcaus::patcaus(
                Mx, My, lib_std, pred_std, ki,
                static_cast<size_t>(std::abs(zero_tolerance)),
                static_cast<size_t>(std::abs(h)),
                dist_metric, relative, weighted, 1);

            result[i][0] = Ei;
            result[i][1] = ki;
            result[i][2] = taui;
            result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
            result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
            result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
            }, threads_sizet);
    }

    
    pc::symdync::PatternCausalityRes res = pc::patcaus::patcaus(
        Mx, My, lib_std, pred_std, 
        static_cast<size_t>(std::abs(num_neighbors)),
        static_cast<size_t>(std::abs(zero_tolerance)),
        static_cast<size_t>(std::abs(h)),
        dist_metric, relative, weighted,
        static_cast<size_t>(std::abs(threads)));

    // --- Create DataFrame for per-sample causality ----------------------------

    const size_t n_samples = res.NoCausality.size();
    Rcpp::LogicalVector real_loop(n_samples, false);
    Rcpp::CharacterVector pattern_labels(n_samples, "no");

    for (size_t rl = 0; rl < res.RealLoop.size(); ++rl) 
    {
        size_t idx = res.RealLoop[rl];
        if (idx < n_samples) 
        {
            // Record validated samples
            real_loop[idx] = true;
            // Map pattern_types (0–3) → descriptive string labels
            switch (res.PatternTypes[rl]) 
            {
                case 0: pattern_labels[idx]  = "no"; break;
                case 1: pattern_labels[idx]  = "positive"; break;
                case 2: pattern_labels[idx]  = "negative"; break;
                case 3: pattern_labels[idx]  = "dark"; break;
                default: pattern_labels[idx] = "unknown"; break;
            }
        }
    }

    Rcpp::DataFrame causality_df = Rcpp::DataFrame::create(
        Rcpp::Named("no") = Rcpp::NumericVector(res.NoCausality.begin(), res.NoCausality.end()),
        Rcpp::Named("positive") = Rcpp::NumericVector(res.PositiveCausality.begin(), res.PositiveCausality.end()),
        Rcpp::Named("negative") = Rcpp::NumericVector(res.NegativeCausality.begin(), res.NegativeCausality.end()),
        Rcpp::Named("dark") = Rcpp::NumericVector(res.DarkCausality.begin(), res.DarkCausality.end()),
        Rcpp::Named("type") = pattern_labels,
        Rcpp::Named("valid") = real_loop
    );

    // --- Create summary DataFrame for causal strengths ------------------------

    Rcpp::CharacterVector causal_type = Rcpp::CharacterVector::create("positive", "negative", "dark");
    Rcpp::NumericVector causal_strength = Rcpp::NumericVector::create(res.TotalPos, res.TotalNeg, res.TotalDark);

    Rcpp::DataFrame summary_df = Rcpp::DataFrame::create(
        Rcpp::Named("type") = causal_type,
        Rcpp::Named("strength") = causal_strength
    );

    // --- Return structured results --------------------------------------------

    return Rcpp::List::create(
        Rcpp::Named("causality") = causality_df,
        Rcpp::Named("summary") = summary_df
    );
}
