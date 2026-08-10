#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include "pc.h"

// Wrapper function to perform delayed mutual information analysis
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppDMI(
    const Rcpp::NumericVector& target,
    const Rcpp::NumericVector& tau,
    const Rcpp::IntegerVector& pred,
    int k = 3,
    int alg = 0,
    double base = 2.0,
    bool normalize = false,
    int threads = 1)
{
    // --- Input Conversion and Validation ---
    std::vector<double> tg = Rcpp::as<std::vector<double>>(target);
    const size_t n_obs = tg.size();

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

    // Construct time delay step tau
    std::vector<size_t> tau_std = Rcpp::as<std::vector<size_t>>(tau);
    size_t max_tau = static_cast<size_t>(*std::max_element(tau_std.begin(), tau_std.end()));

    // ---- sort predict indices ----
    size_t max_lag = (tau == 0) 
            ? (max_E - 1)
            : ((max_E - 1) * static_cast<size_t>(std::abs(tau)));

    pred_std.erase(
        std::remove_if(pred_std.begin(), pred_std.end(), 
            [&](size_t idx){ return idx + 1 < max_tau; }),
        pred_std.end()
    );

    std::sort(pred_std.begin(), pred_std.end());
    pred_std.erase(
        std::unique(pred_std.begin(), pred_std.end()),
        pred_std.end()
    );

    // ---- filter pred (remove NaN in target/source) ----
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

    // --- Perform Delay Mutual Information Analysis ---
    std::vector<double> res = pc::dmi::dmi(
            tg, tau_std, pred_std, 
            static_cast<size_t>(std::abs(k)), 
            static_cast<size_t>(std::abs(alg)), 
            base, normalize,
            static_cast<size_t>(std::abs(threads)));

    // Convert the result back to Rcpp::NumericVector and set names as "E:1", "E:2", ..., "E:n"
    Rcpp::NumericVector result = Rcpp::wrap(res);
    Rcpp::CharacterVector resnames(result.size());
    for (int i = 0; i < result.size(); ++i) {
        resnames[i] = "E:" + std::to_string(i + 1);
    }
    result.names() = resnames;

    // Terminal-friendly hint (one-time, non-intrusive)
    Rcpp::Rcout << "[fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.\n"
                << "[fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.\n"
                << "[fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.\n";

    return result;
}
