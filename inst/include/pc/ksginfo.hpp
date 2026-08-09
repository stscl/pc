/**********************************************************************
 *  File: ksginfo.hpp
 *
 *  Continuous information theoretic measurements
 *  using k-nearest neighbor estimators.
 *
 *  Algorithms:
 *
 *      Kozachenko–Leonenko entropy estimator
 *      Kraskov–Stögbauer–Grassberger MI estimator
 *
 *  Data layout:
 *      Series = std::vector<double>
 *      Matrix = std::vector<std::vector<double>> // mat[var][obs]
 *
 *  Estimator variants (alg parameter):
 *
 *      Entropy / Joint Entropy / Conditional Entropy
 *      ---------------------------------------------
 *      alg = 0
 *            Kozachenko–Leonenko entropy estimator
 *
 *      alg = 1
 *            Bias corrected KL estimator
 *            (adds +1/k correction)
 *
 *      Mutual Information / Conditional Mutual Information
 *      ---------------------------------------------------
 *      alg = 0
 *            Kraskov–Stögbauer–Grassberger estimator I (KSG1)
 *
 *      alg = 1
 *            Kraskov–Stögbauer–Grassberger estimator II (KSG2)
 *
 *  Distance backend:
 *      pc::distance::distance (Chebyshev metric)
 *
 *  Special functions:
 *      pc::numericutils::digamma
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 **********************************************************************/

#ifndef PC_KSGINFO_HPP
#define PC_KSGINFO_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "pc/distance.hpp"
#include "pc/numericutils.hpp"

namespace pc
{

namespace ksginfo
{

    using Series = std::vector<double>;
    using Matrix = std::vector<std::vector<double>>;

    /***********************************************************
     * Utility: subset matrix
     ***********************************************************/
    inline Matrix subset(
        const Matrix& mat,
        const std::vector<size_t>& vars)
    {
        Matrix out;
        out.reserve(vars.size());

        for (size_t v : vars)
            out.push_back(mat[v]);

        return out;
    }

    /***********************************************************
     * Entropy (Kozachenko–Leonenko)
     ***********************************************************/
    inline double entropy(
        const Series& series,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        const size_t n = series.size();

        auto dist = pc::distance::distance(series);

        double avg = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            auto& row = dist[i];
            row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbor count");

            std::nth_element(
                row.begin(),
                row.begin() + (k - 1),
                row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);
            
            avg += (pc::numericutils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                    ? 0.0 : std::log(eps * 2.0);
        }

        avg /= static_cast<double>(n);

        double H = pc::numericutils::digamma(n)
                 - pc::numericutils::digamma(k)
                 + avg;
        
        if (alg == 1)
            H += 1.0 / k;

        if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
            H /= std::log(base);

        return H;
    }

    /***********************************************************
     * Joint Entropy
     ***********************************************************/
    inline double je(
        const Matrix& mat,
        const std::vector<size_t>& vars,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        Matrix sub = subset(mat,vars);

        const size_t d = sub.size();
        const size_t n = sub[0].size();

        auto dist = pc::distance::distance(sub,"maximum",true,false);

        double avg = 0.0;

        for (size_t i = 0; i < n; ++i)
        {   
            auto& row = dist[i];
            row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbor count");

            std::nth_element(
                row.begin(),
                row.begin() + (k - 1),
                row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg += (pc::numericutils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                    ? 0.0 : std::log(eps * 2.0);
        }

        avg /= static_cast<double>(n);

        double H = pc::numericutils::digamma(n)
                 - pc::numericutils::digamma(k)
                 + d * avg;

        if (alg == 1)
            H += 1.0 / k;

        if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
            H /= std::log(base);

        return H;
    }

    /***********************************************************
     * Conditional Entropy
     ***********************************************************/
    inline double ce(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& cond,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        std::vector<size_t> tc = cond;
        tc.insert(tc.end(),target.begin(),target.end());

        return je(mat,tc,k,alg,base) - je(mat,cond,k,alg,base);
    }

    /***********************************************************
     * Mutual Information (KSG estimator)
     ***********************************************************/
    inline double mi(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false)
    {
        std::vector<size_t> xy = target;
        xy.insert(xy.end(), interact.begin(), interact.end());

        auto d_xy = pc::distance::distance(subset(mat,xy),"maximum",true,false);
        auto d_x  = pc::distance::distance(subset(mat,target),"maximum",true,false);
        auto d_y  = pc::distance::distance(subset(mat,interact),"maximum",true,false);

        const size_t n = d_xy.size();
        const size_t d = xy.size();

        double sum = 0.0;
        double avg_log_eps = 0.0;

        for (size_t i = 0; i < n; ++i)
        {   
            auto& row = d_xy[i];
            row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbor count");

            std::nth_element(row.begin(),row.begin()+k-1,row.end());
            
            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg_log_eps += (pc::numericutils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                            ? 0.0 : std::log(eps * 2.0);

            size_t nx = 0, ny = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i == j) continue;

                if (alg == 0)
                {
                    if (!std::isnan(d_x[i][j]) && d_x[i][j] < eps) nx++;
                    if (!std::isnan(d_y[i][j]) && d_y[i][j] < eps) ny++;
                }
                else 
                {
                    if (!std::isnan(d_x[i][j]) && d_x[i][j] <= eps) nx++;
                    if (!std::isnan(d_y[i][j]) && d_y[i][j] <= eps) ny++;
                } 
            }

            if (alg == 0)
                sum += pc::numericutils::digamma(nx+1)
                     + pc::numericutils::digamma(ny+1);
            else
                sum += pc::numericutils::digamma(nx)
                     + pc::numericutils::digamma(ny);
        }

        avg_log_eps /= n;

        double mival = pc::numericutils::digamma(k)
                     + pc::numericutils::digamma(n)
                     - sum / n;

        if (alg == 1) mival -= 1.0 / k;

        mival = std::max(0.0, mival);

        if (!normalize) 
        {
            if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
                mival /= std::log(base);

            return mival;
        } 

        double hxy = pc::numericutils::digamma(n)
                   - pc::numericutils::digamma(k)
                   + d * avg_log_eps;
        if (alg == 1) hxy += 1.0 / k;

        if (hxy <= 0) {
            if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
                mival /= std::log(base);

            return mival;
        } 

        return mival / hxy;
    }

    /***********************************************************
     * Conditional Mutual Information
     ***********************************************************/
    inline double cmi(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        const std::vector<size_t>& conds,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false)
    {
        std::vector<size_t> xyz = conds;
        xyz.insert(xyz.end(), target.begin(), target.end());
        xyz.insert(xyz.end(), interact.begin(), interact.end());

        std::vector<size_t> xy = target;
        xy.insert(xy.end(), interact.begin(), interact.end());

        std::vector<size_t> xz = conds;
        xz.insert(xz.end(), target.begin(), target.end());

        std::vector<size_t> yz = conds;
        yz.insert(yz.end(), interact.begin(), interact.end());

        auto d_xyz = pc::distance::distance(subset(mat,xyz),"maximum",true,false);
        auto d_xz  = pc::distance::distance(subset(mat,xz),"maximum",true,false);
        auto d_yz  = pc::distance::distance(subset(mat,yz),"maximum",true,false);
        auto d_z   = pc::distance::distance(subset(mat,conds),"maximum",true,false);

        const size_t n = d_xyz.size();
        const size_t d = xy.size();

        double sum = 0.0;
        double avg_log_eps = 0.0;

        for (size_t i = 0; i < n; ++i)
        {   
            auto& row = d_xyz[i];
            row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbor count");

            std::nth_element(row.begin(),row.begin()+k-1,row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg_log_eps += (pc::numericutils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                            ? 0.0 : std::log(eps * 2.0);

            size_t nxz = 0, nyz = 0, nz = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i == j) continue;

                if (alg == 0)
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] < eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] < eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  < eps) nz++;
                }
                else
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] <= eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] <= eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  <= eps) nz++;    
                } 
            }

            if (alg == 0)
                sum += pc::numericutils::digamma(nxz+1)
                     + pc::numericutils::digamma(nyz+1)
                     - pc::numericutils::digamma(nz+1);
            else
                sum += pc::numericutils::digamma(nxz)
                     + pc::numericutils::digamma(nyz)
                     - pc::numericutils::digamma(nz);
        }

        avg_log_eps /= n;

        double cmival = pc::numericutils::digamma(k) - sum / n;
        if (alg == 1) cmival -= 1.0 / k;

        if (!normalize)
        {
            if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
                cmival /= std::log(base);

            return cmival;
        } 

        double hxy_z = pc::numericutils::digamma(n)
                     - pc::numericutils::digamma(k)
                     + d * avg_log_eps;
        if (alg == 1) hxy_z += 1.0 / k;

        if (hxy_z <= 0)
        {
            if (!pc::numericutils::doubleNearlyEqual(base,std::exp(1.0)))
                cmival /= std::log(base);

            return cmival;
        } 

        return cmival / hxy_z;
    }

} // namespace ksginfo

}

#endif // PC_KSGINFO_HPP
