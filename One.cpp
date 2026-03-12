#include <vector>
#include <utility>

class Grid
{
public:
    Grid(double xMin, double xMax, int nx,
         double yMin, double yMax, int ny)
    {
        double dx = (xMax - xMin) / nx;
        double dy = (yMax - yMin) / ny;

        for (int i = 0; i <= nx; ++i)
        {
            double x = xMin + i * dx;

            std::vector<std::pair<double, double>> temp_vec;

            for (int j = 0; j <= ny; ++j)
            {
                double y = yMin + j * dy;
                temp_vec.push_back({x, y});
            }

            points_.push_back(temp_vec);
        }
    }

private:
    std::vector<std::vector<std::pair<double, double>>> points_;
};




#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

// ============================================================
//  GLOBAL MODEL INPUTS
//  These are intentionally global / namespace-scoped,
//  because you said you may want to use global parameters.
// ============================================================

namespace Model
{
    // Mean reversion lambda in Cheyette
    inline double lambda = 0.03;

    // Number of Monte Carlo particles
    inline std::size_t numParticles = 20000;

    // Number of Euler substeps per calibration interval
    inline std::size_t subStepsPerInterval = 10;

    // Random seed
    inline unsigned int seed = 42;

    // Initial discount curve P0(T)
    // MUST be set by the user before running.
    inline std::function<double(double)> P0;

    // Initial instantaneous forward f0(t)
    // MUST be set by the user before running.
    inline std::function<double(double)> f0;

    // Underlying swap schedule:
    // swap starts at swapPaymentTimes[0] = T_0
    // and ends at swapPaymentTimes.back() = T_N
    // annuity = sum_{m=0}^{N-1} accruals[m] * P(t, T_{m+1})
    inline std::vector<double> swapPaymentTimes; // size N+1
    inline std::vector<double> accruals;         // size N

    inline void Validate()
    {
        if (!P0)
            throw std::runtime_error("Model::P0 is not set.");
        if (!f0)
            throw std::runtime_error("Model::f0 is not set.");
        if (swapPaymentTimes.size() < 2)
            throw std::runtime_error("Model::swapPaymentTimes must contain at least T0 and T1.");
        if (accruals.size() + 1 != swapPaymentTimes.size())
            throw std::runtime_error("Model::accruals size must be swapPaymentTimes.size() - 1.");
        if (lambda <= 0.0)
            throw std::runtime_error("Model::lambda must be positive.");
        if (numParticles == 0)
            throw std::runtime_error("Model::numParticles must be positive.");
        if (subStepsPerInterval == 0)
            throw std::runtime_error("Model::subStepsPerInterval must be positive.");
    }
}

// ============================================================
//  TYPES
// ============================================================

using GridPoint = std::pair<double, double>;                 // (strike, maturity)
using CalibrationGrid = std::vector<std::vector<GridPoint>>; // one row = one maturity slice
using SurfaceValues = std::vector<std::vector<double>>;      // same shape as grid

struct Particle
{
    double x = 0.0;
    double y = 0.0;
    double intR = 0.0;     // \int_0^t r_s ds
    double discount = 1.0; // exp(-intR)
};

struct SlicePathData
{
    double time = 0.0;
    std::vector<double> strikes;

    std::vector<double> swapRates;  // S_i(t)
    std::vector<double> dSdx;       // \partial_x S_i(t)
    std::vector<double> discounts;  // D_i(t)
};

using LocalVolFunction = std::function<double(double, double)>;
// sigma(t, strike_or_swaprate)

using CombineSliceFunction =
    std::function<std::vector<double>(
        const std::vector<double>& strikes,
        const std::vector<double>& marketPartSlice,
        const SlicePathData& mcData)>;

// ============================================================
//  LINEAR INTERPOLATOR WITH FLAT EXTRAPOLATION
//  This is exactly the smile interpolation piece you mentioned.
// ============================================================

class LinearInterpolator
{
public:
    LinearInterpolator() = default;

    LinearInterpolator(std::vector<double> xs, std::vector<double> ys)
        : xs_(std::move(xs)), ys_(std::move(ys))
    {
        if (xs_.size() != ys_.size())
            throw std::runtime_error("Interpolator: xs and ys size mismatch.");
        if (xs_.empty())
            throw std::runtime_error("Interpolator: empty data.");
        for (std::size_t i = 1; i < xs_.size(); ++i)
        {
            if (!(xs_[i] > xs_[i - 1]))
                throw std::runtime_error("Interpolator: x-grid must be strictly increasing.");
        }
    }

    double operator()(double x) const
    {
        if (x <= xs_.front())
            return ys_.front();
        if (x >= xs_.back())
            return ys_.back();

        auto it = std::lower_bound(xs_.begin(), xs_.end(), x);
        std::size_t right = static_cast<std::size_t>(it - xs_.begin());
        std::size_t left = right - 1;

        const double x0 = xs_[left];
        const double x1 = xs_[right];
        const double y0 = ys_[left];
        const double y1 = ys_[right];

        const double w = (x - x0) / (x1 - x0);
        return y0 + w * (y1 - y0);
    }

private:
    std::vector<double> xs_;
    std::vector<double> ys_;
};

// ============================================================
//  CHEYETTE HELPERS
// ============================================================

double G(double t, double T)
{
    if (T < t)
        throw std::runtime_error("G(t,T): requires T >= t.");

    const double a = Model::lambda;
    return (1.0 - std::exp(-a * (T - t))) / a;
}

double BondPrice(double t, double T, double x, double y)
{
    if (T < t)
        throw std::runtime_error("BondPrice(t,T,...): requires T >= t.");

    const double g = G(t, T);
    const double ratio = Model::P0(T) / Model::P0(t);
    return ratio * std::exp(-g * x - 0.5 * g * g * y);
}

double dBondPrice_dx(double t, double T, double x, double y)
{
    const double p = BondPrice(t, T, x, y);
    return -G(t, T) * p;
}

double ShortRate(double t, double x)
{
    return Model::f0(t) + x;
}

double Annuity(double t, double x, double y)
{
    double A = 0.0;
    const std::size_t n = Model::accruals.size();

    for (std::size_t m = 0; m < n; ++m)
    {
        const double Tpay = Model::swapPaymentTimes[m + 1];
        A += Model::accruals[m] * BondPrice(t, Tpay, x, y);
    }
    return A;
}

double dAnnuity_dx(double t, double x, double y)
{
    double dA = 0.0;
    const std::size_t n = Model::accruals.size();

    for (std::size_t m = 0; m < n; ++m)
    {
        const double Tpay = Model::swapPaymentTimes[m + 1];
        dA += Model::accruals[m] * dBondPrice_dx(t, Tpay, x, y);
    }
    return dA;
}

double SwapRate(double t, double x, double y)
{
    const double T0 = Model::swapPaymentTimes.front();
    const double TN = Model::swapPaymentTimes.back();

    const double Pstart = BondPrice(t, T0, x, y);
    const double Pend   = BondPrice(t, TN, x, y);
    const double A      = Annuity(t, x, y);

    if (std::abs(A) < 1e-14)
        throw std::runtime_error("SwapRate: annuity too close to zero.");

    return (Pstart - Pend) / A;
}

double dSwapRate_dx(double t, double x, double y)
{
    const double T0 = Model::swapPaymentTimes.front();
    const double TN = Model::swapPaymentTimes.back();

    const double Pstart = BondPrice(t, T0, x, y);
    const double Pend   = BondPrice(t, TN, x, y);
    const double A      = Annuity(t, x, y);

    const double dPstart = dBondPrice_dx(t, T0, x, y);
    const double dPend   = dBondPrice_dx(t, TN, x, y);
    const double dA      = dAnnuity_dx(t, x, y);

    if (std::abs(A) < 1e-14)
        throw std::runtime_error("dSwapRate_dx: annuity too close to zero.");

    const double N  = Pstart - Pend;
    const double dN = dPstart - dPend;

    return (dN - (N / A) * dA) / A;
    // same as (dN*A - N*dA)/A^2
    // and since S = N/A, this is (dN - S*dA)/A
}

// ============================================================
//  GRID HELPERS
// ============================================================

std::vector<double> ExtractStrikesFromRow(const std::vector<GridPoint>& row)
{
    if (row.empty())
        throw std::runtime_error("Grid row is empty.");

    std::vector<double> strikes;
    strikes.reserve(row.size());

    const double t = row.front().second;
    for (const auto& pt : row)
    {
        if (std::abs(pt.second - t) > 1e-12)
            throw std::runtime_error("Each inner grid row must have a single common maturity.");
        strikes.push_back(pt.first);
    }

    return strikes;
}

std::vector<double> ExtractMaturities(const CalibrationGrid& grid)
{
    if (grid.empty())
        throw std::runtime_error("Calibration grid is empty.");

    std::vector<double> mats;
    mats.reserve(grid.size());

    for (const auto& row : grid)
    {
        if (row.empty())
            throw std::runtime_error("Calibration grid contains an empty row.");
        mats.push_back(row.front().second);
    }

    for (std::size_t i = 1; i < mats.size(); ++i)
    {
        if (!(mats[i] > mats[i - 1]))
            throw std::runtime_error("Grid maturities must be strictly increasing row by row.");
    }

    return mats;
}

// ============================================================
//  PARTICLE PROPAGATION
// ============================================================

void PropagateParticles(
    std::vector<Particle>& particles,
    double t0,
    double t1,
    const LinearInterpolator& smileAtLeftEndpoint,
    std::mt19937_64& rng)
{
    if (t1 <= t0)
        throw std::runtime_error("PropagateParticles: need t1 > t0.");

    std::normal_distribution<double> normal(0.0, 1.0);

    const std::size_t nSteps = Model::subStepsPerInterval;
    const double dt = (t1 - t0) / static_cast<double>(nSteps);
    const double sqrtDt = std::sqrt(dt);

    for (std::size_t step = 0; step < nSteps; ++step)
    {
        const double t = t0 + static_cast<double>(step) * dt;

        for (auto& p : particles)
        {
            const double s    = SwapRate(t, p.x, p.y);
            const double dsdx = dSwapRate_dx(t, p.x, p.y);

            if (std::abs(dsdx) < 1e-12)
                throw std::runtime_error("Propagation failed: dS/dx too close to zero.");

            const double sigmaLocal = smileAtLeftEndpoint(s);
            const double sigmaHat   = sigmaLocal / dsdx;

            const double z = normal(rng);
            const double dW = sqrtDt * z;

            const double r = ShortRate(t, p.x);

            // Euler step
            const double xOld = p.x;
            const double yOld = p.y;

            p.x = xOld + (yOld - Model::lambda * xOld) * dt + sigmaHat * dW;
            p.y = yOld + (sigmaHat * sigmaHat - 2.0 * Model::lambda * yOld) * dt;

            // Discount accumulation
            p.intR += r * dt;
            p.discount = std::exp(-p.intR);
        }
    }
}

// ============================================================
//  BUILD MONTE CARLO DATA ON A SLICE
// ============================================================

SlicePathData BuildSlicePathData(
    const std::vector<Particle>& particles,
    double t,
    const std::vector<double>& strikes)
{
    SlicePathData out;
    out.time = t;
    out.strikes = strikes;

    out.swapRates.reserve(particles.size());
    out.dSdx.reserve(particles.size());
    out.discounts.reserve(particles.size());

    for (const auto& p : particles)
    {
        const double s    = SwapRate(t, p.x, p.y);
        const double dsdx = dSwapRate_dx(t, p.x, p.y);

        out.swapRates.push_back(s);
        out.dSdx.push_back(dsdx);
        out.discounts.push_back(p.discount);
    }

    return out;
}

// ============================================================
//  MAIN ENGINE
//
//  initialSmileRow0:
//      values sigma(t0, K_j) on the first grid row
//
//  marketPart:
//      same shape as grid
//      already computed externally
//
//  combineSlice:
//      user-supplied function which mixes market part with MC part
//      and returns the calibrated smile at that maturity
// ============================================================

SurfaceValues CalibrateSurfaceMonteCarlo(
    const CalibrationGrid& grid,
    const SurfaceValues& marketPart,
    const std::vector<double>& initialSmileRow0,
    const CombineSliceFunction& combineSlice)
{
    Model::Validate();

    if (grid.size() != marketPart.size())
        throw std::runtime_error("marketPart must have same number of rows as grid.");

    if (grid.empty())
        throw std::runtime_error("grid is empty.");

    const std::vector<double> maturities = ExtractMaturities(grid);
    const std::size_t nSlices = grid.size();

    SurfaceValues calibrated(nSlices);

    // Validate row sizes
    for (std::size_t n = 0; n < nSlices; ++n)
    {
        if (grid[n].size() != marketPart[n].size())
            throw std::runtime_error("marketPart row size mismatch.");
    }

    // First row is taken as already known / initialised
    {
        const std::vector<double> strikes0 = ExtractStrikesFromRow(grid[0]);
        if (initialSmileRow0.size() != strikes0.size())
            throw std::runtime_error("initialSmileRow0 size mismatch with first grid row.");

        calibrated[0] = initialSmileRow0;
    }

    // Initialize particles at x0 = y0 = 0, D0 = 1
    std::vector<Particle> particles(Model::numParticles);

    std::mt19937_64 rng(Model::seed);

    // Sequential calibration in time
    for (std::size_t n = 1; n < nSlices; ++n)
    {
        const double tLeft  = maturities[n - 1];
        const double tRight = maturities[n];

        const std::vector<double> strikesLeft  = ExtractStrikesFromRow(grid[n - 1]);
        const std::vector<double> strikesRight = ExtractStrikesFromRow(grid[n]);

        // Build smile interpolator from already calibrated previous slice
        LinearInterpolator smileLeft(strikesLeft, calibrated[n - 1]);

        // Propagate particles from t_{n-1} to t_n using sigma(t_{n-1}, .)
        PropagateParticles(particles, tLeft, tRight, smileLeft, rng);

        // Build MC data at time t_n
        SlicePathData mcData = BuildSlicePathData(particles, tRight, strikesRight);

        // Combine market part + MC data to obtain sigma(t_n, K_j)
        calibrated[n] = combineSlice(strikesRight, marketPart[n], mcData);

        if (calibrated[n].size() != strikesRight.size())
            throw std::runtime_error("combineSlice returned wrong number of values.");
    }

    return calibrated;
}







std::vector<double> MyCombineSlice(
    const std::vector<double>& strikes,
    const std::vector<double>& marketPartSlice,
    const SlicePathData& mcData)
{
    std::vector<double> sigma(strikes.size(), 0.0);

    // Example only:
    // here you insert the exact formula from the paper / your derivation

    for (std::size_t j = 0; j < strikes.size(); ++j)
    {
        const double K = strikes[j];

        // --------------------------------------------
        // Placeholder Monte Carlo estimator
        // Replace this by your exact MC term formula
        // --------------------------------------------
        double mcTerm = 0.0;

        const double eps = 1e-3; // purely placeholder
        for (std::size_t i = 0; i < mcData.swapRates.size(); ++i)
        {
            if (std::abs(mcData.swapRates[i] - K) < eps)
            {
                mcTerm += mcData.discounts[i];
            }
        }

        mcTerm /= static_cast<double>(mcData.swapRates.size());
        mcTerm = std::max(mcTerm, 1e-12);

        // Example combination:
        sigma[j] = std::sqrt(std::max(marketPartSlice[j] / mcTerm, 0.0));
    }

    return sigma;
}

\\test

int main()
{
    // --------------------------------------------------
    // 1) Fill global model inputs
    // --------------------------------------------------
    Model::lambda = 0.05;
    Model::numParticles = 50000;
    Model::subStepsPerInterval = 20;
    Model::seed = 1234;

    Model::P0 = [](double T)
    {
        // Example flat curve only
        const double r = 0.02;
        return std::exp(-r * T);
    };

    Model::f0 = [](double t)
    {
        (void)t;
        return 0.02; // Example flat instantaneous forward
    };

    // Example swap schedule T0, T1, ..., TN
    Model::swapPaymentTimes = {1.0, 1.5, 2.0, 2.5, 3.0};
    Model::accruals = {0.5, 0.5, 0.5, 0.5};

    // --------------------------------------------------
    // 2) Build grid: one row = one maturity
    //    pair = (strike, maturity)
    // --------------------------------------------------
    CalibrationGrid grid =
    {
        { {0.0100, 0.25}, {0.0150, 0.25}, {0.0200, 0.25} },
        { {0.0100, 0.50}, {0.0150, 0.50}, {0.0200, 0.50} },
        { {0.0100, 1.00}, {0.0150, 1.00}, {0.0200, 1.00} }
    };

    // --------------------------------------------------
    // 3) Market part already computed elsewhere
    //    same shape as grid
    // --------------------------------------------------
    SurfaceValues marketPart =
    {
        {1.0, 1.0, 1.0},
        {1.1, 1.1, 1.1},
        {1.2, 1.2, 1.2}
    };

    // --------------------------------------------------
    // 4) Initial smile on first row
    // --------------------------------------------------
    std::vector<double> initialSmileRow0 = {0.20, 0.19, 0.21};

    // --------------------------------------------------
    // 5) Run
    // --------------------------------------------------
    SurfaceValues calibrated =
        CalibrateSurfaceMonteCarlo(
            grid,
            marketPart,
            initialSmileRow0,
            MyCombineSlice);

    // Print result
    for (std::size_t n = 0; n < calibrated.size(); ++n)
    {
        std::cout << "Slice " << n << ":\n";
        for (double v : calibrated[n])
            std::cout << v << " ";
        std::cout << "\n";
    }

    return 0;
}
