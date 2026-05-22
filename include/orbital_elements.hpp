#pragma once

/**
 * @file orbital_elements.hpp
 * @brief Header-only C++ library for orbital elements conversions.
 *
 * Supported element sets:
 *   - Classical Keplerian       (COE)  : {a, e, i, RAAN, omega, nu}
 *   - Equinoctial               (EQN)  : {a, h, k, p, q, lambda}
 *   - Modified Equinoctial      (MEE)  : {p, f, g, h, k, L}
 *   - Non-singular Keplerian    (NSK)  : {a, u, ex, ey, ix, iy}  (Broucke–Cefola)
 *   - Delaunay                  (DEL)  : {l, g, h, L, G, H}
 *   - Poincaré                  (POI)  : {lambda, L, xi, eta, p, q}
 *   - Cartesian state vector    (CART) : {x, y, z, vx, vy, vz}
 *
 * TLE elements (NORAD) are defined as a data-only struct; their conversion
 * requires the SGP4 propagator which depends on time and atmospheric models
 * and is therefore out of scope.  A helper to_classical_approximate() gives
 * the mean Keplerian elements embedded in a TLE.
 *
 * All angles in RADIANS unless stated otherwise.
 * Gravitational parameter mu = GM  [m³/s²] (SI).  Pass it explicitly so the
 * library works for any central body.
 *
 * Conventions
 * -----------
 *   - Right-hand inertial frame (ECI for Earth, MCI for Mars, …)
 *   - Elliptic orbits only (0 ≤ e < 1) for most conversions.
 *   - Equinoctial and Modified Equinoctial are valid for e < 1, i ≠ 180°.
 *   - Non-singular Keplerian (NSK): valid for all i, all e < 1.
 *   - Delaunay / Poincaré: valid away from e=0, i=0 (same singularities as COE).
 *
 * Usage example
 * -------------
 * @code
 *   #include "orbital_elements.hpp"
 *   using namespace orb;
 *
 *   constexpr double MU_EARTH = 3.986004418e14; // m³/s²
 *
 *   COE<double> coe{7e6, 0.01, 0.5, 1.0, 0.3, 0.0};
 *   auto cart  = coe_to_cart(coe, MU_EARTH);
 *   auto mee   = coe_to_mee(coe);
 *   auto del   = coe_to_delaunay(coe, MU_EARTH);
 *   auto poi   = coe_to_poincare(coe, MU_EARTH);
 *   auto coe2  = mee_to_coe(mee);
 * @endcode
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#define ORBIT_VERSION "0.1.0"
 
namespace orb {

// ============================================================================
// Constants
// ============================================================================

template<typename T> constexpr T PI    = T(3.141592653589793238462643383279502884L);
template<typename T> constexpr T TWO_PI= T(2) * PI<T>;

// ============================================================================
// Utility helpers
// ============================================================================

/// Wrap angle to [0, 2π)
template<typename T>
[[nodiscard]] T wrap_0_2pi(T angle) noexcept {
    angle = std::fmod(angle, TWO_PI<T>);
    if (angle < T(0)) angle += TWO_PI<T>;
    return angle;
}

/// Wrap angle to (-π, π]
template<typename T>
[[nodiscard]] T wrap_neg_pi_pi(T angle) noexcept {
    angle = std::fmod(angle + PI<T>, TWO_PI<T>);
    if (angle < T(0)) angle += TWO_PI<T>;
    return angle - PI<T>;
}

/// Solve Kepler's equation  M = E - e·sin(E)  by Newton–Raphson.
template<typename T>
[[nodiscard]] T mean_to_eccentric_anomaly(T M, T e,
                                           int max_iter = 100,
                                           T   tol      = T(1e-12)) {
    M = wrap_0_2pi(M);
    T E = (e < T(0.8)) ? M : PI<T>;
    for (int i = 0; i < max_iter; ++i) {
        T dE = (M - E + e * std::sin(E)) / (T(1) - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < tol) break;
    }
    return wrap_0_2pi(E);
}

/// Eccentric anomaly → true anomaly
template<typename T>
[[nodiscard]] T eccentric_to_true_anomaly(T E, T e) noexcept {
    return T(2) * std::atan2(std::sqrt(T(1) + e) * std::sin(E / T(2)),
                             std::sqrt(T(1) - e) * std::cos(E / T(2)));
}

/// True anomaly → mean anomaly
template<typename T>
[[nodiscard]] T true_to_mean_anomaly(T nu, T e) noexcept {
    T E = T(2) * std::atan2(std::sqrt(T(1) - e) * std::sin(nu / T(2)),
                             std::sqrt(T(1) + e) * std::cos(nu / T(2)));
    return wrap_0_2pi(E - e * std::sin(E));
}

/// True anomaly → eccentric anomaly
template<typename T>
[[nodiscard]] T true_to_eccentric_anomaly(T nu, T e) noexcept {
    return T(2) * std::atan2(std::sqrt(T(1) - e) * std::sin(nu / T(2)),
                             std::sqrt(T(1) + e) * std::cos(nu / T(2)));
}

// ============================================================================
// Element set structures
// ============================================================================

// ---------------------------------------------------------------------------
/// Classical Orbital Elements  (Keplerian COE)
// ---------------------------------------------------------------------------
template<typename T>
struct COE {
    T a;     ///< Semi-major axis           [m]
    T e;     ///< Eccentricity              [-]
    T i;     ///< Inclination               [rad]
    T RAAN;  ///< Right ascension of the AN [rad]
    T omega; ///< Argument of periapsis     [rad]
    T nu;    ///< True anomaly              [rad]
};

// ---------------------------------------------------------------------------
/// Equinoctial Elements
// ---------------------------------------------------------------------------
template<typename T>
struct EQN {
    T a;      ///< Semi-major axis                [m]
    T h;      ///< e·sin(RAAN+omega)              [-]
    T k;      ///< e·cos(RAAN+omega)              [-]
    T p;      ///< tan(i/2)·sin(RAAN)             [-]
    T q;      ///< tan(i/2)·cos(RAAN)             [-]
    T lambda; ///< Mean longitude  RAAN+omega+M   [rad]
};

// ---------------------------------------------------------------------------
/// Modified Equinoctial Elements  (Walker / MEE)
// ---------------------------------------------------------------------------
template<typename T>
struct MEE {
    T p;  ///< Semi-latus rectum  a(1-e²)         [m]
    T f;  ///< e·cos(RAAN+omega)                  [-]
    T g;  ///< e·sin(RAAN+omega)                  [-]
    T h;  ///< tan(i/2)·cos(RAAN)                 [-]
    T k;  ///< tan(i/2)·sin(RAAN)                 [-]
    T L;  ///< True longitude  RAAN+omega+nu      [rad]
};

// ---------------------------------------------------------------------------
/// Non-singular Keplerian  (Broucke–Cefola)
// ---------------------------------------------------------------------------
template<typename T>
struct NSK {
    T a;   ///< Semi-major axis                   [m]
    T u;   ///< Argument of latitude  omega+nu    [rad]
    T ex;  ///< e·cos(RAAN+omega)                 [-]
    T ey;  ///< e·sin(RAAN+omega)                 [-]
    T ix;  ///< tan(i/2)·cos(RAAN)                [-]
    T iy;  ///< tan(i/2)·sin(RAAN)                [-]
};

// ---------------------------------------------------------------------------
/// Delaunay Elements  (canonical action–angle)
// ---------------------------------------------------------------------------
template<typename T>
struct Delaunay {
    T l;  ///< Mean anomaly  M                    [rad]
    T g;  ///< Argument of periapsis  omega       [rad]
    T h;  ///< RAAN  Ω                            [rad]
    T L;  ///< √(μ·a)                             [m^(3/2)·s^(-1)]  or [m²/s] dep. on μ units
    T G;  ///< √(μ·a·(1-e²))  = L√(1-e²)         [same]
    T H;  ///< G·cos(i)                           [same]
};

// ---------------------------------------------------------------------------
/// Poincaré Elements  (regularised Delaunay)
// ---------------------------------------------------------------------------
template<typename T>
struct Poincare {
    T lambda; ///< Mean longitude  l+g+h           [rad]
    T L;      ///< √(μ·a)                          [m^(3/2)·s^(-1)]
    T xi;     ///< √(2(L−G))·cos(g+h)             [-]
    T eta;    ///< −√(2(L−G))·sin(g+h)            [-]
    T p;      ///< √(2(G−H))·cos(h)               [-]
    T q;      ///< −√(2(G−H))·sin(h)              [-]
};

// ---------------------------------------------------------------------------
/// Cartesian state vector  (position + velocity)
// ---------------------------------------------------------------------------
template<typename T>
struct Cart {
    T x, y, z;        ///< Position  [m]
    T vx, vy, vz;     ///< Velocity  [m/s]
};

// ---------------------------------------------------------------------------
/// TLE (Two-Line Element Set, NORAD format)  — data container only.
///
/// TLE fields are intentionally NOT templated: the format is a fixed NORAD
/// standard with bounded ASCII precision (e.g. eccentricity at 7 sig. figs),
/// and all reference SGP4 implementations (Vallado, Kelso) use double
/// throughout.  Templatising would imply a portability that the format itself
/// does not offer.
///
/// IMPORTANT: these are Brouwer *mean* elements in the TEME frame, valid only
/// as input to SGP4/SDP4.  They are NOT osculating Keplerian elements.
/// tle_to_classical_approx() extracts a rough COE for quick inspection only.
// ---------------------------------------------------------------------------
struct TLE {
    // Line 1
    int    sat_number;       ///< NORAD catalog number
    char   classification;   ///< 'U', 'C', or 'S'
    int    epoch_year;       ///< 2-digit year (e.g. 24 for 2024)
    double epoch_day;        ///< Day-of-year + fractional day
    double ndot;             ///< dn/dt / 2       [rev/day²]
    double nddot;            ///< d²n/dt² / 6     [rev/day³]
    double bstar;            ///< B* drag term     [1/earth_radii]
    // Line 2
    double i_deg;            ///< Inclination      [deg]
    double RAAN_deg;         ///< RAAN             [deg]
    double e;                ///< Eccentricity     [-]
    double omega_deg;        ///< Arg. of periapsis [deg]
    double M_deg;            ///< Mean anomaly      [deg]
    double n;                ///< Mean motion       [rev/day]
    int    rev_number;       ///< Revolution number at epoch
};

// ============================================================================
// COE  ↔  Cartesian
// ============================================================================

/// Classical Keplerian  →  Cartesian state vector
template<typename T>
[[nodiscard]] Cart<T> coe_to_cart(const COE<T>& coe, T mu) {
    const T a = coe.a, e = coe.e, i = coe.i;
    const T W = coe.RAAN, w = coe.omega, nu = coe.nu;

    const T p  = a * (T(1) - e * e);
    const T r  = p / (T(1) + e * std::cos(nu));

    // Perifocal frame
    const T x_pf  = r * std::cos(nu);
    const T y_pf  = r * std::sin(nu);
    const T h_mag = std::sqrt(mu * p);
    const T vx_pf = -(mu / h_mag) * std::sin(nu);
    const T vy_pf =  (mu / h_mag) * (e + std::cos(nu));

    // Rotation matrix  Q = R3(-W)·R1(-i)·R3(-w)
    const T cW = std::cos(W), sW = std::sin(W);
    const T ci = std::cos(i), si = std::sin(i);
    const T cw = std::cos(w), sw = std::sin(w);

    // Columns of Q (row-major: Q[row][col])
    const T Q00 =  cW*cw - sW*sw*ci;
    const T Q10 =  sW*cw + cW*sw*ci;
    const T Q20 =  sw*si;
    const T Q01 = -cW*sw - sW*cw*ci;
    const T Q11 = -sW*sw + cW*cw*ci;
    const T Q21 =  cw*si;

    return {
        Q00*x_pf + Q01*y_pf,
        Q10*x_pf + Q11*y_pf,
        Q20*x_pf + Q21*y_pf,
        Q00*vx_pf + Q01*vy_pf,
        Q10*vx_pf + Q11*vy_pf,
        Q20*vx_pf + Q21*vy_pf
    };
}

/// Cartesian state vector  →  Classical Keplerian
template<typename T>
[[nodiscard]] COE<T> cart_to_coe(const Cart<T>& cart, T mu) {
    const T rx = cart.x,  ry = cart.y,  rz = cart.z;
    const T vx = cart.vx, vy = cart.vy, vz = cart.vz;

    const T r_mag = std::sqrt(rx*rx + ry*ry + rz*rz);
    const T v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);

    // Angular momentum  h = r × v
    const T hx = ry*vz - rz*vy;
    const T hy = rz*vx - rx*vz;
    const T hz = rx*vy - ry*vx;
    const T h_mag = std::sqrt(hx*hx + hy*hy + hz*hz);

    // Node vector  N = K × h
    const T nx = -hy, ny = hx; // nz = 0
    const T n_mag = std::sqrt(nx*nx + ny*ny);

    // Eccentricity vector  e_vec = (v²-μ/r)·r/μ − (r·v)·v/μ
    const T rv = rx*vx + ry*vy + rz*vz;
    const T ex = ((v_mag*v_mag - mu/r_mag)*rx - rv*vx) / mu;
    const T ey = ((v_mag*v_mag - mu/r_mag)*ry - rv*vy) / mu;
    const T ez = ((v_mag*v_mag - mu/r_mag)*rz - rv*vz) / mu;
    const T e  = std::sqrt(ex*ex + ey*ey + ez*ez);

    const T E_mec = v_mag*v_mag / T(2) - mu / r_mag;
    const T a     = -mu / (T(2) * E_mec);
    const T inc   = std::acos(std::clamp(hz / h_mag, T(-1), T(1)));

    // RAAN
    T RAAN = (n_mag > T(1e-12)) ? std::acos(std::clamp(nx / n_mag, T(-1), T(1))) : T(0);
    if (ny < T(0)) RAAN = TWO_PI<T> - RAAN;

    // Argument of periapsis
    T omega = T(0);
    if (n_mag > T(1e-12) && e > T(1e-12)) {
        omega = std::acos(std::clamp((nx*ex + ny*ey) / (n_mag * e), T(-1), T(1)));
        if (ez < T(0)) omega = TWO_PI<T> - omega;
    }

    // True anomaly
    T nu = T(0);
    if (e > T(1e-12)) {
        nu = std::acos(std::clamp((ex*rx + ey*ry + ez*rz) / (e * r_mag), T(-1), T(1)));
        if (rv < T(0)) nu = TWO_PI<T> - nu;
    } else {
        // Circular: use argument of latitude
        if (n_mag > T(1e-12))
            nu = std::acos(std::clamp((nx*rx + ny*ry) / (n_mag * r_mag), T(-1), T(1)));
        if (rz < T(0)) nu = TWO_PI<T> - nu;
    }

    return {a, e, inc, RAAN, omega, nu};
}

// ============================================================================
// COE  ↔  EQN  (Equinoctial)
// ============================================================================

template<typename T>
[[nodiscard]] EQN<T> coe_to_eqn(const COE<T>& coe) {
    const T w_plus_W = coe.omega + coe.RAAN;
    const T M = true_to_mean_anomaly(coe.nu, coe.e);
    return {
        coe.a,
        coe.e * std::sin(w_plus_W),
        coe.e * std::cos(w_plus_W),
        std::tan(coe.i / T(2)) * std::sin(coe.RAAN),
        std::tan(coe.i / T(2)) * std::cos(coe.RAAN),
        wrap_0_2pi(coe.RAAN + coe.omega + M)
    };
}

template<typename T>
[[nodiscard]] COE<T> eqn_to_coe(const EQN<T>& eqn) {
    const T e     = std::sqrt(eqn.h * eqn.h + eqn.k * eqn.k);
    const T i     = T(2) * std::atan2(std::sqrt(eqn.p * eqn.p + eqn.q * eqn.q), T(1));
    const T RAAN  = std::atan2(eqn.p, eqn.q);
    const T omega = wrap_0_2pi(std::atan2(eqn.h, eqn.k) - RAAN);
    const T M     = wrap_0_2pi(eqn.lambda - std::atan2(eqn.h, eqn.k));
    const T E     = mean_to_eccentric_anomaly(M, e);
    const T nu    = eccentric_to_true_anomaly(E, e);
    return {eqn.a, e, i, wrap_0_2pi(RAAN), omega, wrap_0_2pi(nu)};
}

// ============================================================================
// COE  ↔  MEE  (Modified Equinoctial)
// ============================================================================

template<typename T>
[[nodiscard]] MEE<T> coe_to_mee(const COE<T>& coe) {
    const T w_plus_W = coe.omega + coe.RAAN;
    return {
        coe.a * (T(1) - coe.e * coe.e),         // p  = semi-latus rectum
        coe.e * std::cos(w_plus_W),              // f
        coe.e * std::sin(w_plus_W),              // g
        std::tan(coe.i / T(2)) * std::cos(coe.RAAN), // h
        std::tan(coe.i / T(2)) * std::sin(coe.RAAN), // k
        wrap_0_2pi(w_plus_W + coe.nu)            // L = true longitude
    };
}

template<typename T>
[[nodiscard]] COE<T> mee_to_coe(const MEE<T>& mee) {
    const T e     = std::sqrt(mee.f * mee.f + mee.g * mee.g);
    const T i     = T(2) * std::atan(std::sqrt(mee.h * mee.h + mee.k * mee.k));
    const T RAAN  = std::atan2(mee.k, mee.h);
    const T omega = wrap_0_2pi(std::atan2(mee.g, mee.f) - RAAN);
    const T nu    = wrap_0_2pi(mee.L - std::atan2(mee.g, mee.f));
    const T a     = (e < T(1) - T(1e-12)) ? mee.p / (T(1) - e * e)
                                           : throw std::domain_error("MEE: e≥1 not supported");
    return {a, e, i, wrap_0_2pi(RAAN), omega, nu};
}

/// Cartesian  →  MEE  (direct, via perifocal decomposition)
template<typename T>
[[nodiscard]] MEE<T> cart_to_mee(const Cart<T>& cart, T mu) {
    return coe_to_mee(cart_to_coe(cart, mu));
}

template<typename T>
[[nodiscard]] Cart<T> mee_to_cart(const MEE<T>& mee, T mu) {
    return coe_to_cart(mee_to_coe(mee), mu);
}

// ============================================================================
// COE  ↔  NSK  (Non-Singular Keplerian, Broucke–Cefola)
// ============================================================================

template<typename T>
[[nodiscard]] NSK<T> coe_to_nsk(const COE<T>& coe) {
    const T w_plus_W = coe.omega + coe.RAAN;
    return {
        coe.a,
        wrap_0_2pi(coe.omega + coe.nu),            // u = arg. of latitude
        coe.e * std::cos(w_plus_W),                // ex
        coe.e * std::sin(w_plus_W),                // ey
        std::tan(coe.i / T(2)) * std::cos(coe.RAAN), // ix
        std::tan(coe.i / T(2)) * std::sin(coe.RAAN)  // iy
    };
}

template<typename T>
[[nodiscard]] COE<T> nsk_to_coe(const NSK<T>& nsk) {
    const T e     = std::sqrt(nsk.ex * nsk.ex + nsk.ey * nsk.ey);
    const T i     = T(2) * std::atan(std::sqrt(nsk.ix * nsk.ix + nsk.iy * nsk.iy));
    const T RAAN  = std::atan2(nsk.iy, nsk.ix);
    const T omega = wrap_0_2pi(std::atan2(nsk.ey, nsk.ex) - RAAN);
    const T nu    = wrap_0_2pi(nsk.u - omega);
    return {nsk.a, e, i, wrap_0_2pi(RAAN), omega, nu};
}

// ============================================================================
// COE  ↔  Delaunay
// ============================================================================

template<typename T>
[[nodiscard]] Delaunay<T> coe_to_delaunay(const COE<T>& coe, T mu) {
    const T Lval = std::sqrt(mu * coe.a);
    const T Gval = Lval * std::sqrt(T(1) - coe.e * coe.e);
    const T Hval = Gval * std::cos(coe.i);
    const T M    = true_to_mean_anomaly(coe.nu, coe.e);
    return {
        wrap_0_2pi(M),
        wrap_0_2pi(coe.omega),
        wrap_0_2pi(coe.RAAN),
        Lval, Gval, Hval
    };
}

template<typename T>
[[nodiscard]] COE<T> delaunay_to_coe(const Delaunay<T>& del, T mu) {
    const T a  = del.L * del.L / mu;
    const T e  = std::sqrt(std::max(T(0), T(1) - (del.G / del.L) * (del.G / del.L)));
    const T i  = std::acos(std::clamp(del.H / del.G, T(-1), T(1)));
    const T E  = mean_to_eccentric_anomaly(del.l, e);
    const T nu = eccentric_to_true_anomaly(E, e);
    return {a, e, i, wrap_0_2pi(del.h), wrap_0_2pi(del.g), wrap_0_2pi(nu)};
}

// ============================================================================
// COE  ↔  Poincaré
// ============================================================================

template<typename T>
[[nodiscard]] Poincare<T> coe_to_poincare(const COE<T>& coe, T mu) {
    const T Lval  = std::sqrt(mu * coe.a);
    const T Gval  = Lval * std::sqrt(T(1) - coe.e * coe.e);
    const T Hval  = Gval * std::cos(coe.i);
    const T M     = true_to_mean_anomaly(coe.nu, coe.e);
    const T g_h   = coe.omega + coe.RAAN;    // g+h (Delaunay angles)
    const T LmG   = Lval - Gval;
    const T GmH   = Gval - Hval;
    return {
        wrap_0_2pi(M + g_h),                        // lambda = l+g+h
        Lval,
        (LmG >= T(0)) ?  std::sqrt(T(2) * LmG) * std::cos(g_h) : T(0),  // xi
        (LmG >= T(0)) ? -std::sqrt(T(2) * LmG) * std::sin(g_h) : T(0),  // eta
        (GmH >= T(0)) ?  std::sqrt(T(2) * GmH) * std::cos(coe.RAAN) : T(0), // p
        (GmH >= T(0)) ? -std::sqrt(T(2) * GmH) * std::sin(coe.RAAN) : T(0)  // q
    };
}

template<typename T>
[[nodiscard]] COE<T> poincare_to_coe(const Poincare<T>& poi, T mu) {
    const T Lval  = poi.L;
    const T a     = Lval * Lval / mu;
    const T Wxi   = poi.xi  * poi.xi;
    const T Weta  = poi.eta * poi.eta;
    const T Wp    = poi.p   * poi.p;
    const T Wq    = poi.q   * poi.q;
    const T Gval  = Lval - (Wxi + Weta) / T(2);
    const T Hval  = Gval  - (Wp + Wq)   / T(2);
    const T e     = std::sqrt(std::max(T(0), T(1) - (Gval / Lval) * (Gval / Lval)));
    const T i     = std::acos(std::clamp(Hval / Gval, T(-1), T(1)));
    const T g_h   = std::atan2(-poi.eta, poi.xi);   // atan2(-eta, xi) = g+h
    const T RAAN  = std::atan2(-poi.q,   poi.p);    // atan2(-q, p)    = h
    const T omega = wrap_0_2pi(g_h - RAAN);
    const T M     = wrap_0_2pi(poi.lambda - g_h);
    const T E     = mean_to_eccentric_anomaly(M, e);
    const T nu    = eccentric_to_true_anomaly(E, e);
    return {a, e, i, wrap_0_2pi(RAAN), omega, wrap_0_2pi(nu)};
}

// ============================================================================
// EQN  ↔  MEE  (direct)
// ============================================================================

template<typename T>
[[nodiscard]] MEE<T> eqn_to_mee(const EQN<T>& eqn) {
    return coe_to_mee(eqn_to_coe(eqn));
}

template<typename T>
[[nodiscard]] EQN<T> mee_to_eqn(const MEE<T>& mee) {
    return coe_to_eqn(mee_to_coe(mee));
}

// ============================================================================
// NSK  ↔  MEE  (direct)
// ============================================================================

template<typename T>
[[nodiscard]] MEE<T> nsk_to_mee(const NSK<T>& nsk) {
    return coe_to_mee(nsk_to_coe(nsk));
}

template<typename T>
[[nodiscard]] NSK<T> mee_to_nsk(const MEE<T>& mee) {
    return coe_to_nsk(mee_to_coe(mee));
}

// ============================================================================
// Cartesian  ↔  MEE / EQN / NSK / Delaunay / Poincaré
// ============================================================================

template<typename T>
[[nodiscard]] EQN<T> cart_to_eqn(const Cart<T>& cart, T mu) {
    return coe_to_eqn(cart_to_coe(cart, mu));
}

template<typename T>
[[nodiscard]] Cart<T> eqn_to_cart(const EQN<T>& eqn, T mu) {
    return coe_to_cart(eqn_to_coe(eqn), mu);
}

template<typename T>
[[nodiscard]] NSK<T> cart_to_nsk(const Cart<T>& cart, T mu) {
    return coe_to_nsk(cart_to_coe(cart, mu));
}

template<typename T>
[[nodiscard]] Cart<T> nsk_to_cart(const NSK<T>& nsk, T mu) {
    return coe_to_cart(nsk_to_coe(nsk), mu);
}

template<typename T>
[[nodiscard]] Delaunay<T> cart_to_delaunay(const Cart<T>& cart, T mu) {
    return coe_to_delaunay(cart_to_coe(cart, mu), mu);
}

template<typename T>
[[nodiscard]] Cart<T> delaunay_to_cart(const Delaunay<T>& del, T mu) {
    return coe_to_cart(delaunay_to_coe(del, mu), mu);
}

template<typename T>
[[nodiscard]] Poincare<T> cart_to_poincare(const Cart<T>& cart, T mu) {
    return coe_to_poincare(cart_to_coe(cart, mu), mu);
}

template<typename T>
[[nodiscard]] Cart<T> poincare_to_cart(const Poincare<T>& poi, T mu) {
    return coe_to_cart(poincare_to_coe(poi, mu), mu);
}

// ============================================================================
// TLE helpers
// ============================================================================

/// Extract approximate mean Keplerian elements from a TLE (double only).
///
/// WARNING: the result is in the TEME frame and uses Brouwer mean elements.
/// It is NOT an osculating COE.  Use only for quick sanity checks.
/// For accurate propagation, feed the TLE to a proper SGP4 implementation.
///
/// The gravitational parameter is hardcoded to Earth (WGS-72 value used by
/// SGP4) because the TLE format is Earth-specific by definition.
[[nodiscard]] inline COE<double> tle_to_classical_approx(const TLE& tle) {
    constexpr double DEG2RAD            = PI<double> / 180.0;
    constexpr double REV_DAY_TO_RAD_S   = TWO_PI<double> / 86400.0;
    constexpr double MU_EARTH_WGS72     = 3.986008e14;  // m³/s²  (WGS-72, as used by SGP4)

    const double n  = tle.n * REV_DAY_TO_RAD_S;
    const double a  = std::cbrt(MU_EARTH_WGS72 / (n * n));
    const double e  = tle.e;
    const double i  = tle.i_deg     * DEG2RAD;
    const double W  = tle.RAAN_deg  * DEG2RAD;
    const double w  = tle.omega_deg * DEG2RAD;
    const double M  = tle.M_deg     * DEG2RAD;
    const double E  = mean_to_eccentric_anomaly(M, e);
    const double nu = eccentric_to_true_anomaly(E, e);
    return {a, e, i, W, w, nu};
}

/// Build a TLE shell from osculating COE (double only, Earth only).
///
/// NOTE: This is a *rough approximation* — the TLE stores Brouwer mean
/// elements, not osculating elements.  A faithful conversion requires a
/// Brouwer–Lyddane mean/osculating theory.  B* and derivatives are set to 0.
/// Use only for testing or initialisation of a propagator.
[[nodiscard]] inline TLE classical_approx_to_tle(const COE<double>& coe, double mu,
                                                  int sat_number = 0,
                                                  int epoch_year = 0,
                                                  double epoch_day = 0.0) {
    constexpr double RAD2DEG       = 180.0 / PI<double>;
    constexpr double RAD_S_TO_REVD = 86400.0 / TWO_PI<double>;

    const double n = std::sqrt(mu / (coe.a * coe.a * coe.a));
    const double M = true_to_mean_anomaly(coe.nu, coe.e);

    TLE tle{};
    tle.sat_number     = sat_number;
    tle.classification = 'U';
    tle.epoch_year     = epoch_year;
    tle.epoch_day      = epoch_day;
    tle.ndot           = 0.0;
    tle.nddot          = 0.0;
    tle.bstar          = 0.0;
    tle.i_deg          = coe.i     * RAD2DEG;
    tle.RAAN_deg       = coe.RAAN  * RAD2DEG;
    tle.e              = coe.e;
    tle.omega_deg      = coe.omega * RAD2DEG;
    tle.M_deg          = M         * RAD2DEG;
    tle.n              = n * RAD_S_TO_REVD;
    tle.rev_number     = 0;
    return tle;
}

// ============================================================================
// Convenience: array-based generic interface  (useful for ODE integrators)
// ============================================================================

/// Pack / unpack helpers  (COE ↔ std::array<T,6>)
template<typename T>
[[nodiscard]] std::array<T,6> coe_to_array(const COE<T>& c) {
    return {c.a, c.e, c.i, c.RAAN, c.omega, c.nu};
}
template<typename T>
[[nodiscard]] COE<T> array_to_coe(const std::array<T,6>& a) {
    return {a[0], a[1], a[2], a[3], a[4], a[5]};
}

template<typename T>
[[nodiscard]] std::array<T,6> mee_to_array(const MEE<T>& m) {
    return {m.p, m.f, m.g, m.h, m.k, m.L};
}
template<typename T>
[[nodiscard]] MEE<T> array_to_mee(const std::array<T,6>& a) {
    return {a[0], a[1], a[2], a[3], a[4], a[5]};
}

template<typename T>
[[nodiscard]] std::array<T,6> cart_to_array(const Cart<T>& c) {
    return {c.x, c.y, c.z, c.vx, c.vy, c.vz};
}
template<typename T>
[[nodiscard]] Cart<T> array_to_cart(const std::array<T,6>& a) {
    return {a[0], a[1], a[2], a[3], a[4], a[5]};
}

template<typename T>
[[nodiscard]] std::array<T,6> delaunay_to_array(const Delaunay<T>& d) {
    return {d.l, d.g, d.h, d.L, d.G, d.H};
}
template<typename T>
[[nodiscard]] Delaunay<T> array_to_delaunay(const std::array<T,6>& a) {
    return {a[0], a[1], a[2], a[3], a[4], a[5]};
}

} // namespace orb
