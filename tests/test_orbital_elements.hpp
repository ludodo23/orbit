/**
 * @file test_orbital_elements.cpp
 * @brief Catch2 v3 test suite for orbital_elements.hpp
 *
 * Compile:
 *   g++ -std=c++17 -O2 -Wall -Wextra \
 *       test_orbital_elements.cpp -o test_orbital \
 *       $(pkg-config --cflags --libs catch2-with-main)
 *
 * Run:
 *   ./test_orbital            # all tests
 *   ./test_orbital -t "[mee]" # only MEE tests
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "orbital_elements.hpp"

using namespace orb;
using Catch::Approx;

// ============================================================================
// Constants & helpers
// ============================================================================

static constexpr double MU   = 3.986004418e14;  // Earth, WGS-84 [m³/s²]
static constexpr double TOL_R     = 1e-3;        // position tolerance [m]
static constexpr double TOL_V     = 1e-7;        // velocity tolerance [m/s]
static constexpr double TOL_ANGLE = 1e-9;        // angle tolerance [rad]
static constexpr double TOL_A     = 1e-3;        // semi-major axis [m]
static constexpr double TOL_E     = 1e-10;       // eccentricity [-]

// Circular angular distance — handles 0 ↔ 2π aliasing
inline double ang_dist(double a, double b) {
    double d = std::abs(a - b);
    if (d > PI<double>) d = TWO_PI<double> - d;
    return d;
}

// Check COE round-trip field by field
inline void require_coe_eq(const COE<double>& got, const COE<double>& ref,
                            double tol_a     = TOL_A,
                            double tol_angle = TOL_ANGLE) {
    REQUIRE(got.a == Approx(ref.a).margin(tol_a));
    REQUIRE(got.e == Approx(ref.e).margin(TOL_E));
    REQUIRE(got.i == Approx(ref.i).margin(tol_angle));
    REQUIRE(ang_dist(got.RAAN,  ref.RAAN)  < tol_angle);
    REQUIRE(ang_dist(got.omega, ref.omega) < tol_angle);
    REQUIRE(ang_dist(got.nu,    ref.nu)    < tol_angle);
}

// Check Cart round-trip field by field
inline void require_cart_eq(const Cart<double>& got, const Cart<double>& ref) {
    REQUIRE(got.x  == Approx(ref.x ).margin(TOL_R));
    REQUIRE(got.y  == Approx(ref.y ).margin(TOL_R));
    REQUIRE(got.z  == Approx(ref.z ).margin(TOL_R));
    REQUIRE(got.vx == Approx(ref.vx).margin(TOL_V));
    REQUIRE(got.vy == Approx(ref.vy).margin(TOL_V));
    REQUIRE(got.vz == Approx(ref.vz).margin(TOL_V));
}

// ============================================================================
// Reference orbits
// ============================================================================

// ISS-like  (low circular, inclined)
static const COE<double> ISS {
    6.780e6,  1.5e-3,  0.9006, 1.2566, 0.5236, 1.0472
};
// Molniya  (highly eccentric)
static const COE<double> MOLNIYA {
    26560e3,  0.74,    1.1345, 0.0,    PI<double>/2, 0.0
};
// GEO  (near-equatorial, near-circular)
static const COE<double> GEO {
    42164e3,  1e-4,    1e-4,   0.0,    0.0,          0.5
};
// Polar sun-synchronous
static const COE<double> SSO {
    7078e3,   1e-3,    PI<double>/2 + 0.17, 2.5, 1.0, 2.0
};

// ============================================================================
// Section 1 — Kepler equation
// ============================================================================

TEST_CASE("Kepler equation: M → E → ν → M round-trip", "[kepler]") {
    // Tabulated (e, M) pairs covering edge cases
    auto [e, M] = GENERATE(table<double,double>({
        {0.0,  0.0},  {0.0,  3.0},
        {0.1,  0.5},  {0.1,  PI<double>},
        {0.5,  1.5},  {0.5,  5.5},
        {0.8,  0.1},  {0.8,  3.14},
        {0.99, 0.01}, {0.99, 6.0}
    }));

    CAPTURE(e, M);
    const double E  = mean_to_eccentric_anomaly(M, e);
    const double nu = eccentric_to_true_anomaly(E, e);
    const double M2 = true_to_mean_anomaly(nu, e);
    REQUIRE(ang_dist(M2, M) < 1e-11);
}

TEST_CASE("Eccentric → true → eccentric anomaly round-trip", "[kepler]") {
    auto [e, E_in] = GENERATE(table<double,double>({
        {0.0, 0.0}, {0.3, 1.0}, {0.7, 2.5}, {0.95, 0.5}
    }));

    CAPTURE(e, E_in);
    const double nu  = eccentric_to_true_anomaly(E_in, e);
    const double E2  = true_to_eccentric_anomaly(nu, e);
    REQUIRE(ang_dist(E2, wrap_0_2pi(E_in)) < 1e-12);
}

// ============================================================================
// Section 2 — COE ↔ Cartesian
// ============================================================================

TEST_CASE("COE → Cart → COE round-trip", "[coe][cart]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );

    const auto cart = coe_to_cart(coe, MU);
    const auto coe2 = cart_to_coe(cart, MU);
    require_coe_eq(coe2, coe);
}

TEST_CASE("Cart → COE → Cart round-trip", "[coe][cart]") {
    // Build reference cart from ISS COE, then verify the back-and-forth
    const auto cart_ref = coe_to_cart(ISS, MU);
    const auto coe2     = cart_to_coe(cart_ref, MU);
    const auto cart2    = coe_to_cart(coe2, MU);
    require_cart_eq(cart2, cart_ref);
}

TEST_CASE("COE→Cart: vis-viva sanity check", "[coe][cart][physics]") {
    // |v|² = μ(2/r − 1/a)
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO}
    );
    const auto c   = coe_to_cart(coe, MU);
    const double r = std::sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
    const double v2= c.vx*c.vx + c.vy*c.vy + c.vz*c.vz;
    const double v2_expected = MU * (2.0/r - 1.0/coe.a);
    REQUIRE(v2 == Approx(v2_expected).epsilon(1e-9));
}

TEST_CASE("COE→Cart: angular momentum magnitude matches h=√(μp)", "[coe][cart][physics]") {
    const auto c = coe_to_cart(ISS, MU);
    const double hx = c.y*c.vz - c.z*c.vy;
    const double hy = c.z*c.vx - c.x*c.vz;
    const double hz = c.x*c.vy - c.y*c.vx;
    const double h  = std::sqrt(hx*hx + hy*hy + hz*hz);
    const double p  = ISS.a * (1.0 - ISS.e * ISS.e);
    REQUIRE(h == Approx(std::sqrt(MU * p)).epsilon(1e-9));
}

// ============================================================================
// Section 3 — COE ↔ Equinoctial (EQN)
// ============================================================================

TEST_CASE("COE ↔ EQN round-trip", "[coe][eqn]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );
    require_coe_eq(eqn_to_coe(coe_to_eqn(coe)), coe);
}

TEST_CASE("EQN: h²+k² == e²", "[eqn][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto eqn = coe_to_eqn(coe);
    REQUIRE(eqn.h*eqn.h + eqn.k*eqn.k == Approx(coe.e * coe.e).epsilon(1e-12));
}

TEST_CASE("EQN: near-circular orbit (e→0) stays regular", "[eqn][singularity]") {
    COE<double> circ{7000e3, 1e-6, 0.5, 1.0, 0.3, 1.5};
    const auto eqn  = coe_to_eqn(circ);
    const auto coe2 = eqn_to_coe(eqn);
    // We only care that a, e, i are recovered; omega+nu may be degenerate
    REQUIRE(coe2.a == Approx(circ.a).margin(TOL_A));
    REQUIRE(coe2.e == Approx(circ.e).margin(1e-8));
    REQUIRE(coe2.i == Approx(circ.i).margin(TOL_ANGLE));
}

// ============================================================================
// Section 4 — COE ↔ Modified Equinoctial (MEE)
// ============================================================================

TEST_CASE("COE ↔ MEE round-trip", "[coe][mee]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );
    require_coe_eq(mee_to_coe(coe_to_mee(coe)), coe);
}

TEST_CASE("MEE: p == a(1-e²)", "[mee][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA}, COE<double>{GEO});
    const auto mee = coe_to_mee(coe);
    const double p_expected = coe.a * (1.0 - coe.e * coe.e);
    REQUIRE(mee.p == Approx(p_expected).epsilon(1e-12));
}

TEST_CASE("MEE: f²+g² == e²", "[mee][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto mee = coe_to_mee(coe);
    REQUIRE(mee.f*mee.f + mee.g*mee.g == Approx(coe.e * coe.e).epsilon(1e-12));
}

TEST_CASE("Cart ↔ MEE direct path round-trip", "[mee][cart]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA}, COE<double>{GEO});
    const auto cart_ref = coe_to_cart(coe, MU);
    const auto mee      = cart_to_mee(cart_ref, MU);
    const auto cart2    = mee_to_cart(mee, MU);
    require_cart_eq(cart2, cart_ref);
}

// ============================================================================
// Section 5 — COE ↔ Non-Singular Keplerian (NSK)
// ============================================================================

TEST_CASE("COE ↔ NSK round-trip", "[coe][nsk]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );
    require_coe_eq(nsk_to_coe(coe_to_nsk(coe)), coe);
}

TEST_CASE("NSK: ex²+ey² == e²", "[nsk][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto nsk = coe_to_nsk(coe);
    REQUIRE(nsk.ex*nsk.ex + nsk.ey*nsk.ey == Approx(coe.e * coe.e).epsilon(1e-12));
}

TEST_CASE("NSK: ix²+iy² == tan²(i/2)", "[nsk][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{SSO});
    const auto nsk = coe_to_nsk(coe);
    const double ti2 = std::tan(coe.i / 2.0);
    REQUIRE(nsk.ix*nsk.ix + nsk.iy*nsk.iy == Approx(ti2 * ti2).epsilon(1e-12));
}

TEST_CASE("NSK: u = omega + nu (mod 2π)", "[nsk][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto nsk = coe_to_nsk(coe);
    REQUIRE(ang_dist(nsk.u, wrap_0_2pi(coe.omega + coe.nu)) < TOL_ANGLE);
}

// ============================================================================
// Section 6 — COE ↔ Delaunay
// ============================================================================

TEST_CASE("COE ↔ Delaunay round-trip", "[coe][delaunay]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );
    require_coe_eq(delaunay_to_coe(coe_to_delaunay(coe, MU), MU), coe);
}

TEST_CASE("Delaunay: L = √(μa)", "[delaunay][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA}, COE<double>{GEO});
    const auto del = coe_to_delaunay(coe, MU);
    REQUIRE(del.L == Approx(std::sqrt(MU * coe.a)).epsilon(1e-12));
}

TEST_CASE("Delaunay: G = L·√(1-e²)", "[delaunay][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto del = coe_to_delaunay(coe, MU);
    REQUIRE(del.G == Approx(del.L * std::sqrt(1.0 - coe.e*coe.e)).epsilon(1e-12));
}

TEST_CASE("Delaunay: H = G·cos(i)", "[delaunay][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{SSO});
    const auto del = coe_to_delaunay(coe, MU);
    REQUIRE(del.H == Approx(del.G * std::cos(coe.i)).epsilon(1e-12));
}

TEST_CASE("Cart ↔ Delaunay round-trip", "[delaunay][cart]") {
    const auto cart_ref = coe_to_cart(ISS, MU);
    const auto del      = cart_to_delaunay(cart_ref, MU);
    const auto cart2    = delaunay_to_cart(del, MU);
    require_cart_eq(cart2, cart_ref);
}

// ============================================================================
// Section 7 — COE ↔ Poincaré
// ============================================================================

TEST_CASE("COE ↔ Poincaré round-trip", "[coe][poincare]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO},
        COE<double>{SSO}
    );
    require_coe_eq(poincare_to_coe(coe_to_poincare(coe, MU), MU), coe);
}

TEST_CASE("Poincaré: L = √(μa)", "[poincare][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA}, COE<double>{GEO});
    const auto poi = coe_to_poincare(coe, MU);
    REQUIRE(poi.L == Approx(std::sqrt(MU * coe.a)).epsilon(1e-12));
}

TEST_CASE("Poincaré: xi²+eta² = 2(L−G)", "[poincare][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto poi = coe_to_poincare(coe, MU);
    const auto del = coe_to_delaunay(coe, MU);
    REQUIRE(poi.xi*poi.xi + poi.eta*poi.eta == Approx(2.0*(del.L - del.G)).epsilon(1e-10));
}

TEST_CASE("Poincaré: p²+q² = 2(G−H)", "[poincare][invariant]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{SSO});
    const auto poi = coe_to_poincare(coe, MU);
    const auto del = coe_to_delaunay(coe, MU);
    REQUIRE(poi.p*poi.p + poi.q*poi.q == Approx(2.0*(del.G - del.H)).epsilon(1e-10));
}

TEST_CASE("Cart ↔ Poincaré round-trip", "[poincare][cart]") {
    const auto cart_ref = coe_to_cart(MOLNIYA, MU);
    const auto poi      = cart_to_poincare(cart_ref, MU);
    const auto cart2    = poincare_to_cart(poi, MU);
    require_cart_eq(cart2, cart_ref);
}

// ============================================================================
// Section 8 — Cross-conversions  (EQN ↔ MEE, NSK ↔ MEE)
// ============================================================================

TEST_CASE("EQN ↔ MEE cross-conversion", "[eqn][mee]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{MOLNIYA});
    const auto mee  = coe_to_mee(coe);
    const auto eqn2 = mee_to_eqn(mee);
    const auto mee2 = eqn_to_mee(eqn2);
    REQUIRE(mee2.p == Approx(mee.p).epsilon(1e-12));
    REQUIRE(mee2.f == Approx(mee.f).epsilon(1e-12));
    REQUIRE(mee2.g == Approx(mee.g).epsilon(1e-12));
    REQUIRE(mee2.h == Approx(mee.h).epsilon(1e-12));
    REQUIRE(mee2.k == Approx(mee.k).epsilon(1e-12));
    REQUIRE(ang_dist(mee2.L, mee.L) < TOL_ANGLE);
}

TEST_CASE("NSK ↔ MEE cross-conversion", "[nsk][mee]") {
    auto coe = GENERATE(COE<double>{ISS}, COE<double>{SSO});
    const auto mee  = coe_to_mee(coe);
    const auto nsk2 = mee_to_nsk(mee);
    const auto mee2 = nsk_to_mee(nsk2);
    REQUIRE(mee2.p == Approx(mee.p).epsilon(1e-12));
    REQUIRE(mee2.f == Approx(mee.f).epsilon(1e-12));
    REQUIRE(mee2.g == Approx(mee.g).epsilon(1e-12));
}

// ============================================================================
// Section 9 — TLE helpers  (non-templated, double only)
// ============================================================================

TEST_CASE("TLE: classical_approx_to_tle → tle_to_classical_approx round-trip", "[tle]") {
    auto coe = GENERATE(
        COE<double>{ISS},
        COE<double>{MOLNIYA},
        COE<double>{GEO}
    );
    const TLE    tle  = classical_approx_to_tle(coe, MU, 25544, 24, 1.0);
    const auto   coe2 = tle_to_classical_approx(tle);

    // a: Kepler 3rd law from n (rev/day) → some numerical noise but < 1 km
    REQUIRE(coe2.a == Approx(coe.a).margin(1e3));
    REQUIRE(coe2.e == Approx(coe.e).margin(1e-6));
    REQUIRE(coe2.i == Approx(coe.i).margin(1e-7));
    REQUIRE(ang_dist(coe2.RAAN,  coe.RAAN)  < 1e-7);
    REQUIRE(ang_dist(coe2.omega, coe.omega) < 1e-7);
    REQUIRE(ang_dist(coe2.nu,    coe.nu)    < 1e-7);
}

TEST_CASE("TLE: mean motion is positive and in plausible range", "[tle]") {
    // LEO: ~16 rev/day, GEO: ~1 rev/day
    const TLE iss = classical_approx_to_tle(ISS, MU);
    const TLE geo = classical_approx_to_tle(GEO, MU);
    REQUIRE(iss.n > 14.0);
    REQUIRE(iss.n < 18.0);
    REQUIRE(geo.n > 0.9);
    REQUIRE(geo.n < 1.1);
}

TEST_CASE("TLE: struct is not templated (double fields)", "[tle]") {
    // Static type assertions — compile-time checks
    static_assert(std::is_same_v<decltype(TLE::e),          double>, "TLE.e must be double");
    static_assert(std::is_same_v<decltype(TLE::n),          double>, "TLE.n must be double");
    static_assert(std::is_same_v<decltype(TLE::bstar),      double>, "TLE.bstar must be double");
    static_assert(std::is_same_v<decltype(TLE::sat_number), int>,    "TLE.sat_number must be int");
    SUCCEED("TLE is correctly typed as non-template double struct");
}

// ============================================================================
// Section 10 — Array interface
// ============================================================================

TEST_CASE("COE ↔ array pack/unpack", "[array]") {
    const auto arr  = coe_to_array(ISS);
    const auto coe2 = array_to_coe(arr);
    // Exact bit-identical round-trip — check equality directly
    REQUIRE(coe2.a     == ISS.a);
    REQUIRE(coe2.e     == ISS.e);
    REQUIRE(coe2.i     == ISS.i);
    REQUIRE(coe2.RAAN  == ISS.RAAN);
    REQUIRE(coe2.omega == ISS.omega);
    REQUIRE(coe2.nu    == ISS.nu);
}

TEST_CASE("MEE ↔ array pack/unpack", "[array]") {
    const auto mee  = coe_to_mee(ISS);
    const auto arr  = mee_to_array(mee);
    const auto mee2 = array_to_mee(arr);
    REQUIRE(mee2.p == mee.p);
    REQUIRE(mee2.f == mee.f);
    REQUIRE(mee2.g == mee.g);
    REQUIRE(mee2.L == mee.L);
}

TEST_CASE("Delaunay ↔ array pack/unpack", "[array]") {
    const auto del  = coe_to_delaunay(MOLNIYA, MU);
    const auto arr  = delaunay_to_array(del);
    const auto del2 = array_to_delaunay(arr);
    REQUIRE(del2.L == del.L);
    REQUIRE(del2.G == del.G);
    REQUIRE(del2.H == del.H);
    REQUIRE(del2.l == del.l);
}

// ============================================================================
// Section 11 — float precision
// ============================================================================

TEST_CASE("Float: COE ↔ Cart round-trip (single precision)", "[float]") {
    COE<float> coe_f{6.78e6f, 0.01f, 0.9f, 1.2f, 0.5f, 1.0f};
    constexpr float MU_F = 3.986004418e14f;
    const auto cart = coe_to_cart(coe_f, MU_F);
    const auto coe2 = cart_to_coe(cart, MU_F);
    // Float accuracy: ~1e-6 relative → ~10 m on a, ~1e-4 on angles
    REQUIRE(coe2.a == Approx(static_cast<double>(coe_f.a)).margin(1e2));
    REQUIRE(coe2.e == Approx(static_cast<double>(coe_f.e)).margin(1e-4));
    REQUIRE(coe2.i == Approx(static_cast<double>(coe_f.i)).margin(1e-4));
}

TEST_CASE("Float: MEE ↔ COE round-trip (single precision)", "[float]") {
    COE<float> coe_f{7200e3f, 0.05f, 1.1f, 0.8f, 0.4f, 2.0f};
    const auto mee  = coe_to_mee(coe_f);
    const auto coe2 = mee_to_coe(mee);
    REQUIRE(coe2.a == Approx(static_cast<double>(coe_f.a)).margin(1e2));
    REQUIRE(coe2.e == Approx(static_cast<double>(coe_f.e)).margin(1e-4));
}

// ============================================================================
// Section 12 — Angle wrap utilities
// ============================================================================

TEST_CASE("wrap_0_2pi: result always in [0, 2π)", "[utils]") {
    for (double a : {-7.0, -PI<double>, -0.1, 0.0, PI<double>, 6.0, 10.0, 100.0}) {
        const double w = wrap_0_2pi(a);
        REQUIRE(w >= 0.0);
        REQUIRE(w <  TWO_PI<double>);
    }
}

TEST_CASE("wrap_neg_pi_pi: result always in [-π, π)", "[utils]") {
    // The implementation maps to [-pi, pi)  (standard fmod convention)
    for (double a : {-7.0, -PI<double>, 0.0, PI<double>, 6.0, 10.0}) {
        const double w = wrap_neg_pi_pi(a);
        REQUIRE(w >= -PI<double>);
        REQUIRE(w <   PI<double>);
    }
}
