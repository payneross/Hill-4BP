#include "crtbp_cpp/hill4bp.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void check_near(double actual,
                double expected,
                double absolute_tolerance,
                const std::string& message) {
    if (!std::isfinite(actual)
        || std::abs(actual - expected) > absolute_tolerance) {
        ++failures;
        std::cerr << "FAIL: " << message << ": actual=" << actual
                  << " expected=" << expected
                  << " tolerance=" << absolute_tolerance << '\n';
    }
}

void check_scaled(double actual,
                  double expected,
                  double relative_tolerance,
                  const std::string& message) {
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    check_near(actual, expected, relative_tolerance * scale, message);
}

bool finite_vector(const crtbp::Vector& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace

int main() {
    const crtbp::Hill4BPParameters p = crtbp::hill4bp_parameters();

    // Lock the intended current study: highly oblate third body, c1=c2=0.
    check_near(p.c1, 0.0, 0.0, "c1 remains zero");
    check_near(p.c2, 0.0, 0.0, "c2 remains zero");
    check_near(p.c3, -0.097656556510836, 1e-16, "configured c3");
    check_near(p.c, 0.097656556510836, 1e-16, "derived c=-c3");
    check_near(p.lambda1, 0.027254858036835, 1e-16, "configured lambda1");
    check_near(p.lambda2, 2.972745141976870, 1e-15, "configured lambda2");
    check_near(p.A, -0.986372570988435, 1e-15, "derived A");
    check_near(p.B, 0.4863725709815825, 1e-15, "derived B");
    check_near(p.alpha, 3.0, 0.0, "alpha");
    check_near(p.beta, 1.5, 0.0, "beta");
    check_near(p.gamma, 0.4, 1e-16, "gamma");

    const crtbp::Hill4BPEquilibria equilibria = crtbp::hill4bp_equilibria(p);
    check_near(equilibria.E1.cartesian[0], 0.7905547107721489, 2e-14, "E1 x");
    check_near(equilibria.E1.mcgehee[0], 0.5556864759393503, 2e-14, "E1 paper r");
    check_near(equilibria.E1.h, -2.3915368938607138, 2e-14, "E1 energy");
    check_near(equilibria.E3.cartesian[1], 3.3515520023013403, 2e-13, "E3 y");
    check_near(equilibria.E3.h, -0.4540387530293194, 2e-14, "E3 energy");
    check_near(
        std::sqrt(2.0 * p.c), 0.44194243179589804, 2e-15,
        "collision limiting speed");

    for (const double R : {0.05, 0.2, 0.8, 1.3}) {
        const double theta = 0.37;
        const double v = -0.21;
        const double w = 0.48;
        const double paper_r = std::pow(R, p.alpha + 2.0);
        const crtbp::Vector old_state{paper_r, theta, v, w};
        const crtbp::Vector smooth_state{R, theta, v, w};
        const crtbp::Vector old_field =
            crtbp::hill4bp_regularized_field(0.0, old_state, p);
        const crtbp::Vector smooth_field =
            crtbp::hill4bp_sqrt_radius_field(0.0, smooth_state, p);

        check_scaled(
            old_field[0],
            (p.alpha + 2.0) * std::pow(R, p.alpha + 1.0) * smooth_field[0],
            2e-13, "radial push-forward equivalence");
        for (int i = 1; i < 4; ++i) {
            check_scaled(
                old_field[static_cast<std::size_t>(i)],
                smooth_field[static_cast<std::size_t>(i)],
                2e-13, "nonradial field equivalence");
        }
    }

    const crtbp::State4 cartesian{{0.65, 0.12, -0.31, 0.74}};
    const crtbp::State4 old_mcgehee =
        crtbp::hill4bp_cartesian_to_mcgehee(cartesian, p);
    const crtbp::State4 smooth =
        crtbp::hill4bp_cartesian_to_sqrt_radius(cartesian, p);
    check_scaled(
        smooth[0], std::pow(old_mcgehee[0], 1.0 / (p.alpha + 2.0)),
        5e-13, "R=paper_r^(1/(alpha+2))");
    for (int i = 1; i < 4; ++i) {
        check_scaled(
            smooth[static_cast<std::size_t>(i)],
            old_mcgehee[static_cast<std::size_t>(i)],
            5e-13, "coordinate angular/momentum agreement");
    }

    const crtbp::State4 round_trip =
        crtbp::hill4bp_sqrt_radius_to_cartesian(smooth, p);
    for (int i = 0; i < 4; ++i) {
        check_scaled(
            round_trip[static_cast<std::size_t>(i)],
            cartesian[static_cast<std::size_t>(i)],
            5e-13, "Cartesian sqrt-radius round trip");
    }

    const crtbp::Vector near_collision =
        crtbp::hill4bp_sqrt_radius_field(
            0.0, crtbp::Vector{1e-8, 0.7, -0.2, 0.1}, p);
    check(finite_vector(near_collision), "smooth field is finite near collision");

    const double R_section = 0.9;
    const double paper_r_section = std::pow(R_section, p.alpha + 2.0);
    const crtbp::State4 positive_x_up{{
        R_section, 0.0, 0.0, paper_r_section + 0.1}};
    const crtbp::State4 negative_x_up{{
        R_section, std::acos(-1.0), 0.0, paper_r_section - 0.1}};
    check(
        crtbp::hill4bp_sqrt_radius_to_cartesian(positive_x_up, p)[3] > 0.0,
        "positive-x section orientation agrees with positive ydot");
    check(
        crtbp::hill4bp_sqrt_radius_to_cartesian(negative_x_up, p)[3] > 0.0,
        "negative-x section orientation agrees with positive ydot");

    const crtbp::State4 e1_smooth =
        crtbp::hill4bp_cartesian_to_sqrt_radius(equilibria.E1.cartesian, p);
    check_near(e1_smooth[0], 0.8891314361623646, 2e-14, "E1 R");
    const crtbp::Vector e1_field =
        crtbp::hill4bp_sqrt_radius_field(
            0.0, crtbp::to_vector(e1_smooth), p);
    for (double component : e1_field) {
        check_near(component, 0.0, 5e-13, "smooth field vanishes at E1");
    }

    const crtbp::Vector cartesian_e1_field =
        crtbp::hill4bp_vector_field(
            0.0, crtbp::to_vector(equilibria.E1.cartesian), p);
    for (double component : cartesian_e1_field) {
        check_near(component, 0.0, 5e-13, "Cartesian field vanishes at E1");
    }

    check_near(
        crtbp::hill4bp_jacobi_constant(cartesian, p),
        -2.0 * crtbp::hill4bp_hamiltonian(cartesian, p),
        1e-15, "Jacobi constant convention");

    const double closed_h = equilibria.E1.h - 1e-2;
    for (const double root : {0.748680276789277, 0.834670844621541}) {
        const double speed_squared = 2.0 * (
            closed_h + crtbp::hill4bp_effective_potential(
                crtbp::Vec2{{root, 0.0}}, p));
        check_near(speed_squared, 0.0, 2e-14, "closed delta=1e-2 ZVB root");
    }

    // This is the first old sweep point that appeared to hang in CAPD. It is
    // a collision-asymptotic trajectory, not a very late section return.
    crtbp::OdeOptions collision_options;
    collision_options.rel_tol = 1e-11;
    collision_options.abs_tol = 1e-13;
    collision_options.initial_step = 1e-3;
    collision_options.max_step = 0.05;
    collision_options.max_steps = 1000000;
    const crtbp::Vector collision_initial{
        0.806225774829855,
        1.5384615384615385e-12,
        -0.08499307571969406,
        0.5834677490659379,
    };
    const crtbp::OdeResult collision_result = crtbp::integrate_with_events(
        [&p](double tau, const crtbp::Vector& state) {
            return crtbp::hill4bp_sqrt_radius_field(tau, state, p);
        },
        0.0,
        collision_initial,
        60.0,
        collision_options,
        [&p](double tau, const crtbp::Vector& state) {
            return crtbp::hill4bp_sqrt_radius_section_event(tau, state, p);
        },
        -1);
    check(collision_result.success, "known collision approach integrates to tau=60");
    check(
        collision_result.state.size() == 4
            && collision_result.state[0] * collision_result.state[0] < 1e-6
            && collision_result.state[2] < 0.0,
        "known nonreturn reaches the inward collision threshold");
    check(
        collision_result.events.empty(),
        "known collision approach has no intervening downward section crossing");

    if (failures != 0) {
        std::cerr << failures << " regression check(s) failed\n";
        return 1;
    }
    std::cout << "All Hill 4BP regression checks passed\n";
    return 0;
}
