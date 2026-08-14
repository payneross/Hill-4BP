#pragma once

#include "crtbp_cpp/ode.hpp"

#include <array>
#include <utility>

namespace crtbp {

using State4 = std::array<double, 4>;
using Vec2 = std::array<double, 2>;

struct Hill4BPParameters {
    double G = 0.0;
    double mu = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double c = 0.0;
    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double A = 0.0;
    double B = 0.0;
    double nu = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;
};

struct Hill4BPPoincareConfig {
    double h_e1_scale = 0.0;
    int x_count = 0;
    int xdot_count = 0;
    int iterates = 0;
    double x_begin = 0.0;
    double x_end = 0.0;
    double xdot_lower = 0.0;
    double xdot_upper = 0.0;
    double y_section = 0.0;
    double tau_step = 0.0;
    double rel_tol = 0.0;
    double abs_tol = 0.0;
    double initial_step = 0.0;
    double max_step = 0.0;
    int svg_stride = 0;
};

struct Hill4BPEquilibrium {
    State4 mcgehee{};
    State4 cartesian{};
    double h = 0.0;
};

struct Hill4BPEquilibria {
    Hill4BPEquilibrium E1;
    Hill4BPEquilibrium E2;
    Hill4BPEquilibrium E3;
    Hill4BPEquilibrium E4;
    State4 collision_plus{};
    State4 collision_minus{};
};

Hill4BPParameters hill4bp_parameters();
Hill4BPPoincareConfig hill4bp_poincare_config();
Hill4BPEquilibria hill4bp_equilibria(const Hill4BPParameters& params);

Vector hill4bp_vector_field(double t, const Vector& state, const Hill4BPParameters& params);
Vector hill4bp_regularized_field(double tau, const Vector& state, const Hill4BPParameters& params);
// Smooth radial coordinates R=sqrt(Cartesian radius)=r^(1/(alpha+2)).
// For the configured alpha=3 study this removes the fractional r^0.8 term.
Vector hill4bp_sqrt_radius_field(
    double tau, const Vector& state, const Hill4BPParameters& params);

double hill4bp_effective_potential(const Vec2& q, const Hill4BPParameters& params);
double hill4bp_hamiltonian(const State4& state, const Hill4BPParameters& params);
double hill4bp_jacobi_constant(const State4& state, const Hill4BPParameters& params);

std::pair<double, bool> hill4bp_speed(const Vec2& q, double h, const Hill4BPParameters& params);

State4 hill4bp_cartesian_to_mcgehee(const State4& state, const Hill4BPParameters& params);
State4 hill4bp_mcgehee_to_cartesian(const State4& state, const Hill4BPParameters& params);
State4 hill4bp_cartesian_to_sqrt_radius(
    const State4& state, const Hill4BPParameters& params);
State4 hill4bp_sqrt_radius_to_cartesian(
    const State4& state, const Hill4BPParameters& params);

double hill4bp_section_event(double tau, const Vector& state, const Hill4BPParameters& params);
double hill4bp_sqrt_radius_section_event(
    double tau, const Vector& state, const Hill4BPParameters& params);

Vector to_vector(const State4& state);
State4 to_state4(const Vector& values);

}  // namespace crtbp
