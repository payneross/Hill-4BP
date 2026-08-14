#pragma once

#include "crtbp_cpp/ode.hpp"

#include <array>

namespace crtbp {

using State6 = std::array<double, 6>;
using Matrix6 = std::array<std::array<double, 6>, 6>;

struct LibrationPoints {
    std::array<double, 3> L1{};
    std::array<double, 3> L2{};
    std::array<double, 3> L3{};
    std::array<double, 3> L4{};
    std::array<double, 3> L5{};
};

struct CrossingData {
    double x = 0.0;
    double xdot = 0.0;
    double ydot = 0.0;
};

Vector crtbp_vector_field(double t, const Vector& y, double G, double mu);
Vector backward_crtbp_vector_field(double t, const Vector& y, double G, double mu);
Vector crtbp_state_transition_field(double t, const Vector& y, double mu);

double jacobi_constant(const std::array<double, 3>& position,
                       const std::array<double, 3>& velocity,
                       double mu);

double velocity_magnitude(const std::array<double, 3>& position, double C, double mu);

LibrationPoints libration_points(double mu);

Matrix6 state_transition_matrix(double t0,
                                double tf,
                                const State6& state,
                                double mu,
                                const OdeOptions& options = {});

State6 poincare_newton_crossing(double t0,
                                const State6& initial_point,
                                double t1,
                                const State6& final_point,
                                double tolerance,
                                double G,
                                double mu,
                                const OdeOptions& options = {});

CrossingData interpolate_crossing_data(const std::array<double, 4>& p1,
                                       const std::array<double, 4>& p2,
                                       double C,
                                       double mu);

}  // namespace crtbp
