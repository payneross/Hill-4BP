#include "crtbp_cpp/hill4bp.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace crtbp {
namespace {

double equilibrium_radius(double lambda, const Hill4BPParameters& params) {
    const double power = 2.0 - params.gamma * (params.nu + 2.0);
    auto f = [&](double r) {
        return lambda * r * r - params.nu * std::pow(r, power) - params.alpha * params.c;
    };

    double lo = std::numeric_limits<double>::epsilon();
    double hi = std::max(1.0, 2.0 * std::pow(lambda, -1.0 / (2.0 - power)));
    while (f(hi) <= 0.0) {
        hi *= 2.0;
    }

    double flo = f(lo);
    for (int i = 0; i < 160; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = f(mid);
        if (fmid == 0.0 || std::abs(hi - lo) < 1e-15) {
            return mid;
        }
        if (flo * fmid <= 0.0) {
            hi = mid;
        } else {
            lo = mid;
            flo = fmid;
        }
    }

    return 0.5 * (lo + hi);
}

double sqr(double value) {
    return value * value;
}

std::string trim(const std::string& value) {
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) {
        ++begin;
    }

    std::size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return value.substr(begin, end - begin);
}

std::string hill4bp_parameter_file() {
    const std::vector<std::string> candidates{
        "Hill_4BP_parameters.cfg",
        "../Hill_4BP_parameters.cfg",
        "../../Hill_4BP_parameters.cfg",
        "../../../Hill_4BP_parameters.cfg",
    };

    for (const std::string& candidate : candidates) {
        std::ifstream in(candidate);
        if (in) {
            return candidate;
        }
    }

    throw std::runtime_error("could not find Hill_4BP_parameters.cfg");
}

std::map<std::string, double> read_hill4bp_parameter_file() {
    const std::string path = hill4bp_parameter_file();
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open Hill 4BP parameter file: " + path);
    }

    std::map<std::string, double> values;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::runtime_error("invalid Hill 4BP config line " + std::to_string(line_number));
        }

        const std::string key = trim(line.substr(0, equals));
        const std::string raw_value = trim(line.substr(equals + 1));
        if (key.empty() || raw_value.empty()) {
            throw std::runtime_error("invalid Hill 4BP config line " + std::to_string(line_number));
        }

        std::size_t consumed = 0;
        const double value = std::stod(raw_value, &consumed);
        if (consumed != raw_value.size() || !std::isfinite(value)) {
            throw std::runtime_error("invalid numeric value for Hill 4BP config key: " + key);
        }
        if (!values.emplace(key, value).second) {
            throw std::runtime_error("duplicate Hill 4BP config key: " + key);
        }
    }

    return values;
}

const std::map<std::string, double>& hill4bp_config_values() {
    static const std::map<std::string, double> values = read_hill4bp_parameter_file();
    return values;
}

double required_config_value(const std::map<std::string, double>& values, const std::string& key) {
    const auto it = values.find(key);
    if (it == values.end()) {
        throw std::runtime_error("missing Hill 4BP config key: " + key);
    }
    return it->second;
}

double optional_config_value(const std::map<std::string, double>& values,
                             const std::string& key,
                             double default_value) {
    const auto it = values.find(key);
    return it == values.end() ? default_value : it->second;
}

int required_positive_int_config_value(const std::map<std::string, double>& values,
                                       const std::string& key) {
    const double value = required_config_value(values, key);
    const double rounded = std::round(value);
    if (value < 1.0 || std::abs(value - rounded) > std::numeric_limits<double>::epsilon() * value) {
        throw std::runtime_error("Hill 4BP config key must be a positive integer: " + key);
    }
    return static_cast<int>(rounded);
}

}  // namespace

Hill4BPParameters hill4bp_parameters() {
    const auto& values = hill4bp_config_values();
    Hill4BPParameters params;
    params.G = required_config_value(values, "G");
    params.mu = required_config_value(values, "mu");
    params.u1 = required_config_value(values, "u1");
    params.u2 = required_config_value(values, "u2");
    params.c1 = required_config_value(values, "c1");
    params.c2 = required_config_value(values, "c2");
    params.c3 = required_config_value(values, "c3");
    params.nu = required_config_value(values, "nu");
    params.alpha = required_config_value(values, "alpha");

    const double mu = params.mu;
    const double u1 = params.u1;
    const double u2 = params.u2;

    const double delta =
        sqr(mu * std::pow(u1, 3.0) + (1.0 - mu) * std::pow(u2, 3.0))
        - mu * (1.0 - mu) * u1 * u2
            * (-std::pow(u1, 4.0) - std::pow(u2, 4.0) + 2.0 * sqr(u1)
               + 2.0 * sqr(u2) + 2.0 * sqr(u1) * sqr(u2) - 1.0);

    const double base = 2.0 - 2.0 * (1.0 - mu) / std::pow(u1, 5.0)
        - 2.0 * mu / std::pow(u2, 5.0) + 3.0 * (1.0 - mu) / std::pow(u1, 3.0)
        + 3.0 * mu / std::pow(u2, 3.0);

    const double computed_lambda1 =
        0.5 * (base - 3.0 * std::sqrt(delta) / (std::pow(u1, 3.0) * std::pow(u2, 3.0)));
    const double computed_lambda2 =
        0.5 * (base + 3.0 * std::sqrt(delta) / (std::pow(u1, 3.0) * std::pow(u2, 3.0)));

    params.lambda1 = optional_config_value(values, "lambda1", computed_lambda1);
    params.lambda2 = optional_config_value(values, "lambda2", computed_lambda2);
    params.A = (1.0 - params.lambda2) / 2.0;
    params.B = (1.0 - params.lambda1) / 2.0;
    params.beta = params.alpha / 2.0;
    params.gamma = 2.0 / (params.alpha + 2.0);
    params.c = -params.c3;

    return params;
}

Hill4BPPoincareConfig hill4bp_poincare_config() {
    const auto& values = hill4bp_config_values();
    Hill4BPPoincareConfig config;
    config.h_e1_scale = required_config_value(values, "poincare_h_e1_scale");
    config.x_count = required_positive_int_config_value(values, "poincare_k");
    config.xdot_count = required_positive_int_config_value(values, "poincare_l");
    config.iterates = required_positive_int_config_value(values, "poincare_iterates");
    config.x_begin = required_config_value(values, "poincare_x_begin");
    config.x_end = required_config_value(values, "poincare_x_end");
    config.xdot_lower = required_config_value(values, "poincare_xdot_lower");
    config.xdot_upper = required_config_value(values, "poincare_xdot_upper");
    config.y_section = required_config_value(values, "poincare_y_section");
    config.tau_step = required_config_value(values, "poincare_tau_step");
    config.rel_tol = required_config_value(values, "poincare_rel_tol");
    config.abs_tol = required_config_value(values, "poincare_abs_tol");
    config.initial_step = required_config_value(values, "poincare_initial_step");
    config.max_step = required_config_value(values, "poincare_max_step");
    config.svg_stride = required_positive_int_config_value(values, "poincare_svg_stride");
    return config;
}

Hill4BPEquilibria hill4bp_equilibria(const Hill4BPParameters& params) {
    const double r1 = equilibrium_radius(params.lambda2, params);
    const double r2 = equilibrium_radius(params.lambda1, params);
    const double rho1 = std::pow(r1, params.gamma);
    const double rho2 = std::pow(r2, params.gamma);

    Hill4BPEquilibria equilibria;
    equilibria.E1.mcgehee = {{r1, 0.0, 0.0, r1}};
    equilibria.E2.mcgehee = {{r1, std::acos(-1.0), 0.0, r1}};
    equilibria.E3.mcgehee = {{r2, std::acos(-1.0) / 2.0, 0.0, r2}};
    equilibria.E4.mcgehee = {{r2, 3.0 * std::acos(-1.0) / 2.0, 0.0, r2}};

    equilibria.E1.cartesian = {{rho1, 0.0, 0.0, 0.0}};
    equilibria.E2.cartesian = {{-rho1, 0.0, 0.0, 0.0}};
    equilibria.E3.cartesian = {{0.0, rho2, 0.0, 0.0}};
    equilibria.E4.cartesian = {{0.0, -rho2, 0.0, 0.0}};

    equilibria.E1.h = hill4bp_hamiltonian(equilibria.E1.cartesian, params);
    equilibria.E2.h = equilibria.E1.h;
    equilibria.E3.h = hill4bp_hamiltonian(equilibria.E3.cartesian, params);
    equilibria.E4.h = equilibria.E3.h;

    equilibria.collision_plus = {{0.0, std::numeric_limits<double>::quiet_NaN(), std::sqrt(2.0 * params.c), 0.0}};
    equilibria.collision_minus = {{0.0, std::numeric_limits<double>::quiet_NaN(), -std::sqrt(2.0 * params.c), 0.0}};

    return equilibria;
}

Vector hill4bp_vector_field(double, const Vector& state, const Hill4BPParameters& params) {
    if (state.size() != 4) {
        throw std::invalid_argument("Hill 4BP vector field expects a 4-vector");
    }

    const double x = state[0];
    const double y = state[1];
    const double xdot = state[2];
    const double ydot = state[3];
    const double rho = std::sqrt(x * x + y * y);
    if (rho == 0.0) {
        throw std::runtime_error("Hill 4BP Cartesian field is singular at collision");
    }

    const double gravity_factor = 1.0 / std::pow(rho, 3.0) + 3.0 * params.c / std::pow(rho, 5.0);
    const double xddot = 2.0 * ydot + (params.lambda2 - gravity_factor) * x;
    const double yddot = -2.0 * xdot + (params.lambda1 - gravity_factor) * y;
    return {xdot, ydot, xddot, yddot};
}

Vector hill4bp_regularized_field(double, const Vector& state, const Hill4BPParameters& params) {
    if (state.size() != 4) {
        throw std::invalid_argument("regularized Hill 4BP field expects a 4-vector");
    }

    double r = state[0];
    const double theta = state[1];
    const double v = state[2];
    const double w = state[3];

    if (r < 0.0) {
        r = 0.0;
    }

    const double radial_power = 2.0 - params.gamma * (params.nu + 2.0);
    return {
        (params.beta + 1.0) * v * r,
        w - r,
        params.beta * v * v + w * w - params.alpha * params.c
            - params.nu * std::pow(r, radial_power)
            - 2.0 * params.A * r * r * std::cos(theta) * std::cos(theta)
            - 2.0 * params.B * r * r * std::sin(theta) * std::sin(theta),
        (params.beta - 1.0) * v * w
            + 2.0 * (params.A - params.B) * r * r * std::sin(theta) * std::cos(theta),
    };
}

Vector hill4bp_sqrt_radius_field(
    double, const Vector& state, const Hill4BPParameters& params) {
    if (state.size() != 4) {
        throw std::invalid_argument("sqrt-radius Hill 4BP field expects a 4-vector");
    }

    const double R = state[0];
    const double theta = state[1];
    const double v = state[2];
    const double w = state[3];
    if (R < 0.0) {
        throw std::runtime_error("sqrt-radius coordinate R must be nonnegative");
    }

    const double paper_r = std::pow(R, params.alpha + 2.0);
    const double radial_term = params.nu * std::pow(
        R, 2.0 * (params.alpha - params.nu));
    const double anisotropy_scale =
        std::pow(R, 2.0 * (params.alpha + 2.0));
    return {
        0.5 * v * R,
        w - paper_r,
        params.beta * v * v + w * w - params.alpha * params.c
            - radial_term
            - 2.0 * params.A * anisotropy_scale
                * std::cos(theta) * std::cos(theta)
            - 2.0 * params.B * anisotropy_scale
                * std::sin(theta) * std::sin(theta),
        (params.beta - 1.0) * v * w
            + 2.0 * (params.A - params.B) * anisotropy_scale
                * std::sin(theta) * std::cos(theta),
    };
}

double hill4bp_effective_potential(const Vec2& q, const Hill4BPParameters& params) {
    const double rho = std::sqrt(q[0] * q[0] + q[1] * q[1]);
    if (rho == 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return 0.5 * (params.lambda2 * q[0] * q[0] + params.lambda1 * q[1] * q[1])
        + 1.0 / rho + params.c / std::pow(rho, 3.0);
}

double hill4bp_hamiltonian(const State4& state, const Hill4BPParameters& params) {
    const Vec2 q{{state[0], state[1]}};
    return 0.5 * (state[2] * state[2] + state[3] * state[3])
        - hill4bp_effective_potential(q, params);
}

double hill4bp_jacobi_constant(const State4& state, const Hill4BPParameters& params) {
    return -2.0 * hill4bp_hamiltonian(state, params);
}

std::pair<double, bool> hill4bp_speed(const Vec2& q, double h, const Hill4BPParameters& params) {
    const double speed_squared = 2.0 * (h + hill4bp_effective_potential(q, params));
    const bool allowed = speed_squared >= 0.0 && std::isfinite(speed_squared);
    if (!allowed) {
        return {std::numeric_limits<double>::quiet_NaN(), false};
    }
    return {std::sqrt(std::max(speed_squared, 0.0)), true};
}

State4 hill4bp_cartesian_to_mcgehee(const State4& state, const Hill4BPParameters& params) {
    const double q1 = state[0];
    const double q2 = state[1];
    const double xdot = state[2];
    const double ydot = state[3];
    const double rho = std::sqrt(q1 * q1 + q2 * q2);
    if (rho == 0.0) {
        throw std::runtime_error("McGehee coordinates require a nonzero Cartesian position");
    }

    const double theta = std::atan2(q2, q1);
    const double r = std::pow(rho, 1.0 / params.gamma);
    const double p1 = xdot - q2;
    const double p2 = ydot + q1;
    const double scale = std::pow(r, params.gamma * params.beta);

    const double v = scale * (p1 * std::cos(theta) + p2 * std::sin(theta));
    const double w = scale * (-p1 * std::sin(theta) + p2 * std::cos(theta));
    return {{r, theta, v, w}};
}

State4 hill4bp_mcgehee_to_cartesian(const State4& state, const Hill4BPParameters& params) {
    const double r = std::max(state[0], 0.0);
    const double theta = state[1];
    const double v = state[2];
    const double w = state[3];

    const double rho = std::pow(r, params.gamma);
    const double ctheta = std::cos(theta);
    const double stheta = std::sin(theta);

    const double q1 = rho * ctheta;
    const double q2 = rho * stheta;
    const double pscale = r == 0.0
        ? std::numeric_limits<double>::infinity()
        : std::pow(r, -params.gamma * params.beta);

    const double p1 = pscale * (v * ctheta - w * stheta);
    const double p2 = pscale * (v * stheta + w * ctheta);
    const double xdot = p1 + q2;
    const double ydot = p2 - q1;

    return {{q1, q2, xdot, ydot}};
}

State4 hill4bp_cartesian_to_sqrt_radius(
    const State4& state, const Hill4BPParameters& params) {
    const double q1 = state[0];
    const double q2 = state[1];
    const double xdot = state[2];
    const double ydot = state[3];
    const double rho = std::hypot(q1, q2);
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error(
            "sqrt-radius coordinates require a nonzero finite Cartesian position");
    }

    const double R = std::sqrt(rho);
    const double theta = std::atan2(q2, q1);
    const double p1 = xdot - q2;
    const double p2 = ydot + q1;
    const double scale = std::pow(R, params.alpha);
    const double v = scale * (p1 * std::cos(theta) + p2 * std::sin(theta));
    const double w = scale * (-p1 * std::sin(theta) + p2 * std::cos(theta));
    return {{R, theta, v, w}};
}

State4 hill4bp_sqrt_radius_to_cartesian(
    const State4& state, const Hill4BPParameters& params) {
    const double R = state[0];
    if (!(R > 0.0) || !std::isfinite(R)) {
        throw std::runtime_error(
            "sqrt-radius Cartesian conversion requires positive finite R");
    }

    const double theta = state[1];
    const double v = state[2];
    const double w = state[3];
    const double rho = R * R;
    const double ctheta = std::cos(theta);
    const double stheta = std::sin(theta);
    const double pscale = std::pow(R, -params.alpha);
    const double q1 = rho * ctheta;
    const double q2 = rho * stheta;
    const double p1 = pscale * (v * ctheta - w * stheta);
    const double p2 = pscale * (v * stheta + w * ctheta);
    return {{q1, q2, p1 + q2, p2 - q1}};
}

double hill4bp_section_event(double, const Vector& state, const Hill4BPParameters& params) {
    if (state.size() != 4) {
        throw std::invalid_argument("Hill 4BP section event expects a 4-vector");
    }
    const double r = std::max(state[0], 0.0);
    return std::pow(r, params.gamma) * std::sin(state[1]);
}

double hill4bp_sqrt_radius_section_event(
    double, const Vector& state, const Hill4BPParameters&) {
    if (state.size() != 4) {
        throw std::invalid_argument(
            "sqrt-radius Hill 4BP section event expects a 4-vector");
    }
    return std::sin(state[1]);
}

Vector to_vector(const State4& state) {
    return Vector(state.begin(), state.end());
}

State4 to_state4(const Vector& values) {
    if (values.size() != 4) {
        throw std::invalid_argument("expected a 4-vector");
    }
    State4 state{};
    std::copy(values.begin(), values.end(), state.begin());
    return state;
}

}  // namespace crtbp
