// McGehee-regularized planar Hill 4-body problem (oblate bodies) Poincare
// section, built on CAPD. This is the model used by the MATLAB workflow and
// by Belbruno, Gidea & Lam, "Regularization of the Hill four-body problem
// with oblate bodies" (2023), eq. (4.9).
//
// The paper's radial coordinate r has a fractional r^0.8 term for alpha=3.
// That term is only C^0 at collision and is a poor fit for a Taylor solver.
// We therefore integrate R = r^(1/5) = sqrt(rho), where rho is Cartesian
// radius.  The mathematically equivalent, collision-smooth state is
// s = (R, theta, v, w), with dt = R^5 dtau:
//   R'     = 0.5 v R
//   theta' = w - R^5
//   v'     = 1.5 v^2 + w^2 - 3c - R^4
//                 - 2A R^10 cos^2(theta) - 2B R^10 sin^2(theta)
//   w'     = 0.5 v w + 2(A-B) R^10 sin(theta) cos(theta).
//
// Section: Cartesian x2 = 0 with x2dot > 0. Since x2 = R^2 sin(theta),
// sin(theta)=0 is the section away from collision. CAPD is deliberately asked
// for crossings in BOTH directions and this driver filters on physical
// x2dot. This is important: CAPD 6.1's direction-specific "leave section"
// loop does not apply maxReturnTime and can otherwise run forever on a
// collision-asymptotic orbit.
//
// This executable is the dense floating-point reconnaissance driver.  The
// separate hill4bp_validate_capd executable uses CAPD I* types, C0 sets, and
// explicit interval event isolation for proof-producing runs.

#include "capd/capdlib.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using capd::DMap;
using capd::DOdeSolver;
using capd::DNonlinearSection;
using capd::DPoincareMap;
using capd::DVector;

namespace {

// The configured study fixes alpha=3, hence gamma=2/(alpha+2)=2/5.
const double kPaperGamma = 2.0 / 5.0;

// CAPD's parser rejects inline scientific-notation constants, so every numeric
// coefficient is passed as a named parameter instead.
const char* kHill4bpField =
    "par:A,B,ac;"
    "var:R,th,v,w;"
    "fun:"
    "0.5*v*R,"
    "w-R^5,"
    "1.5*v^2+w^2-ac-R^4-2*A*R^10*cos(th)^2-2*B*R^10*sin(th)^2,"
    "0.5*v*w+2*(A-B)*R^10*sin(th)*cos(th);";

struct Params {
    double mu = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double nu = 0.0;
    double alpha = 0.0;
    double lambda1 = 0.0, lambda2 = 0.0, A = 0.0, B = 0.0, c = 0.0;
};

struct PoincareDefaults {
    double h_e1_scale = 0.0;
    double x_begin = 0.0;
    double x_end = 0.0;
    int x_count = 0;
    double xdot_lower = 0.0;
    double xdot_upper = 0.0;
    int xdot_count = 0;
    int iterates = 0;
    double y_section = 0.0;
    double rel_tol = 0.0;
    double abs_tol = 0.0;
    double max_step = 0.0;
};

struct CartesianState {
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    double xdot = std::numeric_limits<double>::quiet_NaN();
    double ydot = std::numeric_limits<double>::quiet_NaN();
};

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

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

std::string normalize_neck(const std::string& value) {
    const std::string normalized = lowercase(value);
    if (normalized == "e1") return "E1";
    if (normalized == "e3") return "E3";
    throw std::runtime_error("neck must be E1 or E3");
}

std::string normalize_side(const std::string& value) {
    const std::string normalized = lowercase(value);
    if (normalized == "closed" || normalized == "open") return normalized;
    throw std::runtime_error("side must be closed or open");
}

bool parse_boolean(const std::string& value, const std::string& key) {
    const std::string normalized = lowercase(value);
    if (normalized == "1" || normalized == "true") return true;
    if (normalized == "0" || normalized == "false") return false;
    throw std::runtime_error(key + " must be 0, 1, false, or true");
}

double parse_double(const std::string& value, const std::string& key) {
    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid numeric value for " + key + ": " + value);
    }
    if (consumed != value.size() || !std::isfinite(parsed)) {
        throw std::runtime_error("invalid finite numeric value for " + key + ": " + value);
    }
    return parsed;
}

int parse_int(const std::string& value, const std::string& key) {
    std::size_t consumed = 0;
    long parsed = 0;
    try {
        parsed = std::stol(value, &consumed, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer value for " + key + ": " + value);
    }
    if (consumed != value.size()
        || parsed < static_cast<long>(std::numeric_limits<int>::min())
        || parsed > static_cast<long>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("invalid integer value for " + key + ": " + value);
    }
    return static_cast<int>(parsed);
}

std::string normalize_velocity_mode(const std::string& value) {
    const std::string normalized = lowercase(value);
    if (normalized == "direction" || normalized == "fraction"
        || normalized == "xdot") {
        return normalized;
    }
    throw std::runtime_error("velocity_mode must be direction, fraction, or xdot");
}

std::string normalize_neck_branch(const std::string& value) {
    const std::string normalized = lowercase(value);
    if (normalized == "positive" || normalized == "negative") return normalized;
    throw std::runtime_error("neck_branch must be positive or negative");
}

void print_help(const char* executable) {
    std::cout
        << "Usage: " << executable << " key=value ...\n\n"
        << "CAPD Poincare sections for the regularized planar Hill 4BP model.\n"
        << "Keys are case-sensitive; unknown or malformed keys are errors.\n\n"
        << "Energy: h, h_e1_scale, h_e1_offset, or positive delta (choose one),\n"
        << "        with neck=E1|E3 and side=closed|open.\n"
        << "Grid:   x_begin, x_end, x_count, xdot_lower, xdot_upper, xdot_count,\n"
        << "        velocity_mode=direction|fraction|xdot, y_section, iterates,\n"
        << "        start_index, end_index.\n"
        << "Neck:   neck_window_sigma (E1 + delta only),\n"
        << "        neck_branch=positive|negative.\n"
        << "CAPD:   solver_order, abs_tol, rel_tol, max_step, max_return_time,\n"
        << "        max_half_crossings, min_transversality, blow_up_norm.\n"
        << "Stops:  collision_radius and outer_radius are Cartesian radii; the\n"
        << "        latter is a study-domain limit, not proof of escape.\n"
        << "Files:  output (two-column plot cloud), returns (full return table),\n"
        << "        outcomes, zvb, svg, metadata; empty disables optional files.\n"
        << "        svg_stride limits preview memory only.\n"
        << "Other:  progress_every, iterate_progress_every, dry_run=0|1.\n\n"
        << "Example:\n  " << executable
        << " neck=E1 side=closed delta=1e-2 neck_window_sigma=6 dry_run=1\n";
}

const char* classify_neck_side(double h, double hcrit) {
    if (h < hcrit) return "closed";
    if (h > hcrit) return "open";
    return "critical_nontransverse";
}

double e1_neck_curvature(const Params& p, double rho) {
    return p.lambda2 + 2.0 / std::pow(rho, 3.0)
        + 12.0 * p.c / std::pow(rho, 5.0);
}

std::string csv_escape(std::string value) {
    bool needs_quotes = false;
    std::string escaped;
    escaped.reserve(value.size() + 8);
    for (char ch : value) {
        if (ch == '"') {
            escaped += "\"\"";
            needs_quotes = true;
        } else if (ch == '\n' || ch == '\r') {
            escaped += ' ';
            needs_quotes = true;
        } else {
            escaped += ch;
            needs_quotes = needs_quotes || ch == ',';
        }
    }
    return needs_quotes ? '"' + escaped + '"' : escaped;
}

bool finite_vector(const DVector& state) {
    if (state.dimension() != 4) return false;
    for (int i = 0; i < 4; ++i) {
        if (!std::isfinite(state[i])) return false;
    }
    return true;
}

std::string parameter_file() {
    const std::vector<std::string> candidates{
        "Hill_4BP_parameters.cfg",
        "../Hill_4BP_parameters.cfg",
        "../../Hill_4BP_parameters.cfg",
        "../../../Hill_4BP_parameters.cfg",
    };

    for (const std::string& candidate : candidates) {
        std::ifstream in(candidate);
        if (in) return candidate;
    }

    throw std::runtime_error("could not find Hill_4BP_parameters.cfg");
}

std::map<std::string, double> read_parameter_file() {
    const std::string path = parameter_file();
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open parameter file: " + path);
    }

    std::map<std::string, double> values;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        line = trim(line);
        if (line.empty()) continue;

        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::runtime_error("invalid parameter line " + std::to_string(line_number));
        }

        const std::string key = trim(line.substr(0, equals));
        const std::string raw_value = trim(line.substr(equals + 1));
        std::size_t consumed = 0;
        const double value = std::stod(raw_value, &consumed);
        if (key.empty() || consumed != raw_value.size()) {
            throw std::runtime_error("invalid parameter value for key: " + key);
        }
        if (!std::isfinite(value)) {
            throw std::runtime_error("parameter value must be finite for key: " + key);
        }
        if (!values.emplace(key, value).second) {
            throw std::runtime_error("duplicate parameter key: " + key);
        }
    }

    return values;
}

double required_value(const std::map<std::string, double>& values, const std::string& key) {
    const auto it = values.find(key);
    if (it == values.end()) {
        throw std::runtime_error("missing parameter key: " + key);
    }
    return it->second;
}

double optional_value(const std::map<std::string, double>& values,
                      const std::string& key,
                      double default_value) {
    const auto it = values.find(key);
    return it == values.end() ? default_value : it->second;
}

int positive_int_value(const std::map<std::string, double>& values, const std::string& key) {
    const double value = required_value(values, key);
    const double rounded = std::round(value);
    if (value < 1.0 || std::abs(value - rounded) > std::numeric_limits<double>::epsilon() * value) {
        throw std::runtime_error("parameter key must be a positive integer: " + key);
    }
    return static_cast<int>(rounded);
}

Params load_params(const std::map<std::string, double>& values) {
    Params p;
    p.mu = required_value(values, "mu");
    p.u1 = required_value(values, "u1");
    p.u2 = required_value(values, "u2");
    p.c1 = required_value(values, "c1");
    p.c2 = required_value(values, "c2");
    p.c3 = required_value(values, "c3");
    p.nu = required_value(values, "nu");
    p.alpha = required_value(values, "alpha");
    return p;
}

PoincareDefaults load_poincare_defaults(const std::map<std::string, double>& values) {
    PoincareDefaults defaults;
    defaults.h_e1_scale = required_value(values, "poincare_h_e1_scale");
    defaults.x_begin = required_value(values, "poincare_x_begin");
    defaults.x_end = required_value(values, "poincare_x_end");
    defaults.x_count = positive_int_value(values, "poincare_k");
    defaults.xdot_lower = required_value(values, "poincare_xdot_lower");
    defaults.xdot_upper = required_value(values, "poincare_xdot_upper");
    defaults.xdot_count = positive_int_value(values, "poincare_l");
    defaults.iterates = positive_int_value(values, "poincare_iterates");
    defaults.y_section = required_value(values, "poincare_y_section");
    defaults.rel_tol = required_value(values, "poincare_rel_tol");
    defaults.abs_tol = required_value(values, "poincare_abs_tol");
    defaults.max_step = required_value(values, "poincare_max_step");
    return defaults;
}

// Eigenvalues of the quadratic part (paper eq. 2.4 / MATLAB hill4bpParameters).
void finalize(Params& p, const std::map<std::string, double>& values) {
    const double mu = p.mu, u1 = p.u1, u2 = p.u2;
    const double Delta = std::pow(mu * u1 * u1 * u1 + (1 - mu) * u2 * u2 * u2, 2)
        - mu * (1 - mu) * u1 * u2
            * (-std::pow(u1, 4) - std::pow(u2, 4) + 2 * u1 * u1 + 2 * u2 * u2
               + 2 * u1 * u1 * u2 * u2 - 1);
    const double base = 2 - 2 * (1 - mu) / std::pow(u1, 5) - 2 * mu / std::pow(u2, 5)
        + 3 * (1 - mu) / std::pow(u1, 3) + 3 * mu / std::pow(u2, 3);
    const double spread = 3 * std::sqrt(Delta) / (std::pow(u1, 3) * std::pow(u2, 3));
    p.lambda1 = optional_value(values, "lambda1", 0.5 * (base - spread));
    p.lambda2 = optional_value(values, "lambda2", 0.5 * (base + spread));
    p.A = (1 - p.lambda2) / 2;
    p.B = (1 - p.lambda1) / 2;
    p.c = -p.c3;   // c = -c3 > 0
}

double equilibrium_radius(double lambda, const Params& p) {
    const double power = 2.0 - kPaperGamma * (p.nu + 2.0);
    auto f = [&](double r) {
        return lambda * r * r - p.nu * std::pow(r, power) - p.alpha * p.c;
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

// Planar effective potential (MATLAB hill4bpEffectivePotential).
double effective_potential(double x, double y, const Params& p) {
    const double rho = std::sqrt(x * x + y * y);
    return 0.5 * (p.lambda2 * x * x + p.lambda1 * y * y)
        + 1.0 / rho + p.c / (rho * rho * rho);
}

// Omega(x,y)-Omega(xe,0), evaluated without subtracting two O(1) values.
// The x-axis expression uses the E1 equilibrium relation and is exact for
// either symmetric branch because the potential is even in x.
double e1_potential_difference(double x, double y, double xe, const Params& p) {
    const double ax = std::abs(x);
    if (!(ax > 0.0) || !(xe > 0.0)) {
        return effective_potential(x, y, p) - effective_potential(xe, 0.0, p);
    }

    const double d = ax - xe;
    const double ax2 = ax * ax;
    const double ax3 = ax2 * ax;
    const double xe2 = xe * xe;
    const double xe3 = xe2 * xe;
    const double xe4 = xe3 * xe;
    const double numerator =
        p.lambda2 * xe3 * ax3
        + 2.0 * p.lambda2 * xe4 * ax2
        + 4.0 * p.c * ax + 2.0 * p.c * xe;
    const double x_axis_difference =
        d * d * numerator / (2.0 * ax3 * xe3);

    if (y == 0.0) return x_axis_difference;
    const double rho = std::hypot(ax, y);
    const double drho = y * y / (rho + ax);
    const double inverse_difference = -drho / (rho * ax);
    const double rho2 = rho * rho;
    const double rho3 = rho2 * rho;
    const double inverse_cube_difference =
        -drho * (rho2 + rho * ax + ax2) / (rho3 * ax3);
    return x_axis_difference + 0.5 * p.lambda1 * y * y
        + inverse_difference + p.c * inverse_cube_difference;
}

double initial_speed_squared(double x,
                             double y,
                             double h,
                             const Params& p,
                             bool stable_e1_delta,
                             double signed_delta,
                             double e1_rho) {
    if (stable_e1_delta) {
        return 2.0 * (
            signed_delta + e1_potential_difference(x, y, e1_rho, p));
    }
    return 2.0 * (h + effective_potential(x, y, p));
}

double equilibrium_position(double lambda, const Params& p) {
    const double r = equilibrium_radius(lambda, p);
    return std::pow(r, kPaperGamma);
}

double equilibrium_energy(double lambda, bool x_axis, const Params& p) {
    const double rho = equilibrium_position(lambda, p);
    return x_axis
        ? -effective_potential(rho, 0.0, p)
        : -effective_potential(0.0, rho, p);
}

double e1_energy(const Params& p) {
    return equilibrium_energy(p.lambda2, true, p);
}

double e3_energy(const Params& p) {
    return equilibrium_energy(p.lambda1, false, p);
}

// Cartesian (x1,x2,x1dot,x2dot) -> smooth McGehee (R,theta,v,w),
// where R=sqrt(Cartesian radius).
DVector to_smooth_mcgehee(double q1, double q2, double xdot, double ydot) {
    const double rho = std::sqrt(q1 * q1 + q2 * q2);
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("smooth McGehee coordinates require nonzero finite radius");
    }
    const double theta = std::atan2(q2, q1);
    const double R = std::sqrt(rho);
    const double p1 = xdot - q2;   // canonical momenta
    const double p2 = ydot + q1;
    const double scale = R * R * R;
    DVector s(4);
    s[0] = R;
    s[1] = theta;
    s[2] = scale * (p1 * std::cos(theta) + p2 * std::sin(theta));
    s[3] = scale * (-p1 * std::sin(theta) + p2 * std::cos(theta));
    return s;
}

CartesianState from_smooth_mcgehee(const DVector& s) {
    CartesianState result;
    if (!finite_vector(s)) return result;

    const double R = s[0];
    if (!(R > 0.0)) return result;

    const double theta = s[1];
    const double v = s[2];
    const double w = s[3];
    const double cth = std::cos(theta);
    const double sth = std::sin(theta);
    const double rho = R * R;
    const double R3 = R * R * R;
    if (!(R3 > 0.0) || !std::isfinite(R3)) return result;

    const double p1 = (v * cth - w * sth) / R3;
    const double p2 = (v * sth + w * cth) / R3;
    result.x = rho * cth;
    result.y = rho * sth;
    result.xdot = p1 + result.y;
    result.ydot = p2 - result.x;
    return result;
}

bool finite_cartesian(const CartesianState& state) {
    return std::isfinite(state.x) && std::isfinite(state.y)
        && std::isfinite(state.xdot) && std::isfinite(state.ydot);
}

double smooth_energy_error(const DVector& s, const Params& p, double h) {
    if (!finite_vector(s) || !(s[0] > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double R = s[0];
    const double theta = s[1];
    const double v = s[2];
    const double w = s[3];
    const double R2 = R * R;
    const double R4 = R2 * R2;
    const double R5 = R4 * R;
    const double R6 = R5 * R;
    const double R10 = R5 * R5;
    const double cth = std::cos(theta);
    const double sth = std::sin(theta);
    const double scaled_residual =
        0.5 * (v * v + w * w) - p.c - w * R5 - R4
        + p.A * R10 * cth * cth + p.B * R10 * sth * sth
        - h * R6;
    return scaled_residual / R6;
}

std::string classify_solver_stop(const DVector& terminal,
                                 const std::string& message,
                                 double collision_radius,
                                 double outer_radius) {
    if (finite_vector(terminal) && terminal[0] >= 0.0) {
        const double rho = terminal[0] * terminal[0];
        if (rho <= collision_radius) {
            return terminal[2] < 0.0 ? "collision" : "inner_radius_limit";
        }
        if (rho >= outer_radius) return "outer_limit";
    }
    const std::string normalized = lowercase(message);
    if (normalized.find("max return time") != std::string::npos) return "timeout";
    if (normalized.find("blow up") != std::string::npos) return "numerical_failure";
    return "numerical_failure";
}

std::string classify_event_surface(const DVector& state,
                                   double collision_R,
                                   double outer_R) {
    if (!finite_vector(state)) return "numerical_failure";
    constexpr double surface_slack = 1e-9;
    if (state[0] <= collision_R * (1.0 + surface_slack)) return "collision";
    if (state[0] >= outer_R * (1.0 - surface_slack)) return "outer_limit";
    const double section_residual = std::abs(std::sin(state[1]));
    const double collision_residual =
        std::abs(state[0] - collision_R) / std::max(1.0, collision_R);
    const double outer_residual =
        std::abs(state[0] - outer_R) / std::max(1.0, outer_R);
    if (collision_residual < section_residual
        && collision_residual <= outer_residual) {
        return "collision";
    }
    if (outer_residual < section_residual
        && outer_residual < collision_residual) {
        return "outer_limit";
    }
    return "section";
}

}  // namespace

int run(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "help" || argument == "--help" || argument == "-h") {
            print_help(argv[0]);
            return 0;
        }
    }

    const std::map<std::string, double> shared_values = read_parameter_file();
    Params p = load_params(shared_values);
    const PoincareDefaults defaults = load_poincare_defaults(shared_values);
    double h = 0.0;
    bool h_overridden = false;
    double h_e1_scale = defaults.h_e1_scale;
    bool h_e1_scale_overridden = false;
    double h_e1_offset = 0.0;
    bool h_e1_offset_overridden = false;
    double delta = 0.0;
    bool delta_overridden = false;
    std::string neck = "E1";
    std::string side = "closed";
    bool dry_run = false;
    double x_begin = defaults.x_begin, x_end = defaults.x_end;
    int x_count = defaults.x_count;
    double xdot_lower = defaults.xdot_lower, xdot_upper = defaults.xdot_upper;
    int xdot_count = defaults.xdot_count;
    int iterates = defaults.iterates;
    int start_index = 0;
    int end_index = -1;
    int progress_every = 1;
    int iterate_progress_every = 100;
    int svg_stride = 1;
    int solver_order = 20;
    double abs_tol = defaults.abs_tol;
    double rel_tol = defaults.rel_tol;
    double max_step = defaults.max_step;
    double max_return_time = 100.0;
    int max_half_crossings = 4;
    double min_transversality = 1e-12;
    double blow_up_norm = 1e10;
    double collision_radius = 1e-6;
    double outer_radius = 5.0;
    double y_section = defaults.y_section;
    std::string velocity_mode = "direction";
    double neck_window_sigma = 0.0;
    std::string neck_branch = "positive";
    bool x_begin_overridden = false;
    bool x_end_overridden = false;
    std::string output = "hill4bp_capd.csv";
    std::string returns_path = "hill4bp_capd_returns.csv";
    std::string outcomes_path = "hill4bp_capd_outcomes.csv";
    std::string zvb_output = "zero_velocity_boundary.csv";
    std::string svg_path = "hill4bp_capd.svg";
    std::string metadata_path = "hill4bp_capd_metadata.txt";

    std::set<std::string> supplied_keys;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const std::size_t eq = a.find('=');
        if (eq == std::string::npos || eq == 0) {
            throw std::runtime_error(
                "arguments must use case-sensitive key=value syntax; got: " + a);
        }
        const std::string k = a.substr(0, eq);
        const std::string val = a.substr(eq + 1);
        if (!supplied_keys.insert(k).second) {
            throw std::runtime_error("duplicate argument key: " + k);
        }
        if (k == "h") {
            h = parse_double(val, k);
            h_overridden = true;
        }
        else if (k == "h_e1_scale") {
            h_e1_scale = parse_double(val, k);
            h_e1_scale_overridden = true;
        }
        else if (k == "h_e1_offset") {
            h_e1_offset = parse_double(val, k);
            h_e1_offset_overridden = true;
        }
        else if (k == "delta") {
            delta = parse_double(val, k);
            delta_overridden = true;
        }
        else if (k == "neck") neck = normalize_neck(val);
        else if (k == "side") side = normalize_side(val);
        else if (k == "dry_run") dry_run = parse_boolean(val, "dry_run");
        else if (k == "x_begin") {
            x_begin = parse_double(val, k);
            x_begin_overridden = true;
        }
        else if (k == "x_end") {
            x_end = parse_double(val, k);
            x_end_overridden = true;
        }
        else if (k == "x_count") x_count = parse_int(val, k);
        else if (k == "xdot_lower") xdot_lower = parse_double(val, k);
        else if (k == "xdot_upper") xdot_upper = parse_double(val, k);
        else if (k == "xdot_count") xdot_count = parse_int(val, k);
        else if (k == "iterates") iterates = parse_int(val, k);
        else if (k == "start_index") start_index = parse_int(val, k);
        else if (k == "end_index") end_index = parse_int(val, k);
        else if (k == "progress_every") progress_every = parse_int(val, k);
        else if (k == "iterate_progress_every") {
            iterate_progress_every = parse_int(val, k);
        }
        else if (k == "svg_stride") svg_stride = parse_int(val, k);
        else if (k == "solver_order") solver_order = parse_int(val, k);
        else if (k == "abs_tol") abs_tol = parse_double(val, k);
        else if (k == "rel_tol") rel_tol = parse_double(val, k);
        else if (k == "max_step") max_step = parse_double(val, k);
        else if (k == "max_return_time") max_return_time = parse_double(val, k);
        else if (k == "max_half_crossings") max_half_crossings = parse_int(val, k);
        else if (k == "min_transversality") {
            min_transversality = parse_double(val, k);
        }
        else if (k == "blow_up_norm") blow_up_norm = parse_double(val, k);
        else if (k == "collision_radius") collision_radius = parse_double(val, k);
        else if (k == "outer_radius") outer_radius = parse_double(val, k);
        else if (k == "y_section") y_section = parse_double(val, k);
        else if (k == "velocity_mode") velocity_mode = normalize_velocity_mode(val);
        else if (k == "neck_window_sigma") neck_window_sigma = parse_double(val, k);
        else if (k == "neck_branch") neck_branch = normalize_neck_branch(val);
        else if (k == "output") output = val;
        else if (k == "returns") returns_path = val;
        else if (k == "outcomes") outcomes_path = val;
        else if (k == "zvb") zvb_output = val;
        else if (k == "svg") svg_path = val;
        else if (k == "metadata") metadata_path = val;
        else {
            throw std::runtime_error(
                "unknown argument key (keys are case-sensitive): " + k);
        }
    }
    if (x_count < 1) {
        throw std::runtime_error("x_count must be positive");
    }
    if (xdot_count < 1) {
        throw std::runtime_error("xdot_count must be positive");
    }
    if (iterates < 1) {
        throw std::runtime_error("iterates must be positive");
    }
    if (end_index < 0) {
        end_index = x_count - 1;
    }
    if (start_index < 0 || end_index < start_index || end_index >= x_count) {
        throw std::runtime_error("invalid start_index/end_index range");
    }
    if (progress_every < 1) {
        progress_every = 1;
    }
    if (iterate_progress_every < 1) {
        iterate_progress_every = iterates + 1;
    }
    if (svg_stride < 1) {
        throw std::runtime_error("svg_stride must be positive");
    }
    if (solver_order < 2) {
        throw std::runtime_error("solver_order must be at least 2");
    }
    if (!(abs_tol > 0.0) || !(rel_tol > 0.0)) {
        throw std::runtime_error("abs_tol and rel_tol must be positive and finite");
    }
    if (!(max_step > 0.0)) {
        throw std::runtime_error("max_step must be positive and finite");
    }
    if (!(max_return_time > 0.0) || !std::isfinite(max_return_time)) {
        throw std::runtime_error("max_return_time must be positive and finite");
    }
    if (max_half_crossings < 2) {
        throw std::runtime_error("max_half_crossings must be at least 2");
    }
    if (min_transversality < 0.0) {
        throw std::runtime_error("min_transversality must be nonnegative and finite");
    }
    if (!(blow_up_norm > 0.0)) {
        throw std::runtime_error("blow_up_norm must be positive and finite");
    }
    if (!(collision_radius > 0.0) || !(outer_radius > collision_radius)) {
        throw std::runtime_error(
            "require 0 < collision_radius < outer_radius, both finite");
    }
    if (neck_window_sigma < 0.0) {
        throw std::runtime_error("neck_window_sigma must be nonnegative and finite");
    }
    if (output.empty()) {
        throw std::runtime_error("output must not be empty");
    }
    const int explicit_energy_choices = static_cast<int>(h_overridden)
        + static_cast<int>(h_e1_scale_overridden)
        + static_cast<int>(h_e1_offset_overridden)
        + static_cast<int>(delta_overridden);
    if (explicit_energy_choices > 1) {
        throw std::runtime_error(
            "specify only one of h, h_e1_scale, h_e1_offset, or delta");
    }
    if (h_overridden && !std::isfinite(h)) {
        throw std::runtime_error("h must be finite");
    }
    if (h_e1_scale_overridden && !std::isfinite(h_e1_scale)) {
        throw std::runtime_error("h_e1_scale must be finite");
    }
    if (h_e1_offset_overridden && !std::isfinite(h_e1_offset)) {
        throw std::runtime_error("h_e1_offset must be finite");
    }
    if (delta_overridden && (!(delta > 0.0) || !std::isfinite(delta))) {
        throw std::runtime_error("delta must be positive, nonzero, and finite");
    }
    finalize(p, shared_values);
    if (std::abs(p.alpha - 3.0) > 1e-14 || std::abs(p.nu - 1.0) > 1e-14) {
        throw std::runtime_error(
            "this CAPD field is specialized to alpha=3 and nu=1");
    }
    if (p.c1 != 0.0 || p.c2 != 0.0) {
        throw std::runtime_error(
            "nonzero c1/c2 require extending the CAPD field coefficients first");
    }
    const double e1_rho = equilibrium_position(p.lambda2, p);
    const double e3_rho = equilibrium_position(p.lambda1, p);
    const double e1_h = e1_energy(p);
    const double e3_h = e3_energy(p);
    const double hcrit = neck == "E3" ? e3_h : e1_h;
    if (delta_overridden) {
        h = side == "open" ? hcrit + delta : hcrit - delta;
        if (h == hcrit) {
            throw std::runtime_error(
                "delta is too small to resolve from the critical energy in double precision");
        }
    } else if (h_e1_offset_overridden) {
        h = e1_h + h_e1_offset;
    } else if (!h_overridden) {
        h = h_e1_scale * e1_h;
    }
    const char* energy_mode = h_overridden ? "direct_h"
        : (delta_overridden ? "neck_delta"
        : (h_e1_offset_overridden ? "h_e1_offset"
        : (h_e1_scale_overridden ? "h_e1_scale" : "default_e1_scale")));
    const char* actual_neck_side = classify_neck_side(h, hcrit);
    const bool stable_e1_delta = delta_overridden && neck == "E1";
    const double signed_delta =
        stable_e1_delta ? (side == "open" ? delta : -delta) : 0.0;
    double neck_length_scale = std::numeric_limits<double>::quiet_NaN();
    if (neck_window_sigma > 0.0) {
        if (!delta_overridden || neck != "E1") {
            throw std::runtime_error(
                "neck_window_sigma requires delta mode with neck=E1");
        }
        if (x_begin_overridden || x_end_overridden) {
            throw std::runtime_error(
                "do not combine neck_window_sigma with x_begin or x_end");
        }
        const double curvature = e1_neck_curvature(p, e1_rho);
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            throw std::runtime_error("could not resolve the E1 neck curvature");
        }
        neck_length_scale = std::sqrt(2.0 * delta / curvature);
        const double center = neck_branch == "negative" ? -e1_rho : e1_rho;
        x_begin = center - neck_window_sigma * neck_length_scale;
        x_end = center + neck_window_sigma * neck_length_scale;
    }

    if (dry_run) {
        std::cerr << std::setprecision(17)
                  << "dry_run=1\n"
                  << "lambda1=" << p.lambda1 << " lambda2=" << p.lambda2
                  << " A=" << p.A << " B=" << p.B << " c=" << p.c
                  << " c1=" << p.c1 << " c2=" << p.c2
                  << " c3=" << p.c3 << '\n'
                  << "E1/E2 |x|=" << e1_rho << " hcrit=" << e1_h
                  << " (Ccrit=" << -2.0 * e1_h << ")\n"
                  << "E3/E4 |y|=" << e3_rho << " hcrit=" << e3_h
                  << " (Ccrit=" << -2.0 * e3_h << ")\n"
                  << "energy_mode=" << energy_mode << " neck=" << neck
                  << " hcrit=" << hcrit
                  << " delta=" << (delta_overridden ? delta : 0.0)
                  << " side=" << (delta_overridden ? side : "not_used")
                  << " actual_neck_side=" << actual_neck_side << '\n'
                  << "h=" << h << " (C=" << -2.0 * h << ")\n"
                  << "x_count=" << x_count << " xdot_count=" << xdot_count
                  << " iterates=" << iterates
                  << " index range=[" << start_index << ',' << end_index << ']'
                  << " velocity_mode=" << velocity_mode
                  << " requested_returns="
                  << static_cast<long long>(end_index - start_index + 1)
                     * static_cast<long long>(xdot_count)
                     * static_cast<long long>(iterates)
                  << " svg_stride=" << svg_stride << '\n'
                  << "x_begin=" << x_begin << " x_end=" << x_end
                  << " neck_window_sigma=" << neck_window_sigma
                  << " neck_length_scale=" << neck_length_scale << '\n'
                  << "solver_order=" << solver_order
                  << " abs_tol=" << abs_tol << " rel_tol=" << rel_tol
                  << " max_step=" << max_step
                  << " max_return_time=" << max_return_time
                  << " max_half_crossings=" << max_half_crossings << '\n'
                  << "collision_radius=" << collision_radius
                  << " outer_radius=" << outer_radius
                  << " min_transversality=" << min_transversality << '\n'
                  << "no output files or CAPD Poincare map were constructed\n";
        return 0;
    }

    // Zero-velocity boundary in the plotted section coordinates.
    // On y=0, H=h gives ydot^2 = 2*(h+Omega(x,0)) - xdot^2, so the
    // allowed region in the (x,xdot) section is bounded by
    // xdot = +/-sqrt(2*(h+Omega(x,0))).
    if (!zvb_output.empty()) {
        std::ofstream zvb(zvb_output);
        if (!zvb) {
            throw std::runtime_error("failed to open zero-velocity output: " + zvb_output);
        }
        zvb.precision(15);

        constexpr int boundary_count = 4000;
        std::vector<std::pair<double, double>> boundary;
        boundary.reserve(2 * boundary_count);
        for (int i = 0; i < boundary_count; ++i) {
            const double t = (boundary_count > 1)
                ? static_cast<double>(i) / static_cast<double>(boundary_count - 1)
                : 0.0;
            const double x = x_begin + t * (x_end - x_begin);
            const double speed2 = initial_speed_squared(
                x, 0.0, h, p, stable_e1_delta, signed_delta, e1_rho);
            if (speed2 >= 0.0 && std::isfinite(speed2)) {
                boundary.emplace_back(x, std::sqrt(speed2));
            }
        }

        for (const auto& q : boundary) {
            zvb << q.first << ',' << q.second << '\n';
        }
        zvb << "NaN,NaN\n";
        for (const auto& q : boundary) {
            zvb << q.first << ',' << -q.second << '\n';
        }
        std::cerr << "zero-velocity boundary -> " << zvb_output << '\n';
    }

    // --- CAPD Poincare map --------------------------------------------------
    DMap field(kHill4bpField, /*map_degree=*/1);
    field.setParameter("A", p.A);
    field.setParameter("B", p.B);
    field.setParameter("ac", p.alpha * p.c);

    DOdeSolver solver(field, solver_order);
    solver.setAbsoluteTolerance(abs_tol);
    solver.setRelativeTolerance(rel_tol);
    solver.setMaxStep(max_step);

    const double collision_R = std::sqrt(collision_radius);
    const double outer_R = std::sqrt(outer_radius);
    DNonlinearSection section(
        "par:Rc,Ro;var:R,th,v,w;fun:sin(th)*(R-Rc)*(Ro-R);");
    section.setParameter("Rc", collision_R);
    section.setParameter("Ro", outer_R);
    DPoincareMap pm(solver, section, capd::poincare::Both);
    pm.setMaxReturnTime(max_return_time);
    pm.setBlowUpMaxNorm(blow_up_norm);

    std::ofstream out(output);
    if (!out) {
        throw std::runtime_error("failed to open Poincare CSV output: " + output);
    }
    out.precision(std::numeric_limits<double>::max_digits10);
    std::ofstream returns;
    if (!returns_path.empty()) {
        returns.open(returns_path);
        if (!returns) {
            throw std::runtime_error(
                "failed to open detailed return CSV output: " + returns_path);
        }
        returns.precision(std::numeric_limits<double>::max_digits10);
        returns
            << "velocity_row,x_index,return_index,x0,initial_xdot,initial_ydot,"
            << "velocity_parameter,velocity_mode,h,hcrit,delta,side,"
            << "segment_return_tau,cumulative_tau,R,theta,v,w,x,y,xdot,ydot,"
            << "section_residual,energy_error\n";
    }
    std::ofstream outcomes;
    if (!outcomes_path.empty()) {
        outcomes.open(outcomes_path);
        if (!outcomes) {
            throw std::runtime_error("failed to open outcome CSV output: " + outcomes_path);
        }
        outcomes.precision(std::numeric_limits<double>::max_digits10);
        outcomes
            << "velocity_row,x_index,x0,velocity_parameter,velocity_mode,"
            << "requested_iterates,completed_iterates,status,last_return_tau,"
            << "max_return_tau,terminal_tau,elapsed_seconds,terminal_R,terminal_rho,"
            << "terminal_theta,terminal_v,terminal_w,message\n";
    }
    std::vector<std::pair<double, double>> pts;

    const double xStep = (x_count > 1) ? (x_end - x_begin) / (x_count - 1) : 0.0;
    const double xdotStep = (xdot_count > 1)
        ? (xdot_upper - xdot_lower) / (xdot_count - 1)
        : 0.0;
    long written = 0, skipped = 0, stopped = 0;
    std::map<std::string, long> outcome_counts;
    long energy_error_samples = 0;
    double max_abs_energy_error = 0.0;
    long return_time_samples = 0;
    double max_observed_return_tau = 0.0;
    double max_section_residual = 0.0;

    std::cerr << "lambda1=" << p.lambda1 << " lambda2=" << p.lambda2
              << " A=" << p.A << " B=" << p.B << " c=" << p.c
              << " c1=" << p.c1 << " c2=" << p.c2 << " c3=" << p.c3
              << " h=" << h << " (C=" << -2 * h << ")\n";
    std::cerr << "E1/E2 neck h=" << e1_h << " (C=" << -2 * e1_h << ")"
              << ", E3/E4 neck h=" << e3_h << " (C=" << -2 * e3_h << ")\n";
    std::cerr << "energy_mode=" << energy_mode << " neck=" << neck
              << " hcrit=" << hcrit
              << " delta=" << (delta_overridden ? delta : 0.0)
              << " side=" << (delta_overridden ? side : "not_used")
              << " actual_neck_side=" << actual_neck_side
              << " delta_h=" << h - hcrit << '\n';
    std::cerr << "x_count=" << x_count << " xdot_count=" << xdot_count
              << " iterates=" << iterates
              << " index range=[" << start_index << ',' << end_index << ']'
              << " velocity_mode=" << velocity_mode
              << " x_range=[" << x_begin << ',' << x_end << ']'
              << " max_return_time=" << max_return_time << '\n';
    std::cerr << "smooth coordinates R=sqrt(rho), CAPD crossing mode=Both"
              << " solver_order=" << solver_order
              << " max_step=" << max_step
              << " abs_tol=" << abs_tol << " rel_tol=" << rel_tol
              << " collision_radius=" << collision_radius
              << " outer_radius=" << outer_radius << '\n';

    const auto start_time = std::chrono::steady_clock::now();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto write_outcome = [&](int velocity_row,
                             int x_index,
                             double x0,
                             double velocity_parameter,
                             int completed_iterates,
                             const std::string& status,
                             double last_return_tau,
                             double maximum_return_tau,
                             double terminal_tau,
                             double elapsed_seconds,
                             const std::array<double, 4>& terminal,
        const std::string& message) {
        ++outcome_counts[status];
        if (!outcomes.is_open()) return;
        const double terminal_rho = std::isfinite(terminal[0])
            ? terminal[0] * terminal[0] : nan;
        outcomes << velocity_row << ',' << x_index << ',' << x0 << ','
                 << velocity_parameter << ',' << velocity_mode << ','
                 << iterates << ',' << completed_iterates << ','
                 << status << ',' << last_return_tau << ','
                 << maximum_return_tau << ',' << terminal_tau << ','
                 << elapsed_seconds << ','
                 << terminal[0] << ',' << terminal_rho << ','
                 << terminal[1] << ',' << terminal[2] << ',' << terminal[3]
                 << ',' << csv_escape(message) << '\n';
        outcomes.flush();
    };

    for (int m = 0; m < xdot_count; ++m) {
        const double velocity_parameter = xdot_lower + m * xdotStep;

        for (int n = start_index; n <= end_index; ++n) {
            const auto ic_start_time = std::chrono::steady_clock::now();
            const double x0 = x_begin + n * xStep;
            const bool show_ic_progress =
                ((n - start_index) % progress_every == 0) || (n == end_index);
            if (show_ic_progress) {
                const auto now = std::chrono::steady_clock::now();
                const double elapsed =
                    std::chrono::duration<double>(now - start_time).count();
                std::cerr << "velocity row " << m << '/' << (xdot_count - 1)
                          << " IC " << n << '/' << (x_count - 1)
                          << " x0=" << x0
                          << " velocity_parameter=" << velocity_parameter
                          << " elapsed=" << elapsed << "s"
                          << " written=" << written << '\n';
            }

            std::array<double, 4> terminal{{nan, nan, nan, nan}};
            const double speed2 = initial_speed_squared(
                x0, y_section, h, p, stable_e1_delta, signed_delta, e1_rho);
            if (!(speed2 >= 0.0) || !std::isfinite(speed2)) {
                ++skipped;
                const double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - ic_start_time).count();
                write_outcome(m, n, x0, velocity_parameter, 0, "forbidden", nan,
                              0.0, nan, elapsed, terminal,
                              "initial position is outside the energy-allowed region");
                if (show_ic_progress) {
                    std::cerr << "  skipped: outside allowed region\n";
                }
                continue;
            }

            const double speed = std::sqrt(speed2);
            const double initial_radius = std::hypot(x0, y_section);
            if (!(initial_radius > collision_radius
                  && initial_radius < outer_radius)) {
                ++skipped;
                const double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - ic_start_time).count();
                write_outcome(m, n, x0, velocity_parameter, 0,
                              "invalid_initial_radius", nan, 0.0, nan, elapsed,
                              terminal,
                              "initial radius must lie strictly between the "
                              "collision and outer limits");
                if (show_ic_progress) {
                    std::cerr << "  skipped: initial radius is outside study domain\n";
                }
                continue;
            }
            double initial_xdot = 0.0;
            double initial_ydot = 0.0;
            if (velocity_mode == "direction") {
                const double direction_norm = std::hypot(velocity_parameter, 1.0);
                initial_xdot = speed * velocity_parameter / direction_norm;
                initial_ydot = speed / direction_norm;
            } else if (velocity_mode == "fraction") {
                if (!(std::abs(velocity_parameter) < 1.0)) {
                    ++skipped;
                    const double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - ic_start_time).count();
                    write_outcome(m, n, x0, velocity_parameter, 0,
                                  "invalid_velocity", nan, 0.0, nan, elapsed, terminal,
                                  "velocity fraction must lie strictly between -1 and 1");
                    if (show_ic_progress) {
                        std::cerr << "  skipped: velocity fraction is not transverse\n";
                    }
                    continue;
                }
                initial_xdot = speed * velocity_parameter;
                initial_ydot = speed * std::sqrt(
                    std::max(0.0, 1.0 - velocity_parameter * velocity_parameter));
            } else {
                const double remaining = speed2 - velocity_parameter * velocity_parameter;
                if (!(remaining > 0.0) || !std::isfinite(remaining)) {
                    ++skipped;
                    const double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - ic_start_time).count();
                    write_outcome(m, n, x0, velocity_parameter, 0,
                                  "invalid_velocity", nan, 0.0, nan, elapsed, terminal,
                                  "requested xdot leaves no positive ydot at this energy");
                    if (show_ic_progress) {
                        std::cerr << "  skipped: xdot is outside the allowed velocity disk\n";
                    }
                    continue;
                }
                initial_xdot = velocity_parameter;
                initial_ydot = std::sqrt(remaining);
            }

            DVector u(4);
            int completed_iterates = 0;
            double last_return_tau = nan;
            double maximum_return_tau = 0.0;
            double terminal_tau = nan;
            double cumulative_tau = 0.0;
            double current_iterate_origin_tau = 0.0;
            std::string status = "completed";
            std::string stop_message;
            try {
                u = to_smooth_mcgehee(
                    x0, y_section, initial_xdot, initial_ydot);
                for (int j = 0; j < 4; ++j) terminal[j] = u[j];

                const CartesianState initial_cartesian = from_smooth_mcgehee(u);
                if (!finite_cartesian(initial_cartesian)
                    || !(initial_cartesian.ydot > min_transversality)) {
                    throw std::runtime_error(
                        "initial condition is not transverse with positive ydot");
                }
                const double initial_energy_error =
                    smooth_energy_error(u, p, h);
                if (std::isfinite(initial_energy_error)) {
                    max_abs_energy_error = std::max(
                        max_abs_energy_error, std::abs(initial_energy_error));
                    ++energy_error_samples;
                }

                for (int it = 0; it < iterates; ++it) {
                    current_iterate_origin_tau = cumulative_tau;
                    if (show_ic_progress && (it % iterate_progress_every == 0)) {
                        std::cerr << "  crossing " << it << '/' << iterates << '\n';
                    }
                    bool accepted = false;
                    double return_tau = 0.0;
                    for (int half = 0; half < max_half_crossings; ++half) {
                        u = pm(u, return_tau);
                        // in_out return_tau is measured from the beginning of
                        // this requested positive-return search and therefore
                        // already accumulates any intervening negative section
                        // crossing.  Offset it only by completed iterations.
                        cumulative_tau = current_iterate_origin_tau + return_tau;
                        terminal_tau = cumulative_tau;
                        for (int j = 0; j < 4; ++j) terminal[j] = u[j];
                        if (return_tau > max_return_time) {
                            status = "timeout";
                            stop_message =
                                "event occurred after the configured return-tau limit";
                            break;
                        }

                        const std::string surface =
                            classify_event_surface(u, collision_R, outer_R);
                        if (surface == "collision" || surface == "outer_limit") {
                            status = surface == "collision" && u[2] >= 0.0
                                ? "inner_radius_limit" : surface;
                            stop_message = "trajectory reached the configured "
                                + surface + " radius";
                            break;
                        }
                        if (surface != "section") {
                            throw std::runtime_error(
                                "CAPD returned a nonfinite event state");
                        }

                        const CartesianState crossing = from_smooth_mcgehee(u);
                        if (!finite_cartesian(crossing)) {
                            throw std::runtime_error(
                                "nonfinite Cartesian state at event");
                        }
                        if (std::abs(crossing.ydot) <= min_transversality) {
                            status = "nontransverse";
                            stop_message =
                                "section event failed the transversality threshold";
                            break;
                        }
                        if (crossing.ydot < 0.0) {
                            continue;
                        }

                        out << crossing.x << ',' << crossing.xdot << '\n';
                        max_section_residual =
                            std::max(max_section_residual, std::abs(std::sin(u[1])));
                        const double crossing_energy_error =
                            smooth_energy_error(u, p, h);
                        if (std::isfinite(crossing_energy_error)) {
                            max_abs_energy_error = std::max(
                                max_abs_energy_error,
                                std::abs(crossing_energy_error));
                            ++energy_error_samples;
                        }
                        if (written % svg_stride == 0) {
                            pts.emplace_back(crossing.x, crossing.xdot);
                        }
                        if (returns.is_open()) {
                            returns << m << ',' << n << ',' << completed_iterates << ','
                                    << x0 << ',' << initial_xdot << ',' << initial_ydot
                                    << ',' << velocity_parameter << ',' << velocity_mode
                                    << ',' << h << ',' << hcrit << ','
                                    << (delta_overridden ? delta : 0.0) << ','
                                    << (delta_overridden ? side : "not_used") << ','
                                    << return_tau << ',' << cumulative_tau << ','
                                    << u[0] << ',' << u[1] << ',' << u[2] << ',' << u[3]
                                    << ',' << crossing.x << ',' << crossing.y << ','
                                    << crossing.xdot << ',' << crossing.ydot << ','
                                    << std::sin(u[1]) << ',' << crossing_energy_error
                                    << '\n';
                            if (written % 1000 == 0) returns.flush();
                        }
                        ++written;
                        ++completed_iterates;
                        last_return_tau = return_tau;
                        maximum_return_tau =
                            std::max(maximum_return_tau, return_tau);
                        max_observed_return_tau =
                            std::max(max_observed_return_tau, return_tau);
                        ++return_time_samples;
                        if (written % 1000 == 0) out.flush();
                        accepted = true;
                        break;
                    }

                    if (status != "completed") break;
                    if (!accepted) {
                        status = "nontransverse";
                        stop_message =
                            "no positive transverse section crossing within "
                            + std::to_string(max_half_crossings)
                            + " event crossings";
                        break;
                    }
                }
            } catch (const capd::dynsys::SolverException<DVector>& e) {
                for (int j = 0; j < 4; ++j) terminal[j] = e.V[j];
                terminal_tau = current_iterate_origin_tau + e.currentTime;
                status = classify_solver_stop(
                    e.V, e.what(), collision_radius, outer_radius);
                stop_message = e.what();
            } catch (const std::exception& e) {
                status = "numerical_failure";
                stop_message = e.what();
            }

            if (status != "completed") {
                ++stopped;
                if (show_ic_progress || stopped <= 10) {
                    std::cerr << "  stopped velocity row " << m << " IC " << n
                              << " x0=" << x0
                              << " velocity_parameter=" << velocity_parameter
                              << " status=" << status
                              << ": " << stop_message << '\n';
                }
            }
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - ic_start_time).count();
            out.flush();
            write_outcome(m, n, x0, velocity_parameter, completed_iterates,
                          status, last_return_tau, maximum_return_tau,
                          terminal_tau, elapsed,
                          terminal, stop_message);
        }
    }

    // --- dependency-free SVG scatter ---------------------------------------
    if (!svg_path.empty() && !pts.empty()) {
        double xmin = std::numeric_limits<double>::max(), xmax = -xmin;
        double ymin = xmin, ymax = -xmin;
        for (const auto& q : pts) {
            xmin = std::min(xmin, q.first);  xmax = std::max(xmax, q.first);
            ymin = std::min(ymin, q.second); ymax = std::max(ymax, q.second);
        }
        const double W = 1000.0, H = 800.0, pad = 40.0;
        const double sx = (xmax > xmin) ? (W - 2 * pad) / (xmax - xmin) : 1.0;
        const double sy = (ymax > ymin) ? (H - 2 * pad) / (ymax - ymin) : 1.0;
        std::ofstream svg(svg_path);
        if (!svg) {
            throw std::runtime_error("failed to open SVG output: " + svg_path);
        }
        svg << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W
            << "' height='" << H << "'>\n<rect width='" << W << "' height='" << H
            << "' fill='white'/>\n";
        for (const auto& q : pts) {
            const double px = pad + (q.first - xmin) * sx;
            const double py = H - pad - (q.second - ymin) * sy;
            svg << "<circle cx='" << px << "' cy='" << py << "' r='0.6' fill='blue'/>\n";
        }
        svg << "<text x='" << W / 2 << "' y='" << H - 8
            << "' font-size='16' text-anchor='middle'>x_1</text>\n"
            << "<text x='14' y='" << H / 2 << "' font-size='16' transform='rotate(-90 14,"
            << H / 2 << ")' text-anchor='middle'>dx_1/dt</text>\n</svg>\n";
        std::cerr << "svg -> " << svg_path << '\n';
    }

    const double total_elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time).count();
    if (!metadata_path.empty()) {
        std::ofstream metadata(metadata_path);
        if (!metadata) {
            throw std::runtime_error("failed to open metadata output: " + metadata_path);
        }
        metadata.precision(17);
        const double h_offset = h - e1_h;
        metadata << "model=regularized_planar_hill4bp_capd\n";
        metadata << "regularized_coordinate=R=sqrt(rho)=paper_r^(1/5)\n";
        metadata << "time_coordinate=regularized_tau\n";
        metadata << "crossing_strategy=both_with_section_collision_outer_surfaces\n";
        metadata << "energy_mode=" << energy_mode << '\n';
        metadata << "h=" << h << '\n';
        metadata << "jacobi_C=" << -2.0 * h << '\n';
        metadata << "neck=" << neck << '\n';
        metadata << "hcrit=" << hcrit << '\n';
        metadata << "jacobi_Ccrit=" << -2.0 * hcrit << '\n';
        metadata << "delta_specified=" << (delta_overridden ? 1 : 0) << '\n';
        metadata << "delta=" << (delta_overridden ? delta : 0.0) << '\n';
        metadata << "side=" << (delta_overridden ? side : "not_used") << '\n';
        metadata << "actual_neck_side=" << actual_neck_side << '\n';
        metadata << "E1_E2_h=" << e1_h << '\n';
        metadata << "E1_E2_jacobi_C=" << -2.0 * e1_h << '\n';
        metadata << "E1_E2_abs_x=" << e1_rho << '\n';
        metadata << "E3_E4_h=" << e3_h << '\n';
        metadata << "E3_E4_jacobi_C=" << -2.0 * e3_h << '\n';
        metadata << "E3_E4_abs_y=" << e3_rho << '\n';
        metadata << "h_e1_offset=" << h_offset << '\n';
        metadata << "h_e1_scale=" << h / e1_h << '\n';
        metadata << "neck_side=" << classify_neck_side(h, e1_h) << '\n';
        metadata << "x_begin=" << x_begin << '\n';
        metadata << "x_end=" << x_end << '\n';
        metadata << "x_count=" << x_count << '\n';
        metadata << "neck_window_sigma=" << neck_window_sigma << '\n';
        metadata << "neck_branch=" << neck_branch << '\n';
        metadata << "neck_length_scale=";
        if (std::isfinite(neck_length_scale)) {
            metadata << neck_length_scale;
        } else {
            metadata << "nan";
        }
        metadata << '\n';
        metadata << "xdot_lower=" << xdot_lower << '\n';
        metadata << "xdot_upper=" << xdot_upper << '\n';
        metadata << "xdot_count=" << xdot_count << '\n';
        metadata << "velocity_mode=" << velocity_mode << '\n';
        metadata << "iterates=" << iterates << '\n';
        metadata << "y_section=" << y_section << '\n';
        metadata << "solver_order=" << solver_order << '\n';
        metadata << "abs_tol=" << abs_tol << '\n';
        metadata << "rel_tol=" << rel_tol << '\n';
        metadata << "max_step=" << max_step << '\n';
        metadata << "max_return_time=" << max_return_time << '\n';
        metadata << "max_half_crossings=" << max_half_crossings << '\n';
        metadata << "min_transversality=" << min_transversality << '\n';
        metadata << "blow_up_norm=" << blow_up_norm << '\n';
        metadata << "collision_radius=" << collision_radius << '\n';
        metadata << "outer_radius=" << outer_radius << '\n';
        metadata << "outcomes_file=" << outcomes_path << '\n';
        metadata << "detailed_returns_file=" << returns_path << '\n';
        metadata << "svg_stride=" << svg_stride << '\n';
        metadata << "svg_preview_points=" << pts.size() << '\n';
        metadata << "elapsed_seconds=" << total_elapsed_seconds << '\n';
        metadata << "point_count=" << written << '\n';
        metadata << "written_points=" << written << '\n';
        metadata << "sampled_initial_conditions="
                 << static_cast<long long>(end_index - start_index + 1)
                    * static_cast<long long>(xdot_count) << '\n';
        metadata << "skipped_initial_conditions=" << skipped << '\n';
        metadata << "stopped_initial_conditions=" << stopped << '\n';
        metadata << "completed_initial_conditions=" << outcome_counts["completed"] << '\n';
        metadata << "forbidden_initial_conditions=" << outcome_counts["forbidden"] << '\n';
        metadata << "invalid_velocity_initial_conditions="
                 << outcome_counts["invalid_velocity"] << '\n';
        metadata << "invalid_initial_radius_initial_conditions="
                 << outcome_counts["invalid_initial_radius"] << '\n';
        metadata << "collision_initial_conditions=" << outcome_counts["collision"] << '\n';
        metadata << "inner_radius_limit_initial_conditions="
                 << outcome_counts["inner_radius_limit"] << '\n';
        metadata << "outer_limit_initial_conditions="
                 << outcome_counts["outer_limit"] << '\n';
        metadata << "timeout_initial_conditions=" << outcome_counts["timeout"] << '\n';
        metadata << "nontransverse_initial_conditions="
                 << outcome_counts["nontransverse"] << '\n';
        metadata << "numerical_failure_initial_conditions="
                 << outcome_counts["numerical_failure"] << '\n';
        metadata << "return_time_samples=" << return_time_samples << '\n';
        metadata << "max_observed_return_tau=";
        if (return_time_samples > 0) {
            metadata << max_observed_return_tau;
        } else {
            metadata << "nan";
        }
        metadata << '\n';
        metadata << "max_section_residual=";
        if (return_time_samples > 0) {
            metadata << max_section_residual;
        } else {
            metadata << "nan";
        }
        metadata << '\n';
        metadata << "energy_error_samples=" << energy_error_samples << '\n';
        metadata << "max_abs_energy_error=";
        if (energy_error_samples > 0) {
            metadata << max_abs_energy_error;
        } else {
            metadata << "nan";
        }
        metadata << '\n';
        metadata << "max_abs_energy_error_over_delta=";
        if (energy_error_samples > 0 && delta_overridden) {
            metadata << max_abs_energy_error / delta;
        } else {
            metadata << "nan";
        }
        metadata << '\n';
        metadata << "lambda1=" << p.lambda1 << '\n';
        metadata << "lambda2=" << p.lambda2 << '\n';
        metadata << "A=" << p.A << '\n';
        metadata << "B=" << p.B << '\n';
        metadata << "c1=" << p.c1 << '\n';
        metadata << "c2=" << p.c2 << '\n';
        metadata << "c3=" << p.c3 << '\n';
        metadata << "c=" << p.c << '\n';
        std::cerr << "metadata -> " << metadata_path << '\n';
    }

    std::cerr << "wrote " << written << " points (" << skipped << " ICs skipped, "
              << stopped << " orbits stopped early, elapsed="
              << total_elapsed_seconds << "s) -> " << output << '\n';
    std::cerr << "outcomes:";
    for (const auto& item : outcome_counts) {
        std::cerr << ' ' << item.first << '=' << item.second;
    }
    std::cerr << '\n';
    if (energy_error_samples > 0) {
        std::cerr << "maximum sampled |H-h|=" << max_abs_energy_error << '\n';
        if (delta_overridden && max_abs_energy_error >= 0.01 * delta) {
            std::cerr << "warning: sampled energy error is at least 1% of delta; "
                      << "treat this level as unresolved until a tighter CAPD run agrees\n";
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
