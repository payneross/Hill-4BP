// Rigorous interval validation for near-E1 Poincare returns in the
// McGehee-regularized planar Hill four-body problem.
//
// This executable deliberately complements, rather than replaces, the dense
// floating-point survey in hill4bp_poincare_capd.cpp.  Every decimal model
// input and every mathematical command-line value is enclosed with
// directed-rounding CAPD intervals (step/tolerance controls are binary64
// algorithm settings).  Initial parameter cells are mapped to boxes containing
// the corresponding fixed-energy graph, and CAPD's rigorous C0 set integrator
// encloses every trajectory in each box.
//
// CAPD 6.1's rigorous PoincareMap path does not honour setMaxReturnTime, and
// ITimeMap overwrites a configured maximum step with its final-time gap.  This
// validator therefore initializes adaptive step control explicitly and moves
// a C0 set by one IOdeSolver step at a time, reasserting and auditing the step
// cap on every move.  On every rigorous dense solution curve it
// examines sin(theta), R-Rc, and Ro-R separately.  Interval Newton plus a
// derivative sign proves a unique root, and disjoint root-time intervals prove
// which event is first.  This catches even numbers of crossings inside one
// solver step.  Ambiguity, wrapping, or a solver exception is an unresolved
// proof obligation, never a dynamical classification.

#include "capd/capdlib.h"
#include "capd_decimal_literal.hpp"

#ifdef __FAST_MATH__
#error "Rigorous interval validation must not be compiled with -ffast-math"
#endif

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using capd::C0HOTripletonSet;
using capd::IMap;
using capd::IOdeSolver;
using capd::IVector;
using capd::interval;

namespace {

static_assert(std::numeric_limits<double>::is_iec559,
              "CAPD NATIVE interval proof requires IEC 60559 binary64");

const char* kValidatedField =
    "par:A,B,ac;"
    "var:R,th,v,w;"
    "fun:"
    "0.5*v*R,"
    "w-R^5,"
    "1.5*v^2+w^2-ac-R^4-2*A*R^10*cos(th)^2-2*B*R^10*sin(th)^2,"
    "0.5*v*w+2*(A-B)*R^10*sin(th)*cos(th);";

struct RawConfig {
    std::map<std::string, std::string> values;
    std::string path;
};

struct ModelIntervals {
    interval lambda1;
    interval lambda2;
    interval A;
    interval B;
    interval c;
    interval ac;
    interval c1;
    interval c2;
    interval c3;
    interval alpha;
    interval nu;
};

struct EquilibriumCertificate {
    interval rho;
    interval h;
    interval derivative;
    bool endpoint_signs = false;
    bool derivative_positive = false;
    int bisections = 0;
};

struct FirstNeckCertificate {
    EquilibriumCertificate e1;
    EquilibriumCertificate e3;
    interval omega_xx;
    interval omega_yy;
    interval curvature;
    bool e1_index_one = false;
    bool e1_before_e3 = false;
    bool positive_lambda_c = false;
    bool global_axis_root_uniqueness = false;
    bool no_mixed_axis_equilibria = false;
    bool passed = false;
};

struct Settings {
    std::string side = "closed";
    std::string delta_text;
    interval delta;
    std::string branch = "positive";
    int branch_sign = 1;
    bool dry_run = false;
    bool require_complete = true;
    bool require_returns = false;
    std::string grid_mode = "points";

    bool explicit_window = false;
    interval x_begin;
    interval x_end;
    interval neck_window_sigma = interval("4", "4");
    int x_count = 41;
    int start_index = 0;
    int end_index = -1;
    interval x_radius = interval(0.0);

    interval fraction_lower = interval("-0.9", "-0.9");
    interval fraction_upper = interval("0.9", "0.9");
    int fraction_count = 21;
    interval fraction_radius = interval(0.0);

    int iterates = 1;
    int solver_order = 30;
    double abs_tol = 1e-14;
    double rel_tol = 1e-14;
    double max_step = 0.005;
    int max_step_retries = 24;
    interval tau_cap = interval("100", "100");
    int max_half_crossings = 4;
    int event_time_subdivision_depth = 30;
    long long max_event_scan_nodes = 65536;
    int max_subdivision_depth = 10;
    long long max_leaf_boxes = 200000;
    interval inner_radius = interval("1e-6", "1e-6");
    interval outer_radius = interval("5", "5");
    interval min_transversality = interval("1e-12", "1e-12");
    int progress_every = 25;

    std::string output = "hill4bp_interval_enclosures.csv";
    std::string metadata = "hill4bp_interval_proof_metadata.txt";
};

struct ParameterCell {
    interval x;
    interval fraction;
    int x_index = 0;
    int fraction_index = 0;
    int depth = 0;
    std::string path = "r";
};

struct EventRecord {
    std::string type;
    int event_index = -1;
    int positive_return_index = -1;
    interval tau = interval(0.0);
    interval clock = interval(0.0);
    interval x = interval(0.0);
    interval xdot = interval(0.0);
    interval ydot = interval(0.0);
    interval transversality = interval(0.0);
    interval R = interval(0.0);
    interval energy_residual = interval(0.0);
};

struct AttemptResult {
    bool resolved = false;
    bool admissible = false;
    bool forbidden = false;
    int completed_returns = 0;
    std::string status;
    std::string message;
    interval speed_squared = interval(0.0);
    interval energy_graph_residual = interval(0.0);
    long long event_scan_nodes = 0;
    long long solver_step_retries = 0;
    std::vector<EventRecord> events;
};

struct LeafResult {
    ParameterCell cell;
    AttemptResult attempt;
};

struct Bounds {
    double lo = std::numeric_limits<double>::quiet_NaN();
    double hi = std::numeric_limits<double>::quiet_NaN();
};

std::string trim(const std::string& value) {
    std::size_t begin = 0;
    while (begin < value.size()
           && std::isspace(static_cast<unsigned char>(value[begin]))) {
        ++begin;
    }
    std::size_t end = value.size();
    while (end > begin
           && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }
    return value.substr(begin, end - begin);
}

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) {
                       return static_cast<char>(std::tolower(ch));
                   });
    return value;
}

std::string csv_escape(const std::string& value) {
    bool quote = false;
    std::string result;
    result.reserve(value.size() + 8);
    for (char ch : value) {
        if (ch == '"') {
            result += "\"\"";
            quote = true;
        } else if (ch == '\n' || ch == '\r') {
            result += ' ';
            quote = true;
        } else {
            result += ch;
            quote = quote || ch == ',';
        }
    }
    return quote ? '"' + result + '"' : result;
}

void validate_decimal(const std::string& text, const std::string& key) {
    (void)hill4bp_capd_input::canonical_decimal(text, key);
}

interval decimal_interval(const std::string& text, const std::string& key) {
    // CAPD NATIVE deliberately applies predecessor/successor to string
    // endpoints, even when the mathematical literal is exactly zero.  That
    // turns "0" into [-min_subnormal,+min_subnormal], which is a valid outer
    // enclosure but not a valid nonnegative radius.  Zero is exactly
    // representable, so preserve it as the rigorous singleton [0,0].
    if (hill4bp_capd_input::exact_decimal_zero(text, key)) {
        return interval(0.0);
    }
    return interval(text, text);
}

interval required_exact_constant(const std::string& text,
                                 const std::string& key,
                                 const std::string& expected_text,
                                 double expected_value) {
    if (!hill4bp_capd_input::exact_decimal_equal(
            text, expected_text, key)) {
        throw std::runtime_error(
            "rigorous field requires configured " + key + " exactly equal to "
            + expected_text);
    }
    // The only callers use 0, 1, and 3, all exactly representable in binary64.
    return interval(expected_value);
}

double decimal_double(const std::string& text, const std::string& key) {
    validate_decimal(text, key);
    return std::stod(text);
}

int parse_int(const std::string& text, const std::string& key) {
    std::size_t consumed = 0;
    long value = 0;
    try {
        value = std::stol(text, &consumed, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer for " + key + ": " + text);
    }
    if (consumed != text.size()
        || value < std::numeric_limits<int>::min()
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("invalid integer for " + key + ": " + text);
    }
    return static_cast<int>(value);
}

long long parse_long_long(const std::string& text, const std::string& key) {
    std::size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &consumed, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer for " + key + ": " + text);
    }
    if (consumed != text.size()) {
        throw std::runtime_error("invalid integer for " + key + ": " + text);
    }
    return value;
}

bool parse_bool(const std::string& text, const std::string& key) {
    const std::string value = lowercase(text);
    if (value == "1" || value == "true") return true;
    if (value == "0" || value == "false") return false;
    throw std::runtime_error(key + " must be 0, 1, false, or true");
}

bool contains_zero(const interval& value) {
    return value.leftBound() <= 0.0 && value.rightBound() >= 0.0;
}

bool strictly_positive(const interval& value) {
    return value.leftBound() > 0.0;
}

bool strictly_negative(const interval& value) {
    return value.rightBound() < 0.0;
}

double width(const interval& value) {
    return value.rightBound() - value.leftBound();
}

interval hull(const interval& a, const interval& b) {
    return interval(std::min(a.leftBound(), b.leftBound()),
                    std::max(a.rightBound(), b.rightBound()));
}

interval point_from_bound(double value) {
    // A binary64 value is exactly representable as a degenerate native CAPD
    // interval.  This is used only after the enclosing arithmetic produced it.
    return interval(value);
}

Bounds bounds(const interval& value) {
    return Bounds{value.leftBound(), value.rightBound()};
}

std::string find_parameter_file() {
    const std::vector<std::string> candidates{
        "Hill_4BP_parameters.cfg",
        "../Hill_4BP_parameters.cfg",
        "../../Hill_4BP_parameters.cfg",
        "../../../Hill_4BP_parameters.cfg",
    };
    for (const std::string& candidate : candidates) {
        std::ifstream input(candidate);
        if (input) return candidate;
    }
    throw std::runtime_error("could not find Hill_4BP_parameters.cfg");
}

RawConfig read_raw_config() {
    RawConfig result;
    result.path = find_parameter_file();
    std::ifstream input(result.path);
    if (!input) throw std::runtime_error("failed to open " + result.path);

    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        line = trim(line);
        if (line.empty()) continue;
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos || equals == 0) {
            throw std::runtime_error(
                "invalid parameter line " + std::to_string(line_number));
        }
        const std::string key = trim(line.substr(0, equals));
        const std::string value = trim(line.substr(equals + 1));
        validate_decimal(value, key);
        if (!result.values.emplace(key, value).second) {
            throw std::runtime_error("duplicate parameter key: " + key);
        }
    }
    return result;
}

const std::string& required_raw(const RawConfig& config,
                                const std::string& key) {
    const auto found = config.values.find(key);
    if (found == config.values.end()) {
        throw std::runtime_error("missing parameter key: " + key);
    }
    return found->second;
}

ModelIntervals load_model(const RawConfig& config) {
    ModelIntervals p;
    p.lambda1 = decimal_interval(required_raw(config, "lambda1"), "lambda1");
    p.lambda2 = decimal_interval(required_raw(config, "lambda2"), "lambda2");
    p.c1 = required_exact_constant(
        required_raw(config, "c1"), "c1", "0", 0.0);
    p.c2 = required_exact_constant(
        required_raw(config, "c2"), "c2", "0", 0.0);
    p.c3 = decimal_interval(required_raw(config, "c3"), "c3");
    p.alpha = required_exact_constant(
        required_raw(config, "alpha"), "alpha", "3", 3.0);
    p.nu = required_exact_constant(
        required_raw(config, "nu"), "nu", "1", 1.0);
    if (!strictly_negative(p.c3)) {
        throw std::runtime_error("require c3<0 so c=-c3>0");
    }
    p.c = -p.c3;
    if (!strictly_positive(p.lambda1) || !strictly_positive(p.lambda2)) {
        throw std::runtime_error(
            "global first-neck certificate requires lambda1>0 and lambda2>0");
    }
    p.A = (interval(1.0) - p.lambda2) / interval(2.0);
    p.B = (interval(1.0) - p.lambda1) / interval(2.0);
    p.ac = p.alpha * p.c;
    return p;
}

interval pow3(const interval& x) { return x * x * x; }
interval pow4(const interval& x) { const interval x2 = x * x; return x2 * x2; }
interval pow5(const interval& x) { return pow4(x) * x; }
interval pow6(const interval& x) { return pow3(x) * pow3(x); }
interval pow10(const interval& x) { const interval x5 = pow5(x); return x5 * x5; }

// R^6 times the physical Hamiltonian residual H-h in the regularized
// coordinates.  Evaluating this polynomial form avoids divisions near the
// inner study limit and is therefore the useful quantity to audit in the
// interval event table.
interval scaled_energy_residual(const IVector& state,
                                const interval& h,
                                const ModelIntervals& p) {
    const interval& R = state[0];
    const interval& th = state[1];
    const interval& v = state[2];
    const interval& w = state[3];
    return interval(0.5) * (v * v + w * w) - p.c
        - w * pow5(R) - pow4(R)
        + p.A * pow10(R) * sqr(cos(th))
        + p.B * pow10(R) * sqr(sin(th))
        - h * pow6(R);
}

interval equilibrium_function(const interval& rho,
                              const interval& lambda,
                              const interval& c) {
    return lambda * rho - interval(1.0) / (rho * rho)
        - interval(3.0) * c / pow4(rho);
}

interval equilibrium_derivative(const interval& rho,
                                const interval& lambda,
                                const interval& c) {
    return lambda + interval(2.0) / pow3(rho)
        + interval(12.0) * c / pow5(rho);
}

interval equilibrium_energy(const interval& rho,
                            const interval& lambda,
                            const interval& c) {
    return -(interval(0.5) * lambda * rho * rho
             + interval(1.0) / rho + c / pow3(rho));
}

EquilibriumCertificate certify_equilibrium(const interval& lambda,
                                            const interval& c,
                                            double initial_lo,
                                            double initial_hi,
                                            const std::string& name) {
    const interval initial(initial_lo, initial_hi);
    const interval derivative = equilibrium_derivative(initial, lambda, c);
    const interval f_lo = equilibrium_function(
        point_from_bound(initial_lo), lambda, c);
    const interval f_hi = equilibrium_function(
        point_from_bound(initial_hi), lambda, c);
    if (!(strictly_negative(f_lo) && strictly_positive(f_hi)
          && strictly_positive(derivative))) {
        throw std::runtime_error(
            "failed to prove a unique " + name + " root in its initial bracket");
    }

    double lo = initial_lo;
    double hi = initial_hi;
    int iterations = 0;
    for (; iterations < 256; ++iterations) {
        const double mid = lo + (hi - lo) * 0.5;
        if (!(mid > lo && mid < hi)) break;
        const interval f_mid = equilibrium_function(
            point_from_bound(mid), lambda, c);
        if (strictly_negative(f_mid)) {
            lo = mid;
        } else if (strictly_positive(f_mid)) {
            hi = mid;
        } else {
            // Parameter/literal rounding has reached its irreducible width.
            break;
        }
    }

    EquilibriumCertificate certificate;
    certificate.rho = interval(lo, hi);
    certificate.h = equilibrium_energy(certificate.rho, lambda, c);
    certificate.derivative = equilibrium_derivative(
        certificate.rho, lambda, c);
    certificate.endpoint_signs = true;
    certificate.derivative_positive = strictly_positive(certificate.derivative);
    certificate.bisections = iterations;
    if (!certificate.derivative_positive) {
        throw std::runtime_error(name + " derivative proof was lost after isolation");
    }
    return certificate;
}

FirstNeckCertificate certify_first_neck(const ModelIntervals& p) {
    FirstNeckCertificate result;
    result.e1 = certify_equilibrium(p.lambda2, p.c, 0.5, 1.5, "E1/E2");
    result.e3 = certify_equilibrium(p.lambda1, p.c, 1.0, 8.0, "E3/E4");
    result.omega_xx = equilibrium_derivative(
        result.e1.rho, p.lambda2, p.c);
    result.curvature = result.omega_xx;
    result.omega_yy = p.lambda1 - interval(1.0) / pow3(result.e1.rho)
        - interval(3.0) * p.c / pow5(result.e1.rho);
    result.e1_index_one = strictly_positive(result.omega_xx)
        && strictly_negative(result.omega_yy);
    result.e1_before_e3 = result.e1.h.rightBound()
        < result.e3.h.leftBound();
    result.positive_lambda_c = strictly_positive(p.lambda1)
        && strictly_positive(p.lambda2) && strictly_positive(p.c);
    // For lambda,c>0, g_lambda(rho)=lambda*rho-rho^-2-3c*rho^-4
    // increases strictly from -infinity to +infinity on rho>0.  The bracketed
    // root is therefore the unique positive root on its complete axis.
    result.global_axis_root_uniqueness = result.positive_lambda_c
        && result.e1.endpoint_signs && result.e3.endpoint_signs;
    // A mixed-axis critical point would require the same radial factor to
    // equal both lambda1 and lambda2.  Strict separation rules it out.
    result.no_mixed_axis_equilibria =
        p.lambda1.rightBound() < p.lambda2.leftBound()
        || p.lambda2.rightBound() < p.lambda1.leftBound();
    result.passed = result.e1.endpoint_signs
        && result.e1.derivative_positive
        && result.e3.endpoint_signs
        && result.e3.derivative_positive
        && result.e1_index_one
        && result.e1_before_e3
        && result.global_axis_root_uniqueness
        && result.no_mixed_axis_equilibria;
    if (!result.passed) {
        throw std::runtime_error(
            "could not certify E1/E2 as the first index-one neck pair");
    }
    return result;
}

interval e1_potential_difference(const interval& x,
                                 const interval& xe,
                                 const ModelIntervals& p,
                                 int branch_sign) {
    const interval ax = branch_sign > 0 ? x : -x;
    if (!(ax.leftBound() > 0.0 && xe.leftBound() > 0.0)) {
        throw std::runtime_error("E1 potential difference requires one x branch");
    }
    const interval d = ax - xe;
    const interval ax2 = ax * ax;
    const interval ax3 = ax2 * ax;
    const interval xe2 = xe * xe;
    const interval xe3 = xe2 * xe;
    const interval xe4 = xe3 * xe;
    const interval numerator =
        p.lambda2 * xe3 * ax3
        + interval(2.0) * p.lambda2 * xe4 * ax2
        + interval(4.0) * p.c * ax
        + interval(2.0) * p.c * xe;
    return sqr(d) * numerator / (interval(2.0) * ax3 * xe3);
}

interval stable_speed_squared(const interval& x,
                              const Settings& settings,
                              const ModelIntervals& p,
                              const FirstNeckCertificate& neck) {
    const interval signed_delta = settings.side == "open"
        ? settings.delta : -settings.delta;
    return interval(2.0) * (
        signed_delta
        + e1_potential_difference(
            x, neck.e1.rho, p, settings.branch_sign));
}

interval affine_grid_value(const interval& begin,
                           const interval& end,
                           int numerator,
                           int denominator) {
    if (denominator == 0) return (begin + end) / interval(2.0);
    const interval t = interval(static_cast<double>(numerator))
        / interval(static_cast<double>(denominator));
    return (interval(1.0) - t) * begin + t * end;
}

interval grid_cell(const interval& begin,
                   const interval& end,
                   int index,
                   int count,
                   const std::string& mode,
                   const interval& radius) {
    if (mode == "cover") {
        const interval left = affine_grid_value(begin, end, index, count);
        const interval right = affine_grid_value(begin, end, index + 1, count);
        return hull(left, right);
    }
    const interval center = count == 1
        ? (begin + end) / interval(2.0)
        : affine_grid_value(begin, end, index, count - 1);
    return center + interval(-1.0, 1.0) * radius;
}

std::pair<interval, interval> split_interval(const interval& value) {
    const double midpoint = value.leftBound()
        + (value.rightBound() - value.leftBound()) * 0.5;
    if (!(midpoint > value.leftBound() && midpoint < value.rightBound())) {
        return std::make_pair(value, value);
    }
    return std::make_pair(interval(value.leftBound(), midpoint),
                          interval(midpoint, value.rightBound()));
}

IVector initial_box(const ParameterCell& cell,
                    const Settings& settings,
                    const ModelIntervals& p,
                    const FirstNeckCertificate& neck,
                    interval& speed_squared,
                    interval& graph_residual) {
    speed_squared = stable_speed_squared(cell.x, settings, p, neck);
    if (!strictly_positive(speed_squared)) {
        return IVector(4);
    }
    if (!(cell.fraction.leftBound() > -1.0
          && cell.fraction.rightBound() < 1.0)) {
        throw std::runtime_error("fraction cell must be strictly inside (-1,1)");
    }

    const interval speed = sqrt(speed_squared);
    const interval xdot = speed * cell.fraction;
    const interval one_minus_fraction_squared =
        interval(1.0) - sqr(cell.fraction);
    if (!strictly_positive(one_minus_fraction_squared)) {
        throw std::runtime_error("fraction cell is not uniformly transverse");
    }
    const interval ydot = speed * sqrt(one_minus_fraction_squared);
    const interval rho = settings.branch_sign > 0 ? cell.x : -cell.x;
    if (!strictly_positive(rho)) {
        throw std::runtime_error("x cell crosses the origin or wrong neck branch");
    }
    const interval R = sqrt(rho);
    const interval R3 = pow3(R);

    IVector state(4);
    state[0] = R;
    state[1] = interval(0.0);  // angle relative to the selected E1 branch
    state[2] = interval(static_cast<double>(settings.branch_sign)) * R3 * xdot;
    state[3] = interval(static_cast<double>(settings.branch_sign))
        * R3 * (ydot + cell.x);
    // This residual is evaluated in the stable, cancellation-free energy
    // graph coordinates.  It must enclose zero; width is dependency wrapping,
    // not a floating-point energy drift.
    const interval signed_delta = settings.side == "open"
        ? settings.delta : -settings.delta;
    const interval potential_difference = e1_potential_difference(
        cell.x, neck.e1.rho, p, settings.branch_sign);
    graph_residual = interval(0.5) * (xdot * xdot + ydot * ydot)
        - (signed_delta + potential_difference);
    if (!contains_zero(graph_residual)) {
        throw std::runtime_error(
            "internal error: interval energy graph residual excludes zero");
    }
    return state;
}

void cartesian_enclosure(const IVector& state,
                         int branch_sign,
                         interval& x,
                         interval& xdot,
                         interval& ydot) {
    const interval& R = state[0];
    const interval& th = state[1];
    const interval& v = state[2];
    const interval& w = state[3];
    if (!strictly_positive(R)) {
        throw std::runtime_error("validated event enclosure does not prove R>0");
    }
    const interval sigma(static_cast<double>(branch_sign));
    const interval cth = cos(th);
    const interval sth = sin(th);
    const interval R2 = R * R;
    const interval R3 = R2 * R;
    x = sigma * R2 * cth;
    const interval y = sigma * R2 * sth;
    const interval p1 = sigma * (v * cth - w * sth) / R3;
    const interval p2 = sigma * (v * sth + w * cth) / R3;
    xdot = p1 + y;
    ydot = p2 - x;
}

enum class SurfaceKind { section, inner, outer };

struct RootCandidate {
    SurfaceKind surface = SurfaceKind::section;
    interval relative_time = interval(0.0);
};

struct SurfaceScan {
    std::vector<RootCandidate> candidates;
    std::vector<interval> unresolved_times;
    bool node_budget_exhausted = false;
};

struct StepScan {
    bool has_event = false;
    bool unresolved = false;
    bool node_budget_exhausted = false;
    std::string message;
    RootCandidate event;
    long long nodes_visited = 0;
};

const char* surface_name(SurfaceKind surface) {
    switch (surface) {
        case SurfaceKind::section: return "section";
        case SurfaceKind::inner: return "inner_radius";
        case SurfaceKind::outer: return "outer_radius";
    }
    return "unknown";
}

interval surface_value(SurfaceKind surface,
                       const IVector& state,
                       const interval& Rc,
                       const interval& Ro) {
    switch (surface) {
        case SurfaceKind::section:
            return sin(state[1]);
        case SurfaceKind::inner:
            return state[0] - Rc;
        case SurfaceKind::outer:
            return Ro - state[0];
    }
    throw std::logic_error("unknown event surface");
}

interval surface_derivative(SurfaceKind surface,
                            const IVector& state) {
    switch (surface) {
        case SurfaceKind::section:
            return cos(state[1]) * (state[3] - pow5(state[0]));
        case SurfaceKind::inner:
            return interval(0.5) * state[2] * state[0];
        case SurfaceKind::outer:
            return -interval(0.5) * state[2] * state[0];
    }
    throw std::logic_error("unknown event surface");
}

void scan_surface_interval(const IOdeSolver::SolutionCurve& curve,
                           SurfaceKind surface,
                           const interval& domain,
                           const interval& Rc,
                           const interval& Ro,
                           bool ignore_initial_section_root,
                           int depth,
                           int maximum_depth,
                           long long& nodes_remaining,
                           SurfaceScan& scan) {
    if (nodes_remaining <= 0) {
        scan.unresolved_times.push_back(domain);
        scan.node_budget_exhausted = true;
        return;
    }
    --nodes_remaining;
    const IVector state = curve(domain);
    const interval value = surface_value(surface, state, Rc, Ro);
    if (!contains_zero(value)) return;

    const interval derivative = surface_derivative(surface, state);

    // Every restarted section chart has theta=0 exactly.  If sin(theta) is
    // strictly monotone on a prefix beginning at t=0 and has the expected
    // strict endpoint sign, the only zero in that prefix is the deliberately
    // ignored initial one.
    if (surface == SurfaceKind::section
        && ignore_initial_section_root
        && domain.leftBound() == curve.getLeftDomain()
        && !contains_zero(derivative)) {
        const interval left_value = surface_value(
            surface, curve(domain.left()), Rc, Ro);
        const interval right_value = surface_value(
            surface, curve(domain.right()), Rc, Ro);
        const bool departed = contains_zero(left_value)
            && ((strictly_positive(derivative) && strictly_positive(right_value))
                || (strictly_negative(derivative) && strictly_negative(right_value)));
        if (departed) return;
    }

    if (!contains_zero(derivative)) {
        interval midpoint = domain.mid();
        const interval value_at_midpoint = surface_value(
            surface, curve(midpoint), Rc, Ro);
        interval newton = midpoint - value_at_midpoint / derivative;
        interval contracted;
        if (!intersection(domain, newton, contracted)) {
            // Interval Newton exclusion: no member of the parameterized flow
            // can have a zero in this time cell.
            return;
        }

        const interval left_value = surface_value(
            surface, curve(domain.left()), Rc, Ro);
        const interval right_value = surface_value(
            surface, curve(domain.right()), Rc, Ro);
        const bool endpoint_existence =
            (strictly_positive(derivative)
             && strictly_negative(left_value)
             && strictly_positive(right_value))
            || (strictly_negative(derivative)
                && strictly_positive(left_value)
                && strictly_negative(right_value));
        const bool parameterized_newton_existence =
            newton.subsetInterior(domain);
        // A uniform strict endpoint sign change or a strict parameterized
        // interval-Newton inclusion proves that every trajectory in the set
        // has one root here.  The derivative interval excluding zero gives
        // uniqueness across the complete time cell.
        if (endpoint_existence || parameterized_newton_existence) {
            interval root_time = contracted;
            for (int iteration = 0; iteration < 12; ++iteration) {
                const IVector root_state = curve(root_time);
                const interval root_derivative = surface_derivative(
                    surface, root_state);
                if (contains_zero(root_derivative)) break;
                const interval root_midpoint = root_time.mid();
                const interval root_value = surface_value(
                    surface, curve(root_midpoint), Rc, Ro);
                const interval next_newton =
                    root_midpoint - root_value / root_derivative;
                interval next_root;
                if (!intersection(root_time, next_newton, next_root)
                    || next_root == root_time) {
                    break;
                }
                root_time = next_root;
            }
            scan.candidates.push_back(RootCandidate{surface, root_time});
            return;
        }
    }

    if (depth >= maximum_depth) {
        scan.unresolved_times.push_back(domain);
        return;
    }
    const std::pair<interval, interval> pieces = split_interval(domain);
    if (pieces.first == pieces.second) {
        scan.unresolved_times.push_back(domain);
        return;
    }
    scan_surface_interval(curve, surface, pieces.first, Rc, Ro,
                          ignore_initial_section_root,
                          depth + 1, maximum_depth, nodes_remaining, scan);
    scan_surface_interval(curve, surface, pieces.second, Rc, Ro,
                          false, depth + 1, maximum_depth,
                          nodes_remaining, scan);
}

StepScan scan_step(const IOdeSolver::SolutionCurve& curve,
                   const interval& domain,
                   const interval& Rc,
                   const interval& Ro,
                   bool ignore_initial_section_root,
                   int maximum_depth,
                   long long maximum_nodes) {
    std::vector<RootCandidate> candidates;
    std::vector<std::pair<SurfaceKind, interval>> unresolved;
    StepScan result;
    long long nodes_remaining = maximum_nodes;
    for (SurfaceKind surface : {
             SurfaceKind::section, SurfaceKind::inner, SurfaceKind::outer}) {
        SurfaceScan surface_scan;
        scan_surface_interval(curve, surface, domain, Rc, Ro,
                              surface == SurfaceKind::section
                                  && ignore_initial_section_root,
                              0, maximum_depth, nodes_remaining, surface_scan);
        candidates.insert(candidates.end(),
                          surface_scan.candidates.begin(),
                          surface_scan.candidates.end());
        for (const interval& time : surface_scan.unresolved_times) {
            unresolved.emplace_back(surface, time);
        }
        result.node_budget_exhausted = result.node_budget_exhausted
            || surface_scan.node_budget_exhausted;
    }

    result.nodes_visited = maximum_nodes - nodes_remaining;
    if (candidates.empty()) {
        if (!unresolved.empty()) {
            result.unresolved = true;
            result.message = "possible event could not be excluded or isolated on "
                + std::string(surface_name(unresolved.front().first));
        }
        return result;
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const RootCandidate& a, const RootCandidate& b) {
                  return a.relative_time.leftBound()
                      < b.relative_time.leftBound();
              });
    const RootCandidate& earliest = candidates.front();
    const double earliest_upper = earliest.relative_time.rightBound();
    for (std::size_t i = 1; i < candidates.size(); ++i) {
        if (!(earliest_upper
              < candidates[i].relative_time.leftBound())) {
            result.unresolved = true;
            result.message = "root-time enclosures do not prove event ordering between "
                + std::string(surface_name(earliest.surface)) + " and "
                + surface_name(candidates[i].surface);
            return result;
        }
    }
    for (const auto& item : unresolved) {
        if (!(earliest_upper < item.second.leftBound())) {
            result.unresolved = true;
            result.message = "an unresolved " + std::string(surface_name(item.first))
                + " interval may precede the candidate event";
            return result;
        }
    }
    result.has_event = true;
    result.event = earliest;
    return result;
}

bool intersect_interval(const interval& a,
                        const interval& b,
                        interval& result) {
    return intersection(a, b, result);
}

void section_cartesian_enclosure(const IVector& reduced_state,
                                 int physical_axis_sign,
                                 interval& x,
                                 interval& xdot,
                                 interval& ydot) {
    const interval& R = reduced_state[0];
    const interval& v = reduced_state[2];
    const interval& w = reduced_state[3];
    if (!strictly_positive(R)) {
        throw std::runtime_error("section enclosure does not prove R>0");
    }
    const interval sign(static_cast<double>(physical_axis_sign));
    const interval R2 = R * R;
    const interval R3 = R2 * R;
    const interval R5 = pow5(R);
    x = sign * R2;
    xdot = sign * v / R3;
    ydot = sign * (w - R5) / R3;
}

AttemptResult validate_once(const ParameterCell& cell,
                            const Settings& settings,
                            const ModelIntervals& p,
                            const FirstNeckCertificate& neck) {
    AttemptResult result;
    result.speed_squared = stable_speed_squared(cell.x, settings, p, neck);
    if (strictly_negative(result.speed_squared)) {
        result.resolved = true;
        result.forbidden = true;
        result.status = "rigorously_forbidden";
        result.message = "speed_squared is strictly negative on the parameter cell";
        return result;
    }
    if (!strictly_positive(result.speed_squared)) {
        result.status = "unresolved_zero_velocity_boundary";
        result.message = "speed_squared contains zero; subdivide the x cell";
        return result;
    }
    result.admissible = true;

    IVector initial(4);
    try {
        initial = initial_box(cell, settings, p, neck,
                              result.speed_squared,
                              result.energy_graph_residual);
    } catch (const std::exception& error) {
        result.status = "unresolved_initial_graph";
        result.message = error.what();
        return result;
    }

    const interval initial_section_derivative =
        cos(initial[1]) * (initial[3] - pow5(initial[0]));
    const interval physical_initial_transversality =
        interval(static_cast<double>(settings.branch_sign))
        * initial_section_derivative / pow3(initial[0]);
    if (!(physical_initial_transversality.leftBound()
          > settings.min_transversality.rightBound())) {
        result.status = "unresolved_initial_transversality";
        result.message = "initial ydot is not uniformly above the proof threshold";
        return result;
    }
    if (!(initial[0].leftBound() > sqrt(settings.inner_radius).rightBound()
          && initial[0].rightBound() < sqrt(settings.outer_radius).leftBound())) {
        result.status = "unresolved_initial_radius";
        result.message = "initial box is not strictly inside both radius guards";
        return result;
    }

    try {
        const interval Rc = sqrt(settings.inner_radius);
        const interval Ro = sqrt(settings.outer_radius);
        const interval signed_delta = settings.side == "open"
            ? settings.delta : -settings.delta;
        const interval h = neck.e1.h + signed_delta;

        IVector restart_state = initial;
        interval restart_time(0.0);
        int physical_axis_sign = settings.branch_sign;
        int event_index = 0;
        int half_crossings = 0;
        const int event_limit = std::max(
            settings.max_half_crossings,
            settings.iterates * settings.max_half_crossings);

        while (result.completed_returns < settings.iterates
               && event_index < event_limit) {
            if (!(restart_time.rightBound() < settings.tau_cap.leftBound())) {
                if (restart_time.leftBound() >= settings.tau_cap.rightBound()) {
                    result.resolved = true;
                    result.status = "validated_tau_cap";
                    result.message =
                        "no further event occurs strictly before the validated tau cap";
                } else {
                    result.status = "unresolved_tau_order";
                    result.message =
                        "restart-time enclosure overlaps the validated tau cap";
                }
                return result;
            }

            IMap field(kValidatedField, 1);
            field.setParameter("A", p.A);
            field.setParameter("B", p.B);
            field.setParameter("ac", p.ac);
            IOdeSolver solver(field, settings.solver_order);
            solver.setAbsoluteTolerance(settings.abs_tol);
            solver.setRelativeTolerance(settings.rel_tol);
            C0HOTripletonSet set(restart_state, restart_time);
            const interval configured_step_cap(settings.max_step);
            solver.setMaxStep(configured_step_cap);
            solver.getStepControl().init(
                solver, set.getCurrentTime(), set);

            bool found_event = false;
            bool first_step_after_restart = true;
            while (set.getCurrentTime().leftBound()
                   < settings.tau_cap.rightBound()) {
                const interval time_before_step = set.getCurrentTime();
                interval retry_step_cap = configured_step_cap;
                bool step_accepted = false;
                std::string last_step_error;
                for (int retry = 0;
                     retry <= settings.max_step_retries; ++retry) {
                    try {
                        solver.setMaxStep(retry_step_cap);
                        solver(set);
                        step_accepted = true;
                        break;
                    } catch (const std::exception& error) {
                        last_step_error = error.what();
                        if (retry == settings.max_step_retries) break;
                        retry_step_cap /= interval(1.5);
                        solver.clearCoefficients();
                        ++result.solver_step_retries;
                        if (retry_step_cap.rightBound()
                            < solver.getStepControl().getMinStepAllowed()) {
                            break;
                        }
                    }
                }
                if (!step_accepted) {
                    result.status = "unresolved_solver_step_retries";
                    result.message =
                        "CAPD could not validate a capped step after retries: "
                        + last_step_error;
                    return result;
                }
                const interval accepted_step = solver.getStep();
                if (!(strictly_positive(accepted_step)
                      && accepted_step.rightBound()
                         <= configured_step_cap.rightBound())) {
                    result.status = "unresolved_solver_step_contract";
                    result.message =
                        "CAPD accepted a nonpositive step or violated max_step";
                    return result;
                }
                const IOdeSolver::SolutionCurve& curve = solver.getCurve();
                const interval step_domain(
                    curve.getLeftDomain(), curve.getRightDomain());
                const StepScan scan = scan_step(
                    curve, step_domain, Rc, Ro,
                    first_step_after_restart,
                    settings.event_time_subdivision_depth,
                    settings.max_event_scan_nodes);
                result.event_scan_nodes += scan.nodes_visited;
                first_step_after_restart = false;
                if (scan.unresolved) {
                    result.status = scan.node_budget_exhausted
                        ? "unresolved_event_scan_node_budget"
                        : "unresolved_dense_event_scan";
                    result.message = scan.node_budget_exhausted
                        ? scan.message + "; per-step event scan node budget exhausted"
                        : scan.message;
                    return result;
                }
                if (!scan.has_event) continue;

                found_event = true;
                const interval event_time =
                    time_before_step + scan.event.relative_time;
                // A non-degenerate restart-time enclosure can make the final
                // dense curve extend slightly beyond the requested decimal
                // cap.  Only classify an event if its complete time enclosure
                // is before the cap; otherwise prove the cap is first or leave
                // the ordering as an explicit proof obligation.
                if (event_time.leftBound() >= settings.tau_cap.rightBound()) {
                    result.resolved = true;
                    result.status = "validated_tau_cap";
                    result.message =
                        "the earliest separately isolated event is not before the requested tau cap";
                    return result;
                }
                if (!(event_time.rightBound() < settings.tau_cap.leftBound())) {
                    result.status = "unresolved_tau_order";
                    result.message =
                        "earliest event-time enclosure overlaps the requested tau cap";
                    return result;
                }
                IVector image = curve(scan.event.relative_time);
                EventRecord event;
                event.event_index = event_index;
                event.tau = event_time;
                event.clock = event_time;
                event.R = image[0];
                event.energy_residual = scaled_energy_residual(image, h, p);
                if (!contains_zero(event.energy_residual)) {
                    result.status = "unresolved_event_energy_consistency";
                    result.message =
                        "event enclosure fails the scaled Hamiltonian consistency check";
                    return result;
                }

                if (scan.event.surface == SurfaceKind::section) {
                    const interval local_cosine = cos(image[1]);
                    int local_axis_sign = 0;
                    if (strictly_positive(local_cosine)) local_axis_sign = 1;
                    else if (strictly_negative(local_cosine)) local_axis_sign = -1;
                    else {
                        result.status = "unresolved_section_chart";
                        result.message =
                            "cos(theta) is not separated from zero at the section root";
                        return result;
                    }
                    const interval section_derivative = surface_derivative(
                        SurfaceKind::section, image);
                    const interval physical_derivative =
                        interval(static_cast<double>(physical_axis_sign))
                        * section_derivative / pow3(image[0]);
                    const int next_physical_axis_sign =
                        physical_axis_sign * local_axis_sign;
                    const interval R4 = pow4(image[0]);
                    const interval R5 = R4 * image[0];
                    const interval Q = interval(2.0) * p.c
                        + interval(2.0) * R4
                        + interval(2.0) * h * pow6(image[0])
                        + p.lambda2 * pow10(image[0])
                        - image[2] * image[2];
                    if (!strictly_positive(Q)) {
                        result.status = "unresolved_energy_reconditioning";
                        result.message =
                            "section energy radicand is not strictly positive";
                        return result;
                    }
                    const interval w_minus_R5 = image[3] - R5;
                    interval energy_w;
                    if (strictly_positive(w_minus_R5)) {
                        energy_w = R5 + sqrt(Q);
                    } else if (strictly_negative(w_minus_R5)) {
                        energy_w = R5 - sqrt(Q);
                    } else {
                        result.status = "unresolved_energy_branch";
                        result.message = "w-R^5 is not separated from zero";
                        return result;
                    }
                    interval tightened_w;
                    if (!intersect_interval(image[3], energy_w, tightened_w)) {
                        result.status = "unresolved_energy_intersection";
                        result.message =
                            "flow and fixed-energy section enclosures do not intersect";
                        return result;
                    }

                    IVector reduced(4);
                    reduced[0] = image[0];
                    reduced[1] = interval(0.0);
                    reduced[2] = image[2];
                    reduced[3] = tightened_w;
                    section_cartesian_enclosure(
                        reduced, next_physical_axis_sign,
                        event.x, event.xdot, event.ydot);
                    interval section_ydot;
                    if (!intersect_interval(
                            event.ydot, physical_derivative, section_ydot)) {
                        result.status = "unresolved_section_velocity_consistency";
                        result.message =
                            "section chart and event derivative ydot enclosures do not intersect";
                        return result;
                    }
                    event.ydot = section_ydot;
                    event.transversality = section_ydot;

                    if (section_ydot.leftBound()
                        > settings.min_transversality.rightBound()) {
                        ++result.completed_returns;
                        half_crossings = 0;
                        event.type = "validated_positive_return";
                        event.positive_return_index = result.completed_returns - 1;
                    } else if (section_ydot.rightBound()
                               < -settings.min_transversality.rightBound()) {
                        ++half_crossings;
                        event.type = "validated_negative_half_crossing";
                        event.positive_return_index = result.completed_returns;
                        if (half_crossings >= settings.max_half_crossings) {
                            result.events.push_back(event);
                            result.status = "unresolved_half_crossing_limit";
                            result.message =
                                "positive return not found within the section-event limit";
                            return result;
                        }
                    } else {
                        result.status = "unresolved_section_orientation";
                        result.message =
                            "physical section derivative is not separated from zero";
                        return result;
                    }
                    result.events.push_back(event);
                    restart_state = reduced;
                    restart_time = event_time;
                    physical_axis_sign = next_physical_axis_sign;
                } else if (scan.event.surface == SurfaceKind::inner) {
                    const interval radial_derivative =
                        interval(0.5) * image[2] * image[0];
                    event.transversality = radial_derivative;
                    if (strictly_negative(radial_derivative)) {
                        event.type = "validated_inner_limit_inward";
                    } else if (strictly_positive(radial_derivative)) {
                        event.type = "validated_inner_limit_outward";
                    } else {
                        result.status = "unresolved_inner_orientation";
                        result.message = "inner-radius derivative contains zero";
                        return result;
                    }
                    cartesian_enclosure(image, physical_axis_sign,
                                        event.x, event.xdot, event.ydot);
                    result.events.push_back(event);
                    result.resolved = true;
                    result.status = event.type;
                    result.message =
                        "first event is the inner study limit; this is not a collision proof";
                    return result;
                } else {
                    const interval radial_derivative =
                        interval(0.5) * image[2] * image[0];
                    event.transversality = radial_derivative;
                    if (!strictly_positive(radial_derivative)) {
                        result.status = "unresolved_outer_orientation";
                        result.message = "outer-radius crossing is not proved outward";
                        return result;
                    }
                    event.type = "validated_outer_limit_outward";
                    cartesian_enclosure(image, physical_axis_sign,
                                        event.x, event.xdot, event.ydot);
                    result.events.push_back(event);
                    result.resolved = true;
                    result.status = event.type;
                    result.message =
                        "trajectory box reached the outer study limit; this is not an escape proof";
                    return result;
                }
                break;
            }

            if (!found_event) {
                result.resolved = true;
                result.status = "validated_tau_cap";
                result.message =
                    "separate capped-step enclosures exclude all events strictly before the tau cap";
                return result;
            }
            ++event_index;
        }

        if (result.completed_returns == settings.iterates) {
            result.resolved = true;
            result.status = "validated_requested_returns";
            result.message = "all requested positive returns were rigorously enclosed";
        } else {
            result.status = "unresolved_event_limit";
            result.message = "event-call limit reached before completing the proof";
        }
    } catch (const capd::dynsys::SolverException<IVector>& error) {
        result.status = "unresolved_solver_exception";
        result.message = error.what();
    } catch (const std::exception& error) {
        result.status = "unresolved_capd_exception";
        result.message = error.what();
    }
    return result;
}

class AdaptiveValidator {
public:
    AdaptiveValidator(const Settings& settings,
                      const ModelIntervals& model,
                      const FirstNeckCertificate& neck)
        : m_settings(settings), m_model(model), m_neck(neck) {}

    void validate(const ParameterCell& cell) {
        if (static_cast<long long>(m_leaves.size()) >= m_settings.max_leaf_boxes) {
            AttemptResult attempt;
            attempt.status = "unresolved_leaf_limit";
            attempt.message = "max_leaf_boxes reached";
            m_leaves.push_back(LeafResult{cell, attempt});
            return;
        }
        AttemptResult attempt = validate_once(cell, m_settings, m_model, m_neck);
        if (attempt.resolved || cell.depth >= m_settings.max_subdivision_depth) {
            if (!attempt.resolved && cell.depth >= m_settings.max_subdivision_depth) {
                attempt.message += "; maximum subdivision depth reached";
            }
            m_leaves.push_back(LeafResult{cell, attempt});
            return;
        }

        const bool prefer_x = attempt.status == "unresolved_zero_velocity_boundary"
            || width(cell.x) >= width(cell.fraction);
        const std::pair<interval, interval> pieces = prefer_x
            ? split_interval(cell.x) : split_interval(cell.fraction);
        if (pieces.first == pieces.second) {
            attempt.message += "; interval cannot be bisected further";
            m_leaves.push_back(LeafResult{cell, attempt});
            return;
        }

        ParameterCell left = cell;
        ParameterCell right = cell;
        ++left.depth;
        ++right.depth;
        left.path += prefer_x ? "x0" : "f0";
        right.path += prefer_x ? "x1" : "f1";
        if (prefer_x) {
            left.x = pieces.first;
            right.x = pieces.second;
        } else {
            left.fraction = pieces.first;
            right.fraction = pieces.second;
        }
        validate(left);
        validate(right);
    }

    const std::vector<LeafResult>& leaves() const { return m_leaves; }

private:
    const Settings& m_settings;
    const ModelIntervals& m_model;
    const FirstNeckCertificate& m_neck;
    std::vector<LeafResult> m_leaves;
};

void print_help(const char* executable) {
    std::cout
        << "Usage: " << executable << " key=value ...\n\n"
        << "Rigorous CAPD interval validation near the first Hill-4BP neck.\n"
        << "Required: delta=<positive decimal>.  Only neck=E1 is supported.\n\n"
        << "Energy/window: side=closed|open, neck_branch=positive|negative,\n"
        << "  neck_window_sigma, or explicit x_begin and x_end.\n"
        << "Grid: grid_mode=points|cover, x_count, start_index, end_index,\n"
        << "  x_radius, fraction_lower, fraction_upper, fraction_count,\n"
        << "  fraction_radius.  cover partitions the full parameter rectangle.\n"
        << "Proof: iterates, solver_order, abs_tol, rel_tol, max_step,\n"
        << "  max_return_time (validated total tau cap), max_half_crossings,\n"
        << "  max_step_retries,\n"
        << "  event_time_subdivision_depth, max_subdivision_depth,\n"
        << "  max_event_scan_nodes (per accepted solver step),\n"
        << "  max_leaf_boxes, collision_radius,\n"
        << "  outer_radius, min_transversality, require_complete,\n"
        << "  require_returns, output, metadata, dry_run.\n\n"
        << "Every numeric literal is parsed as an outward-rounded interval.\n";
}

Settings parse_settings(int argc, char** argv) {
    Settings settings;
    bool x_begin_supplied = false;
    bool x_end_supplied = false;
    bool delta_supplied = false;
    std::set<std::string> supplied;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "help" || argument == "--help" || argument == "-h") {
            print_help(argv[0]);
            std::exit(0);
        }
        const std::size_t equals = argument.find('=');
        if (equals == std::string::npos || equals == 0) {
            throw std::runtime_error("arguments require key=value syntax: " + argument);
        }
        const std::string key = argument.substr(0, equals);
        const std::string value = argument.substr(equals + 1);
        if (!supplied.insert(key).second) {
            throw std::runtime_error("duplicate argument key: " + key);
        }

        if (key == "neck") {
            if (lowercase(value) != "e1") {
                throw std::runtime_error("rigorous validator currently supports neck=E1 only");
            }
        } else if (key == "side") {
            settings.side = lowercase(value);
        } else if (key == "delta") {
            settings.delta_text = value;
            settings.delta = decimal_interval(value, key);
            delta_supplied = true;
        } else if (key == "neck_branch") {
            settings.branch = lowercase(value);
        } else if (key == "neck_window_sigma") {
            settings.neck_window_sigma = decimal_interval(value, key);
        } else if (key == "x_begin") {
            settings.x_begin = decimal_interval(value, key);
            x_begin_supplied = true;
        } else if (key == "x_end") {
            settings.x_end = decimal_interval(value, key);
            x_end_supplied = true;
        } else if (key == "grid_mode") {
            settings.grid_mode = lowercase(value);
        } else if (key == "x_count") {
            settings.x_count = parse_int(value, key);
        } else if (key == "start_index") {
            settings.start_index = parse_int(value, key);
        } else if (key == "end_index") {
            settings.end_index = parse_int(value, key);
        } else if (key == "x_radius") {
            settings.x_radius = decimal_interval(value, key);
        } else if (key == "fraction_lower") {
            settings.fraction_lower = decimal_interval(value, key);
        } else if (key == "fraction_upper") {
            settings.fraction_upper = decimal_interval(value, key);
        } else if (key == "fraction_count") {
            settings.fraction_count = parse_int(value, key);
        } else if (key == "fraction_radius") {
            settings.fraction_radius = decimal_interval(value, key);
        } else if (key == "iterates") {
            settings.iterates = parse_int(value, key);
        } else if (key == "solver_order") {
            settings.solver_order = parse_int(value, key);
        } else if (key == "abs_tol") {
            settings.abs_tol = decimal_double(value, key);
        } else if (key == "rel_tol") {
            settings.rel_tol = decimal_double(value, key);
        } else if (key == "max_step") {
            settings.max_step = decimal_double(value, key);
        } else if (key == "max_step_retries") {
            settings.max_step_retries = parse_int(value, key);
        } else if (key == "max_return_time" || key == "tau_cap") {
            settings.tau_cap = decimal_interval(value, key);
        } else if (key == "max_half_crossings") {
            settings.max_half_crossings = parse_int(value, key);
        } else if (key == "event_time_subdivision_depth") {
            settings.event_time_subdivision_depth = parse_int(value, key);
        } else if (key == "max_event_scan_nodes") {
            settings.max_event_scan_nodes = parse_long_long(value, key);
        } else if (key == "max_subdivision_depth") {
            settings.max_subdivision_depth = parse_int(value, key);
        } else if (key == "max_leaf_boxes") {
            settings.max_leaf_boxes = parse_long_long(value, key);
        } else if (key == "collision_radius" || key == "inner_radius") {
            settings.inner_radius = decimal_interval(value, key);
        } else if (key == "outer_radius") {
            settings.outer_radius = decimal_interval(value, key);
        } else if (key == "min_transversality") {
            settings.min_transversality = decimal_interval(value, key);
        } else if (key == "require_complete") {
            settings.require_complete = parse_bool(value, key);
        } else if (key == "require_returns") {
            settings.require_returns = parse_bool(value, key);
        } else if (key == "progress_every") {
            settings.progress_every = parse_int(value, key);
        } else if (key == "output") {
            settings.output = value;
        } else if (key == "metadata") {
            settings.metadata = value;
        } else if (key == "dry_run") {
            settings.dry_run = parse_bool(value, key);
        } else {
            throw std::runtime_error("unknown case-sensitive argument key: " + key);
        }
    }

    if (!delta_supplied || !strictly_positive(settings.delta)) {
        throw std::runtime_error("delta=<positive nonzero decimal> is required");
    }
    if (settings.side != "closed" && settings.side != "open") {
        throw std::runtime_error("side must be closed or open");
    }
    if (settings.branch == "positive") settings.branch_sign = 1;
    else if (settings.branch == "negative") settings.branch_sign = -1;
    else throw std::runtime_error("neck_branch must be positive or negative");
    if (settings.grid_mode != "points" && settings.grid_mode != "cover") {
        throw std::runtime_error("grid_mode must be points or cover");
    }
    if (x_begin_supplied != x_end_supplied) {
        throw std::runtime_error("set both x_begin and x_end, or neither");
    }
    settings.explicit_window = x_begin_supplied;
    if (settings.x_count < 1 || settings.fraction_count < 1
        || settings.iterates < 1 || settings.solver_order < 3
        || settings.solver_order > 64
        || settings.max_half_crossings < 1
        || settings.max_step_retries < 0
        || settings.event_time_subdivision_depth < 1
        || settings.max_event_scan_nodes < 1
        || settings.max_subdivision_depth < 0
        || settings.max_leaf_boxes < 1) {
        throw std::runtime_error("counts/orders/limits are outside their valid ranges");
    }
    if (settings.end_index < 0) settings.end_index = settings.x_count - 1;
    if (settings.start_index < 0 || settings.end_index < settings.start_index
        || settings.end_index >= settings.x_count) {
        throw std::runtime_error("invalid start_index/end_index");
    }
    if (!(settings.abs_tol > 0.0 && settings.rel_tol > 0.0
          && settings.max_step > 0.0)) {
        throw std::runtime_error("solver tolerances and max_step must be positive");
    }
    if (!strictly_positive(settings.tau_cap)) {
        throw std::runtime_error("max_return_time/tau_cap must be strictly positive");
    }
    if (!strictly_positive(settings.inner_radius)) {
        throw std::runtime_error("collision_radius/inner_radius must be strictly positive");
    }
    if (!(settings.outer_radius.leftBound()
          > settings.inner_radius.rightBound())) {
        throw std::runtime_error(
            "outer_radius must be rigorously greater than inner_radius");
    }
    if (settings.min_transversality.leftBound() < 0.0) {
        throw std::runtime_error("min_transversality must be nonnegative");
    }
    if (settings.x_radius.leftBound() < 0.0) {
        throw std::runtime_error("x_radius must be nonnegative");
    }
    if (settings.fraction_radius.leftBound() < 0.0) {
        throw std::runtime_error("fraction_radius must be nonnegative");
    }
    if (settings.neck_window_sigma.leftBound() <= 0.0) {
        throw std::runtime_error("neck_window_sigma must be strictly positive");
    }
    if (!(settings.fraction_lower.rightBound()
          < settings.fraction_upper.leftBound())
        || settings.fraction_lower.leftBound() <= -1.0
        || settings.fraction_upper.rightBound() >= 1.0) {
        throw std::runtime_error("require -1 < fraction_lower < fraction_upper < 1");
    }
    if (settings.output.empty() || settings.metadata.empty()) {
        throw std::runtime_error("output and metadata paths must not be empty");
    }
    if (settings.progress_every < 1) settings.progress_every = 1;
    return settings;
}

void write_interval_pair(std::ostream& output, const Bounds& value) {
    output << value.lo << ',' << value.hi;
}

void write_interval_pair(std::ostream& output, const interval& value) {
    write_interval_pair(output, bounds(value));
}

void write_nan_pair(std::ostream& output) {
    output << "nan,nan";
}

int run(int argc, char** argv) {
    const Settings settings = parse_settings(argc, argv);
    const RawConfig raw_config = read_raw_config();
    const ModelIntervals model = load_model(raw_config);
    const FirstNeckCertificate neck = certify_first_neck(model);

    interval x_begin;
    interval x_end;
    interval neck_length_scale = sqrt(
        interval(2.0) * settings.delta / neck.curvature);
    if (settings.explicit_window) {
        x_begin = settings.x_begin;
        x_end = settings.x_end;
    } else {
        const interval center = interval(static_cast<double>(settings.branch_sign))
            * neck.e1.rho;
        x_begin = center - settings.neck_window_sigma * neck_length_scale;
        x_end = center + settings.neck_window_sigma * neck_length_scale;
    }
    if (!(x_begin.rightBound() < x_end.leftBound())) {
        throw std::runtime_error("x window endpoints are not rigorously ordered");
    }
    if ((settings.branch_sign > 0 && x_begin.leftBound() <= 0.0)
        || (settings.branch_sign < 0 && x_end.rightBound() >= 0.0)) {
        throw std::runtime_error("x window crosses zero or disagrees with neck_branch");
    }

    const interval energy_offset = settings.side == "open"
        ? settings.delta : -settings.delta;
    const interval h = neck.e1.h + energy_offset;
    const bool side_certified = settings.side == "open"
        ? strictly_positive(energy_offset) : strictly_negative(energy_offset);
    if (!side_certified) {
        throw std::runtime_error("failed to certify the requested side of E1");
    }

    std::cerr << std::setprecision(17)
              << "rigorous_interval_mode=1\n"
              << "parameter_file=" << raw_config.path << '\n'
              << "E1_rho=" << neck.e1.rho << " E1_hcrit=" << neck.e1.h << '\n'
              << "E3_rho=" << neck.e3.rho << " E3_hcrit=" << neck.e3.h << '\n'
              << "Omega_xx(E1)=" << neck.omega_xx
              << " Omega_yy(E1)=" << neck.omega_yy << '\n'
              << "first_neck_certified=" << neck.passed
              << " global_axis_roots=" << neck.global_axis_root_uniqueness
              << " no_mixed_axis=" << neck.no_mixed_axis_equilibria
              << " side=" << settings.side
              << " delta=" << settings.delta
              << " h=" << h << '\n'
              << "x_window=" << hull(x_begin, x_end)
              << " neck_length_scale=" << neck_length_scale
              << " grid_mode=" << settings.grid_mode << '\n'
              << "grid=" << settings.x_count << 'x' << settings.fraction_count
              << " shard=[" << settings.start_index << ',' << settings.end_index << ']'
              << " iterates=" << settings.iterates
              << " tau_cap=" << settings.tau_cap
              << " max_step_retries=" << settings.max_step_retries
              << " event_scan_node_budget=" << settings.max_event_scan_nodes
              << '\n';

    if (settings.dry_run) {
        std::cerr << "dry_run=1; no CAPD flow or output files constructed\n";
        return 0;
    }

    const auto started = std::chrono::steady_clock::now();
    AdaptiveValidator validator(settings, model, neck);
    long long root_cells = 0;
    for (int x_index = settings.start_index;
         x_index <= settings.end_index; ++x_index) {
        const interval x_cell = grid_cell(
            x_begin, x_end, x_index, settings.x_count,
            settings.grid_mode, settings.x_radius);
        for (int fraction_index = 0;
             fraction_index < settings.fraction_count; ++fraction_index) {
            const interval fraction_cell = grid_cell(
                settings.fraction_lower, settings.fraction_upper,
                fraction_index, settings.fraction_count,
                settings.grid_mode, settings.fraction_radius);
            ParameterCell cell;
            cell.x = x_cell;
            cell.fraction = fraction_cell;
            cell.x_index = x_index;
            cell.fraction_index = fraction_index;
            validator.validate(cell);
            ++root_cells;
            if (root_cells % settings.progress_every == 0) {
                std::cerr << "validated root cells=" << root_cells
                          << " leaves=" << validator.leaves().size() << '\n';
            }
        }
    }

    std::ofstream output(settings.output);
    if (!output) throw std::runtime_error("failed to open " + settings.output);
    output << std::setprecision(17);
    output
        << "x_index,fraction_index,leaf_path,subdivision_depth,"
        << "x0_lo,x0_hi,fraction_lo,fraction_hi,admissible,forbidden,resolved,"
        << "completed_returns,requested_returns,status,event_type,event_index,"
        << "positive_return_index,tau_lo,tau_hi,clock_lo,clock_hi,"
        << "x_lo,x_hi,xdot_lo,xdot_hi,ydot_lo,ydot_hi,"
        << "transversality_lo,transversality_hi,R_lo,R_hi,"
        << "h_lo,h_hi,event_energy_residual_lo,event_energy_residual_hi,"
        << "speed_squared_lo,speed_squared_hi,"
        << "energy_graph_residual_lo,energy_graph_residual_hi,"
        << "event_scan_nodes,solver_step_retries,message\n";

    long long resolved_leaves = 0;
    long long unresolved_leaves = 0;
    long long admissible_leaves = 0;
    long long forbidden_leaves = 0;
    long long requested_return_leaves = 0;
    long long event_scan_nodes_total = 0;
    long long solver_step_retries_total = 0;
    std::map<std::string, long long> status_counts;

    for (const LeafResult& leaf : validator.leaves()) {
        const AttemptResult& attempt = leaf.attempt;
        ++status_counts[attempt.status];
        resolved_leaves += attempt.resolved ? 1 : 0;
        unresolved_leaves += attempt.resolved ? 0 : 1;
        admissible_leaves += attempt.admissible ? 1 : 0;
        forbidden_leaves += attempt.forbidden ? 1 : 0;
        requested_return_leaves +=
            attempt.completed_returns == settings.iterates ? 1 : 0;
        event_scan_nodes_total += attempt.event_scan_nodes;
        solver_step_retries_total += attempt.solver_step_retries;

        const std::size_t rows = std::max<std::size_t>(1, attempt.events.size());
        for (std::size_t row = 0; row < rows; ++row) {
            output << leaf.cell.x_index << ',' << leaf.cell.fraction_index << ','
                   << leaf.cell.path << ',' << leaf.cell.depth << ',';
            write_interval_pair(output, leaf.cell.x); output << ',';
            write_interval_pair(output, leaf.cell.fraction); output << ',';
            output << attempt.admissible << ',' << attempt.forbidden << ','
                   << attempt.resolved << ',' << attempt.completed_returns << ','
                   << settings.iterates << ',' << attempt.status << ',';
            if (attempt.events.empty()) {
                output << "none,-1,-1,";
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_nan_pair(output); output << ',';
                write_interval_pair(output, h); output << ',';
                write_nan_pair(output); output << ',';
            } else {
                const EventRecord& event = attempt.events[row];
                output << event.type << ',' << event.event_index << ','
                       << event.positive_return_index << ',';
                write_interval_pair(output, event.tau); output << ',';
                write_interval_pair(output, event.clock); output << ',';
                write_interval_pair(output, event.x); output << ',';
                write_interval_pair(output, event.xdot); output << ',';
                write_interval_pair(output, event.ydot); output << ',';
                write_interval_pair(output, event.transversality); output << ',';
                write_interval_pair(output, event.R); output << ',';
                write_interval_pair(output, h); output << ',';
                write_interval_pair(output, event.energy_residual); output << ',';
            }
            write_interval_pair(output, attempt.speed_squared); output << ',';
            write_interval_pair(output, attempt.energy_graph_residual); output << ',';
            output << attempt.event_scan_nodes << ','
                   << attempt.solver_step_retries << ',';
            output << csv_escape(attempt.events.empty()
                                     ? attempt.message
                                     : attempt.events[row].type + ": " + attempt.message)
                   << '\n';
        }
    }

    const bool all_cells_resolved = unresolved_leaves == 0;
    const bool all_admissible_returned = admissible_leaves == requested_return_leaves;
    const bool proof_complete = neck.passed && side_certified && all_cells_resolved
        && (!settings.require_returns || all_admissible_returned);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::ofstream metadata(settings.metadata);
    if (!metadata) throw std::runtime_error("failed to open " + settings.metadata);
    metadata << std::setprecision(17);
    metadata << "format=hill4bp_capd_interval_proof_v1\n";
    metadata << "rigorous_interval_arithmetic=1\n";
#if defined(__USE_NATIVE__)
    metadata << "capd_interval_backend=NATIVE\n";
#elif defined(__USE_FILIB__)
    metadata << "capd_interval_backend=FILIB\n";
#elif defined(__USE_CXSC__)
    metadata << "capd_interval_backend=CXSC\n";
#else
    metadata << "capd_interval_backend=unknown\n";
#endif
    metadata << "compiler=" << __VERSION__ << '\n';
    metadata << "cplusplus=" << __cplusplus << '\n';
    metadata << "fast_math=0\n";
    metadata << "time_coordinate=regularized_tau\n";
    metadata << "clock_column=absolute_regularized_tau\n";
    metadata << "proof_claim=outer_enclosures_of_fixed_energy_graph_returns_and_separately_ordered_events\n";
    metadata << "proof_limitations=no_claim_of_collision_escape_invariant_curve_chaos_or_delta_to_zero_limit\n";
    metadata << "parameter_file=" << raw_config.path << '\n';
    metadata << "first_neck_certified=" << neck.passed << '\n';
    metadata << "positive_lambda_c_certified="
             << neck.positive_lambda_c << '\n';
    metadata << "global_axis_root_uniqueness_certified="
             << neck.global_axis_root_uniqueness << '\n';
    metadata << "no_mixed_axis_equilibria_certified="
             << neck.no_mixed_axis_equilibria << '\n';
    metadata << "E1_unique_root_certified="
             << neck.global_axis_root_uniqueness << '\n';
    metadata << "E3_unique_root_certified="
             << neck.global_axis_root_uniqueness << '\n';
    metadata << "E1_index_one_certified=" << neck.e1_index_one << '\n';
    metadata << "E1_before_E3_certified=" << neck.e1_before_e3 << '\n';
    metadata << "E1_rho_lo=" << neck.e1.rho.leftBound() << '\n';
    metadata << "E1_rho_hi=" << neck.e1.rho.rightBound() << '\n';
    metadata << "E1_hcrit_lo=" << neck.e1.h.leftBound() << '\n';
    metadata << "E1_hcrit_hi=" << neck.e1.h.rightBound() << '\n';
    metadata << "E3_hcrit_lo=" << neck.e3.h.leftBound() << '\n';
    metadata << "E3_hcrit_hi=" << neck.e3.h.rightBound() << '\n';
    metadata << "Omega_xx_E1_lo=" << neck.omega_xx.leftBound() << '\n';
    metadata << "Omega_xx_E1_hi=" << neck.omega_xx.rightBound() << '\n';
    metadata << "Omega_yy_E1_lo=" << neck.omega_yy.leftBound() << '\n';
    metadata << "Omega_yy_E1_hi=" << neck.omega_yy.rightBound() << '\n';
    metadata << "side=" << settings.side << '\n';
    metadata << "side_certified=" << side_certified << '\n';
    metadata << "delta_literal=" << settings.delta_text << '\n';
    metadata << "delta_lo=" << settings.delta.leftBound() << '\n';
    metadata << "delta_hi=" << settings.delta.rightBound() << '\n';
    metadata << "h_lo=" << h.leftBound() << '\n';
    metadata << "h_hi=" << h.rightBound() << '\n';
    metadata << "neck_branch=" << settings.branch << '\n';
    metadata << "grid_mode=" << settings.grid_mode << '\n';
    metadata << "x_begin_lo=" << x_begin.leftBound() << '\n';
    metadata << "x_begin_hi=" << x_begin.rightBound() << '\n';
    metadata << "x_end_lo=" << x_end.leftBound() << '\n';
    metadata << "x_end_hi=" << x_end.rightBound() << '\n';
    metadata << "x_count=" << settings.x_count << '\n';
    metadata << "start_index=" << settings.start_index << '\n';
    metadata << "end_index=" << settings.end_index << '\n';
    metadata << "fraction_count=" << settings.fraction_count << '\n';
    metadata << "iterates=" << settings.iterates << '\n';
    metadata << "solver_order=" << settings.solver_order << '\n';
    metadata << "abs_tol=" << settings.abs_tol << '\n';
    metadata << "rel_tol=" << settings.rel_tol << '\n';
    metadata << "max_step=" << settings.max_step << '\n';
    metadata << "max_step_retries=" << settings.max_step_retries << '\n';
    metadata << "tau_cap_lo=" << settings.tau_cap.leftBound() << '\n';
    metadata << "tau_cap_hi=" << settings.tau_cap.rightBound() << '\n';
    metadata << "tau_horizon=half_open_[0,tau_cap)\n";
    metadata << "tau_cap_proof_stop=all_trajectory_times_at_or_above_requested_upper_endpoint\n";
    metadata << "configured_max_step_enforced_and_checked=1\n";
    metadata << "integration_driver=explicit_adaptive_IOdeSolver_steps\n";
    metadata << "event_detection=separate_capped_IOdeSolver_dense_curve_scan\n";
    metadata << "event_surfaces=sin(th),R-Rc,Ro-R\n";
    metadata << "event_derivatives=direct_interval_vector_field_formulas\n";
    metadata << "section_chart_reset=pi_periodicity_with_cosine_parity\n";
    metadata << "event_existence=strict_endpoint_bracket_or_strict_parameterized_interval_newton_inclusion\n";
    metadata << "event_localization=interval_newton\n";
    metadata << "event_time_subdivision_depth="
             << settings.event_time_subdivision_depth << '\n';
    metadata << "event_scan_node_budget_per_step="
             << settings.max_event_scan_nodes << '\n';
    metadata << "event_scan_nodes_total=" << event_scan_nodes_total << '\n';
    metadata << "solver_step_retries_total="
             << solver_step_retries_total << '\n';
    metadata << "set_representation=C0HOTripletonSet\n";
    metadata << "root_cells=" << root_cells << '\n';
    metadata << "leaf_cells=" << validator.leaves().size() << '\n';
    metadata << "resolved_leaf_cells=" << resolved_leaves << '\n';
    metadata << "unresolved_leaf_cells=" << unresolved_leaves << '\n';
    metadata << "admissible_leaf_cells=" << admissible_leaves << '\n';
    metadata << "forbidden_leaf_cells=" << forbidden_leaves << '\n';
    metadata << "requested_return_leaf_cells=" << requested_return_leaves << '\n';
    metadata << "all_cells_resolved=" << all_cells_resolved << '\n';
    metadata << "all_admissible_requested_returns=" << all_admissible_returned << '\n';
    metadata << "require_complete=" << settings.require_complete << '\n';
    metadata << "require_returns=" << settings.require_returns << '\n';
    metadata << "proof_complete=" << proof_complete << '\n';
    metadata << "elapsed_seconds=" << elapsed << '\n';
    metadata << "lambda1=" << model.lambda1 << '\n';
    metadata << "lambda2=" << model.lambda2 << '\n';
    metadata << "A=" << model.A << '\n';
    metadata << "B=" << model.B << '\n';
    metadata << "c=" << model.c << '\n';
    for (const auto& item : status_counts) {
        metadata << "status_count_" << item.first << '=' << item.second << '\n';
    }

    std::cerr << "interval output -> " << settings.output << '\n'
              << "proof metadata -> " << settings.metadata << '\n'
              << "proof_complete=" << proof_complete
              << " resolved=" << resolved_leaves
              << " unresolved=" << unresolved_leaves
              << " elapsed=" << elapsed << "s\n";

    if (settings.require_complete && !proof_complete) {
        std::cerr << "proof incomplete: inspect unresolved rows before making a claim\n";
        return 2;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
