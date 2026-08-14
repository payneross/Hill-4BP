#include "crtbp_cpp/hill4bp.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Point {
    double x = 0.0;
    double xdot = 0.0;
};

struct Config {
    double h = 0.0;
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
    std::string output_csv = "poincare_grid.csv";
    std::string output_svg = "poincare_section.svg";
    std::string output_metadata = "poincare_metadata.txt";
    bool write_svg = true;
    bool quiet = false;
    int svg_stride = 1;
};

void print_help(const char* executable) {
    std::cout
        << "Usage: " << executable << " [options]\n\n"
        << "Generates the regularized planar Hill 4BP Poincare section from the MATLAB\n"
        << "Hill_4BPpoincareGrid.m workflow. Full-resolution points are written to CSV.\n\n"
        << "Options:\n"
        << "  --h VALUE                Hamiltonian energy (default from shared config)\n"
        << "  --x-count N              Number of x samples (default from shared config)\n"
        << "  --xdot-count N           Number of xdot-direction rows (default from shared config)\n"
        << "  --iterates N             Crossings per initial condition (default from shared config)\n"
        << "  --x-begin VALUE          First x value (default from shared config)\n"
        << "  --x-end VALUE            Last x value (default from shared config)\n"
        << "  --xdot-lower VALUE       First xdot direction value (default from shared config)\n"
        << "  --xdot-upper VALUE       Last xdot direction value (default from shared config)\n"
        << "  --y-section VALUE        Initial y offset (default from shared config)\n"
        << "  --tau-step VALUE         Integration segment length (default from shared config)\n"
        << "  --rel-tol VALUE          Relative ODE tolerance (default from shared config)\n"
        << "  --abs-tol VALUE          Absolute ODE tolerance (default from shared config)\n"
        << "  --max-step VALUE         Maximum adaptive step (default from shared config)\n"
        << "  --initial-step VALUE     Initial adaptive step (default from shared config)\n"
        << "  --output PATH            CSV output path (default poincare_grid.csv)\n"
        << "  --svg PATH               SVG output path (default poincare_section.svg)\n"
        << "  --metadata PATH          Metadata output path (default poincare_metadata.txt)\n"
        << "  --svg-stride N           Plot every Nth point in SVG (default 1)\n"
        << "  --no-svg                 Skip SVG preview generation\n"
        << "  --quiet                  Reduce progress logging\n"
        << "  --quick                  Verification preset: 8 x nodes, 12 crossings\n"
        << "  --help                   Show this help\n";
}

double parse_double(const std::string& value, const std::string& option) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0') {
        throw std::invalid_argument("invalid numeric value for " + option + ": " + value);
    }
    return parsed;
}

int parse_int(const std::string& value, const std::string& option) {
    char* end = nullptr;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' || parsed < 1
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("invalid integer value for " + option + ": " + value);
    }
    return static_cast<int>(parsed);
}

std::string require_value(int& i, int argc, char** argv) {
    if (i + 1 >= argc) {
        throw std::invalid_argument(std::string("missing value for ") + argv[i]);
    }
    ++i;
    return argv[i];
}

Config parse_args(int argc, char** argv, Config config) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            print_help(argv[0]);
            std::exit(0);
        } else if (arg == "--h") {
            config.h = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--x-count") {
            config.x_count = parse_int(require_value(i, argc, argv), arg);
        } else if (arg == "--xdot-count") {
            config.xdot_count = parse_int(require_value(i, argc, argv), arg);
        } else if (arg == "--iterates") {
            config.iterates = parse_int(require_value(i, argc, argv), arg);
        } else if (arg == "--x-begin") {
            config.x_begin = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--x-end") {
            config.x_end = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--xdot-lower") {
            config.xdot_lower = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--xdot-upper") {
            config.xdot_upper = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--y-section") {
            config.y_section = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--tau-step") {
            config.tau_step = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--rel-tol") {
            config.rel_tol = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--abs-tol") {
            config.abs_tol = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--max-step") {
            config.max_step = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--initial-step") {
            config.initial_step = parse_double(require_value(i, argc, argv), arg);
        } else if (arg == "--output") {
            config.output_csv = require_value(i, argc, argv);
        } else if (arg == "--svg") {
            config.output_svg = require_value(i, argc, argv);
        } else if (arg == "--metadata") {
            config.output_metadata = require_value(i, argc, argv);
        } else if (arg == "--svg-stride") {
            config.svg_stride = parse_int(require_value(i, argc, argv), arg);
        } else if (arg == "--no-svg") {
            config.write_svg = false;
        } else if (arg == "--quiet") {
            config.quiet = true;
        } else if (arg == "--quick") {
            config.x_count = 8;
            config.iterates = 12;
            config.max_step = 0.1;
            config.initial_step = 1e-3;
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    if (!(config.tau_step > 0.0)) {
        throw std::invalid_argument("--tau-step must be positive");
    }
    if (!(config.max_step > 0.0) || !(config.initial_step > 0.0)) {
        throw std::invalid_argument("--max-step and --initial-step must be positive");
    }
    return config;
}

Config config_from_shared(const crtbp::Hill4BPPoincareConfig& shared,
                          const crtbp::Hill4BPEquilibria& equilibria) {
    Config config;
    config.h = shared.h_e1_scale * equilibria.E1.h;
    config.x_count = shared.x_count;
    config.xdot_count = shared.xdot_count;
    config.iterates = shared.iterates;
    config.x_begin = shared.x_begin;
    config.x_end = shared.x_end;
    config.xdot_lower = shared.xdot_lower;
    config.xdot_upper = shared.xdot_upper;
    config.y_section = shared.y_section;
    config.tau_step = shared.tau_step;
    config.rel_tol = shared.rel_tol;
    config.abs_tol = shared.abs_tol;
    config.initial_step = shared.initial_step;
    config.max_step = shared.max_step;
    config.svg_stride = shared.svg_stride;
    return config;
}

bool finite_state(const crtbp::State4& state) {
    return std::all_of(state.begin(), state.end(), [](double value) {
        return std::isfinite(value);
    });
}

std::string fmt(double value, int precision = 6) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();
}

void write_csv(const std::string& path, const std::vector<Point>& points) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("failed to open CSV output: " + path);
    }

    out << "x,xdot\n";
    out << std::setprecision(17);
    for (const Point& point : points) {
        out << point.x << ',' << point.xdot << '\n';
    }
}

void write_metadata(const std::string& path,
                    const Config& config,
                    const crtbp::Hill4BPParameters& params,
                    const crtbp::Hill4BPEquilibria& equilibria,
                    std::size_t point_count,
                    int skipped,
                    int truncated) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("failed to open metadata output: " + path);
    }

    out << std::setprecision(17);
    out << "model=regularized_planar_hill4bp\n";
    out << "regularized_coordinate=R=sqrt(rho)=paper_r^(1/(alpha+2))\n";
    out << "h=" << config.h << '\n';
    out << "jacobi_C=" << -2.0 * config.h << '\n';
    out << "lambda1=" << params.lambda1 << '\n';
    out << "lambda2=" << params.lambda2 << '\n';
    out << "A=" << params.A << '\n';
    out << "B=" << params.B << '\n';
    out << "c=" << params.c << '\n';
    out << "E1_h=" << equilibria.E1.h << '\n';
    out << "E3_h=" << equilibria.E3.h << '\n';
    out << "x_begin=" << config.x_begin << '\n';
    out << "x_end=" << config.x_end << '\n';
    out << "x_count=" << config.x_count << '\n';
    out << "xdot_lower=" << config.xdot_lower << '\n';
    out << "xdot_upper=" << config.xdot_upper << '\n';
    out << "xdot_count=" << config.xdot_count << '\n';
    out << "iterates=" << config.iterates << '\n';
    out << "tau_step=" << config.tau_step << '\n';
    out << "rel_tol=" << config.rel_tol << '\n';
    out << "abs_tol=" << config.abs_tol << '\n';
    out << "max_step=" << config.max_step << '\n';
    out << "point_count=" << point_count << '\n';
    out << "skipped_initial_conditions=" << skipped << '\n';
    out << "truncated_initial_conditions=" << truncated << '\n';
}

void write_svg(const std::string& path,
               const std::vector<Point>& points,
               const Config& config) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("failed to open SVG output: " + path);
    }

    constexpr double width = 1200.0;
    constexpr double height = 820.0;
    constexpr double left = 86.0;
    constexpr double right = 34.0;
    constexpr double top = 54.0;
    constexpr double bottom = 78.0;

    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const Point& point : points) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.xdot);
        max_y = std::max(max_y, point.xdot);
    }

    if (points.empty()) {
        min_x = config.x_begin;
        max_x = config.x_end;
        min_y = -1.0;
        max_y = 1.0;
    }

    if (min_x == max_x) {
        min_x -= 0.5;
        max_x += 0.5;
    }
    if (min_y == max_y) {
        min_y -= 0.5;
        max_y += 0.5;
    }

    const double pad_x = 0.04 * (max_x - min_x);
    const double pad_y = 0.06 * (max_y - min_y);
    min_x -= pad_x;
    max_x += pad_x;
    min_y -= pad_y;
    max_y += pad_y;

    auto map_x = [&](double x) {
        return left + (x - min_x) / (max_x - min_x) * (width - left - right);
    };
    auto map_y = [&](double y) {
        return height - bottom - (y - min_y) / (max_y - min_y) * (height - top - bottom);
    };

    out << std::fixed << std::setprecision(3);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height
        << "\" width=\"" << width << "\" height=\"" << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    out << "<text x=\"" << width / 2.0
        << "\" y=\"28\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" "
        << "font-size=\"20\" fill=\"#172033\">Planar regularized Hill 4BP Poincare section, h = "
        << fmt(config.h, 6) << "</text>\n";

    out << "<g stroke=\"#d8dee8\" stroke-width=\"1\">\n";
    for (int i = 0; i <= 10; ++i) {
        const double gx = left + i * (width - left - right) / 10.0;
        const double gy = top + i * (height - top - bottom) / 10.0;
        out << "<line x1=\"" << gx << "\" y1=\"" << top << "\" x2=\"" << gx << "\" y2=\""
            << height - bottom << "\"/>\n";
        out << "<line x1=\"" << left << "\" y1=\"" << gy << "\" x2=\"" << width - right
            << "\" y2=\"" << gy << "\"/>\n";
    }
    out << "</g>\n";

    out << "<rect x=\"" << left << "\" y=\"" << top << "\" width=\"" << width - left - right
        << "\" height=\"" << height - top - bottom
        << "\" fill=\"none\" stroke=\"#263445\" stroke-width=\"1.5\"/>\n";

    out << "<g font-family=\"Helvetica,Arial,sans-serif\" font-size=\"13\" fill=\"#334155\">\n";
    for (int i = 0; i <= 5; ++i) {
        const double alpha = static_cast<double>(i) / 5.0;
        const double x_value = min_x + alpha * (max_x - min_x);
        const double y_value = min_y + alpha * (max_y - min_y);
        out << "<text x=\"" << map_x(x_value) << "\" y=\"" << height - bottom + 24
            << "\" text-anchor=\"middle\">" << fmt(x_value, 4) << "</text>\n";
        out << "<text x=\"" << left - 12 << "\" y=\"" << map_y(y_value) + 4
            << "\" text-anchor=\"end\">" << fmt(y_value, 4) << "</text>\n";
    }
    out << "<text x=\"" << width / 2.0 << "\" y=\"" << height - 22
        << "\" text-anchor=\"middle\" font-size=\"15\">x1</text>\n";
    out << "<text x=\"22\" y=\"" << height / 2.0
        << "\" text-anchor=\"middle\" transform=\"rotate(-90 22 " << height / 2.0
        << ")\" font-size=\"15\">dx1/dt</text>\n";
    out << "</g>\n";

    out << "<g fill=\"#1f5fbf\" fill-opacity=\"0.68\">\n";
    for (std::size_t i = 0; i < points.size(); i += static_cast<std::size_t>(config.svg_stride)) {
        out << "<circle cx=\"" << map_x(points[i].x) << "\" cy=\"" << map_y(points[i].xdot)
            << "\" r=\"1.15\"/>\n";
    }
    out << "</g>\n";
    out << "</svg>\n";
}

std::vector<Point> generate_poincare(const Config& config,
                                     const crtbp::Hill4BPParameters& params,
                                     int& skipped_initial_conditions,
                                     int& truncated_initial_conditions) {
    crtbp::OdeOptions ode_options;
    ode_options.rel_tol = config.rel_tol;
    ode_options.abs_tol = config.abs_tol;
    ode_options.initial_step = config.initial_step;
    ode_options.max_step = config.max_step;
    ode_options.min_step = 1e-13;
    ode_options.max_steps = 10000000;

    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(config.x_count)
                   * static_cast<std::size_t>(config.xdot_count)
                   * static_cast<std::size_t>(config.iterates));

    const double x_step = config.x_count > 1
        ? (config.x_end - config.x_begin) / static_cast<double>(config.x_count - 1)
        : 0.0;
    const double xdot_step = config.xdot_count > 1
        ? (config.xdot_upper - config.xdot_lower) / static_cast<double>(config.xdot_count - 1)
        : 0.0;

    skipped_initial_conditions = 0;
    truncated_initial_conditions = 0;
    const int max_segments_per_initial_condition = 5 * config.iterates + 100;

    for (int m = 0; m < config.xdot_count; ++m) {
        for (int n = 0; n < config.x_count; ++n) {
            if (!config.quiet) {
                std::cout << "velocity row " << (m + 1) << "/" << config.xdot_count
                          << ", x node " << (n + 1) << "/" << config.x_count
                          << ", crossings " << points.size() << '\n';
            }

            const double x0 = config.x_begin + static_cast<double>(n) * x_step;
            const double xdot_direction = config.xdot_lower + static_cast<double>(m) * xdot_step;
            const crtbp::Vec2 q0{{x0, config.y_section}};
            const double direction_norm = std::hypot(xdot_direction, 1.0);
            const double dir_x = xdot_direction / direction_norm;
            const double dir_y = 1.0 / direction_norm;

            const auto speed = crtbp::hill4bp_speed(q0, config.h, params);
            if (!speed.second) {
                ++skipped_initial_conditions;
                continue;
            }

            const crtbp::State4 initial_cartesian{{q0[0], q0[1], speed.first * dir_x, speed.first * dir_y}};
            crtbp::State4 initial_regularized =
                crtbp::hill4bp_cartesian_to_sqrt_radius(initial_cartesian, params);

            int local_iterates = 0;
            int segment_count = 0;
            while (local_iterates < config.iterates
                   && segment_count < max_segments_per_initial_condition) {
                ++segment_count;

                crtbp::OdeResult result = crtbp::integrate_with_events(
                    [&params](double tau, const crtbp::Vector& state) {
                        return crtbp::hill4bp_sqrt_radius_field(tau, state, params);
                    },
                    0.0,
                    crtbp::to_vector(initial_regularized),
                    config.tau_step,
                    ode_options,
                    [&params](double tau, const crtbp::Vector& state) {
                        return crtbp::hill4bp_sqrt_radius_section_event(
                            tau, state, params);
                    },
                    1);

                if (!result.success || !crtbp::all_finite(result.state)) {
                    break;
                }

                for (const crtbp::OdeEvent& event : result.events) {
                    const crtbp::State4 event_state = crtbp::to_state4(event.state);
                    const crtbp::State4 event_cartesian =
                        crtbp::hill4bp_sqrt_radius_to_cartesian(event_state, params);
                    if (!finite_state(event_cartesian) || event_cartesian[3] <= 0.0) {
                        continue;
                    }

                    points.push_back({event_cartesian[0], event_cartesian[2]});
                    ++local_iterates;
                    if (local_iterates >= config.iterates) {
                        break;
                    }
                }

                initial_regularized = crtbp::to_state4(result.state);
            }

            if (local_iterates < config.iterates) {
                ++truncated_initial_conditions;
            }
        }
    }

    return points;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help") {
                print_help(argv[0]);
                return 0;
            }
        }

        const crtbp::Hill4BPParameters params = crtbp::hill4bp_parameters();
        const crtbp::Hill4BPEquilibria equilibria = crtbp::hill4bp_equilibria(params);
        const crtbp::Hill4BPPoincareConfig shared = crtbp::hill4bp_poincare_config();
        const Config config = parse_args(argc, argv, config_from_shared(shared, equilibria));

        if (!config.quiet) {
            std::cout << std::setprecision(15);
            std::cout << "Hill 4BP Hamiltonian h = " << config.h
                      << ", Jacobi C = " << -2.0 * config.h << '\n';
            std::cout << "E1/E2 h = " << equilibria.E1.h
                      << ", E3/E4 h = " << equilibria.E3.h << '\n';
        }

        int skipped = 0;
        int truncated = 0;
        const std::vector<Point> points = generate_poincare(config, params, skipped, truncated);

        write_csv(config.output_csv, points);
        write_metadata(config.output_metadata, config, params, equilibria, points.size(), skipped, truncated);
        if (config.write_svg) {
            write_svg(config.output_svg, points, config);
        }

        std::cout << "Wrote " << points.size() << " Poincare points to " << config.output_csv
                  << '\n';
        std::cout << "Metadata: " << config.output_metadata << '\n';
        if (config.write_svg) {
            std::cout << "SVG preview: " << config.output_svg << '\n';
        }
        if (skipped > 0 || truncated > 0) {
            std::cout << "Skipped initial conditions: " << skipped
                      << ", truncated initial conditions: " << truncated << '\n';
        }
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
