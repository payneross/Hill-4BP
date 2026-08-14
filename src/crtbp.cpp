#include "crtbp_cpp/crtbp.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace crtbp {
namespace {

State6 to_state6(const Vector& values) {
    if (values.size() != 6) {
        throw std::invalid_argument("expected a 6-vector");
    }
    State6 state{};
    std::copy(values.begin(), values.end(), state.begin());
    return state;
}

Vector to_vector(const State6& state) {
    return Vector(state.begin(), state.end());
}

double cube(double value) {
    return value * value * value;
}

double crtbp_collinear_equation(double x, double mu) {
    const double l = 1.0 - mu;
    const double dx1 = x + mu;
    const double dx2 = x - l;
    const double r1 = std::abs(dx1);
    const double r2 = std::abs(dx2);
    return x - (1.0 - mu) * dx1 / cube(r1) - mu * dx2 / cube(r2);
}

double bisect_root(double lo, double hi, double mu) {
    double flo = crtbp_collinear_equation(lo, mu);
    double fhi = crtbp_collinear_equation(hi, mu);
    if (!std::isfinite(flo) || !std::isfinite(fhi) || flo * fhi > 0.0) {
        throw std::runtime_error("invalid collinear root bracket");
    }

    for (int i = 0; i < 120; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = crtbp_collinear_equation(mid, mu);
        if (fmid == 0.0 || std::abs(hi - lo) < 1e-15) {
            return mid;
        }
        if (flo * fmid <= 0.0) {
            hi = mid;
            fhi = fmid;
        } else {
            lo = mid;
            flo = fmid;
        }
        (void)fhi;
    }
    return 0.5 * (lo + hi);
}

double find_root_by_scan(double lo, double hi, double mu, int samples) {
    double previous_x = lo;
    double previous_f = crtbp_collinear_equation(previous_x, mu);
    for (int i = 1; i <= samples; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(samples);
        const double x = lo + alpha * (hi - lo);
        const double fx = crtbp_collinear_equation(x, mu);
        if (std::isfinite(previous_f) && std::isfinite(fx) && previous_f * fx <= 0.0) {
            return bisect_root(previous_x, x, mu);
        }
        previous_x = x;
        previous_f = fx;
    }
    throw std::runtime_error("failed to bracket collinear CRTBP root");
}

std::array<std::array<double, 3>, 3> crtbp_position_jacobian(
    const std::array<double, 3>& q,
    double mu) {
    std::array<std::array<double, 3>, 3> g{};
    g[0][0] = 1.0;
    g[1][1] = 1.0;

    const std::array<std::array<double, 3>, 2> centers{{
        {{-mu, 0.0, 0.0}},
        {{1.0 - mu, 0.0, 0.0}},
    }};
    const std::array<double, 2> masses{{1.0 - mu, mu}};

    for (std::size_t body = 0; body < centers.size(); ++body) {
        std::array<double, 3> d{};
        for (std::size_t i = 0; i < 3; ++i) {
            d[i] = q[i] - centers[body][i];
        }
        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        const double r = std::sqrt(r2);
        const double r3 = r2 * r;
        const double r5 = r3 * r2;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double delta = (i == j) ? 1.0 : 0.0;
                g[i][j] += masses[body] * (3.0 * d[i] * d[j] / r5 - delta / r3);
            }
        }
    }

    return g;
}

}  // namespace

Vector crtbp_vector_field(double, const Vector& y, double G, double mu) {
    if (y.size() != 6) {
        throw std::invalid_argument("CRTBP vector field expects a 6-vector");
    }

    const double r1 = std::sqrt((mu + y[0]) * (mu + y[0]) + y[1] * y[1] + y[2] * y[2]);
    const double r2 =
        std::sqrt((1.0 - mu - y[0]) * (1.0 - mu - y[0]) + y[1] * y[1] + y[2] * y[2]);
    const double m1 = 1.0 - mu;
    const double m2 = mu;

    return {
        y[3],
        y[4],
        y[5],
        y[0] + 2.0 * y[4] + G * m1 * (-mu - y[0]) / cube(r1)
            + G * m2 * (1.0 - mu - y[0]) / cube(r2),
        y[1] - 2.0 * y[3] - G * m1 * y[1] / cube(r1) - G * m2 * y[1] / cube(r2),
        -G * m1 * y[2] / cube(r1) - G * m2 * y[2] / cube(r2),
    };
}

Vector backward_crtbp_vector_field(double t, const Vector& y, double G, double mu) {
    Vector forward = crtbp_vector_field(t, y, G, mu);
    for (double& value : forward) {
        value = -value;
    }
    return forward;
}

Vector crtbp_state_transition_field(double, const Vector& y, double mu) {
    if (y.size() != 42) {
        throw std::invalid_argument("CRTBP state-transition field expects a 42-vector");
    }

    std::array<double, 3> q{{y[36], y[37], y[38]}};
    auto lower_left = crtbp_position_jacobian(q, mu);

    Matrix6 df{};
    for (int i = 0; i < 3; ++i) {
        df[i][i + 3] = 1.0;
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            df[i + 3][j] = lower_left[i][j];
        }
    }
    df[3][4] = 2.0;
    df[4][3] = -2.0;

    Matrix6 a{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            a[i][j] = y[6 * i + j];
        }
    }

    Vector out(42, 0.0);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double value = 0.0;
            for (int k = 0; k < 6; ++k) {
                value += df[i][k] * a[k][j];
            }
            out[6 * i + j] = value;
        }
    }

    Vector state(y.begin() + 36, y.end());
    Vector field = crtbp_vector_field(0.0, state, 1.0, mu);
    std::copy(field.begin(), field.end(), out.begin() + 36);
    return out;
}

double jacobi_constant(const std::array<double, 3>& position,
                       const std::array<double, 3>& velocity,
                       double mu) {
    const double r1 = std::sqrt((mu + position[0]) * (mu + position[0])
                                + position[1] * position[1]
                                + position[2] * position[2]);
    const double r2 = std::sqrt((position[0] - (1.0 - mu)) * (position[0] - (1.0 - mu))
                                + position[1] * position[1]
                                + position[2] * position[2]);
    const double v2 = velocity[0] * velocity[0] + velocity[1] * velocity[1]
        + velocity[2] * velocity[2];
    return -0.5 * v2
        + 2.0
            * (0.5 * (position[0] * position[0] + position[1] * position[1])
               + (1.0 - mu) / r1 + mu / r2);
}

double velocity_magnitude(const std::array<double, 3>& position, double C, double mu) {
    const double r1 = std::sqrt((mu + position[0]) * (mu + position[0])
                                + position[1] * position[1]
                                + position[2] * position[2]);
    const double r2 = std::sqrt((1.0 - mu - position[0]) * (1.0 - mu - position[0])
                                + position[1] * position[1]
                                + position[2] * position[2]);
    const double potential = 0.5 * (position[0] * position[0] + position[1] * position[1])
        + (1.0 - mu) / r1 + mu / r2;
    return std::sqrt(2.0 * potential - C);
}

LibrationPoints libration_points(double mu) {
    if (!(mu > 0.0 && mu < 0.5)) {
        throw std::invalid_argument("mu must satisfy 0 < mu < 0.5");
    }

    const double l = 1.0 - mu;
    const double eps = 1e-10;
    const double l1 = find_root_by_scan(-mu + eps, l - eps, mu, 10000);

    double right = 2.0;
    while (crtbp_collinear_equation(right, mu) < 0.0) {
        right *= 2.0;
    }
    const double l2 = find_root_by_scan(l + eps, right, mu, 10000);

    double left = -2.0;
    while (crtbp_collinear_equation(left, mu) > 0.0) {
        left *= 2.0;
    }
    const double l3 = find_root_by_scan(left, -mu - eps, mu, 10000);

    LibrationPoints points;
    points.L1 = {{l1, 0.0, 0.0}};
    points.L2 = {{l2, 0.0, 0.0}};
    points.L3 = {{l3, 0.0, 0.0}};
    points.L4 = {{0.5 - mu, std::sqrt(3.0) / 2.0, 0.0}};
    points.L5 = {{0.5 - mu, -std::sqrt(3.0) / 2.0, 0.0}};
    return points;
}

Matrix6 state_transition_matrix(double t0,
                                double tf,
                                const State6& state,
                                double mu,
                                const OdeOptions& options) {
    Vector initial(42, 0.0);
    for (int i = 0; i < 6; ++i) {
        initial[6 * i + i] = 1.0;
    }
    std::copy(state.begin(), state.end(), initial.begin() + 36);

    OdeResult result = integrate(
        [mu](double t, const Vector& y) { return crtbp_state_transition_field(t, y, mu); },
        t0,
        initial,
        tf,
        options);
    if (!result.success) {
        throw std::runtime_error(result.message);
    }

    Matrix6 a{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            a[i][j] = result.state[6 * i + j];
        }
    }
    return a;
}

State6 poincare_newton_crossing(double t0,
                                const State6& initial_point,
                                double t1,
                                const State6& final_point,
                                double tolerance,
                                double G,
                                double mu,
                                const OdeOptions& options) {
    constexpr int max_iterations = 16;
    constexpr double max_time = 4.0;
    State6 newton_failed{{-mu, 0.0, 0.0, 0.0, 0.0, 0.0}};
    double t_n = t1 - t0;
    State6 crossing = final_point;

    auto integrate_endpoint = [&](const State6& start, double dt, bool backward) -> State6 {
        const Vector initial = to_vector(start);
        OdeResult result = integrate(
            [G, mu, backward](double t, const Vector& y) {
                return backward ? backward_crtbp_vector_field(t, y, G, mu)
                                : crtbp_vector_field(t, y, G, mu);
            },
            0.0,
            initial,
            dt,
            options);
        if (!result.success) {
            return newton_failed;
        }
        return to_state6(result.state);
    };

    if (std::abs(final_point[1]) <= std::abs(initial_point[1])) {
        crossing = final_point;
        double y_n = final_point[1];
        t_n = t1 - t0;
        int iterations = 0;

        while (std::abs(y_n) >= tolerance) {
            if (iterations >= max_iterations || t_n >= max_time) {
                return newton_failed;
            }
            crossing = integrate_endpoint(initial_point, t_n, false);
            const double dfx = crossing[4];
            if (std::abs(dfx) < std::numeric_limits<double>::epsilon()) {
                return newton_failed;
            }
            t_n -= crossing[1] / dfx;
            crossing = integrate_endpoint(initial_point, t_n, false);
            y_n = crossing[1];
            ++iterations;
        }
        return crossing;
    }

    crossing = initial_point;
    double y_n = initial_point[1];
    t_n = t1 - t0;
    int iterations = 0;

    while (std::abs(y_n) >= tolerance) {
        if (iterations >= max_iterations || t_n >= max_time) {
            return newton_failed;
        }
        crossing = integrate_endpoint(final_point, t_n, true);
        const double dfx = -crossing[4];
        if (std::abs(dfx) < std::numeric_limits<double>::epsilon()) {
            return newton_failed;
        }
        t_n -= crossing[1] / dfx;
        crossing = integrate_endpoint(final_point, t_n, true);
        y_n = crossing[1];
        ++iterations;
    }
    return crossing;
}

CrossingData interpolate_crossing_data(const std::array<double, 4>& p1,
                                       const std::array<double, 4>& p2,
                                       double C,
                                       double mu) {
    const double x1 = p1[1];
    const double x2 = p2[1];
    const double yx1 = p1[0];
    const double s1 = -p1[3];
    const double s2 = -p2[3];

    const double a = (s2 - s1) / (x2 - x1);
    const double b = s1 - 2.0 * a * x1;
    const double c = yx1 - a * x1 * x1 - b * x1;

    const double theta = std::atan(-1.0 / b);
    const double speed = velocity_magnitude({{c, 0.0, 0.0}}, C, mu);

    return {c, speed * std::cos(theta), speed * std::sin(theta)};
}

}  // namespace crtbp
