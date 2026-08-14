#include "crtbp_cpp/ode.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace crtbp {
namespace {

struct StepResult {
    Vector y5;
    Vector y4;
    Vector f0;
    Vector f1;
};

void require_same_size(const Vector& a, const Vector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("ODE vectors have inconsistent sizes");
    }
}

Vector add_scaled(const Vector& y,
                  double h,
                  const std::vector<const Vector*>& ks,
                  const std::vector<double>& coeffs) {
    Vector out = y;
    for (std::size_t i = 0; i < ks.size(); ++i) {
        const Vector& k = *ks[i];
        require_same_size(out, k);
        for (std::size_t j = 0; j < out.size(); ++j) {
            out[j] += h * coeffs[i] * k[j];
        }
    }
    return out;
}

StepResult dormand_prince_step(const OdeFunction& f,
                               double t,
                               const Vector& y,
                               double h) {
    const Vector k1 = f(t, y);
    const Vector y2 = add_scaled(y, h, {&k1}, {1.0 / 5.0});
    const Vector k2 = f(t + h * 1.0 / 5.0, y2);

    const Vector y3 = add_scaled(y, h, {&k1, &k2}, {3.0 / 40.0, 9.0 / 40.0});
    const Vector k3 = f(t + h * 3.0 / 10.0, y3);

    const Vector y4 = add_scaled(y, h, {&k1, &k2, &k3},
                                 {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0});
    const Vector k4 = f(t + h * 4.0 / 5.0, y4);

    const Vector y5stage = add_scaled(
        y,
        h,
        {&k1, &k2, &k3, &k4},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0});
    const Vector k5 = f(t + h * 8.0 / 9.0, y5stage);

    const Vector y6 = add_scaled(
        y,
        h,
        {&k1, &k2, &k3, &k4, &k5},
        {9017.0 / 3168.0,
         -355.0 / 33.0,
         46732.0 / 5247.0,
         49.0 / 176.0,
         -5103.0 / 18656.0});
    const Vector k6 = f(t + h, y6);

    Vector fifth = add_scaled(y,
                              h,
                              {&k1, &k3, &k4, &k5, &k6},
                              {35.0 / 384.0,
                               500.0 / 1113.0,
                               125.0 / 192.0,
                               -2187.0 / 6784.0,
                               11.0 / 84.0});
    const Vector k7 = f(t + h, fifth);

    Vector fourth = add_scaled(y,
                               h,
                               {&k1, &k3, &k4, &k5, &k6, &k7},
                               {5179.0 / 57600.0,
                                7571.0 / 16695.0,
                                393.0 / 640.0,
                                -92097.0 / 339200.0,
                                187.0 / 2100.0,
                                1.0 / 40.0});

    return {fifth, fourth, k1, k7};
}

double error_norm(const Vector& y,
                  const Vector& y_next,
                  const Vector& error,
                  const OdeOptions& options) {
    require_same_size(y, y_next);
    require_same_size(y, error);

    double sum = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i) {
        const double scale = options.abs_tol
            + options.rel_tol * std::max(std::abs(y[i]), std::abs(y_next[i]));
        const double normalized = error[i] / scale;
        sum += normalized * normalized;
    }

    return std::sqrt(sum / static_cast<double>(std::max<std::size_t>(1, y.size())));
}

Vector hermite_state(const Vector& y0,
                     const Vector& f0,
                     const Vector& y1,
                     const Vector& f1,
                     double h,
                     double alpha) {
    require_same_size(y0, y1);
    require_same_size(y0, f0);
    require_same_size(y0, f1);

    const double a2 = alpha * alpha;
    const double a3 = a2 * alpha;
    const double h00 = 2.0 * a3 - 3.0 * a2 + 1.0;
    const double h10 = a3 - 2.0 * a2 + alpha;
    const double h01 = -2.0 * a3 + 3.0 * a2;
    const double h11 = a3 - a2;

    Vector out(y0.size(), 0.0);
    for (std::size_t i = 0; i < y0.size(); ++i) {
        out[i] = h00 * y0[i] + h10 * h * f0[i] + h01 * y1[i] + h11 * h * f1[i];
    }
    return out;
}

bool crossed(double g0, double g1, int direction) {
    if (direction > 0) {
        return g0 < 0.0 && g1 >= 0.0;
    }
    if (direction < 0) {
        return g0 > 0.0 && g1 <= 0.0;
    }
    return (g0 <= 0.0 && g1 >= 0.0) || (g0 >= 0.0 && g1 <= 0.0);
}

OdeEvent refine_event(const EventFunction& event,
                      double t0,
                      const Vector& y0,
                      const Vector& f0,
                      double t1,
                      const Vector& y1,
                      const Vector& f1,
                      double g0,
                      int direction) {
    const double h = t1 - t0;
    double lo = 0.0;
    double hi = 1.0;
    double glo = g0;

    for (int iter = 0; iter < 70; ++iter) {
        const double mid = 0.5 * (lo + hi);
        Vector y_mid = hermite_state(y0, f0, y1, f1, h, mid);
        const double gmid = event(t0 + mid * h, y_mid);

        if (crossed(glo, gmid, direction)) {
            hi = mid;
        } else {
            lo = mid;
            glo = gmid;
        }

        if (std::abs(hi - lo) * std::abs(h) < 1e-13) {
            break;
        }
    }

    const double alpha = 0.5 * (lo + hi);
    return {t0 + alpha * h, hermite_state(y0, f0, y1, f1, h, alpha)};
}

OdeResult integrate_impl(const OdeFunction& f,
                         double t0,
                         const Vector& y0,
                         double t1,
                         const OdeOptions& options,
                         const EventFunction* event,
                         int event_direction) {
    OdeResult result;
    result.t = t0;
    result.state = y0;

    if (y0.empty()) {
        result.success = false;
        result.message = "empty initial state";
        return result;
    }
    if (t0 == t1) {
        return result;
    }

    const double direction = (t1 > t0) ? 1.0 : -1.0;
    double h = std::abs(options.initial_step) > 0.0 ? std::abs(options.initial_step)
                                                    : std::abs(t1 - t0) / 100.0;
    h = std::min(h, std::abs(options.max_step));
    h = std::max(h, std::abs(options.min_step));
    h *= direction;

    double previous_event_t = std::numeric_limits<double>::quiet_NaN();
    double g_previous = event ? (*event)(result.t, result.state) : 0.0;

    while ((t1 - result.t) * direction > 0.0) {
        if (result.accepted_steps + result.rejected_steps > options.max_steps) {
            result.success = false;
            result.message = "maximum ODE step count exceeded";
            return result;
        }

        if (std::abs(h) > std::abs(t1 - result.t)) {
            h = t1 - result.t;
        }
        if (std::abs(h) < std::abs(options.min_step)) {
            result.success = false;
            result.message = "ODE step size fell below min_step";
            return result;
        }

        StepResult step;
        try {
            step = dormand_prince_step(f, result.t, result.state, h);
        } catch (const std::exception& ex) {
            result.success = false;
            result.message = ex.what();
            return result;
        }

        if (!all_finite(step.y5) || !all_finite(step.y4)) {
            result.success = false;
            result.message = "non-finite ODE state";
            return result;
        }

        Vector local_error(step.y5.size(), 0.0);
        for (std::size_t i = 0; i < local_error.size(); ++i) {
            local_error[i] = step.y5[i] - step.y4[i];
        }

        const double err = error_norm(result.state, step.y5, local_error, options);
        const bool accept = err <= 1.0 || std::abs(h) <= std::abs(options.min_step);

        if (accept) {
            const double old_t = result.t;
            const Vector old_state = result.state;
            const double old_g = g_previous;

            result.t += h;
            result.state = step.y5;
            ++result.accepted_steps;

            if (event) {
                const double new_g = (*event)(result.t, result.state);
                if (crossed(old_g, new_g, event_direction)) {
                    OdeEvent ode_event = refine_event(*event,
                                                      old_t,
                                                      old_state,
                                                      step.f0,
                                                      result.t,
                                                      result.state,
                                                      step.f1,
                                                      old_g,
                                                      event_direction);
                    if (!std::isfinite(previous_event_t)
                        || std::abs(ode_event.t - previous_event_t)
                            > options.event_time_epsilon) {
                        result.events.push_back(std::move(ode_event));
                        previous_event_t = result.events.back().t;
                    }
                }
                g_previous = new_g;
            }
        } else {
            ++result.rejected_steps;
        }

        const double safety = 0.9;
        const double min_factor = 0.1;
        const double max_factor = 5.0;
        double factor = max_factor;
        if (err > 0.0 && std::isfinite(err)) {
            factor = safety * std::pow(err, -0.2);
            factor = std::min(max_factor, std::max(min_factor, factor));
        }
        if (!accept) {
            factor = std::min(1.0, factor);
        }

        h *= factor;
        if (std::abs(h) > std::abs(options.max_step)) {
            h = direction * std::abs(options.max_step);
        }
    }

    return result;
}

}  // namespace

OdeResult integrate(const OdeFunction& f,
                    double t0,
                    const Vector& y0,
                    double t1,
                    const OdeOptions& options) {
    return integrate_impl(f, t0, y0, t1, options, nullptr, 0);
}

OdeResult integrate_with_events(const OdeFunction& f,
                                double t0,
                                const Vector& y0,
                                double t1,
                                const OdeOptions& options,
                                const EventFunction& event,
                                int direction) {
    return integrate_impl(f, t0, y0, t1, options, &event, direction);
}

bool all_finite(const Vector& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace crtbp
