#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace crtbp {

using Vector = std::vector<double>;
using OdeFunction = std::function<Vector(double, const Vector&)>;
using EventFunction = std::function<double(double, const Vector&)>;

struct OdeOptions {
    double rel_tol = 1e-10;
    double abs_tol = 1e-12;
    double initial_step = 1e-3;
    double max_step = 0.05;
    double min_step = 1e-13;
    std::size_t max_steps = 10000000;
    double event_time_epsilon = 1e-11;
};

struct OdeEvent {
    double t = 0.0;
    Vector state;
};

struct OdeResult {
    double t = 0.0;
    Vector state;
    std::vector<OdeEvent> events;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    bool success = true;
    std::string message;
};

OdeResult integrate(const OdeFunction& f,
                    double t0,
                    const Vector& y0,
                    double t1,
                    const OdeOptions& options = {});

OdeResult integrate_with_events(const OdeFunction& f,
                                double t0,
                                const Vector& y0,
                                double t1,
                                const OdeOptions& options,
                                const EventFunction& event,
                                int direction);

bool all_finite(const Vector& values);

}  // namespace crtbp
