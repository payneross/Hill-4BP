// CRTBP Poincare section generator built on CAPD (non-rigorous "D" types)

// The CRTBP vector field is handed to CAPD as a formula string.
// CAPD's DPoincareMap integrates straight to each crossing.

// To get validated section points, swap the D* types for I* types
// (DMap->IMap, DOdeSolver->IOdeSolver, DPoincareMap->IPoincareMap, DVector->
// IVector) and feed interval initial boxes.

#include "capd/capdlib.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using capd::DMap;
using capd::DOdeSolver;
using capd::DCoordinateSection;
using capd::DPoincareMap;
using capd::DVector;

namespace {

// (G = 1, primaries at x = -mu and x = 1 - mu)
// Omega = 0.5*(x^2 + y^2) + (1-mu)/r1 + mu/r2, Jacobi C = 2*Omega - v^2
const char* kCrtbpFormula =
    "par:mu;"
    "var:x,y,z,dx,dy,dz;"
    "fun:"
    "dx,"
    "dy,"
    "dz,"
    "x+2*dy-(1-mu)*(x+mu)/sqrt((x+mu)^2+y^2+z^2)^3"
    "-mu*(x-1+mu)/sqrt((x-1+mu)^2+y^2+z^2)^3,"
    "y-2*dx-(1-mu)*y/sqrt((x+mu)^2+y^2+z^2)^3"
    "-mu*y/sqrt((x-1+mu)^2+y^2+z^2)^3,"
    "-(1-mu)*z/sqrt((x+mu)^2+y^2+z^2)^3"
    "-mu*z/sqrt((x-1+mu)^2+y^2+z^2)^3;";

// Planar in-plane speed^2 available at (x, y=0, z=0) for a given Jacobi C
double speed_squared(double x, double mu, double C) {
    const double r1 = std::sqrt((x + mu) * (x + mu));
    const double r2 = std::sqrt((x - 1.0 + mu) * (x - 1.0 + mu));
    const double Omega = 0.5 * x * x + (1.0 - mu) / r1 + mu / r2;
    return 2.0 * Omega - C;
}

}  

int main(int argc, char** argv) {
    
    double mu = 0.0121505856;  // Earth-Moon-ish; set to your system
    double C = 3.18; // Jacobi constant
    double x_begin = 0.085, x_end = 0.85;
    int x_count = 300;
    int iterates = 1000;
    std::string output = "crtbp_capd.csv";
    std::string svg_path = "crtbp_capd.svg";  // set svg= to "" to disable

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto eq = a.find('=');
        if (eq == std::string::npos) continue;
        std::string k = a.substr(0, eq), v = a.substr(eq + 1);
        if (k == "mu") mu = std::stod(v);
        else if (k == "C") C = std::stod(v);
        else if (k == "x_begin") x_begin = std::stod(v);
        else if (k == "x_end") x_end = std::stod(v);
        else if (k == "x_count") x_count = std::stoi(v);
        else if (k == "iterates") iterates = std::stoi(v);
        else if (k == "output") output = v;
        else if (k == "svg") svg_path = v;
    }

    // build the CAPD Poincare map 
    // Two *different* orders:
    //   1.  map degree: order of spatial (phase-variable) derivatives the map can
    //     produce. For a Poincare map we need none, so degree 1. Large
    //     values allocate a huge jet (C(n+d,d) coeffs) and blow up memory.
    //   2.  solver order: the time-Taylor truncation order of the integrator.
    const int map_degree = 1;
    const int solver_order = 20;
    DMap vf(kCrtbpFormula, map_degree);
    vf.setParameter("mu", mu);

    DOdeSolver solver(vf, solver_order);
    solver.setAbsoluteTolerance(1e-14);
    solver.setRelativeTolerance(1e-12);

    // Section y = 0, record dy/dt > 0.
    DCoordinateSection section(6, 1);
    DPoincareMap pm(solver, section, capd::poincare::MinusPlus);

    std::ofstream out(output);
    out.precision(15);

    const double xStep = (x_count > 1) ? (x_end - x_begin) / (x_count - 1) : 0.0;
    long written = 0, skipped = 0;
    std::vector<std::pair<double, double>> pts; 

    for (int n = 0; n < x_count; ++n) {
        const double x0 = x_begin + n * xStep;

        // Place the IC on the section with dy > 0. Here we put all in-plane
        // speed into ydot (xdot = 0)
        const double v2 = speed_squared(x0, mu, C);
        if (v2 <= 0.0) { ++skipped; continue; }
        const double ydot0 = std::sqrt(v2);

        DVector u(6);
        u[0] = x0; u[1] = 0.0; u[2] = 0.0;
        u[3] = 0.0; u[4] = ydot0; u[5] = 0.0;

        try {
            for (int it = 0; it < iterates; ++it) {
                u = pm(u);  // integrate to the next crossing
                out << u[0] << ',' << u[3] << '\n';  // (x, xdot)
                pts.emplace_back(u[0], u[3]);
                ++written;
            }
        } catch (const std::exception& e) {
            // Orbit escaped / never returns to the section: move on.
            std::cerr << "x0=" << x0 << " stopped: " << e.what() << '\n';
        }
    }

    // dependency-free SVG scatter so the section can be eyeballed 
    if (!svg_path.empty() && !pts.empty()) {
        double xmin = std::numeric_limits<double>::max(), xmax = -xmin;
        double ymin = xmin, ymax = -xmin;
        for (const auto& p : pts) {
            xmin = std::min(xmin, p.first);  xmax = std::max(xmax, p.first);
            ymin = std::min(ymin, p.second); ymax = std::max(ymax, p.second);
        }
        const double W = 1000.0, H = 800.0, pad = 40.0;
        const double sx = (xmax > xmin) ? (W - 2 * pad) / (xmax - xmin) : 1.0;
        const double sy = (ymax > ymin) ? (H - 2 * pad) / (ymax - ymin) : 1.0;

        std::ofstream svg(svg_path);
        svg << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W
            << "' height='" << H << "'>\n"
            << "<rect width='" << W << "' height='" << H << "' fill='white'/>\n";
        for (const auto& p : pts) {
            const double px = pad + (p.first - xmin) * sx;
            const double py = H - pad - (p.second - ymin) * sy; // flip y
            svg << "<circle cx='" << px << "' cy='" << py
                << "' r='0.6' fill='blue'/>\n";
        }
        svg << "<text x='" << W / 2 << "' y='" << H - 8
            << "' font-size='16' text-anchor='middle'>x</text>\n"
            << "<text x='14' y='" << H / 2
            << "' font-size='16' transform='rotate(-90 14," << H / 2
            << ")' text-anchor='middle'>xdot</text>\n</svg>\n";
        std::cerr << "svg -> " << svg_path << '\n';
    }

    std::cerr << "wrote " << written << " points (" << skipped
              << " ICs skipped), mu=" << mu << " C=" << C << " -> " << output << '\n';
    return 0;
}
