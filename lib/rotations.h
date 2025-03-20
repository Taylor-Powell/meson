#ifndef __rotations_h__
#define __rotations_h__

#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <numbers>
#include "generic_funcs.h"

namespace {
    typedef std::complex<double> cd;

    template <typename T>
    using vec2D = std::vector<std::vector<T>>;

    // Defined for convenience...
    double pi = std::numbers::pi;
    double pi2 = std::numbers::pi / 2.0;
    double pi4 = std::numbers::pi / 4.0;
    double twopi = 2.0 * std::numbers::pi;
}
namespace rotations {
    // Lists of symmetries, initial angles, allowed momenta for each symmetry, and chosen angles for each momentum.
    const std::vector<std::string> symList = {"OhD", "Dic4", "Dic4", "Dic2", "Dic3", "C40mn", "C4nnm"};
    const std::vector<std::vector<double>> symInit_angles = {
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, pi2, 0}, {0, pi4, 0.955316618124509}, {0, pi2, 0.4636476090008061}, {0, -3.0 * pi4, 0.6154797086703874}
    };
    const std::vector<std::vector<std::vector<int>>> sym_moms = {
        {{0, 0, 0}}, 
        {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, 
        {{0, 0, 2}, {2, 0, 0}, {0, 2, 0}, {-2, 0, 0}, {0, -2, 0}, {0, 0, -2}},    
        {{0, 1, 1}, {1, 1, 0}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}},
        {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}},
        {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}},
        {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}}
    };
    const std::vector<std::vector<std::vector<double>>> sym_angles = {
        {{0, 0, 0}},
        {{0, 0, 0}, {0, pi2, 0}, {pi2, pi2, 0}, {pi, pi2, 0}, {-pi2, pi2, 0}, {0, pi, 0}},
        {{0, 0, 0}, {0, pi2, 0}, {pi2, pi2, 0}, {pi, pi2, 0}, {-pi2, pi2, 0}, {0, pi, 0}},
        {{0, 0, 0}, {0, pi2, 0}, {twopi, 0, 3.0 * pi2}, {-pi2, pi2, 0}, {0, -pi, 0}, {twopi, 0, pi2}, {0, -pi2, 0}, {twopi, 0, pi}, {twopi, pi2, 3.0 * pi2}, {7.0 * pi2, pi2, pi}, {twopi, -pi2, pi2}, {twopi, -pi2, pi2}},
        {{0, 0, 0}, {0, 0, pi2}, {0, 0, -pi2}, {-pi2, -pi2, pi2}, {0, 0, -pi}, {0, -pi, -pi}, {0, -pi, 0}, {0, -pi, -pi2}},
        {{0, 0, 0}, {3.0 * pi2, -pi2, 0}, {0, pi2, pi2}, {3.0 * pi2, -pi2, 3.0 * pi2}, {0, pi2, 0}, {0, 0, 3.0 * pi2}, {0, -pi, 0}, {3.0 * pi2, pi2, 0}, {0, -pi2, 3.0 * pi2}, {3.0 * pi2, pi2, pi2}, {0, -pi2, 0}, {0, -pi, pi2}, {0, 0, pi}, {3.0 * pi2, -pi2, pi}, {0, pi2, 3.0 * pi2}, {3.0 * pi2, -pi2, pi2}, {0, pi2, pi}, {0, 0, pi2}, {0, -pi, pi}, {3.0 * pi2, pi2, pi}, {0, -pi2, pi2}, {3.0 * pi2, pi2, 3.0 * pi2}, {0, -pi2, pi}, {0, -pi, 3.0 * pi2}},
        {{0, 0, 0}, {3.0 * pi2, -pi2, 0}, {0, pi2, pi2}, {0, 0, pi2}, {3.0 * pi2, -pi2, 3.0 * pi2}, {0, 0, 3.0 * pi2}, {3.0 * pi2, -pi2, pi2}, {0, pi2, pi}, {0, pi2, 0}, {0, -pi, pi2}, {3.0 * pi2, pi2, pi2}, {0, -pi2, 0}, {0, 0, pi}, {3.0 * pi2, -pi2, pi}, {0, pi2, 3.0 * pi2}, {0, -pi, 0}, {3.0 * pi2, pi2, pi}, {0, -pi, pi}, {3.0 * pi2, pi2, 0}, {0, -pi2, 3.0 * pi2}, {0, -pi2, pi2}, {0, -pi, 3.0 * pi2}, {3.0 * pi2, pi2, 3.0 * pi2}, {0, -pi2, pi}}
    };


    /////////////////// Forward declarations ///////////////////

    // Get the rotation angles for a given symmetry and momentum
    std::vector<double> getRotAngles(const std::string sym, const std::vector<int> mom3_i);

    // Initialize the polarization vector for a given symmetry from z-axis
    void rotPolVec_init(std::vector<cd>& polVec, const std::string sym);

    // Rotate the polarization vector from angle given by rotPolVec_init for a given momentum
    void rotPolVec(std::vector<cd>& polVec, const std::string sym, const std::vector<int> mom3_i);

    // Get the polarization vector for a given spin along the z-axis
    std::vector<cd> getPolz4(double E, double mom_sq, int Jz, bool current);

    // Get the polarization vector for a given Jz, rotated from the z-axis
    std::vector<cd> getPol4_Jz(double E, double mom_sq, const std::vector<int> mom3_i, int Jz, const std::string sym, bool current);

    // Get the polarization vector for a given helicity, rotated from the z-axis
    std::vector<cd> getPol4_hel(double E, double mom_sq, std::vector<int> mom3_i, int helicity, const std::string sym, bool current);
}
#endif