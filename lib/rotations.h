#ifndef __rotations_h__
#define __rotations_h__

#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "generic_funcs.h"

namespace {
    typedef std::complex<double> cd;
    cd j1 = cd(0,1);

    template <typename T>
    using vec2D = std::vector<std::vector<T>>;

    // Defined for convenience...
    double pi = std::numbers::pi;
    double pi2 = std::numbers::pi / 2.0;
    double pi4 = std::numbers::pi / 4.0;
    double twopi = 2.0 * std::numbers::pi;
}
namespace rotations {
    const std::vector<std::string> symList = {"Dic4", "Dic2", "Dic3", "C40mn", "C4nnm"};
    const std::vector<std::vector<double>> symInit_angles = {
        {0, 0, 0}, {0, pi2, 0}, {0, pi4, 0.955316618124509}, {0, pi2, 0.4636476090008061}, {0, -3.0 * pi4, 0.6154797086703874}
    };
    const std::vector<std::vector<std::vector<int>>> sym_moms = {
        {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}},    
        {{0, 1, 1}, {1, 1, 0}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}},
        {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}},
        {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}},
        {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}}
    };
    const std::vector<std::vector<std::vector<double>>> sym_angles = {
        {{0, 0, 0}, {0, pi2, 0}, {pi2, pi2, 0}, {pi, pi2, 0}, {-pi2, pi2, 0}, {0, pi, 0}},
        {{0, 0, 0}, {0, pi2, 0}, {twopi, 0, 3.0 * pi2}, {-pi2, pi2, 0}, {0, -pi, 0}, {twopi, 0, pi2}, {0, -pi2, 0}, {twopi, 0, pi}, {twopi, pi2, 3.0 * pi2}, {7.0 * pi2, pi2, pi}, {twopi, -pi2, pi2}, {twopi, -pi2, pi2}},
        {{0, 0, 0}, {0, 0, pi2}, {0, 0, -pi2}, {-pi2, -pi2, pi2}, {0, 0, -pi}, {0, -pi, -pi}, {0, -pi, 0}, {0, -pi, -pi2}},
        {{0, 0, 0}, {3.0 * pi2, -pi2, 0}, {0, pi2, pi2}, {3.0 * pi2, -pi2, 3.0 * pi2}, {0, pi2, 0}, {0, 0, 3.0 * pi2}, {0, -pi, 0}, {3.0 * pi2, pi2, 0}, {0, -pi2, 3.0 * pi2}, {3.0 * pi2, pi2, pi2}, {0, -pi2, 0}, {0, -pi, pi2}, {0, 0, pi}, {3.0 * pi2, -pi2, pi}, {0, pi2, 3.0 * pi2}, {3.0 * pi2, -pi2, pi2}, {0, pi2, pi}, {0, 0, pi2}, {0, -pi, pi}, {3.0 * pi2, pi2, pi}, {0, -pi2, pi2}, {3.0 * pi2, pi2, 3.0 * pi2}, {0, -pi2, pi}, {0, -pi, 3.0 * pi2}},
        {{0, 0, 0}, {3.0 * pi2, -pi2, 0}, {0, pi2, pi2}, {0, 0, pi2}, {3.0 * pi2, -pi2, 3.0 * pi2}, {0, 0, 3.0 * pi2}, {3.0 * pi2, -pi2, pi2}, {0, pi2, pi}, {0, pi2, 0}, {0, -pi, pi2}, {3.0 * pi2, pi2, pi2}, {0, -pi2, 0}, {0, 0, pi}, {3.0 * pi2, -pi2, pi}, {0, pi2, 3.0 * pi2}, {0, -pi, 0}, {3.0 * pi2, pi2, pi}, {0, -pi, pi}, {3.0 * pi2, pi2, 0}, {0, -pi2, 3.0 * pi2}, {0, -pi2, pi2}, {0, -pi, 3.0 * pi2}, {3.0 * pi2, pi2, 3.0 * pi2}, {0, -pi2, pi}}
    };

    void rotPolVec_init(vec2D<cd>& polVec, std::string sym) {
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                polVec = basics::rotVec(polVec, symInit_angles[i][0], symInit_angles[i][1], symInit_angles[i][2]);
                return;
            }
        }
    }

    void rotPolVec(vec2D<cd>& polVec, std::string sym, std::vector<int> mom) {
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                for (int j = 0; j < sym_moms[i].size(); j++) {
                    if (mom == sym_moms[i][j]) {
                        polVec = basics::rotVec(polVec, sym_angles[i][j][0], sym_angles[i][j][1], sym_angles[i][j][2]);
                        return;
                    }
                }
                throw std::string("Momentum " + std::to_string(mom[0]) + std::to_string(mom[1]) + std::to_string(mom[2]) + " for symmetry " + sym + " not recognized in rotations::rotPolVec.\n");
            }
            throw std::string("Symmetry " + sym + " not recognized in rotations::rotPolVec.\n");
        }
    }
}
#endif