#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "generic_funcs.h"
#include "rotations.h"

namespace rotations {
    std::string getSym(std::string momstr) {
        if ((momstr == "001") || (momstr == "002")) return "Dic4";
        else if (momstr == "011") return "Dic2";
        else if (momstr == "111") return "Dic3";
        else if (momstr == "012") return "C40mn";
        else if (momstr == "112") return "C4nnm";
        else throw std::string("Momentum " + momstr + " not recognized in rotations::getSym().\n");
    }
    
    std::vector<double> getRotAngles(std::string sym, std::vector<int> mom) {
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                for (int j = 0; j < sym_moms[i].size(); j++) {
                    if (mom == sym_moms[i][j]) {
                        return sym_angles[i][j];
                    }
                }
                throw std::string("Momentum " + std::to_string(mom[0]) + std::to_string(mom[1]) + std::to_string(mom[2]) + " for symmetry " + sym + " not recognized in rotations::getSymAngles.\n");
            }            
        }
        throw std::string("Symmetry " + sym + " not recognized in rotations::getSymAngles.\n");
        return {0, 0, 0};
    }
    
    void rotPolVec_init(std::vector<cd>& polVec, std::string sym) {
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                polVec = basics::rotVec(polVec, symInit_angles[i][0], symInit_angles[i][1], symInit_angles[i][2]);
                return;
            }
        }
    }

    void rotPolVec(std::vector<cd>& polVec, std::string sym, std::vector<int> mom) {
        std::vector<double> angles = getRotAngles(sym, mom);
        polVec = basics::rotVec(polVec, angles[0], angles[1], angles[2]);
    }    

    std::vector<cd> getPolz4(double E, double mom_sq, int helicity, bool current) {
        std::vector<cd> polz(4, cd(0,0));
        if (helicity == 0) {
            if (current) polz = {cd(0,0), cd(0,0), cd(0,0), cd(1,0)};
            else {
                double m = std::sqrt(std::pow(E, 2) - mom_sq);
                double p = std::sqrt(mom_sq);
                polz = {cd(p/m, 0), cd(0, 0), cd(0, 0), cd(E/m, 0)};
            }
        }
        else if (helicity == 1) {
            double f = 1.0 / std::sqrt(2.0);
            polz = {cd(0,0), cd(-f,0), cd(0,-f), cd(0,0)};
        }
        else if (helicity == -1) {
            double f = 1.0 / std::sqrt(2.0);
            polz = {cd(0,0), cd(f,0), cd(0,-f), cd(0,0)};
        }
        else throw std::string("Invalid helicity in basics::getPolz4().");
        return polz;
    }

    std::vector<cd> getPol4(double E, double mom_sq, std::vector<int> mom3, int helicity, std::string sym, bool current) {
        std::vector<cd> polz = getPolz4(E, mom_sq, helicity, current);
        rotPolVec_init(polz, sym);
        rotPolVec(polz, sym, mom3);
        return polz;
    }
}