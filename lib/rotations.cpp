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
        else if (momstr == "000") return "full";
        else throw std::string("Momentum " + momstr + " not recognized in rotations::getSym().\n");
    }
    
    std::vector<double> getRotAngles(std::string sym, std::vector<int> mom) {
        if (sym == "full") return {0, 0, 0};
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
        throw std::string("Symmetry " + sym + " not recognized in rotations::getRotAngles.\n");
        return {0, 0, 0};
    }
    
    void rotPolVec_init(std::vector<cd>& polVec, std::string sym) {
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                std::vector<cd> polVec_Spatial(polVec.begin() + 1, polVec.end());
                polVec_Spatial = basics::rotVec(polVec_Spatial, symInit_angles[i][0], symInit_angles[i][1], symInit_angles[i][2]);
                // Replace last three components of polVec with polVec_Spatial
                for (int j = 0; j < 3; j++) {
                    polVec[j+1] = polVec_Spatial[j];
                }
                return;
            }
        }
    }

    void rotPolVec(std::vector<cd>& polVec, std::string sym, std::vector<int> mom) {
        std::vector<double> angles = getRotAngles(sym, mom);
        std::vector<cd> polVec_Spatial(polVec.begin() + 1, polVec.end());
        polVec_Spatial = basics::rotVec(polVec_Spatial, angles[0], angles[1], angles[2]);
        // Replace last three components of polVec with polVec_Spatial
        for (int j = 0; j < 3; j++) {
            polVec[j+1] = polVec_Spatial[j];
        }
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