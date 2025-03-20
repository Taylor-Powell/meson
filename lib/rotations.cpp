#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "generic_funcs.h"
#include "rotations.h"

namespace rotations {
    std::vector<int> getMom3_i(const std::string momStr) {
        std::vector<int> mom3_i;
        for (int i = 0; i < 3; i++) {
            mom3_i.push_back(momStr[i] - '0');
        }
        return mom3_i;
    }

    std::string getSym(const std::string momType) {
        if (momType == "00n") return "Dic4";
        else if (momType == "0nn") return "Dic2";
        else if (momType == "nnn") return "Dic3";
        else if (momType == "0mn") return "C40mn";
        else if (momType == "nnm") return "C4nnm";
        else if (momType == "000") return "OhD";
        else throw std::string("Momentum " + momType + " not recognized in rotations::getSym().\n");
    }

    std::string getMomType(const std::vector<int> mom3_i) {
        std::string sym;
        int mom_sq = basics::dot(mom3_i, mom3_i);
        if (mom_sq == 0) return "000";
        else if (mom_sq == 1) return "00n";
        else if (mom_sq == 2) return "0nn";
        else if (mom_sq == 3) return "nnn";
        else if (mom_sq == 4) return "00n"; 
        else if (mom_sq == 5) return "0mn";
        else if (mom_sq == 6) return "nnm";
        else throw std::string("Momentum " + std::to_string(mom3_i[0]) + std::to_string(mom3_i[1]) + std::to_string(mom3_i[2]) + " not recognized in rotations::getSym.\n");


    }
    
    std::vector<double> getRotAngles(const std::string sym, const std::vector<int> mom3_i) {
        if (sym == "OhD") return {0, 0, 0};
        for (int i = 0; i < symList.size(); i++) {
            if (symList[i] == sym) {
                for (int j = 0; j < sym_moms[i].size(); j++) {
                    if (mom3_i == sym_moms[i][j]) {
                        return sym_angles[i][j];
                    }
                }
            }          
        }
        throw std::string("Symmetry " + sym + " not recognized in rotations::getRotAngles.\n");
        return {0, 0, 0};
    }
    
    void rotPolVec_init(std::vector<cd>& polVec, const std::string sym) {
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
        if (sym == "OhD") return;
        else throw std::string("Symmetry " + sym + " not recognized in rotations::rotPolVec_init.\n");
    }

    void rotPolVec(std::vector<cd>& polVec, const std::string sym, const std::vector<int> mom3_i) {
        //std::cout << "rotPolVec called with sym = " << sym << ", mom = " << mom[0] << mom[1] << mom[2] << std::endl;
        std::vector<double> angles = getRotAngles(sym, mom3_i);
        std::vector<cd> polVec_Spatial(polVec.begin() + 1, polVec.end());
        polVec_Spatial = basics::rotVec(polVec_Spatial, angles[0], angles[1], angles[2]);
        // Replace last three components of polVec with polVec_Spatial
        for (int j = 0; j < 3; j++) {
            polVec[j+1] = polVec_Spatial[j];
        }
    }    

    std::vector<cd> getPolz4(double E, double mom3_sq, int helicity, bool current) {
        std::vector<cd> polz(4, cd(0,0));
        if (helicity == 0) {
            if (current) polz = {cd(0,0), cd(0,0), cd(0,0), cd(1,0)};
            else {
                double m = std::sqrt(std::pow(E, 2) - mom3_sq);
                double p3 = std::sqrt(mom3_sq);
                polz = {cd(p3/m, 0), cd(0, 0), cd(0, 0), cd(E/m, 0)};
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

    std::vector<cd> getPol4_hel(double E, double mom3_sq, const std::vector<int> mom3_i, int helicity, const std::string sym, bool current) {
        std::vector<cd> polz = getPolz4(E, mom3_sq, helicity, current);
        rotPolVec_init(polz, sym);
        rotPolVec(polz, sym, mom3_i);
        return polz;
    }

    std::vector<cd> getPol4_Jz(double E, double mom3_sq, const std::vector<int> mom3_i, int Jz, const std::string sym, bool current) {
        std::vector<cd> polz = getPolz4(E, mom3_sq, Jz, current);
        rotPolVec_init(polz, sym);
        rotPolVec(polz, sym, mom3_i);
        return polz;
    }
}