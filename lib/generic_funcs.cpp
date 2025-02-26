#include <cmath>
#include <vector>
#include "generic_funcs.h"

namespace basics {

    /* Set of momentum lists with allowed cubic rotations and their dims
     * Includes {000,001,011,111,002,012,112} */
    vec2D<int> momList_001 {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    vec2D<int> momList_002 {{2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2}};
    vec2D<int> momList_011 {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}};
    vec2D<int> momList_111 {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};
    vec2D<int> momList_012 {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}};
    vec2D<int> momList_112 {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}};
    
    /* Given a momentum string in ascending order, returns the set of 
     * allowed permutations as a 2D vector of integers */
    vec2D<int> getMomPerms(std::string mom) {
        if (mom == "000") return vec2D<int>{{0,0,0}};
        else if (mom == "001") return momList_001;
        else if (mom == "011") return momList_011;
        else if (mom == "111") return momList_111;
        else if (mom == "002") return momList_002;
        else if (mom == "012") return momList_012;
        else if (mom == "112") return momList_112;

        // If not in allowed set, throw error.
        else throw std::string("Momentum " + mom + " not in allowed set {000,001,011,111,002,012,112}.");        
    }

    std::vector<std::string> getIrreps(int etaTilde, std::string mom, int helicity) {
        std::vector<std::string> irrepList;
        if ((helicity == 0) && (etaTilde == 1)) irrepList.push_back("A1");
        else if (helicity == 0) irrepList.push_back("A2");
        else if (mom == "000") {
            if (etaTilde == 1) irrepList.push_back("T1p");
            else irrepList.push_back("T2m");
        }
        else if ((mom == "001") || (mom == "002")) irrepList.push_back("E2");
        else if (mom == "011") {
            irrepList.push_back("B1");
            irrepList.push_back("B2");
        }
        else if (mom == "111") irrepList.push_back("E2");        
        else if ((mom == "210") || (mom == "211")) {
            irrepList.push_back("A1");
            irrepList.push_back("A2");
        }
        return irrepList;
    }

    // Only written up to helicity 1 for current need.
    double subductHelicity(int etaTilde, std::string irrep, std::string mom, int helicity, int irrepRow) {
        if (helicity == 0) {
            if ((etaTilde == 1) && (irrep == "A1")) return 1.0;
            else if ((etaTilde == -1) && (irrep == "A2")) return 1.0;
            else return 0.0;
        }
        if (std::abs(helicity) == 1) {
            int s = 0;
            if ((irrep == "E2") && ((mom == "001") || (mom == "111") || (mom == "002"))) {
                s = (irrepRow == 1) ? 1 : -1;
            }
            else if (((irrep == "B1") || (irrep == "B2")) && (mom == "011")) {
                s = (irrep == "B1") ? 1 : -1;
            }
            else if ((irrep == "A1") || (irrep == "A2")) {
                s = (irrep == "A1") ? -1 : 1;
            }
            else return 0.0;
            return (kDelta(helicity, 1) + s * etaTilde * kDelta(helicity, -1)) / std::sqrt(2.0);
        }
        return 0.0;
    }
}