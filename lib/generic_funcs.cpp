#include <cmath>
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

    std::string getIrrep(int etaTilde, std::string mom, int spin, int helicity) {
        if ((helicity == 0) && (etaTilde == 1)) return "A1";
        else if (helicity ==0) return "A2";
        if (mom == "000") {
            if ((helicity == 1) && (etaTilde == 1)) return "T1p";
            else if (helicity == 1) return "T2m";
        }
        if ((mom == "001") || (mom == "002")) {
            if (helicity == 1) return "E2";
        }
        if (mom == "011") {
            if (helicity == 1) return "B1B2";
        }
        if (mom == "111") {
            if (helicity == 1) return "E2";
        }
        if ((mom == "210") || (mom == "211")) {
            if (helicity == 1) return "A1A2";
        }
    }
}