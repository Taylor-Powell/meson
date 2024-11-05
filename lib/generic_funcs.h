#ifndef __generic_funcs_h__
#define __generic_funcs_h__

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <tuple>
#include <string>
#include <cmath>

namespace basics {
    template <typename T>
    using vec2D = std::vector<std::vector<T>>;
    
    template <typename T>
    T dot(std::vector<T>& v1, std::vector<T>& v2) {
        T val = T();
        if (v1.size() != v2.size()) 
            throw std::string("Attempting to dot vectors of unequal size.\n");
        for (int i = 0; i < v1.size(); i++) val += v1[i] * v2[i];
        return val;
    }

    /*
    * Function to search through vec2D list for a matching pair.
    *   - Returns false is a match is found
    *   - Default indices are {0,1}
    *   - Default precision is 1.0e-10
    */
    bool isUniquePair_vec2D(vec2D<double>& list, double v1, double v2, int i1=0, int i2=1, double eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return false;
        return true;
    }
    /*
    * Function to search through vec2D list for a matching pair. 
    *   - v1 and v2 are the values to find matches for
    *   - i1 and i2 are the column indices to search over
    *   - Returns -1 if no match is found
    *   - otherwise returns row index of FIRST match
    *   - Default indices are {0,1}
    *   - Default precision is 1.0e-10
    */
    int findMatchingPair_vec2D(vec2D<double>& list, double v1, double v2, int i1=0, int i2=1, double eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return i;
        return -1;
    }
    
    /*
    * Set of momentum lists with allowed cubic rotations and their dims
    * Includes {000,001,011,111,002,012,112}
    */
    vec2D<int> momList_001 {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    vec2D<int> momList_002 {{2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2}};
    vec2D<int> momList_011 {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}};
    vec2D<int> momList_111 {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};
    vec2D<int> momList_012 {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}};
    vec2D<int> momList_112 {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}};

    /*
    * Given a momentum string in ascending order, returns the set of 
    * allowed permutations as a 2D vector of integers
    */
    vec2D<int> momPerms(std::string mom) {
        if (mom == "000") return vec2D<int>{{0,0,0}};
        else if (mom == "001") return momList_001;
        else if (mom == "011") return momList_011;
        else if (mom == "111") return momList_111;
        else if (mom == "002") return momList_002;
        else if (mom == "012") return momList_012;
        else if (mom == "112") return momList_112;

        // If not in allowed set, throw error.
        else throw std::string("Momentum " + mom + " not in allowed set {000,001,011,111,002,012,112}.\n");        
    }
}


#endif