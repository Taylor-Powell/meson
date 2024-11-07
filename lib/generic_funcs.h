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
   template <typename T>
    bool isUniquePair_vec2D(vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
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
   template <typename T>
    int findMatchingPair_vec2D(vec2D<T>& list, T v1, T v2, int i1=0, int i2=1, T eps = 1.0e-10) {
        for (int i = 0; i < list.size(); i++)
            if ((std::abs(list[i][i1] - v1) < eps) && 
                (std::abs(list[i][i2] - v2) < eps))
                return i;
        return -1; 
    }

    /*
    * Given a momentum string in ascending order, returns the set of 
    * allowed permutations as a 2D vector of integers
    */
    vec2D<int> getMomPerms(std::string mom);
}


#endif