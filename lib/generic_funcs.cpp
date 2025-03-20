#include <cmath>
#include <vector>
#include <algorithm>
#include "generic_funcs.h"

namespace basics {
    // Flexible function to get all permutations of a 3-digit momentum string
    // Allows for string to be unordered
    vec2D<int> getMomPerms(const std::string momstr) {
        if (momstr.size() != 3) {
            throw std::string("Momentum " + momstr + " is not 3 digits in basics::getMomPerms.\n");
        }
        // Initialize variables
        vec2D<int> perms;
        std::vector<int> mom3_i;
        
        // Convert string to vector of integers
        for (char c : momstr) {
            if (!std::isdigit(c)) {
                throw std::string("Momentum " + momstr + " is not a set of integers in basics::getMomPerms.\n");
            }
            mom3_i.push_back(c - '0');
        }
        
        std::sort(mom3_i.begin(), mom3_i.end());
        do { // Generate all permutations
            int numFlips = 8; // 2^3 possible sign flips
            for (int i = 0; i < numFlips; ++i) { // Loop over all possible sign flips
                std::vector<int> flippedVec = mom3_i;
                bool skip = false;
                for (int j = 0; j < 3; ++j) {
                    if (i & (1 << j)) { // If jth bit is set
                        flippedVec[j] = -flippedVec[j];
                    }
                }
                perms.push_back(flippedVec);
            }
        } while (std::next_permutation(mom3_i.begin(), mom3_i.end()));

        // Remove duplicate permutations
        std::sort(perms.begin(), perms.end());
        perms.erase(std::unique(perms.begin(), perms.end()), perms.end());
        return perms;
    }

    vec2D<int> getMomPerms(const std::vector<int> mom3_i) {
        if (mom3_i.size() != 3) {
            throw std::string("Momentum " + std::to_string(mom3_i[0]) + std::to_string(mom3_i[1]) + std::to_string(mom3_i[2]) + " is not 3 digits in basics::getMomPerms.\n");
        }
        // Initialize variables
        vec2D<int> perms;
        
        std::sort(mom3_i.begin(), mom3_i.end());
        do { // Generate all permutations
            int numFlips = 8; // 2^3 possible sign flips
            for (int i = 0; i < numFlips; ++i) { // Loop over all possible sign flips
                std::vector<int> flippedVec = mom3_i;
                bool skip = false;
                for (int j = 0; j < 3; ++j) {
                    if (i & (1 << j)) { // If jth bit is set
                        flippedVec[j] = -flippedVec[j];
                    }
                }
                perms.push_back(flippedVec);
            }
        } while (std::next_permutation(mom3_i.begin(), mom3_i.end()));

        // Remove duplicate permutations
        std::sort(perms.begin(), perms.end());
        perms.erase(std::unique(perms.begin(), perms.end()), perms.end());
        return perms;
    }

    // std::vector<std::string> getIrreps(int etaTilde, const std::string mom, int helicity) {
    //     std::vector<std::string> irrepList;
    //     if ((helicity == 0) && (etaTilde == 1)) irrepList.push_back("A1");
    //     else if (helicity == 0) irrepList.push_back("A2");
    //     else if (mom == "000") {
    //         if (etaTilde == 1) irrepList.push_back("T1p");
    //         else irrepList.push_back("T2m");
    //     }
    //     else if ((mom == "001") || (mom == "002")) irrepList.push_back("E2");
    //     else if (mom == "011") {
    //         irrepList.push_back("B1");
    //         irrepList.push_back("B2");
    //     }
    //     else if (mom == "111") irrepList.push_back("E2");        
    //     else if ((mom == "210") || (mom == "211")) {
    //         irrepList.push_back("A1");
    //         irrepList.push_back("A2");
    //     }
    //     return irrepList;
    // }

    std::vector<std::string> getIrreps(int etaTilde, const std::string momType, int helicity) {
        std::vector<std::string> irrepList;
        if ((helicity == 0) && (etaTilde == 1)) irrepList.push_back("A1");
        else if (helicity == 0) irrepList.push_back("A2");
        else if (momType == "000") {
            if (etaTilde == 1) irrepList.push_back("T1p");
            else irrepList.push_back("T1m");
        }
        else if (momType == "00n") irrepList.push_back("E2");
        else if (momType == "0nn") {
            irrepList.push_back("B1");
            irrepList.push_back("B2");
        }
        else if (momType == "nnn") irrepList.push_back("E2");        
        else if ((momType == "0mn") || (momType == "mnn")) {
            irrepList.push_back("A1");
            irrepList.push_back("A2");
        }
        return irrepList;
    }

    // // Only written up to helicity 1 for current need.
    // double subductHelicity(int etaTilde, const std::string irrep, const std::string momstr, int helicity, int irrepRow) {
    //     if (helicity == 0) {
    //         if ((etaTilde == 1) && (irrep == "A1")) return 1.0;
    //         else if ((etaTilde == -1) && (irrep == "A2")) return 1.0;
    //         else return 0.0;
    //     }
    //     if (std::abs(helicity) == 1) {
    //         int s = 0;
    //         if ((irrep == "E2") && ((momstr == "001") || (momstr == "111") || (momstr == "002"))) {
    //             s = (irrepRow == 1) ? 1 : -1;
    //         }
    //         else if (((irrep == "B1") || (irrep == "B2")) && (momstr == "011")) {
    //             s = (irrep == "B1") ? 1 : -1;
    //         }
    //         else if ((irrep == "A1") || (irrep == "A2")) {
    //             s = (irrep == "A1") ? -1 : 1;
    //         }
    //         else return 0.0;
    //         return (kDelta(helicity, 1) + s * etaTilde * kDelta(helicity, -1)) / std::sqrt(2.0);
    //     }
    //     return 0.0;
    // }

    // Only written up to helicity 1 for current need.
    double subductHelicity(int etaTilde, const std::string irrep, const std::string momType, int helicity, int irrepRow) {
        if (helicity == 0) {
            if ((etaTilde == 1) && (irrep == "A1")) return 1.0;
            else if ((etaTilde == -1) && (irrep == "A2")) return 1.0;
            else return 0.0;
        }
        if (std::abs(helicity) == 1) {
            int s = 0;
            if ((irrep == "E2") && ((momType == "00n") || (momType == "nnn"))) {
                s = (irrepRow == 1) ? 1 : -1;
            }
            else if (((irrep == "B1") || (irrep == "B2")) && (momType == "0nn")) {
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

    /** Simple function to check if momentum is <= 211 in all permutations */
    bool check3Mom(const std::vector<int> mom3_i) {
        if (mom3_i.size() != 3) {
            throw std::string("Momentum vector is not 3 digits in basics::check3Mom.\n");
        }
        for (int i = 0; i < 3; i++) {
            if ((abs(mom3_i[i]) < 3) && (abs(mom3_i[(i+1)%3]) < 2) && (abs(mom3_i[(i+2)%3]) < 2))
                return true;
        }
        return false;
    }
}