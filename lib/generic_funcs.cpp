#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cctype>
#include "generic_funcs.h"

namespace basics {
    // Flexible function to get all permutations of a 3-digit momentum string
    // Allows for string to be unordered
    vec2D<int> getMomPerms(const std::string momstr) {
        if (momstr.size() != 3) {
            throw std::invalid_argument("Momentum " + momstr + " is not 3 digits in basics::getMomPerms.\n");
        }
        // Initialize variables
        vec2D<int> perms;
        std::vector<int> mom3_i;
        
        // Convert string to vector of integers
        for (char c : momstr) {
            if (!std::isdigit(c)) {
                throw std::invalid_argument("Momentum " + momstr + " is not a set of integers in basics::getMomPerms.\n");
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

    vec2D<int> getMomPerms(const std::vector<int>& mom3_i) {
        if (mom3_i.size() != 3) {
            throw std::invalid_argument("Momentum " + std::to_string(mom3_i[0]) + std::to_string(mom3_i[1]) + std::to_string(mom3_i[2]) + " is not 3 digits in basics::getMomPerms.\n");
        }
        // Initialize variables
        vec2D<int> perms;
        std::vector<int> sortedMom3_i = mom3_i;
        
        std::sort(sortedMom3_i.begin(), sortedMom3_i.end());
        do { // Generate all permutations
            int numFlips = 8; // 2^3 possible sign flips
            for (int i = 0; i < numFlips; ++i) { // Loop over all possible sign flips
                std::vector<int> flippedVec = sortedMom3_i;
                for (int j = 0; j < 3; ++j) {
                    if (i & (1 << j)) { // If jth bit is set
                        flippedVec[j] = -flippedVec[j];
                    }
                }
                perms.push_back(flippedVec);
            }
        } while (std::next_permutation(sortedMom3_i.begin(), sortedMom3_i.end()));

        // Remove duplicate permutations
        std::sort(perms.begin(), perms.end());
        perms.erase(std::unique(perms.begin(), perms.end()), perms.end());
        return perms;
    }

    std::vector<std::string> getIrreps(const std::vector<int> mom3_i, int parity, int spin) {
        // Sanity checks
        if (spin > 1) throw std::invalid_argument("Spin > 1 in Data::getIrreps()\n");
        else if (spin < 0) throw std::invalid_argument("Spin < 0 in Data::getIrreps()\n");
        else if (parity != 1 && parity != -1) throw std::invalid_argument("Invalid parity in Data::getIrreps()\n");

        std::vector<std::string> irreps;
        int etaTilde = parity * std::pow(-1, spin);
        std::string momType = getMomType(mom3_i);
        if (momType == "000") {
            if (spin == 1) {
                if (parity == 1) irreps.push_back("T1p");
                else if (parity == -1) irreps.push_back("T1m");
            }
            else if (spin == 0) irreps.push_back("A1");
        }
        else { // Catch-all for hel=0 subductions
            if (etaTilde == 1) irreps.push_back("A1");
            else if (etaTilde == -1) irreps.push_back("A2");
        }
        if (momType == "00n") {
            irreps.push_back("E2");
        }
        else if (momType == "0nn") {
            irreps.push_back("B1");
            irreps.push_back("B2");
        }
        else if (momType == "nnn") {
            irreps.push_back("E2");
        }
        else if (momType == "0mn") { // Covers the hel=1 subductions for "0mn" and "nnm"
            if (etaTilde == 1) irreps.push_back("A2");
            else if (etaTilde == -1) irreps.push_back("A1");
        }
        return irreps;
    }

    std::vector<std::string> getIrreps(const std::string momType, int parity, int spin) {
        // Sanity checks
        if (spin > 1) throw std::invalid_argument("Spin > 1 in basics::getIrreps()\n");
        else if (spin < 0) throw std::invalid_argument("Spin < 0 in basics::getIrreps()\n");
        else if (parity != 1 && parity != -1) throw std::invalid_argument("Invalid parity in basics::getIrreps()\n");

        std::vector<std::string> irreps;
        int etaTilde = parity * std::pow(-1, spin);
        if (momType == "000") {
            if (spin == 1) {
                if (parity == 1) irreps.push_back("T1p");
                else if (parity == -1) irreps.push_back("T1m");
            }
            else if (spin == 0) irreps.push_back("A1");
        }
        else { // Catch-all for hel=0 subductions
            if (etaTilde == 1) irreps.push_back("A1");
            else if (etaTilde == -1) irreps.push_back("A2");
        }
        if (momType == "00n") {
            irreps.push_back("E2");
        }
        else if (momType == "0nn") {
            irreps.push_back("B1");
            irreps.push_back("B2");
        }
        else if (momType == "nnn") {
            irreps.push_back("E2");
        }
        else if (momType == "0mn") { // Covers the hel=1 subductions for "0mn" and "nnm"
            if (etaTilde == 1) irreps.push_back("A2");
            else if (etaTilde == -1) irreps.push_back("A1");
        }
        return irreps;
    }

    // Only written up to helicity 1 for current need.
    double subductHelicity(int etaTilde, const std::string irrep, const std::string momType, int helicity, int irrepRow) {

        // Should split this into subductions for OhD and helicity
        if (momType == "000" ) {
            if (irrep == "A1") return 1.0;
            else if (etaTilde == 1) {
                if (irrep == "T1m") {
                    if (irrepRow - 2 == helicity) return 1.0;
                    else return 0.0;
                }
            }
            else if (etaTilde == -1) {
                if (irrep == "T1p") {
                    if (irrepRow - 2 == helicity) return 1.0;
                    else return 0.0;
                }
            }
            else return 0.0;
        }
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
            else if (((irrep == "A1") || (irrep == "A2")) && ((momType == "0mn") || (momType == "nnm"))) {
                s = (irrep == "A1") ? -1 : 1;
            }
            else return 0.0;
            return (kDelta(helicity, 1) + s * etaTilde * kDelta(helicity, -1)) / std::sqrt(2.0);
        }
        return 0.0;
    }

    std::vector<int> getMom3_i(const std::string momStr) {
        std::vector<int> mom3_i;
        for (int i = 0; i < 3; i++) {
            // if not a digit, throw an error
            if (!std::isdigit(momStr[i])) {
                throw std::invalid_argument("Momentum " + momStr + " is not a set of integers in basics::getMom3_i.\n");
            }
            mom3_i.push_back(momStr[i] - '0');
        }
        return mom3_i;
    }

    /** Convert from mom3_i to string for momentum type. 
     * Function is safe against all permutations of mom3_i
     * Allowed set is {000, 001, 011, 111, 002, 012, 112, 022, 003, 122}
     */
    std::string getMomType(const std::vector<int> mom3_i) {
        if (mom3_i.size() != 3) {
            throw std::invalid_argument("Momentum " + std::to_string(mom3_i[0]) + std::to_string(mom3_i[1]) + std::to_string(mom3_i[2]) + " is not 3 digits in basics::getMomType.\n");
        }
        std::string sym;
        int mom_sq = basics::dot(mom3_i, mom3_i);
        if (mom_sq == 0) return "000";
        else if (mom_sq == 1) return "00n"; // 001
        else if (mom_sq == 2) return "0nn"; // 011
        else if (mom_sq == 3) return "nnn"; // 111
        else if (mom_sq == 4) return "00n"; // 002
        else if (mom_sq == 5) return "0mn"; // 012
        else if (mom_sq == 6) return "nnm"; // 112

        // Start of "extra terms" > 112
        else if (mom_sq == 8) return "0nn"; // 022
        else if (mom_sq == 9) { // 003 or 122
            if (mom3_i[0] == 3 || mom3_i[1] == 3 || mom3_i[2] == 3) 
                return "00n"; // 003
            else return "nnm"; // 112
        }
        else throw std::invalid_argument("Momentum " + std::to_string(mom3_i[0]) + std::to_string(mom3_i[1]) + std::to_string(mom3_i[2]) + " not recognized in basics::getMomType.\n");
    }
}