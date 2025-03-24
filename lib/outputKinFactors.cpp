#include <iomanip>
#include <complex>
#include "outputKinFactors.h"
#include "matelem.h"

namespace kinFactors {
    void Data::readData(const std::string filename) {
        std::ifstream file (filename);
        inState s;
        if (!file) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        else {
            std::string ignore, var, val;

            while(std::getline(file, var)) {
                if (var[0] != 'V') {
                    std::stringstream ss(var);
                    std::getline(ss, var, '=');
                    if (var == "anis") ss >> anis;
                    else if (var == "at_mpi") ss >> at_mpi;
                    else if (var == "at_inv") ss >> at_inv;
                    else if (var == "P") ss >> parity;
                    else if (var == "J") ss >> spin;
                    else if (var == "inAbsHelicity") ss >> inAbsHelicity;
                    else if (var == "outAbsHelicity") ss >> outAbsHelicity;
                    else continue;
                }
                else if ((var[0] == 'V') && (var[1] != ' ')) {
                    var.erase(0, 1); // Remove "V" from string
                    std::stringstream ss(var);
                    std::string momstr, tempFile;
                    ss >> s.V >> momstr >> s.irrep >> s.E >> s.Eerr;
                    s.mom3_i = basics::getMom3_i(momstr);
                    s.momType = basics::getMomType(s.mom3_i);
                    s.params = var;
                    tempFile = "data/kinFactors_V_" + std::to_string(s.V);
                    tempFile += "_b1mom_" + momstr;
                    tempFile += "_irrep_" + s.irrep;
                    tempFile+= "_E_" + std::to_string(s.E).substr(0, std::to_string(s.E).find('.') + 6);
                    double twopi_chiL = 2.0 * std::numbers::pi / (anis * s.V);
                    s.irrepRow = 1;
                    s.outfile = tempFile + "_irrepRow_" + std::to_string(s.irrepRow) + ".dat";
                    inStates.push_back(s);
                    // Rigging it simply for now... b1 states are either A2 or T1p
                    // If generalizing, CHANGE THIS
                    if (s.irrep == "T1p") {
                        s.irrepRow = 2;
                        s.outfile = tempFile + "_irrepRow_" + std::to_string(s.irrepRow) + ".dat";
                        inStates.push_back(s);
                        s.irrepRow = 3;
                        s.outfile = tempFile + "_irrepRow_" + std::to_string(s.irrepRow) + ".dat";
                        inStates.push_back(s);
                    }
                }
                else continue;
            }
            numLvls = inStates.size();
            etaTilde = parity * std::pow(-1, spin);
            #if 0
            printParams();
            #endif
            
            if ((anis == 0.0) || (at_mpi == 0.0) || (numLvls == 0)) {
                std::string errormsg = "One or more variables not initialized in ";
                errormsg += __func__;
                throw errormsg;
            }
        }
    }

    void Data::printParams() {
        std::cout << "anis = " << anis << std::endl;
        std::cout << "at_mpi = " << at_mpi << std::endl;
        std::cout << "at_inv = " << at_inv << std::endl;
        std::cout << "P" << parity << std::endl;
        std::cout << "J" << spin << std::endl;
        std::cout << "numLvls = " << numLvls << std::endl;
        std::cout << std::left;
        for (int i = 0; i < numLvls; i++) {
            std::cout << "V = " << std::setw(5) << inStates[i].V 
                      << "E = " << std::setw(10) << inStates[i].E 
                      << "+/- " << std::setw(11) << inStates[i].Eerr
                      << "mom = " << std::setw(6) << inStates[i].momstr
                      << "irrep = " << std::setw(6) << inStates[i].irrep
                      << std::endl;
        }
    }

    void Data::outputKinematics(const basics::vec2D<int> qMomList) {
        // Fix helicity of b1 to 0
        int inHelicity = 0;

        // Define a small epsilon
        double epsilon = 1.0e-10;

        // Iterate over the qMomList to get the qTuples
        // These are fixed for each qMom, so only need to be calculated once 
        std::vector<stateTuple> qTuples;
        for (int i = 0; i < qMomList.size(); i++) {
            // Use version allowing for absHel = {0, 1}
            std::vector<stateTuple> temp = getTuples(qMomList[i], parity, spin);
            qTuples.insert(qTuples.end(), temp.begin(), temp.end());
        }

        // Iterate over lines of inState txt file (b1 states)
        for (int i = 0; i < inStates.size(); i++) {
            inState s = inStates[i];
            double coeff = basics::subductHelicity(etaTilde, s.irrep, s.momType, inHelicity, s.irrepRow);
            if (std::abs(coeff) < epsilon) continue;

            // New file for each level
            std::ofstream fout(s.outfile);
            double twopi_chiL = 2.0 * std::numbers::pi / (anis * s.V);

            // Make the 3-momentum base vector, check if it's <211, and get the momentum permutations
            // The check here should be redundant since the input file should only
            //      contain valid 3-momenta
            if (!basics::check3Mom(s.mom3_i)) continue;
            basics::vec2D<int> pMomList = basics::getMomPerms(s.mom3_i);

            // Iterate over the momentum permutations of inState
            for (int k = 0; k < pMomList.size(); k++) {
                matelem::state in(s.V, s.irrep, s.E, s.Eerr, pMomList[k], spin, parity, s.irrepRow, inHelicity, twopi_chiL, false);

                // Iterate over all qTuples (mom3, irrep, irrepRow)
                for(int qdx = 0; qdx < qTuples.size(); qdx++) {
                    // Fix the pion momentum and irrep, then ensure piMom < 211
                    std::vector<int> piMom3_i;
                    std::string pi_irrep = "A1"; // Inflexible, CHANGE THIS if generalizing
                    for (int ii = 0; ii < 3; ii++) {
                        piMom3_i.push_back(pMomList[k][ii] - qTuples[qdx].mom3_i[ii]);
                        if (piMom3_i[ii] != 0) pi_irrep = "A2";
                    }
                    if (!basics::check3Mom(piMom3_i)) continue;

                    // Calculate Epi and create the states for the pion and current
                    double Epi = std::sqrt(std::pow(at_mpi, 2) + std::pow(twopi_chiL, 2) * basics::dot(piMom3_i, piMom3_i));
                    matelem::state out(s.V, pi_irrep, Epi, 0.0, piMom3_i, 0, -1, 1, 0, twopi_chiL, false);

                    // Create the current state
                    matelem::state cur(s.V, qTuples[qdx].irrep, s.E - Epi, 0.0, qTuples[qdx].mom3_i, 1, -1, qTuples[qdx].irrepRow, qTuples[qdx].absHelicity, twopi_chiL, true);

                    // Create the matrix element container, then calculate the kinematic factors and append them to the output container
                    matelem::matelem m(in, cur, out, s.V, anis, twopi_chiL);
                    m.subductAll(true);
                    bool anythingUseful = m.calcKinFactors();

                    if (anythingUseful) {
                        // Append m.kFactors to kFactors using m.getKinFactors()
                        std::vector<std::vector<cd>> kF = m.getKinFactors();
                        kFactors.insert(kFactors.end(), kF.begin(), kF.end());
                        // Append m.outstring to outStrings using m.getOutStrings()
                        std::vector<std::string> oS = m.getOutStrings();
                        outStrings.insert(outStrings.end(), oS.begin(), oS.end());
                    }
                }
            }           

            
            // sort kFactors and outStrings by Q_sq
            std::vector<std::pair<double, int>> QsqIndex;
            std::cout << "There are " << kFactors.size() << " kinematic factors." << std::endl;
            for (int i = 0; i < kFactors.size(); i++) {
                QsqIndex.push_back(std::make_pair(kFactors[i][0].real(), i));
            }
            std::sort(QsqIndex.begin(), QsqIndex.end());
            for (int i = 0; i < QsqIndex.size(); i++) {
                int index = QsqIndex[i].second;
                fout << outStrings[index] << std::endl;
                std::cout << "idx=" << index << ": " << outStrings[index] << std::endl;
            }
            fout.close();
        }
    }

    /** Overloaded function to get the tuples {mom3_i, irrep, irrepRow} over which to iterate. */
    std::vector<stateTuple> Data::getTuples(const std::vector<int> mom3_i, int parity, int spin) {
        // Initialize the output vector and the stateTuple struct
        std::vector<stateTuple> tuples;
        stateTuple s;

        // Also have to iterate over |\lambda|
        for (int idx = 0; idx <= spin; idx++) {
            s.absHelicity = idx;

            // Set the momType and get the irreps for the given mom3_i
            s.momType = basics::getMomType(mom3_i);
            std::vector<std::string> irreps = basics::getIrreps(mom3_i, parity, spin);

            // Get the permutations of the mom3_i
            basics::vec2D<int> perms = basics::getMomPerms(mom3_i);

            // Iterate over the irreps and permutations to get the tuples
            for (int i = 0; i < irreps.size(); i++) {            
                s.irrep = irreps[i];
                s.irrepRow = 1;
                for (int j = 0; j < perms.size(); j++) {
                    s.mom3_i = perms[j];
                    tuples.push_back(s); // Always push back the irrepRow 1 tuple
                    if (s.irrep == "E2") { // 2-dimensional irreps
                        s.irrepRow = 2;
                        tuples.push_back(s);
                    }
                    else if (s.irrep == "T1p" || s.irrep == "T1m") { // 3-dimensional irreps
                        for (int row = 2; row <= 3; row++) {
                            s.irrepRow = row;
                            tuples.push_back(s);
                        }
                    }    
                }        
            }   
        }
        return tuples;
    }

    /** Overloaded function to get the tuples {mom3_i, irrep, irrepRow} over which to iterate. */
    std::vector<stateTuple> Data::getTuples(const std::vector<int> mom3_i, int parity, int spin, int targetAbsHel) {
        // Initialize the output vector and the stateTuple struct
        std::vector<stateTuple> tuples;
        stateTuple s;

        // Check that absHelMax is valid. If must be >= 0 and <= spin
        if (targetAbsHel < 0 || targetAbsHel > spin) {
            throw std::string("Invalid targetAbsHel in Data::getMoreTuples()\n");
        }
        s.absHelicity = targetAbsHel;

        // Set the momType and get the irreps for the given mom3_i
        s.momType = basics::getMomType(mom3_i);
        std::vector<std::string> irreps = basics::getIrreps(mom3_i, parity, spin);

        // Get the permutations of the mom3_i
        basics::vec2D<int> perms = basics::getMomPerms(mom3_i);

        // Iterate over the irreps and permutations to get the tuples
        for (int i = 0; i < irreps.size(); i++) {            
            s.irrep = irreps[i];
            s.irrepRow = 1;
            for (int j = 0; j < perms.size(); j++) {
                s.mom3_i = perms[j];
                tuples.push_back(s); // Always push back the irrepRow 1 tuple
                if (s.irrep == "E2") { // 2-dimensional irreps
                    s.irrepRow = 2;
                    tuples.push_back(s);
                }
                else if (s.irrep == "T1p" || s.irrep == "T1m") { // 3-dimensional irreps
                    for (int row = 2; row <= 3; row++) {
                        s.irrepRow = row;
                        tuples.push_back(s);
                    }
                }    
            }        
        }   
        return tuples;
    }

    // Get subduced helicity states with coeff from basics::subductHelicity()
    std::vector<std::pair<cd, int>> Data::getHelCoeffs(int etaTilde, std::string momType, std::string irrep, int irrepRow, int absHel) {
        // Sanity checks
        if (absHel > 1) std::string("abs(helicity) > 1 in Data::getHelStates()\n ");
        else if (absHel < 0) std::string("abs(helicity) < 0 in Data::getHelStates()\n ");

        // Initialize the output vector and the pair<cd, int> struct
        std::vector<std::pair<cd, int>> vals;
        vals.push_back(std::pair(basics::subductHelicity(etaTilde, irrep, momType, absHel, irrepRow), absHel));
        if (absHel != 0) { // If helicity is not 0, add the helicity -absHel state
            vals.push_back(std::pair(basics::subductHelicity(etaTilde, irrep, momType, -absHel, irrepRow), -absHel));
        }
        return vals;
    }

    // Another helper function to declutter Data::outputKinematics()
    std::string Data::getPiParamString(const std::vector<int> piMom, double anis, double at_mpi, double twopi_chiL, int V) {
        std::string pi_param = std::to_string(V) + " " + basics::getMomType(piMom);
        if (piMom[0] == 0 && piMom[1] == 0 && piMom[2] == 0) {
            pi_param += " A1 ";
        }
        else {
            pi_param += " A2 ";
        }
        double Epi = std::sqrt(std::pow(at_mpi,2) + std::pow(twopi_chiL, 2) * basics::dot(piMom, piMom));
        pi_param += std::to_string(Epi) + " 0.0";
        return pi_param;
    }
}