#include <iomanip>
#include <complex>
#include "outputKinFactors.h"
#include "matelem.h"

namespace kinFactors {
    void Data::readData(const std::string filename) {
        std::ifstream file (filename);
        outState s;
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
                    else continue;
                }
                else if ((var[0] == 'V') && (var[1] != ' ')) {
                    var.erase(0, 1); // Remove "V" from string
                    std::stringstream ss(var);
                    std::string momstr;
                    ss >> s.V >> momstr >> s.irrep >> s.E >> s.Eerr;
                    s.mom3_i = basics::getMom3_i(momstr);
                    s.momType = basics::getMomType(s.mom3_i);
                    s.params = var;
                    s.outfile = "data/kinFactors_V_" + std::to_string(s.V);
                    s.outfile += "_E_" + std::to_string(s.E).substr(0, std::to_string(s.E).find('.') + 7);
                    s.outfile += "_b1mom_" + s.momstr + ".dat";
                    outStates.push_back(s);
                }
                else continue;
            }
            numLvls = outStates.size();
            etaTilde = parity * std::pow(-1, spin);
            #if 0
            printParams();
            #endif
            
            if ((anis == 0.0) || (at_mpi == 0.0) || (numLvls == 0)) {
                std::string errormsg = "One or more variables not initialized in ";
                errormsg += __func__;
                throw errormsg;
            }

            // Test outputKinematics with first line of .txt file
            #if 0

            // Verbose output
            bool printSteps = false;
            #if 0
            printSteps = true;
            #endif

            basics::vec2D<int> qMomList = {{0,0,1}};
            outState s = outStates[0];
            std::ofstream fout(s.outfile);
            double twopi_chiL = 2.0 * std::numbers::pi / (anis * s.V);
            std::cout << "Param string = " << s.params << std::endl;
            std::vector<int> mom3;
            for (int i = 0; i < s.momstr.size(); i++) {
                mom3.push_back(s.momstr[i] - '0');
            }
            if (printSteps) std::cout << "Making in state..." << std::endl;
            matelem::state in(s.V, s.irrep, s.E, s.Eerr, mom3, spin, parity, 0, 0, twopi_chiL, false);


            if (printSteps) std::cout << "Making out state..." << std::endl;
            std::vector<qTuple> qTuples = getqTuples(qMomList);
            std::vector<int> piMom;
            std::string pi_irrep = "A1";
            for (int i = 0; i < 3; i++) {
                piMom.push_back(mom3[i] - qTuples[0].qMom3[i]);
                if (piMom[i] != 0) pi_irrep = "A2";
            }
            double Epi = std::sqrt(std::pow(at_mpi, 2) + std::pow(twopi_chiL, 2) * basics::dot(piMom, piMom));
            matelem::state out(s.V, pi_irrep, Epi, 0.0, piMom, 0, -1, 0, 0, twopi_chiL, false);

            
            if (printSteps) std::cout << "Making current state..." << std::endl;
            for (int i = 0; i < qTuples.size(); i++) {
                matelem::state cur(s.V, qTuples[i].irrep, s.E - Epi, 0.0, qTuples[i].qMom3, 1, -1, qTuples[i].irrepRow, qTuples[i].helicity, twopi_chiL, true);
                if (printSteps) std::cout << "Making matelem..." << std::endl;
                matelem::matelem m(in, cur, out, s.V, anis, twopi_chiL);
                if (printSteps) std::cout << "Calculating kinematic factors..." << std::endl;
                m.calcKinFactors();
                if (printSteps) std::cout << "Writing kinematic factors..." << std::endl;
                m.writeKinFactors(fout);
            }
            
            fout.close();



            #endif
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
            std::cout << "V = " << std::setw(5) << outStates[i].V 
                      << "E = " << std::setw(10) << outStates[i].E 
                      << "+/- " << std::setw(11) << outStates[i].Eerr
                      << "mom = " << std::setw(6) << outStates[i].momstr
                      << "irrep = " << std::setw(6) << outStates[i].irrep
                      << std::endl;
        }
    }

    void Data::outputKinematics(const basics::vec2D<int> qMomList) {
        for (int i = 0; i < outStates.size(); i++) {
            outState s = outStates[i];
            std::ofstream fout(s.outfile);
            double twopi_chiL = 2.0 * std::numbers::pi / (anis * s.V);

            std::vector<stateTuple> qTuples, inTuples;
            for (int i = 0; i < qMomList.size(); i++) {
                std::vector<stateTuple> temp = getTuples(qMomList[i], parity, spin);
                qTuples.insert(qTuples.end(), temp.begin(), temp.end());
            }

            // Make the 3-momentum base vector, check if it's <211, and get the momentum permutations
            if (!basics::check3Mom(s.mom3_i)) continue;
            basics::vec2D<int> pMomList = basics::getMomPerms(s.mom3_i);

            // Iterate over the momentum permutations
            for (int j = 0; j < pMomList.size(); j++) {
                matelem::state in(s.V, s.irrep, s.E, s.Eerr, pMomList[j], spin, parity, 0, 0, twopi_chiL, false);                
                // Iterate over all qTuples (mom3, irrep, irrepRow)
                for(int k = 0; k < qTuples.size(); k++) {
                    // Fix the pion momentum and irrep, then ensure piMom < 211
                    std::vector<int> piMom;
                    std::string pi_irrep = "A1";
                    for (int ii = 0; ii < 3; ii++) {
                        piMom.push_back(pMomList[j][ii] - qTuples[k].mom3_i[ii]);
                        if (piMom[ii] != 0) pi_irrep = "A2";
                    }
                    if (!basics::check3Mom(piMom)) continue;

                    // Calculate Epi and create the states for the pion and current
                    double Epi = std::sqrt(std::pow(at_mpi, 2) + std::pow(twopi_chiL, 2) * basics::dot(piMom, piMom));
                    matelem::state out(s.V, pi_irrep, Epi, 0.0, piMom, 0, -1, 0, 0, twopi_chiL, false);
                    
                    #if 0 // Gotta fix the current state creation
                    matelem::state cur(s.V, qTuples[k].irrep, s.E - Epi, 0.0, qTuples[k].qMom3, 1, -1, qTuples[k].irrepRow, qTuples[k].helicity, twopi_chiL, true);

                    // Create the matrix element container, then calculate the kinematic factors and append them to the output container
                    matelem::matelem m(in, cur, out, s.V, anis, twopi_chiL);
                    m.calcKinFactors();
                    // Append m.kFactors to kFactors using m.getKinFactors()
                    std::vector<std::vector<cd>> kF = m.getKinFactors();
                    kFactors.insert(kFactors.end(), kF.begin(), kF.end());
                    // Append m.outstring to outStrings using m.getOutStrings()
                    std::vector<std::string> oS = m.getOutStrings();
                    outStrings.insert(outStrings.end(), oS.begin(), oS.end());
                    #endif
                }
            }
            // sort kFactors and outStrings by Q_sq
            std::vector<std::pair<double, int>> QsqIndex;
            for (int i = 0; i < kFactors.size(); i++) {
                QsqIndex.push_back(std::make_pair(kFactors[i][0].real(), i));
            }
            std::sort(QsqIndex.begin(), QsqIndex.end());
            for (int i = 0; i < QsqIndex.size(); i++) {
                int index = QsqIndex[i].second;
                fout << outStrings[index] << std::endl;
            }
            fout.close();
        }
    }

    /** Function to get the tuples {mom3_i, irrep, irrepRow} over which to iterate. */
    std::vector<stateTuple> Data::getTuples(const std::vector<int> mom3_i, int parity, int spin) {
        // Initialize the output vector and the stateTuple struct
        std::vector<stateTuple> tuples;
        stateTuple s;

        // Set the momType and get the irreps for the given mom3_i
        s.momType = basics::getMomType(mom3_i);
        std::vector<std::string> irreps = basics::getIrreps(mom3_i, parity, spin);

        // Get the permutations of the mom3_i
        basics::vec2D<int> perms = basics::getMomPerms(mom3_i);

        // Iterate over the irreps and permutations to get the tuples
        for (int i = 0; i < irreps.size(); i++) {            
            s.irrep = irreps[i];
            s.irrepRow = 0;
            for (int j = 0; j < perms.size(); j++) {
                s.mom3_i = perms[j];
                tuples.push_back(s); // Always push back the irrepRow 0 tuple
                if (s.irrep == "E2") { // 2-dimensional irreps
                    s.irrepRow = 1;
                    tuples.push_back(s);
                }
                else if (s.irrep == "T1p" || s.irrep == "T1m") { // 3-dimensional irreps
                    for (int row = 1; row < 3; row++) {
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