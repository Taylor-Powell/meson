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
            std::vector<qTuple> qTuples = getqTuples(qMomList);

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

    // For a given current 3-momentum, expand over the possible irreps and rows
    // "000" -> T1m rows 0-2
    // "00n" -> E2 rows 0-1
    // "0nn" -> B1 row 0, B2 row 0
    // "nnn" -> E2 rows 0-1
    // "0mn" -> A1 row 0, A2 row 0
    std::vector<qTuple> Data::getqTuples(const basics::vec2D<int> qMomList) {
        std::vector<qTuple> qTuples;
        for (int k = 0; k < qMomList.size(); k++) {
            qTuple q;
            q.mom3_i = qMomList[k];
            q.momType = basics::getMomType(q.mom3_i);
            if (q.momType == "000") {
                q.irrep = "T1m";
                for (int row = 0; row < 3; row++) {
                    q.irrepRow = row;
                    qTuples.push_back(q);
                }
            }
            else { // All other cases get an A2 contribution for hel=0
                q.irrep = "A2";
                q.irrepRow = 0;
                qTuples.push_back(q);
            }
            if (q.momType == "00n") {
                q.irrep = "E2";
                for (int row = 0; row < 2; row++) {
                    q.irrepRow = row;
                    qTuples.push_back(q);
                }
            }
            else if (q.momType == "0mn") {
                q.irrepRow = 0;
                q.irrep = "A1";
                qTuples.push_back(q);
                q.irrep = "A2";
                qTuples.push_back(q);
            }
            
            if (q.momType == "0nn") {
                q.irrepRow = 0;
                q.irrep = "B1";
                qTuples.push_back(q);
                q.irrep = "B2";
                qTuples.push_back(q);
            }
            else { // "nnn" case
                q.irrep = "E2";
                for (int row = 0; row < 2; row++) {
                    q.irrepRow = row;
                    qTuples.push_back(q);
                }
            }
        }
        return qTuples;
    }

    

    // Given a "qTuple", expand into superposition of helicity states with coeff basics::subductHelicity
    std::vector<std::pair<cd, int>> Data::getHelStates(int etaTilde, std::string momType, std::string irrep, int irrepRow, int spin) {
        std::vector<std::pair<cd, int>> coeffANDhelicities;
        if (spin == 0) {
            coeffANDhelicities.push_back(std::pair(1.0, 0));
        }
        else if (spin == 1) {
            for (int hel = -1; hel < 2; hel++) { // Iterate over possible helicities, calculate subduction coeffs, add to stack.
                double val = basics::subductHelicity(etaTilde, irrep, momType, hel, irrepRow);
                coeffANDhelicities.push_back(std::pair(val, hel));
            }
        }
        else throw std::string("Spin > 1 in Data::getHelStates()\n ");
        return coeffANDhelicities;
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