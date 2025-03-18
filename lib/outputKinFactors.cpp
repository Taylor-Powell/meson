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
                    ss >> s.V >> s.momstr >> s.irrep >> s.E >> s.Eerr;
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
            
            if ((anis == 0.0) || (at_mpi == 0.0) || (numLvls == 0)) {
                //std::cout << "Printing file parameters as read from file:\n";
                //printParams();
                std::string errormsg = "One or more variables not initialized in ";
                errormsg += __func__;
                throw errormsg;
            }

            // Test outputKinematics with first line of .txt file
            #if 1

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
            outState s = outStates[i];
            std::cout << "V = " << std::setw(5) << s.V 
                      << "E = " << std::setw(10) << s.E 
                      << "+/- " << std::setw(11) << s.Eerr
                      << "mom = " << std::setw(6) << s.momstr
                      << "irrep = " << std::setw(6) << s.irrep
                      << std::endl;
        }
    }

    void Data::outputKinematics(const std::string outfile, const basics::vec2D<int> qMomList) {
        std::ofstream out(outfile);
        // matelem::state s;
        // basics::vec2D<int> pMomList;
        // std::vector<qTuple> qTuples;
        // std::vector<int> piMom;
        // int index;

        // if (!out) {
        //     std::string errormsg = "Failed to open file in ";
        //     errormsg += __func__;
        //     throw errormsg;
        // }
        // // Loop over lines in input .txt file for b1
        // for (int i = 0; i < outStates.size(); i++) {
        //     s = outStates[i];

        //     pMomList = basics::getMomPerms(s.momstr);
        //     // Loop over pMomList for b1 meson
        //     for (int j = 0; j < pMomList.size(); j++) {
        //         // Initialize the state for the b1 meson
        //         // Hardcoding Jz=0 and the irrep is 1-dimensional
        //         matelem::state in(s.params, pMomList[j], anis, spin, parity, 0, 0, false);

        //         // Loop over qMomList to create unique qTuples
        //         qTuples = getqTuples(qMomList);
        //         for (int k = 0; k < qTuples.size(); k++) {
        //             // Create the parameter string for the current. 
        //             // I don't care about E and Eerr for now, so its hard-coded to 0.0
        //             std::string cur_param = std::to_string(s.V) + " " + qTuples[k].momStr + " " + qTuples[k].irrep + " 0.0 0.0";

        //             // Initialize the state for the current
        //             // Hardcoding mass=0 and J^P=1^-
        //             matelem::state cur(cur_param, qTuples[k].qMom3, anis, 1, -1, qTuples[k].irrepRow, qTuples[k].helicity, true);

        //             // Fix the pion state for the current qMom
        //             for (int idx = 0; idx < 3; idx++) {
        //                 piMom.push_back(pMomList[j][idx] - qTuples[k].qMom3[idx]);
        //             }

        //             // Create the parameter string for the pion with helper function
        //             std::string pi_param = getPiParamString(piMom, anis, at_mpi, s.twopi_chiL, s.V);
                    
        //             // Initialize the state for the pion
        //             // Hardcoding JP=0-, helicity=0 and the irrep is 1-dimensional
        //             matelem::state out(pi_param, piMom, anis, 0, -1, 0, 0, false);

        //             // Create the matelem object and push it to the queue
        //             matelem::matelem m(in, cur, out);
        //             m.subductAll(false);
        //             m.calcKinFactors();
        //             std::complex<double> Qsq = m.getQsq(in, out);

        //             // Clear piMom for next iteration
        //             piMom.clear();
        //         }
        //     }
        // }        
        out.close();
    }

    // Kind of a silly function... But it gets the job done.
    std::vector<qTuple> Data::getqTuples(const basics::vec2D<int> qMomList) {
        std::vector<qTuple> qTuples;
        for (int k = 0; k < qMomList.size(); k++) {
            qTuple q;
            q.qMom3 = qMomList[k];
            q.momStr = std::to_string(qMomList[k][0]) + std::to_string(qMomList[k][1]) + std::to_string(qMomList[k][2]);
            for (int hel = -1; hel <= 1; hel++) {
                q.helicity = hel;
                if (q.momStr == "000") {
                    q.irrep = "T1m";
                    for (int row = 0; row < 3; row++) {
                        q.irrepRow = row;
                        qTuples.push_back(q);
                    }
                }
                else if (hel == 0) {
                    q.irrep = "A2";
                    q.irrepRow = 0;
                    qTuples.push_back(q);
                }
                else if (q.momStr == "011") {
                    q.irrepRow = 0;
                    q.irrep = "B1";
                    qTuples.push_back(q);
                    q.irrep = "B2";
                    qTuples.push_back(q);
                }
                else {
                    q.irrep = "E2";
                    for (int row = 0; row < 2; row++) {
                        q.irrepRow = row;
                        qTuples.push_back(q);
                    }
                }
            }
        }
        return qTuples;
    }

    // Another helper function to declutter Data::outputKinematics()
    std::string Data::getPiParamString(const std::vector<int> piMom, double anis, double at_mpi, double twopi_chiL, int V) {
        std::string pi_param = std::to_string(V) + " " + rotations::getMomStr(piMom);
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