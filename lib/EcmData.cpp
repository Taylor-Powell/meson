#include <numbers>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "EcmData.h"

namespace ecm {
    /* Function to read the data from file for a particular energy.
    * Format of the file should be as:
    * 
    * HEADER
    * Anistropy at_mass_pion
    * 
    * HEADER
    * V momemtum_string irrep Evalue error
    * V momemtum_string irrep Evalue error
    * ... etc
    */
    void EcmData::readData(std::string filename) {
        std::ifstream file (filename);
        energyLevel E;
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
                    else if (var == "P") ss >> parity;
                    else if (var == "J") ss >> spin;
                    else continue;
                }
                else if ((var[0] == 'V') && (var[1] != ' ')) {
                    var.erase(0, 1); // Remove "V" from string
                    std::stringstream ss(var);
                    ss >> E.V >> E.mom >> E.irrep >> E.E >> E.err;
                    E.at_mpi = at_mpi;
                    E.twopi_chiL = 2.0 * std::numbers::pi / (anis * E.V);
                    Evals.push(E);
                }
                else continue;
            }
            numLvls = Evals.size();
            etaTilde = parity * std::pow(-1, spin);
            
            if ((anis == 0.0) || (at_mpi == 0.0) || (numLvls == 0)) {
                std::cout << "Printing file parameters as read from file:\n";
                printParams();
                std::string errormsg = "One or more variables not initialized in ";
                errormsg += __func__;
                throw errormsg;
            }
        }
    }

    /* Function to pop off the top level of the Evals queue
    *   - removes next item in queue  */
    energyLevel EcmData::popLevelInfo() {
        energyLevel l = Evals.front();
        Evals.pop();
        return l;
    }

    /**  Function to print the values stored in the EcmData class */
    void EcmData::printParams() {
        std::cout << "anis = " << anis << std::endl;
        std::cout << "at_mpi = " << at_mpi << std::endl;
        std::cout << "numLvls = " << numLvls << std::endl;
        std::cout << std::left;
        for (int i = 0; i < numLvls; i++) {
            energyLevel l = Evals.front();
            std::cout << "V = " << std::setw(5) <<  l.V 
                      << "E = " << std::setw(10) <<  l.E 
                      << "+/- " << std::setw(11) << l.err
                      << "mom = " << std::setw(6) << l.mom
                      << "irrep = " << std::setw(6) << l.irrep
                      << std::endl;
            Evals.pop();
            Evals.push(l);
        }
    }

    double getQsq(ecm::energyLevel& l, std::vector<int> q, std::vector<int> mom) {
        double Eval = std::sqrt(pow(l.at_mpi,2) + std::pow(l.twopi_chiL, 2) * basics::dot(mom, mom));
        return -std::pow(l.E - Eval, 2) + std::pow(l.twopi_chiL, 2) * basics::dot(q, q);
    }

    /** Function to check if momentum is <= 211 in all permutations */
    bool check3Mom(std::vector<int> mom3) {
        if (((abs(mom3[0]) < 3) && (abs(mom3[0]) < 2) && (abs(mom3[0]) < 2))
            || ((abs(mom3[0]) < 2) && (abs(mom3[0]) < 3) && (abs(mom3[0]) < 2)) 
            || ((abs(mom3[0]) < 2) && (abs(mom3[0]) < 2) && (abs(mom3[0]) < 3)))
            return true;
        return false;
    }

    void EcmData::outputEvsQsq(std::string outfile, basics::vec2D<int> qMomList) {
        std::ofstream out(outfile);
        double val;
        energyLevel l;
        basics::vec2D<int> pMomList;
        std::vector<int> qmom, pmom, mom3;
        basics::vec2D<double> outvals;
        std::vector<outVals> output;
        int index;

        if (!out) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        while (!isQueueEmpty()) {
            l = popLevelInfo();
            pMomList = basics::getMomPerms(l.mom);
            for (int i = 0; i < qMomList.size(); i++) {
                qmom = qMomList[i];
                for (int j = 0; j < pMomList.size(); j++) {
                    pmom = pMomList[j];
                    mom3 = { pmom[0]-qmom[0], pmom[1]-qmom[1], pmom[2]-qmom[2]};
                    if (!check3Mom(mom3)) continue;
                    val = getQsq(l, qmom, mom3);
                    
                    index = basics::findMatchingPair_vec2D(outvals, l.E, val);
                    if (index == -1) {
                        outvals.push_back({l.E, val});
                        ecm::outVals ov({qmom}, {pmom}, {mom3}, val, l.E);
                        output.push_back(ov);
                    }
                    else {
                        output[index].count += 1;
                        output[index].qmoms.push_back(qmom);
                        output[index].Pmoms.push_back(pmom);
                        output[index].pmoms.push_back(mom3);
                    }
                }
            }
        }
        out << "Ecm qsq count {Pmom pmom qmom}" << std::endl;
        for (int i = 0; i < output.size(); i++) {
            out << output[i].Ecm << " " << output[i].qsq << " " << output[i].count;
            for (int j = 0; j < output[i].qmoms.size(); j++) {
                out << " {" << output[i].Pmoms[j][0] << output[i].Pmoms[j][1] << output[i].Pmoms[j][2];
                out << " " << output[i].pmoms[j][0] << output[i].pmoms[j][1] << output[i].pmoms[j][2];
                out << " " << output[i].qmoms[j][0] << output[i].qmoms[j][1] << output[i].qmoms[j][2] << "}";
            }
            if (i != outvals.size() - 1) out << std::endl;
        }
        out.close();
    }
    
}