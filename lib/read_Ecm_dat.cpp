#include <numbers>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "read_Ecm_dat.h"

namespace ecm {
    /*
    * Function to read the data from file for a particular energy.
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
            
            if ((anis == 0.0) || (at_mpi == 0.0) || (numLvls == 0)) {
                std::cout << "Printing file parameters as read from file:\n";
                printParams();
                std::string errormsg = "One or more variables not initialized in ";
                errormsg += __func__;
                throw errormsg;
            }
        }
    }

    /*
    * Function to pop off the top level of the Evals queue
    *   - removes next item in queue
    */
    energyLevel EcmData::popLevelInfo() {
        energyLevel l = Evals.front();
        Evals.pop();
        return l;
    }

    /** 
    * Function to print the values stored in the EcmData class
    */
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

    
}