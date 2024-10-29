#include "read_Ecm_dat.h"

namespace readEcm {

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
            std::string ignore;
            getline(file, ignore); // Opening header "anis at_mpi"
            file >> E.anis >> E.at_mpi >> ignore;
            getline(file, ignore); // Secondary header "V mom irrep E err"

            while (file >> E.V >> E.mom >> E.irrep >> E.E >> E.err) {
                Evals.push(E);
            }
        }
        numLvls = Evals.size();
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
    
}