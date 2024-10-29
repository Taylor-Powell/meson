#include "read_Ecm_dat.h"
#include "generic_funcs.h"

namespace readEcm {

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
            getline(file, ignore); // Secondary header "V px py pz irrep E err"

            while (file >> E.V >> E.mom >> E.irrep >> E.E >> E.err) {
                Evals.push(E);
            }
        }
        numLvls = Evals.size();
    }

    energyLevel EcmData::popLevelInfo() {
        energyLevel l = Evals.front();
        Evals.pop();
        return l;
    }
    
}