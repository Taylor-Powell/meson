// -*- C++ -*-
/* File: read_Ecm_dat.h
 *
 * Description: 
 *
 * This file contains the EcmData class to hold information read 
 * in from a file for a set of energy levels.
 *
 */

#ifndef __read_Ecm_dat_h__
#define __read_Ecm_dat_h__

#include <iostream>
#include <fstream>
#include <string>
#include <queue>

namespace readEcm {
    struct energyLevel {
        int V;
        std::string irrep, mom;
        double E, err, anis, at_mpi;
    };
    class EcmData {
        public:
            // Constructors and Destructor
            EcmData() {}
            EcmData(std::string filename) { readData(filename); }
            ~EcmData() {}

            // Inline functions
            inline int getNumLevels() { return numLvls; };
            inline energyLevel frontLevelInfo() { return Evals.front(); }
            inline void pushBackLevelInfo(energyLevel l) {Evals.push(l);}
            bool isQueueEmpty() { return Evals.empty(); }

            // Functions defined in read_Ecm_dat.cpp
            void readData(std::string filename);
            energyLevel popLevelInfo();

        private:
            std::queue<energyLevel> Evals;
            int numLvls;
    };

    
}




#endif