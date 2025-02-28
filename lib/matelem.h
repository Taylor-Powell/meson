// -*- C++ -*-
/* File: matelem.h
 *
 * Description: 
 *
 * This file contains the matelem class
 */

 #ifndef __matelem_h__
 #define __matelem_h__
 
 #include <iostream>
 #include <fstream>
 #include <string>
 #include "generic_funcs.h"
 #include "rotations.h"

namespace {
    typedef std::complex<double> cd;
    const cd j1 = cd(0,1);
}
namespace matelem {
    struct state {
        state(std::string p, std::vector<int> threemom, double anis, double mass, int J, int P, int row, int hel, bool current);
        state() {}
        std::vector<int> mom;
        std::vector<cd> fourMom, polVec;
        std::string params, momstr, irrep, sym;
        int irrepRow, helicity, spin, spinZ, parity, etaTilde, V;
        double E, Eerr, twopi_chiL, mState;
        cd coeff;
    };

    class matelem {
        public:
            matelem() {}
            matelem(state in, state J, state out) {
                init.push_back(in);
                cur.push_back(J);
                fin.push_back(out);
            }
            // Functions
            void projectAll();
            void calcKinFactors();
            ~matelem() {}
        private:
            // Variables
            std::vector<state> init, cur, fin;
            std::vector<std::vector<cd>> kFactors;

            // Forward Declarations
            void ExpandHelOps(std::vector<state>& s);
            void subductState(std::vector<state>& s);
            cd getOmegaVal(state& in, state& out);
            std::vector<cd> kinFactors(state& in, state& cur, state& out);
    };

}

 #endif