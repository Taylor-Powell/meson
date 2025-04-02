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
        /**
        Constructor for state struct
        @param V: Integer length of the spatial volume
        @param momstr: Momentum string in ascending order
        @param irrep: Irrep
        @param E: Energy
        @param Eerr: Error in the energy
        @param threemom: Integer vector of the three-momentum
        @param J: Spin
        @param P: Parity
        @param row: Row of the irrep
        @param absHel: Absolute value of target helicity 
        @param twopi_chiL: 2 * pi / (anis * V)
        @param current: Whether the state is a current
        */ 
        state(int V, std::string irrep, double E, double Eerr, std::vector<int> mom3_i, int J, int P, int row, int absHel, double twopi_chiL, bool current);
        state() {}
        std::vector<int> mom3_i;
        std::vector<cd> mom4, polVec;
        std::string params, momType, irrep, sym;
        int V, irrepRow, helicity, spin, absHel, parity, etaTilde;
        double E, Eerr, mState;
        double coeff;
    };

    struct outputStruct {
        std::string outstring;
        cd Qsq;
        std::vector<std::vector<cd>> kFactors;
        std::vector<int> helicity;
    };

    class matelem {
        public:
            /** Constructor for matelem
             * @param in: Initial state
             * @param J: Current state
             * @param out: Final state
             * @param V: Integer length of the spatial volume
             * @param anis: Anisotropy of the lattice
             * @param twopi_chiL: 2 * pi / (anis * V)
             */
            matelem(state in, state J, state out, int V, double anis, double twopi_chiL) : V(V), anis(anis), twopi_chiL(twopi_chiL) {
                init.push_back(in);
                cur.push_back(J);
                fin.push_back(out);
                Qsq = getQsq(in, out);
                outStruct.Qsq = Qsq;
            }
            /** Default constructor */
            matelem() {}

            // Functions
            void subductAll(bool isHelState = true);
            cd getQsq(const state& in, const state& out);
            bool calcKinFactors();

            // Inline functions
            inline std::vector<std::vector<cd>> getKinFactors() { return outStruct.kFactors; }
            inline std::string getOutString() { return outStruct.outstring; }
            inline void writeKinFactors(std::ofstream& fout) { fout << outStruct.outstring << std::endl; }
            inline outputStruct getOutStruct() { return outStruct; }
            // Destructor
            ~matelem() {}

            // Variables
            double anis, twopi_chiL;
            int V;            

        private:
            // Variables
            std::vector<state> init, cur, fin;
            cd Qsq;
            outputStruct outStruct;
            

            // Forward Declarations
            void subductHelicityState(std::vector<state>& s);
            cd getOmegaVal(const state& in, const state& out);
            std::vector<cd> kinFactors(const state& in, const state& cur, const state& out);
            std::string formatOutString(const std::vector<double>& coeffs, const std::vector<std::vector<int>>& indices);
    };

}

 #endif