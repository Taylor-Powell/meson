#ifndef __rotations_h__
#define __rotations_h__

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <tuple>
#include <string>
#include <cmath>
#include <complex>
#include "generic_funcs.h"

namespace rotations {
    typedef std::complex<double> cd;
    cd j1 = cd(0,1);

    template <typename T>
    using vec2D = std::vector<std::vector<T>>;

    vec2D<cd> polRot_Dic4 = {{1,0,0},{0,1,0},{0,0,1}};
    vec2D<cd> polRot_Dic2 = {{0.853553,0,0},{0,0.57735,0},{0,0,0.853553}};
    vec2D<cd> polRot_Dic3 = {{0.670146-0.415829 * j1,0,0},{0,0,0.57735},{0,0,0.670146+0.415829 * j1}};
    vec2D<cd> polRot_C40mn = {{0.420578-0.848721*j1,0,0},{0,0,0.894427},{0,0,0.420578+0.848721*j1}}; 
    vec2D<cd> polRot_C4nnm = {{-0.0864113+0.904128*j1,0,0},{0,0,0.816497},{0,0,-0.0864113-0.904128*j1}};

    void rotPolVec_init(vec2D<cd>& polVec, std::string sym) {
        if (sym == "Dic4") polVec = basics::matMult(polVec, polRot_Dic4);
        else if (sym == "Dic2") polVec = basics::matMult(polVec, polRot_Dic2);
        else if (sym == "Dic3") polVec = basics::matMult(polVec, polRot_Dic3);
        else if (sym == "C40mn") polVec = basics::matMult(polVec, polRot_C40mn);
        else if (sym == "C4nnm") polVec = basics::matMult(polVec, polRot_C4nnm);
        else throw std::string("Symmetry " + sym + " not recognized in basics::rotPolVec.\n");
    }

    // Dic4
    vec2D<cd> rot_0_0_0 = {{1,0,0},{0,1,0},{0,0,1}};
    vec2D<cd> rot_0_pi2_0 = {{0,0,1},{0,1,0},{-1,0,0}};
    vec2D<cd> rot_pi2_pi2_0 = {{0,-1,0},{0,0,1},{-1,0,0}};
    vec2D<cd> rot_pi_pi2_0 = {{0,0,-1},{0,-1,0},{-1,0,0}};
    vec2D<cd> rot_npi2_pi2_0 = {{0,1,0},{0,0,-1},{-1,0,0}};
    vec2D<cd> rot_0_pi_0 = {{-1,0,0},{0,1,0},{0,0,-1}};
    vec2D<cd> rot_2pi_0_3pi2 = {{0,1,0},{0,0,-1},{-1,0,0}};
    vec2D<cd> rot_0_npi_0 = {{-1,0,0},{0,1,0},{0,0,-1}};

    void rotPolVec(vec2D<cd>& polVec, std::string sym, std::string mom) {
        if (sym == "Dic4") {
            if (mom == "001") polVec = basics::matMult(polVec, rot_0_0_0);
            else if (mom == "100") polVec = basics::matMult(polVec, rot_0_pi2_0);
            else if (mom == "010") polVec = basics::matMult(polVec, rot_pi2_pi2_0);
            else if (mom == "-100") polVec = basics::matMult(polVec, rot_pi_pi2_0);
            else if (mom == "0-10") polVec = basics::matMult(polVec, rot_npi2_pi2_0);
            else if (mom == "00-1") polVec = basics::matMult(polVec, rot_0_pi_0);
            else throw std::string("Momentum " + mom + " not recognized in basics::rotPolVec.\n");
        }
        else if (sym == "Dic2") {
            if (mom == "011") polVec = basics::matMult(polVec, rot_0_0_0);
            else if (mom == "110") polVec = basics::matMult(polVec, rot_0_pi2_0);
            else if (mom == "101") polVec = basics::matMult(polVec, rot_2pi_0_3pi2);
            else if (mom == "1-10") polVec = basics::matMult(polVec, rot_npi2_pi2_0);
            else if (mom == "01-1") polVec = basics::matMult(polVec, rot_0_npi_0);
            // Add rest or generalize
        }
    }
}
#endif