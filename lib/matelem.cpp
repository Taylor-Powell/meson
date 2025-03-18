#include <numbers>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "matelem.h"
#include "Wigner.h"


namespace matelem {
    state::state(int V, std::string irrep, double E, double Eerr, std::vector<int> mom3_i, int J, int P, int row, int Jz, double twopi_chiL, bool current) 
        : V(V), irrep(irrep), E(E), Eerr(Eerr), mom3_i(mom3_i), spin(J), spinZ(Jz), parity(P), irrepRow(row) {        
        momstr = rotations::getMomStr(mom3_i);
        std::vector<double> mom3 = {mom3_i[0] * twopi_chiL, mom3_i[1] * twopi_chiL, mom3_i[2] * twopi_chiL};
        mom4 = {E, mom3[0], mom3[1], mom3[2]};
        double mom3_sq = basics::dot(mom3, mom3);
        if (current) mState = 0.0;
        else mState = std::sqrt((basics::fourDot(mom4, mom4)).real()); 
        coeff = 1.0;
        etaTilde = parity * std::pow(-1, spin);
        sym = rotations::getSym(momstr);
        polVec = rotations::getPol4_Jz(E, mom3_sq, mom3_i, Jz, sym, current);
    }

    void matelem::expandAllHelOps() {
        ExpandHelOps(init);
        ExpandHelOps(cur);
        ExpandHelOps(fin);
    }

    void matelem::subductAll(bool isHelState) {
        if (isHelState) {
            subductHelicityState(init);
            subductHelicityState(cur);
            subductHelicityState(fin);
        }
        else {
            subductJzState(init);
            subductJzState(cur);
            subductJzState(fin);
        }
    }
    
    void matelem::calcKinFactors() {
        std::cout << "Calculating kinematic factors..." << std::endl;
        std::vector<cd> kin;
        for (int i = 0; i < init.size(); i++) {
            for (int j = 0; j < cur.size(); j++) {
                for (int k = 0; k < fin.size(); k++) {
                    kin = kinFactors(init[i], cur[j], fin[k]);
                    #if 1
                    std::cout << "Kinematic factors for state " << i << " " << j << " " << k << ": ";
                    std::cout << "E1 = " << kin[0];
                    std::cout << ", C1 = " << kin[1] << std::endl;
                    #endif
                    kFactors.push_back(kin);

                    std::string out  = "p_b1=";
                    out += std::to_string(init[i].mom3_i[0]) + std::to_string(init[i].mom3_i[1]) + std::to_string(init[i].mom3_i[2]) + " q=";
                    out += std::to_string(cur[j].mom3_i[0]) + std::to_string(cur[j].mom3_i[1]) + std::to_string(cur[j].mom3_i[2]) + " irrep=";
                    out += cur[j].irrep + " row=" + std::to_string(cur[j].irrepRow) + " ";
                    cd Qsq = getQsq(init[i], fin[k]);
                    out += "Qsq=" + std::to_string(Qsq.real()) + " ";
                    out += " E1=(" + std::to_string(kin[0].real()) + ", " + std::to_string(kin[0].imag()) + ") ";
                    out += " C1=(" + std::to_string(kin[1].real()) + ", " + std::to_string(kin[1].imag()) + ")";
                    outstring.push_back(out);
                }
            }
        }
    }

    cd matelem::getQsq(state& in, state& out) {
        bool print = false;
        #if 0
        print = true;
        #endif
        
        std::vector<cd> qMom;
        if (print) std::cout << "qMom = ";
        for (int i = 0; i < 4; i++) {
            qMom.push_back(out.mom4[i] - in.mom4[i]);
            if (print) std::cout << out.mom4[i] << " - " << in.mom4[i] << " = " << qMom[i] << ", ";
        }
        if (print) std::cout << "qMom size = " << qMom.size() << std::endl;
        return (-1.0) * basics::fourDot(qMom, qMom);
    }

    void matelem::subductHelicityState(std::vector<state>& s) {
        if (s.size() != 1) throw std::string("subductState called with vector of size != 1.\n");
        if (s[0].etaTilde == 0) s[0].etaTilde = s[0].parity * std::pow(-1, s[0].spin);
        if (s[0].helicity == 0) s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momstr, s[0].helicity, s[0].irrepRow);
        else {
            s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momstr, s[0].helicity, s[0].irrepRow);
            s.push_back(s[0]);
            s[1].helicity = -s[1].helicity;
            s[1].coeff *= basics::subductHelicity(s[1].etaTilde, s[1].irrep, s[1].momstr, s[1].helicity, s[1].irrepRow);
        }
    }

    void matelem::subductJzState(std::vector<state>& s) {
        if (s.size() != 1) throw std::string("subductState called with vector of size != 1.\n");
        
    }
    
    void matelem::ExpandHelOps(std::vector<state>& s) {
        std::vector<state> s2;
        for (int i = 0; i < s.size(); i++) {
            state sTemp;
            std::vector<double> angles;
            std::vector<state> newStates;
            for (int i = 0; i <= 2 * s[i].spin + 1; i++) {
                sTemp = s[i];
                angles = rotations::getRotAngles(sTemp.momstr, sTemp.mom3_i);
                sTemp.spinZ = -s[i].spin + i;
                sTemp.coeff *= std::conj(WignerD::Wigner_D(s[i].spin, s[i].spinZ, s[i].helicity, angles[0], angles[1], angles[2]));
                newStates.push_back(sTemp);
            }
            s2.insert(s2.end(), newStates.begin(), newStates.end());
        }
        s = s2;
    }

    cd matelem::getOmegaVal(state& in, state& out) { // Checked
        cd val = std::pow(basics::fourDot(in.mom4, out.mom4), 2);
        val -= std::pow(in.mState, 2) * std::pow(out.mState, 2);        
        return val;
    }

    std::vector<cd> matelem::kinFactors(state& in, state& cur, state& out) { // Checked
        std::vector<cd> kin;
        cd Omega = getOmegaVal(in, out);
        cd tempVal;

        //////////////////////////////////////////////////////////////////////////////////////////
        #if 0
        // Output input values for each state for debugging
        std::cout << "in.formMom = (";
        for (int i = 0; i < 4; i++) {
            std::cout << in.fourMom[i] << ", ";
        }
        std::cout << ")" << std::endl;
        std::cout << "in.polVec = (";
        for (int i = 0; i < 4; i++) {
            std::cout << in.polVec[i] << ", ";
        }
        std::cout <<  ")" << std::endl;
        std::cout << "cur.fourMom = (";
        for (int i = 0; i < 4; i++) {
            std::cout << cur.fourMom[i] << ", ";
        }
        std::cout <<  ")" << std::endl;
        std::cout << "cur.polVec = (";
        for (int i = 0; i < 4; i++) {
            std::cout << cur.polVec[i] << ", ";
        }
        std::cout <<  ")" << std::endl;
        std::cout << "out.fourMom = (";
        for (int i = 0; i < 4; i++) {
            std::cout << out.fourMom[i] << ", ";
        }
        std::cout <<  ")" << std::endl;
        std::cout << "out.polVec = (";
        for (int i = 0; i < 4; i++) {
            std::cout << out.polVec[i] << ", ";
        }
        std::cout <<  ")\n" << std::endl;
        #endif
        //////////////////////////////////////////////////////////////////////////////////////////
        // Calculate Qsq
        kin.push_back(getQsq(in, out));


        // Calculate the E1 coefficient
        std::vector<cd> Ecoeff;
        for (int i = 0; i < 4; i++) { // Checked
            tempVal = basics::fourDot(in.mom4, out.mom4) * in.mom4[i];
            tempVal -= std::pow(in.mState, 2) * out.mom4[i];
            tempVal *= basics::fourDot(in.polVec, out.mom4) / Omega;
            tempVal = out.polVec[i] - tempVal; 
            Ecoeff.push_back(tempVal);
        }
        kin.push_back(basics::fourDot(Ecoeff, cur.polVec));


        // Calculate the C1 coefficient
        std::vector<cd> Ccoeff;
        for (int i = 0; i < 4; i++) { // Checked
            tempVal = basics::fourDot(in.mom4, out.mom4);
            tempVal *= in.mom4[i] * out.mom4[i];
            tempVal -= std::pow(in.mState, 2) * out.mom4[i];
            tempVal -= std::pow(out.mState, 2) * in.mom4[i];
            tempVal *= basics::fourDot(in.polVec, out.mom4) * in.mState / Omega;
            tempVal /= std::sqrt(basics::fourDot(cur.mom4, cur.mom4));
            Ccoeff.push_back(tempVal);
        }
        kin.push_back(basics::fourDot(Ccoeff, cur.polVec));

        return kin;
    }
}