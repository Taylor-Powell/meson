#include <numbers>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "matelem.h"
#include "Wigner.h"


namespace matelem {
    state::state(int V, std::string irrep, double E, double Eerr, std::vector<int> mom3_i, int J, int P, int row, int helicity, double twopi_chiL, bool current) 
        : V(V), irrep(irrep), E(E), Eerr(Eerr), mom3_i(mom3_i), spin(J), helicity(helicity), parity(P), irrepRow(row) {        
        momType = basics::getMomType(mom3_i);
        std::vector<double> mom3 = {mom3_i[0] * twopi_chiL, mom3_i[1] * twopi_chiL, mom3_i[2] * twopi_chiL};
        mom4 = {E, mom3[0], mom3[1], mom3[2]};
        double mom3_sq = basics::dot(mom3, mom3);
        if (current) mState = 0.0;
        else mState = std::sqrt((basics::fourDot(mom4, mom4)).real()); 
        coeff = 1.0;
        etaTilde = parity * std::pow(-1, spin);
        sym = basics::getSym(momType);
        polVec = rotations::getPol4_hel(E, mom3_sq, mom3_i, helicity, sym, current);
    }

    void matelem::expandAllHelOps() {
        ExpandHelOps(init);
        ExpandHelOps(cur);
        ExpandHelOps(fin);
    }

    void matelem::subductAll(bool isHelState) {
        if (!isHelState) {
            throw std::string("subductAll called with isHelState = false.\n");
        }
        double eps = 1e-10;
        bool zeroCheck = false;

        #if 0
        zeroCheck = true;
        #endif

        subductHelicityState(init);
        if (zeroCheck) for (int i = 0; i < init.size(); i++) {
            if (std::abs(init[i].coeff) < eps) {
                std::cout << "Zero for b1(" << std::to_string(init[i].mom3_i[0]) << std::to_string(init[i].mom3_i[1]) << std::to_string(init[i].mom3_i[2]) << ") ";
                std::cout << init[i].irrep << "(" << init[i].irrepRow << ") etaTilde=" << init[i].etaTilde;
                std::cout << ", hel=" << init[i].helicity << std::endl;
            }
        }
        subductHelicityState(cur);
        if (zeroCheck) for (int i = 0; i < cur.size(); i++) {
            if (std::abs(cur[i].coeff) < eps) {
                std::cout << "Zero for current(" << std::to_string(cur[i].mom3_i[0]) << std::to_string(cur[i].mom3_i[1]) << std::to_string(cur[i].mom3_i[2]) << ") ";
                std::cout << cur[i].irrep << "(" << cur[i].irrepRow << ") etaTilde=" << cur[i].etaTilde;
                std::cout << ", hel=" << cur[i].helicity << std::endl;
            }
        }
        subductHelicityState(fin);
        if (zeroCheck) for (int i = 0; i < fin.size(); i++) {
            if (std::abs(fin[i].coeff) < eps) {
                std::cout << "Zero for pion(" << std::to_string(fin[i].mom3_i[0]) << std::to_string(fin[i].mom3_i[1]) << std::to_string(fin[i].mom3_i[2]) << ") ";
                std::cout << fin[i].irrep << "(" << fin[i].irrepRow << ") etaTilde=" << fin[i].etaTilde;
                std::cout << ", hel=" << fin[i].helicity << std::endl;
            }
        }
    }
    
    bool matelem::calcKinFactors() {
        // std::cout << "Calculating kinematic factors..." << std::endl;
        std::vector<cd> kin;
        std::vector<std::vector<cd>> kins;
        std::vector<std::vector<int>> indices;
        std::vector<double> coeffs;
        cd Qsq = 0.0;
        std::string out;
        double epsilon = 1e-10;
        int skipCount = 0, totalCount = 0;
        bool zeroCheck = false;

        #if 0
        zeroCheck = true;
        #endif

        for (int i = 0; i < init.size(); i++) {
            if (abs(init[i].coeff) < epsilon) {
                if (zeroCheck) std::cout << "Skipping b1 with irrep=" << init[i].irrep << "(" << init[i].irrepRow << ") and momType=" << init[i].momType << std::endl;
                skipCount++;
                continue;
            }
            for (int j = 0; j < cur.size(); j++) {
                if (abs(cur[j].coeff) < epsilon) {
                    if (zeroCheck) std::cout << "Skipping current with irrep=" << cur[j].irrep << "(" << cur[j].irrepRow << ") and momType=" << cur[j].momType << std::endl;
                    skipCount++;
                    continue;
                }
                for (int k = 0; k < fin.size(); k++) {
                    if (abs(fin[k].coeff) < epsilon) {
                        if (zeroCheck) std::cout << "Skipping pion with irrep=" << fin[k].irrep << "(" << fin[k].irrepRow << ") and momType=" << fin[k].momType << std::endl;
                        skipCount++;
                        continue;
                    }
                    Qsq = getQsq(init[i], fin[k]);
                    coeffs.push_back(init[i].coeff * cur[j].coeff * fin[k].coeff);
                    kin.clear();
                    kin = kinFactors(init[i], cur[j], fin[k]);
                    kins.push_back(kin);
                    #if 0
                    std::cout << "Kinematic factors for state " << i << " " << j << " " << k << ": ";
                    std::cout << "E1 = " << kin[0];
                    std::cout << ", C1 = " << kin[1] << std::endl;
                    #endif
                    /** ///////////////////////////////////////////////
                     * Need to fix HERE 
                     * kFactors should have same indices as outsttrings
                     * *///////////////////////////////////////////////

                    kFactors.push_back(kin);
                    indices.push_back({i, j, k});                    
                    totalCount++;
                }
            }
        }
        #if 0
        if (zeroCheck) std::cout << "Skipped " << skipCount << " states with zero coefficient.    ";
        if (zeroCheck) std::cout << "Calculated " << totalCount << " kinematic factors." << std::endl;
        #endif

        if (totalCount == 0) return false;

        out = "Qsq=" + basics::formatValue(Qsq.real()) + ", ";
        out += "b1(";
        out += std::to_string(init[0].mom3_i[0]) + std::to_string(init[0].mom3_i[1]) + std::to_string(init[0].mom3_i[2]) + ") ";
        out += init[0].irrep + "(" + std::to_string(init[0].irrepRow) + "), ";
        out += "q(";                    
        out += std::to_string(cur[0].mom3_i[0]) + std::to_string(cur[0].mom3_i[1]) + std::to_string(cur[0].mom3_i[2]) + ") ";
        out += cur[0].irrep + "(" + std::to_string(cur[0].irrepRow) + "), ";
        out += "pion(";
        out += std::to_string(fin[0].mom3_i[0]) + std::to_string(fin[0].mom3_i[1]) + std::to_string(fin[0].mom3_i[2]) + ") ";
        out += fin[0].irrep + "(" + std::to_string(fin[0].irrepRow) + ")\n";
        for (int i = 0; i < totalCount; i++) {
            out += "     b1_hel=" + std::to_string(init[indices[i][0]].helicity) + ", ";
            out += "q_hel=" + std::to_string(cur[indices[i][1]].helicity) + ", ";
            out += "coeff=" + basics::formatValue(coeffs[i]) + ", ";
            out += "E1=(" + basics::formatValue(kins[i][0].real()) + ", " + basics::formatValue(kins[i][0].imag()) + "), ";
            out += "C1=(" + basics::formatValue(kins[i][1].real()) + ", " + basics::formatValue(kins[i][1].imag()) + ")";
            if (i != totalCount - 1) out += "\n";
        }
        outstring.push_back(out);
        return true;
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
        if (s[0].helicity == 0) s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momType, s[0].helicity, s[0].irrepRow);
        else {
            s[0].coeff *= basics::subductHelicity(s[0].etaTilde, s[0].irrep, s[0].momType, s[0].helicity, s[0].irrepRow);
            state sTemp = s[0];
            sTemp.helicity = -sTemp.helicity;
            sTemp.coeff *= basics::subductHelicity(sTemp.etaTilde, sTemp.irrep, sTemp.momType, sTemp.helicity, sTemp.irrepRow);
            s.push_back(sTemp);
        }
    }
    
    void matelem::ExpandHelOps(std::vector<state>& s) {
        // std::vector<state> s2;
        // for (int i = 0; i < s.size(); i++) {
        //     state sTemp;
        //     std::vector<double> angles;
        //     std::vector<state> newStates;
        //     for (int i = 0; i <= 2 * s[i].spin + 1; i++) {
        //         sTemp = s[i];
        //         angles = rotations::getRotAngles(sTemp.sym, sTemp.mom3_i);
        //         sTemp.spinZ = -s[i].spin + i;
        //         sTemp.coeff *= std::conj(WignerD::Wigner_D(s[i].spin, s[i].spinZ, s[i].helicity, angles[0], angles[1], angles[2]));
        //         newStates.push_back(sTemp);
        //     }
        //     s2.insert(s2.end(), newStates.begin(), newStates.end());
        // }
        // s = s2;
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