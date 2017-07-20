//
//  PCOS
//
//  This code is provided as a supplement to the paper titled
//  "Optimal Stopping with a Probabilistic Constraint"
//  by Aaron Zeff Palmer and Alexander Vladimirsky, submitted to
//  the Journal of Optimization Theory and Applications in 2017. 
//  A link will be provided to the journal upon acceptance.
//
//  All the code was written by Aaron Zeff Palmer. Any questions on the
//  content of the paper may be addressed to the corresponding author,
//  Aaron Zeff Palmer, at azp6@cornell.edu. Please do not contact with
//  requests to extend this code or any solicitations.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  PCOS.hpp
//
//  The class PCOS contains the algorithms required to solve a PCOS problem
//

#ifndef PCOS_hpp
#define PCOS_hpp

#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Policy.hpp"
#include "PI.hpp"

//Either define FILEPATH here or in 'Makefile'.
//#define FILEPATH "/Users/aaronpalmer/Documents/Examples/"

using namespace std;


class PCOS {
    
private:
    
    int T_1;
    int *T_0; // the discontinuity of \chi(x,t)
    
    PI *Pi; //Problem parameters
    
    double *U; //unconstrained value function
    
    double *R_f;
    double *R_s; //conditional constrained values at time 0.
    
    double E_0;
    double P_0;
    
    double E_m;
    double P_m;
    
    double lambda_f;
    double lambda_s;
    double lambda_flat;
    double lambda_sharp; //Lagrange multipliers
    
    double P_f;
    double P_s;
    double P_flat;
    double P_sharp; //constrained values
 
    double E_f;
    double E_s;
    double E_flat;
    double E_sharp; //expected costs
    
    Policy *A_f;
    Policy *A_s;
    Policy *A_flat;
    Policy *A_sharp; //Policies are in piecewise-monotonic form


public:
    
    PCOS(PI* Pi);
    ~PCOS();
    
    void output(string name); //outputs R_f,R_s, Phi_0, U, A_f, A_s, lambda_s, lambda_f in a format that is readable by Matlab
    
    void solve_unconstrained(double threshold); //See line 2 of Algorithm 1
    
    void solve_minimal_constrained(); //Solves for E_m and P_m of the minimal constrained value policy A^m. See obs 1.
    
    void solve_constrained(double Delta); //Algorithm 1
    
    void resolve_degeneracies_forward(); //Algorithm 2
    void resolve_degeneracies_backward(); //Algorithm 3
    
    double get_lambda_f();
    double get_lambda_s();
    double get_E_f();
    double get_E_s();
    
    
};

#endif /* PCOS_hpp */
