//
//  PCOS
//
//  This code is provided as a supplement to the paper titled
//  "Optimal Stopping-Time Problems with a Probabilistic Constraint"
//  by Aaron Zeff Palmer and Alexander Vladimirsky, submitted to
//  JOTA in 2017. A link will be provided to the journal upon acceptance.
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
//  PI.hpp
//
//  The class PI has all the parameters of the PCOS problem.
//

#ifndef PI_hpp
#define PI_hpp

class PI {
private:
    
    int ** neighbors; //The neighboring nodes that are not in X_0, i.e. N(x)\X_0
    
    int * neighbors_size; //size of the neighbors[x]
    
    int * N; //number of neighbering nodes in X, i.e. |N(x)|
    
    int size; //size of X_1
    
    double *p; //coefficient for the Markov process 
    
    double *psi; //terminal cost
    
    double *phi_0;
    
    double k; //running cost
    
    double pi; //probabilistic constraint cost threshold
    
    double epsilon; //probabilistic constraint bound
    

public:
    
    PI(int size);
    ~PI();
    
    //The difference operator for p.  W should be size big. We always assume W is 0 on X_0
    double M(double *W, int x);
    
    //expectation of W.
    double expectation(double *W);
    
    double chi(int x, int t);
    
    // initializes neighbors and p constant for an interval with reflecting boundary conditions at one end.
    void initialize_half_interval(double p);
    
    // initializes neighbors and p constant for an interval with absorbing boundary conditions on both ends.
    void initialize_interval(double p);
    
    // initializes uniform starting distribution;
    void initialize_phi0_uniform();
    
    //initializes starting distribution proportional to |N(x)|
    void initialize_phi0_N();//Use this for half-interval
    
    //initializes a point-mass at x.
    void initialize_phi0_delta(int x);
    
    void initialize_psi_constant(double psi);
    
    
    void set_neighbors(int i, int *n);
    int * get_neighbors(int i);
    
    void set_neighbors_size(int x, int n);
    int get_neighbors_size(int x);
    
    void set_N(int x, int n);
    int get_N(int x);
    
    void set_p(int x, double p);
    double get_p(int x);
    
    void set_psi(int x, double psi);
    double get_psi(int x);
    
    void set_phi_0(int x, double phi);
    double get_phi_0(int x);
    
    void set_k(double k);
    double get_k();
    
    void set_pi(double pi);
    double get_pi();
    
    void set_epsilon(double epsilon);
    double get_epsilon();
    
    int get_size();
    
};


#endif /* PI_hpp */
