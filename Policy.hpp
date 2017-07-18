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
//  Policy.hpp
//
//  The class Policy stores a piecewise-monotonic policy.
//

#ifndef Policy_hpp
#define Policy_hpp


class Policy{

private:
    
    int *S_0; //first switching time
    int *S_1; //second switching time

    //values at switching times
    double *A_0;
    double *A_1;


public:
    
    Policy(int size, int *T_0, int T_1);
    ~Policy();
    
    void set_S(bool a, int x, int t);
    int get_S(bool a, int x);
    
    void set_A(bool a, int x, double A);
    double get_A(bool a, int x);
    
    double get_A(bool a, int x, int t);
    
    //requires the policies of the same size.
    void set_equal(int size, Policy *A);
};

#endif /* Policy_hpp */
