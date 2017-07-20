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
//  Policy.cpp
//

#include "Policy.hpp"

Policy::Policy(int size, int *T_0, int T_1){
    
    S_0=new int[size];
    S_1=new int[size];
    A_0=new double[size];
    A_1=new double[size];


    for (int x=0; x<size; x++) {
        S_0[x]=T_0[x]+1;
        A_0[x]=1;
        S_1[x]=T_1+1;
        A_1[x]=1;
    }
}

Policy::~Policy(){
    delete [] S_0;
    delete [] S_1;
    delete [] A_0;
    delete [] A_1;

}

void Policy::set_S(bool a, int x, int t){
    if (a) {
        S_1[x]=t;
    } else {
        S_0[x]=t;
    }
}

int Policy::get_S(bool a, int x){
    if (a) {
        return S_1[x];
    }
    return S_0[x];
}

void Policy::set_A(bool a, int x, double A){
    if (a) {
        A_1[x]=A;
    } else {
        A_0[x]=A;
    }
}

double Policy::get_A(bool a, int x){
    if (a) {
        return A_1[x];
    }
    return A_0[x];
}

double Policy::get_A(bool a, int x, int t){
    if (a) {
        if (t<S_1[x]) {
            return 0;
        }
        if (t==S_1[x]) {
            return A_1[x];
        }
        return 1;
    }
    if (t<S_0[x]) {
        return 0;
    }
    if (t==S_0[x]) {
        return A_0[x];
    }
    return 1;
}

void Policy::set_equal(int size, Policy *A){
    for (int x=0; x<size; x++) {
        S_1[x]=A->get_S(1, x);
        S_0[x]=A->get_S(0, x);
        A_1[x]=A->get_A(1, x);
        A_0[x]=A->get_A(0, x);
    }
}
