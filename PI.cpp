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
//  PI.cpp
//

#include "PI.hpp"

PI::PI(int size){
    neighbors=new int*[size];
    
    neighbors_size= new int[size];
    N = new int[size];
    
    p=new double[size];
    
    psi = new double [size];
    
    phi_0 = new double[size];
    
    this->size=size;
}

PI::~PI(){
    for (int i=0; i<size; i++) {
        delete []neighbors[i];
    }
    delete []neighbors;
    delete []neighbors_size;
    delete []N;
    delete []p;
    delete []psi;
    delete []phi_0;
}

double PI::M(double *W, int x){

    double val=(1-p[x])*W[x];
    
    for (int i=0; i<neighbors_size[x]; i++) {
        val+=p[x]/N[x]*W[neighbors[x][i]];
    }
    return val;
}

double PI::expectation(double *W){
    double E=0;
    for (int x=0; x<size; x++) {
        E+=W[x]*phi_0[x];
    }
    return E;
}

double PI::chi(int x, int t){
    if (k*t+psi[x]<=pi) {
        return 0;
    }
    return 1;
}

void PI::initialize_half_interval(double p){
    
    if (size>1) {
        neighbors[0]=new int[1];
        neighbors[0][0]=1;
        this->p[0]=p;
        neighbors_size[0]=1;
        N[0]=2;
        
        for (int x=1; x<size-1; x++) {
            neighbors[x]=new int[2];
            neighbors[x][0]=x-1;
            neighbors[x][1]=x+1;
            this->p[x]=p;
            N[x]=2;
            neighbors_size[x]=2;
        }
        
        neighbors[size-1]=new int[1];
        neighbors[size-1][0]=size-2;
        this->p[size-1]=p;
        neighbors_size[size-1]=1;
        N[size-1]=1;
    } else {
        if (size==1) {
            neighbors[0]=new int[0];
            this->p[0]=p;
            N[0]=2;
            neighbors_size[0]=0;
        }
    }
    
}

void PI::initialize_interval(double p){
    
    if (size>1) {
        neighbors[0]=new int[1];
        neighbors[0][0]=1;
        this->p[0]=p;
        N[0]=2;
        neighbors_size[0]=1;
        
        for (int x=1; x<size-1; x++) {
            neighbors[x]=new int[2];
            neighbors[x][0]=x-1;
            neighbors[x][1]=x+1;
            this->p[x]=p;
            N[x]=2;
            neighbors_size[x]=2;
        }
        
        neighbors[size-1]=new int[1];
        neighbors[size-1][0]=size-2;
        this->p[size-1]=p;
        N[size-1]=2;
        neighbors_size[size-1]=1;
    } else {
        if (size==1) {
            neighbors[0]=new int[0];
            this->p[0]=p;
            N[0]=2;
            neighbors_size[0]=0;
        }
    }
    
}


void PI::initialize_psi_constant(double psi){
    for (int x=0; x<size; x++) {
        this->psi[x]=psi;
    }
}


void PI::initialize_phi0_uniform(){
    double phi=1./size;
    for (int x=0; x<size; x++) {
        this->phi_0[x]=phi;
    }
}

void PI::initialize_phi0_N(){
    double sum=0;
    for (int x=0; x<size; x++) {
        phi_0[x]=N[x];
        sum+=phi_0[x];
    }
    for (int x=0; x<size; x++) {
        phi_0[x]=phi_0[x]/sum;
    }
}

void PI::initialize_phi0_delta(int x){
    for (int i=0; i<size; i++) {
        phi_0[i]=0;
    }
    phi_0[x]=1;
}


void PI::set_neighbors(int i, int *n){
    neighbors[i]=n;
}

int * PI::get_neighbors(int i){
    return neighbors[i];
}

void PI::set_neighbors_size(int x, int n){
    neighbors_size[x]=n;
}

int PI::get_neighbors_size(int x){
    return neighbors_size[x];
}

void PI::set_N(int x, int n){
    N[x]=n;
}

int PI::get_N(int x){
    return N[x];
}
void PI:: set_p(int i, double p) {
    this->p[i]=p;
}

double PI::get_p(int i){
    return p[i];
}

void PI::set_psi(int x, double psi){
    this->psi[x]=psi;
}

double PI::get_psi(int x){
    return psi[x];
}

void PI::set_phi_0(int x, double phi){
    this->phi_0[x]=phi;
}

double PI::get_phi_0(int x){
    return phi_0[x];
}

void PI::set_k(double k){
    this->k=k;
}

double PI::get_k(){
    return k;
}

void PI::set_pi(double pi){
    this->pi=pi;
}

double PI::get_pi(){
    return pi;
}

void PI::set_epsilon(double epsilon){
    this->epsilon=epsilon;
}

double PI::get_epsilon(){
    return epsilon;
}

int PI::get_size(){
    return size;
}


