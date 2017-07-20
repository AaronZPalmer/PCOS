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
//  PCOS.cpp
//
//  This file contains the algorithms used to solve a PCOS problem.
//

#include "PCOS.hpp"


PCOS::PCOS(PI* Pi){
    
    this->Pi=Pi;
    
    T_1=floor(Pi->get_pi()/Pi->get_k());
    
    int size=Pi->get_size();
    
    T_0=new int[size];
    for(int x=0; x<size; x++){
        T_0[x]=floor((Pi->get_pi()-Pi->get_psi(x))/Pi->get_k());
        
        //This is required because floor of an int might subtract 1 due to round off error.
        while(Pi->chi(x,T_0[x]+1)==0){
            T_0[x]++;
        }
    }
    
    U=new double[size];;

    R_f = new double[size];
    R_s = new double[size];
    
    A_f = new Policy(size,T_0,T_1);
    A_s = new Policy(size,T_0,T_1);
    A_flat = new Policy(size,T_0,T_1);
    A_sharp = new Policy(size,T_0,T_1);

}

PCOS::~PCOS(){

    delete []T_0;

    delete [] U;
    
    delete [] R_s;
    delete [] R_f;
    
    delete A_f;
    delete A_s;
    delete A_flat;
    delete A_sharp;
}


void PCOS::solve_unconstrained(double threshold){
    double error=threshold+1;
    double temp_val;
    int counter=0;
    int size=Pi->get_size();
    for (int x=0; x<size; x++) {
        U[x]=Pi->get_psi(x);
    }
    
    while (error>threshold) {
        error=0;
        for (int x=0; x<size; x++) {
            temp_val=U[x];
            
            U[x]=min(Pi->get_k()+Pi->M(U, x),Pi->get_psi(x));

            error=max(error,temp_val-U[x]);
        }
        counter++;
    }
    
    
    cout<<"Unconstrained problem solved in "<<counter<<" iterations."<<endl;
    E_0=Pi->expectation(U);
    cout<<"Unconstrained value is "<<E_0<<endl;

}

void PCOS::solve_minimal_constrained(){
    
    int size=Pi->get_size();
    
    double val=0;
    
    bool minimal_val=true;
    
    for(int x=0; x<size; x++){
        val+=Pi->get_phi_0(x)*Pi->get_psi(x);
        
        if(T_0[x]<0 && Pi->get_phi_0(x)>0){
            minimal_val=false;
        }
        
        A_f->set_S(0, x, 0);
        A_f->set_S(1,x,T_1+1);
    }
    
    if(minimal_val){
        E_m=val;
        P_m =0;
    } else {
        
        
        double *Z_c=new double[size];
        double *Z_n=new double[size];
        
        double *R_c=new double[size];
        double *R_n=new double[size];
        
        for(int t=T_1; t>=0; t--){
            for(int x=size-1;x>=0; x--){
                if(t<T_1){
                    bool a=Pi->chi(x,t);
                    double A_temp = A_f->get_A(a, x, t);
        
                    R_n[x]=(1-A_temp)*Pi->M(R_c,x)+A_temp*Pi->chi(x,t);
        
                    Z_n[x]=(1-A_temp)*(Pi->M(Z_c,x)+Pi->get_k())+A_temp*Pi->get_psi(x);
                } else {
                    Z_n[x]=U[x];
                    R_n[x]=1;
                }
        
           }
            for(int x=size-1; x>=0;x--){
                R_c[x]=R_n[x];
                Z_c[x]=Z_n[x];
            }
        }
    
        P_m=Pi->expectation(R_n);
        E_m=Pi->expectation(Z_n);
    
    
        delete []Z_c;
        delete []Z_n;
        delete []R_c;
        delete []R_n;
    }
    
    
}



void PCOS::solve_constrained(double Delta){
    int size=Pi->get_size();
    
    double *V_c=new double[size];
    double *V_n=new double[size];
    
    double *R_c=new double[size];
    double *R_n=new double[size];
    
    Policy *A = new Policy(size,T_0, T_1);
    
    solve_minimal_constrained();
    
    if(P_m>=Pi->get_epsilon()){
        cout<<"Minimal constrained value is P_m = "<<P_m<<" >= "<<Pi->get_epsilon()<<" = epsilon "<<endl;
    }
    
    double lambda=0;
    lambda_s=0;
    lambda_f=(E_m-E_0)/(Pi->get_epsilon()-P_m);
    
    int counter=0;
    
    double P;
    double E;
    
    //Fix this
    E_f=Pi->get_psi(0);
    P_f=0;
    
    while (lambda_f-lambda_s>Delta) {
        
        //See line 4 of Algorithm 1
        for (int x=0; x<size; x++) {
            V_c[x]=U[x]+lambda;
            V_n[x]=U[x]+lambda;
            R_c[x]=1;
            R_n[x]=1;
            
            if( U[x]>Pi->get_psi(x)){
                A->set_S(1, x, T_1);
            } else {
                A->set_S(1, x, T_1+1);
            }
            A->set_S(0, x, T_0[x]+1);
            A->set_A(0, x, 1);
            A->set_A(1, x, 1);
            
        }
        
        for (int t=T_1-1; t>=0; t--) {
            
            for (int x=0; x<size; x++) {
                double C=Pi->chi(x,t);
                bool a=C;
                V_n[x]=Pi->get_k()+Pi->M(V_c,x);
                
                if (V_n[x]>Pi->get_psi(x)+lambda*C) {
                    V_n[x]=Pi->get_psi(x)+lambda*C;
                    R_n[x]=C;
                    A->set_S(a, x, t);
                    
                } else {
                    if (V_n[x]==Pi->get_psi(x)+lambda*C ) {
                        if(a) {
                            R_n[x]=Pi->M(R_c,x);
                        } else {
                            A->set_S(a, x, t);
                            R_n[x]=C;
                        }
                        
                    } else {
                        R_n[x]=Pi->M(R_c,x);
                    }
                }
                
                
            }
            for (int x=0; x<size; x++) {
                V_c[x]=V_n[x];
                R_c[x]=R_n[x];
            }
        }
        
        //Lines 6 and 7 of Algorithm 1
        P=Pi->expectation(R_n);
        E=Pi->expectation(V_n)-lambda*P;
        
        if(lambda==0){
            P_0=P;   
        }
        
        std::cout<<"Iteration: "<<counter<<std::endl;
        std::cout<<"Constraint is: "<<P<<std::endl;
        std::cout<<"Expected cost is: "<<E<<std::endl;
        std::cout<<"lambda is: "<<lambda<<std::endl;
         
        //Lines 8 - 11
        if (P>Pi->get_epsilon()) {
            lambda_s=lambda;
            P_s=P;
            E_s=E;
            A_s->set_equal(size, A);
            
            for (int x=0; x<size; x++) {
                R_s[x]=R_n[x];
            }
            
        } else {
            lambda_f=lambda;
            P_f=P;
            E_f=E;
            A_f->set_equal(size, A);
            
            for (int x=0; x<size; x++) {
                R_f[x]=R_n[x];
            }
            
        }
        lambda=(lambda_s+lambda_f)/2; //Line 12
        counter++;
    }
    
    //Lines 14-21
    if( P_f < Pi->get_epsilon() && lambda_f>0){
        
        this->resolve_degeneracies_forward();
        
        if(P_flat < Pi->get_epsilon()){
            this->resolve_degeneracies_backward();
        } else {
            
            P_sharp=P_flat;
            E_sharp=E_flat;
            A_sharp->set_equal(size, A_flat);
            lambda_sharp=lambda_flat;
        }
        
    } else{
        //In this case P_s,E_s,A_s,lambda_s have not yet been initialized them so we set them as the optimal solution.
        P_s=P_f;
        E_s=E_f;
        A_s->set_equal(size,A_f);
        lambda_s=lambda_f;
        for (int x=0; x<size; x++) {
            R_s[x]=R_f[x];
        }
        
        P_flat=P_f;
        E_flat=E_f;
        A_flat->set_equal(size,A_f);
        lambda_flat=lambda_f;
        
        P_sharp=P_f;
        E_sharp=E_f;
        A_sharp->set_equal(size,A_f);
        lambda_sharp=lambda_f;
    }
    
    //Performs the optimality check
    for (int t=T_1; t>=0; t--) {
        
        for (int x=0; x<size; x++) {
            if(t<T_1){
 
                V_n[x]=min(Pi->get_k()+Pi->M(V_c,x), Pi->get_psi(x)+lambda_sharp*Pi->chi(x,t));
            } else {
                V_n[x]=U[x]+lambda_sharp;
            }
            
        }
        for (int x=0; x<size; x++) {
            V_c[x]=V_n[x];
        }
    }
    
    double D=Pi->expectation(V_n);
    
    cout<<"E_sharp = "<<E_sharp<<", P_sharp  = "<<P_sharp<<", E_* >= "<<D-lambda_sharp*Pi->get_epsilon()<<endl;
    
    delete [] V_c;
    delete [] V_n;
    
    delete [] R_c;
    delete [] R_n;
    
    delete A;
    
}

void PCOS::resolve_degeneracies_forward(){
    int size=Pi->get_size();
    
    //Step 1: compute \tilde{D}_0 and the values of R^f and Z^f
    
    int D_size_0=0;
    
    double *V_c=new double[size];
    double *V_n=new double[size];
    
    double *R_c=new double[size];
    double *R_n=new double[size];
    
    //Line 1
    for( int x=0; x<size; x++){
        D_size_0 += A_s->get_S(0,x)-A_f->get_S(0, x);
        
        R_c[x]=1;
        R_n[x]=1;
        
        V_c[x]=U[x]+lambda_f;
        V_n[x]=U[x]+lambda_f;
    }
    
    double *MR_f = new double[D_size_0];
    double *MZ_f = new double[D_size_0];
    
    //Line 2
    Policy *A = new Policy(size, T_0, T_1);
    
    A->set_equal(size, A_f);
    
    int d_count_0=D_size_0-1;
    
    //Line 3
    //We store the degenerate points in the reverse order they are processed
    for( int t=T_1; t>=0; t--) {
        for (int x=size-1; x>=0; x--) {
            double C=Pi->chi(x,t);
            V_n[x]=Pi->get_k()+Pi->M(V_c,x);
            
            if (V_n[x]>Pi->get_psi(x)+lambda_f*C) {
                V_n[x]=Pi->get_psi(x)+lambda_f*C;
                R_n[x]=C;
                
            } else {
                if (V_n[x]==Pi->get_psi(x)+lambda_f*C) {
                    R_n[x]=min(Pi->M(R_c,x),C);
                    
                } else {
                    R_n[x]=Pi->M(R_c,x);
                }
            }
            
            if ( t>= A_f->get_S(0,x) && t < A_s->get_S(0,x) ) {
                MR_f[d_count_0] = Pi->M(R_c,x);
                MZ_f[d_count_0] = Pi->M(V_c,x)-lambda_f*MR_f[d_count_0];
                d_count_0--;
            }
            
        }
        for (int x=size-1; x>=0; x--) {
            V_c[x]=V_n[x];
            R_c[x]=R_n[x];
        }
    }
    
    delete [] V_c;
    delete [] V_n;
    
    delete [] R_c;
    delete [] R_n;
    
    //step 2 solve for phi forward in time, resolving the degenerate points when t<=T_0
    
    double *phi_c=new double[size];

    double *phi_n=new double[size];
    
    double P=P_f;
    
    double E=E_f;
    
    double A_temp;
    
    
    
    d_count_0=0;
    
    bool resolved=false;
    int t=0;
    
    while(t<T_1 && !resolved) {
        
        int x=0;
        
        while(x<size && !resolved) {

            bool a=Pi->chi(x,t-1);
            
            //Lines 8-9
            if (t>0) {
                phi_n[x]=(1-A->get_A(a, x, t-1))*phi_c[x]*(1-Pi->get_p(x));
                int *n=Pi->get_neighbors(x);
                int size_n=Pi->get_neighbors_size(x);
                for (int i=0; i<size_n; i++) {
                    phi_n[x]+=(1-A->get_A(a, n[i], t-1))*phi_c[n[i]]*Pi->get_p(n[i])/Pi->get_N(n[i]);
                }
            } else {
                phi_n[x]=Pi->get_phi_0(x);
                phi_c[x]=Pi->get_phi_0(x);
            }
            
            a=Pi->chi(x,t);
            
            if (a==0 && t < A_s->get_S(0,x) && A->get_S(0,x)==t) {

                //Find A_temp \in [0,1] to maximize A_temp*(1-R_temp) if a=1 or A_temp*R_temp if a=0 subject to P+phi[x]*A_temp*(R_temp-1)<=epsilon
                
                //cout<<"Resolving at x="<<x<< "and t="<<t<<endl;
                
                if(phi_n[x]>0){

                    if (MZ_f[d_count_0]+Pi->get_k() <= Pi->get_psi(x)) {
                        double PMR=phi_n[x]*MR_f[d_count_0];
                        double PMZ=phi_n[x]*MZ_f[d_count_0];
                        
                        if (PMR>0) {
                            A_temp=max((P-Pi->get_epsilon())/PMR+1,0.);

                            P=P+(1-A_temp)*PMR;
                    
                        } else {
                            A_temp=0;
                        }
                        
                        E=E+(1-A_temp)*(PMZ+phi_n[x]*(Pi->get_k()-Pi->get_psi(x)));
    
                        
                        if (P>=Pi->get_epsilon()) {
                            A->set_A(0, x, A_temp);
                            
                            cout<<"Resolved at x = "<<x+1<<" and t = "<<t<<" with A = "<<A_temp<<endl;
                    
                            resolved=true;
                        } else {
                            A->set_S(0, x, t+1); A->set_A(0, x, 1);                        }
                    }
                    
                } else {
                    A->set_S(0,x,t+1); A->set_A(0, x, 1);
                }
                d_count_0++;
            }
            x++;
        }
        for (int x=0; x<size; x++) {
            phi_c[x]=phi_n[x];
        }
        t++;
    }
    
    A_flat->set_equal(size, A);
    E_flat=E;
    P_flat=P;

    lambda_flat=(E_f+lambda_f*P_f-E_flat)/P_flat;

    
    cout<<"E_flat = "<<E_flat<<", P_flat = "<<P_flat<<", lambda_flat = "<<lambda_flat<<endl;
    
    delete []MR_f;
    delete []MZ_f;
    
    delete A;
    delete [] phi_c;
    delete [] phi_n;
}


void PCOS::resolve_degeneracies_backward(){
    int size=Pi->get_size();

    double *Z_c=new double[size];
    double *Z_n= new double [size];
    
    double *R_c=new double[size];
    double *R_n= new double [size];
    
    double *phi_c=new double[size];
    double *phi_n= new double[size];
    
    int D_size_1=0;
    
    //Line 1
    for(int x=0; x<size; x++){
        D_size_1 += A_flat->get_S(1,x)-A_s->get_S(1, x);
        
        phi_c[x]=Pi->get_phi_0(x);
        phi_n[x]=phi_c[x];
    }
    
    //Line 2
    Policy *A = new Policy(size,T_0,T_1);
    A->set_equal(size,A_flat);
    
    double *phiD=new double[D_size_1];
    
    int d_count_1=D_size_1-1;
    
    //Line 3
    for(int t=0; t<=T_1; t++){
        for(int x=0; x<size; x++){
            
            bool a;
            
            if (t>0) {
                a=Pi->chi(x,t-1);

                phi_n[x]=(1-A->get_A(a, x, t-1))*phi_c[x]*(1-Pi->get_p(x));
                int *n=Pi->get_neighbors(x);
                int size_n=Pi->get_neighbors_size(x);
                for (int i=0; i<size_n; i++) {
                    phi_n[x]+=(1-A->get_A(a, n[i], t-1))*phi_c[n[i]]*Pi->get_p(n[i])/Pi->get_N(n[i]);
                }
            }
            
            a=Pi->chi(x,t);
            
            if( t>=A_s->get_S(1,x) && t<A->get_S(1,x)){
                phiD[d_count_1]=phi_n[x];
                d_count_1--;
            }
        }
        
        for(int x=0; x<size; x++){
            phi_c[x]=phi_n[x];
        }
    }
    
    double P=P_flat; double E=E_flat;
    
    int t=T_1-1;
    d_count_1=0;
    
    bool resolved = false;
    
    while(t>=0 && !resolved) {
        
        int x=size-1;
        
        while(x>=0 && !resolved) {
            
            if ( t>= A_s->get_S(1,x) && A->get_S(1,x)==t+1 ) {
                
                //cout<<"Resolving A_1 at x = "<<x<<" and t = "<<t<<endl;
                
                double phi=phiD[d_count_1];
                
                
                
                if(phi>0){
                    
                    if (Pi->M(Z_c,x)+Pi->get_k() - Pi->get_psi(x)>= 0) {
                        A->set_S(1,x,t);
                        double A_temp=1;
                        if(phi*(1-Pi->M(R_c,x))>0) {
                            A_temp=min((Pi->get_epsilon()-P)/(phi*(1-Pi->M(R_c,x))),1.);
                        }
                        
                        P=P+A_temp*phi*(1-Pi->M(R_c,x));
                            
                        E=E+A_temp*phi*(Pi->get_psi(x)-Pi->M(Z_c,x)-Pi->get_k());
                            
                        A->set_A(1,x,A_temp);
                        
                        if (P>=Pi->get_epsilon()) {
                            
                            cout<<"Resolved at x = "<<x+1<<" and t = "<<t<<" with A = "<<A_temp<<endl;

                            
                            resolved=true;
                           
                        }
                        
                    }
                    
                } else {
                    A->set_S(1,x,t);
                    A->set_A(1,x,1);
                }
                
                d_count_1++;
            }
            
            bool a=Pi->chi(x,t);
            double A_temp = A->get_A(a, x, t);
            
            R_n[x]=(1-A_temp)*Pi->M(R_c,x)+A_temp*Pi->chi(x,t);
            
            Z_n[x]=(1-A_temp)*(Pi->M(Z_c,x)+Pi->get_k())+A_temp*Pi->get_psi(x);
            
            x--;
        }
        for(x=0; x<size;x++){
            R_c[x]=R_n[x];
            Z_c[x]=Z_n[x];
        }
        t--;
    }
    
    A_sharp->set_equal(size, A);
    E_sharp=E;
    P_sharp=P;

    lambda_sharp=(E_f+lambda_f*P_f-E_sharp)/P_sharp;
    
    
    delete []Z_c;
    delete []Z_n;
    
    delete []R_c;
    delete [] R_n;
    
    delete [] phiD;
    
    delete [] phi_c;
    delete [] phi_n;
    
    delete A;
}

double PCOS::get_lambda_f(){
    return lambda_f;
}
double PCOS::get_lambda_s(){
    return lambda_s;
}
double PCOS::get_E_f(){
    return E_f;
}
double PCOS::get_E_s(){
    return E_s;
}

void PCOS::output(string name){
    
    int size=Pi->get_size();
    
    ofstream file;
    stringstream str;
    str<<FILEPATH<<"PCOS_"<<name<<".m";
    string s=str.str();
    
    char const *fileName=s.c_str();
    file.precision(15);
    file.open(fileName, ios::out);
    
    file<<"lambdaf="<<lambda_f<<";"<<endl;
    file<<"lambdas="<<lambda_s<<";"<<endl;
    file<<"Rf=zeros("<<size<<",1);"<<endl;
    file<<"Rs=zeros("<<size<<",1);"<<endl;
    file<<"Phi0=zeros("<<size<<",1);"<<endl;
    file<<"U=zeros("<<size<<",1);"<<endl;
    file<<"S0f=zeros("<<size<<",1);"<<endl;
    file<<"S1f=zeros("<<size<<",1);"<<endl;
    file<<"S0s=zeros("<<size<<",1);"<<endl;
    file<<"S1s=zeros("<<size<<",1);"<<endl;

    for (int x=0; x<size; x++) {
        file<<"Rf("<<x+1<<")="<<R_f[x]<<"; Rs("<<x+1<<")="<<R_s[x]<<";"<<endl;
        file<<"Phi0("<<x+1<<")="<<Pi->get_phi_0(x)<<"; U("<<x+1<<")="<<U[x]<<";"<<endl;
        file<<"S0f("<<x+1<<")="<<A_f->get_S(0, x)<<"; S1f("<<x+1<<")="<<A_f->get_S(1, x)<<";"<<endl;
        file<<"S0s("<<x+1<<")="<<A_s->get_S(0, x)<<"; S1s("<<x+1<<")="<<A_s->get_S(1, x)<<";"<<endl;

    }
}


