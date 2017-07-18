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
//  main.cpp
//
//  This file contains the primary examples whose data appears in our paper.
//

#include "PI.hpp"
#include "PCOS.hpp"

void lambdaEpsilonTest(int numTests, double maxEpsilon, int N, int time_steps);

int main(int argc, const char * argv[]) {

    if (argc == 2) {
        
        cout<<FILEPATH<<endl;
        
        double d;
        int N=200;
        int size=2*N-1;
        int time_steps;
        double psi=0.9;
        double pi=1;
        double epsilon=0.02;
        PI *Pi=new PI(size);;
        
        if( string(argv[1])=="ex1") {
            
            cout<<"Initializing Example 1"<<endl;
            d=0.25;
            time_steps=100000;
            
            Pi->initialize_phi0_delta(N-1);
            
        } else if( string(argv[1])=="ex2") {
            
            cout<<"Initializing Example 2"<<endl;
            
            d=0.05;
            time_steps=20000;
            
            Pi->initialize_phi0_uniform();
            
        } else {
            cerr << "Usage: " << argv[0] << " e1 OR e2 OR lambdaepsilon" << endl;
            
            delete Pi;
            return 1;
        }
        
        double k=1./(double)time_steps;
        double p=2*d*(2*N)*(2*N)/(double)time_steps;
        
        Pi->initialize_interval(p);
        Pi->initialize_psi_constant(psi);
        
        Pi->set_k(k);
        Pi->set_epsilon(epsilon);
        Pi->set_pi(pi);
        
        cout<<Pi->get_k()<<" "<<Pi->get_psi(0)<<" "<<Pi->get_p(0)<<endl;
        
        PCOS *example=new PCOS(Pi);
    
        example->solve_unconstrained(pow(10.,-10));
        
        example->solve_constrained(pow(10.,-6));
        
        example->output(string(argv[1]));
        
        delete Pi;
        delete example;
        
        return 0;
    }


     if (argc == 3) {
         
         if( string(argv[1])=="lambdaepsilon") {
             
             stringstream str(argv[2]);
             //int N=stoi(argv[2]);
             int N;
             str>>N;
             
             lambdaEpsilonTest(111, 0.12, N, N*N*50);
             return 0;
         }
     }
    
    //Use one of the three commands to run examples 1 or 2 or the calculation of lambda and epsilon
    cerr << "Usage: " << argv[0] << " e1 OR e2 OR lambdaepsilon" << endl;
    
    return 1;
}


void lambdaEpsilonTest(int numTests, double maxEpsilon, int N, int time_steps){
    
    double d=0.25;
    
    int size=2*N-1;
    //int size=N; //Use for half-interval (equivalent problem)
    
    double p=2*d*(2*N)*(2*N)/time_steps;
    
    double psi=0.9;
    double pi=1;
    double k=1./time_steps;
    
    ofstream file;
    stringstream str;
    str<<FILEPATH<<"PCOS_lambdaepsilon_"<<N<<".m";
    string s=str.str();
    
    char const *fileName=s.c_str();
    file.precision(15);
    file.open(fileName, ios::out);
    file<<"epsilon=zeros("<<numTests<<",1);"<<endl;
    file<<"lambda=zeros("<<numTests<<",2);"<<endl;
    file<<"E=zeros("<<numTests<<",2);"<<endl;
    
    
    for (int i=0; i<numTests; i++) {
        double epsilon = (numTests-i)*maxEpsilon/(numTests);
        PI *Pi=new PI(size);
        
        Pi->initialize_phi0_delta(N-1);
        Pi->initialize_interval(p);
        Pi->initialize_psi_constant(psi);
        
        Pi->set_k(k);
        Pi->set_pi(pi);
        Pi->set_epsilon(epsilon);
        
        //Pi->initialize_phi0_N();
        //Pi->initialize_half_interval(p);

        
        PCOS *example=new PCOS(Pi);
        
        example->solve_unconstrained(pow(10.,-8));
        
        example->solve_constrained(pow(10.,-6));
        cout<<"epsilon: "<<epsilon<<", lambda: "<<example->get_lambda_f()<<endl;
        
        file<<"epsilon("<<i+1<<")="<<epsilon<<";"<<endl;
        file<<"lambda("<<i+1<<",1)="<<example->get_lambda_f()<<";"<<endl;
        file<<"lambda("<<i+1<<",2)="<<example->get_lambda_s()<<";"<<endl;
        file<<"E("<<i+1<<",1)="<<example->get_E_f()<<";"<<endl;
        file<<"E("<<i+1<<",2)="<<example->get_E_s()<<";"<<endl;
        
        delete Pi;
        delete example;
    }
    
    
}

