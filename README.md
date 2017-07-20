# PCOS
Numerical Experiments for Probabilistically Constrained Optimal Stopping Problems

This code is provided as a supplement to the paper titled "Optimal Stopping with a Probabilistic Constraint" by Aaron Zeff Palmer and Alexander Vladimirsky, submitted to the Journal of Optimization Theory and Applications in 2017.  A link will be provided to the journal upon acceptance.

All the code was written by Aaron Zeff Palmer. Any questions on the content of the paper may be addressed to the corresponding author, Aaron Zeff Palmer, at azp6@cornell.edu. Please do not contact with requests to extend this code or any solicitations.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Abstract
We present an efficient method for optimal stopping problems with a probabilistic constraint. The goal is to optimize the expected cumulative cost, but with the prescribed upper bound on the probability that this cost exceeds  a threshold value.  This probabilistic constraint causes optimal policies to be time-dependent and randomized, however,  we show that an optimal policy can always be selected with "piecewise-monotonic" time-dependence and "nearly-deterministic" randomization.  We prove these properties using the Bellman optimality equations for a Lagrangian relaxation of the original problem. Our algorithm exploits these properties for computational efficiency.  We illustrate its performance and the structure of optimal policies on two numerical examples. 

# To Compile
To compile the code, edit the first line of 'Makefile' to define a path to a directory where the output is saved.  Alternatively, FILEPATH can defined at the top of PCOS.hpp

On some machines it will be sufficient to type 'make' in the Unix terminal at the directory with the source files and 'Makefile'.  This code was developed using Xcode 8.3.3 on Mac OS Sierra 10.12.5 and was also tested on Linux Ubuntu.   There are no dependencies besides basic input/output (\<iostream>, \<fstream>,\<sstream>) and \<cmath>, which are included in PCOS.hpp.

# To Run
The program is executed in Unix using either 

'./PCOS.out ex1' 

to run example 1, 

'./PCOS.out ex2' 

to run example 2, or 

'./PCOS.out lambdaepsilon N'

to run example 1 for a series of epsilon values, where N is an integer that is half the size of the problem. (N=5 and N=25 are included in the paper).

The output are files labeled 'PCOS_ex1.m', 'PCOS_ex2.m', or 'PCOS_lambdaepsilon_N.m' that can be read by Matlab. The Matlab functions are provided:

'[lambdaf, lambdas,Rf,Rs,Phi0,U,S0f,S1f,S0s,S1s,x,t] = plot_example(num)'

and

'[epsilon,lambda,E] = plot_lambdaepsilon(N).'

These functions can be run after the directory containing the files 'plot_example.m', 'plot_lambdaepsilon.m', and the output of PCOS is added to the path for Matlab. They reproduce the plots appearing in the paper.

# To Experiment
The parameters of the examples can be altered in main.cpp. For any other experimentation with the algorithm, the essential code is in PCOS.cpp.
