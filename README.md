# Generalized-Monogamy-of-Entanglement-and-Quantification-of-Bipartite-Entanglement
This repository contains all the code needed to generate the figures in the article Generalized Monogamy of Entanglement and Quantification of Bipartite Entanglement. 
The file Cbounds3qubits.m samples the random three-qubit pure states used to generate the data points in Figure 2.
The file NegativityBounds.m samples the random 2+N states and classifies them depending on whether any of the eigenvalues lambda_i is zero. The data generated with this code is used in Figure 3.
The file BoundariesNeg.m is the numerical implementation of the minimization/maximization of N_{A_1A_2}^{max} described in detail in Appendix A. The data generated with this code is used in Figure 3.
The file 4qubit.m is used to generate Figures 5 and 6. This file simulates the dynamics for a system of 4 qubits (2 in A and 2 in B). Inside the file, the Hamiltonian H_{AB} and H_{A} can be changed to generate the different plots.
The file 8qubit.m is used to generate Figure 7. This file simulates the dynamics for a system of 8 qubits (4 in A and 4 in B). Inside the file, the Hamiltonian H_{AB} and H_{A} can be changed to generate the different plots. In this same code, if the data of the dynamics of subsystem A is chosen to be stored, Figures 8 and 9 can be generated with it.
