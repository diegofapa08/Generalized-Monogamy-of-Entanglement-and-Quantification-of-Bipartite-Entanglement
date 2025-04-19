function rho_t = UnitaryEvolutionRho(rho0, t, V, D)
%Here rho0 is the initial state density matrix, V is the diagonalization transformation of
%the Hamiltonian and D is a string with the corresponding eigenvalues. t is
%the elapsed time
    expD = diag(exp(-1i * D * t)); 
    U_t = V * expD * V'; 
    rho_t = U_t * rho0 * U_t';
end
