function yt = UnitaryEvolution(y0,t,V,D)
%Here y0 is the initial state, V is the diagonalization transformation of
%the Hamiltonian and D is a string with the corresponding eigenvalues. t is
%the elapsed time
expD = diag(exp(-1i *D*t)); 
U_t = V*expD*V';
yt = U_t*y0;
end