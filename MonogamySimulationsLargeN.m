%This code can be used to simulate larger number of qubits exploiting the
%symmetries of the system. It is important that we will assume that we
%always initialize the system with all qubits pointing down or all qubits
%pointing up, that way we can reduce the size of the Hilbert space
%considerably for the Hamiltonians of interest.

%Any questions of the code can be sent to Diego Fallas Paddill at difa1788@colorado.edu
%To reference this code use the doi in the repository where you downloaded
%it from.

%First we define the number of total qubits, we assume half the qubits are
%in subsystem A and half in B

N=100; %Number of qubits
dim = N+1; %Dimension of the Hilbert space under OAT evolution
dimA =N/2 +1;  %Dimension of subsystem A

[T]= CGTransformationDicke(N); %This defines a transformation matrix that maps states in the |J M_J> representation of total spin
% to states in the |JA M_JA JB M_JB> representation.

%We define the initial state of all evolutions which is all atoms pointing
%down, this correspond to the Dicke state |N/2,-N/2>

state0 = sparse(dim,1);
state0(dim,1)=1;

%Now we define the Hamiltonian in this basis

diagSm = zeros(dim-1,1);
for n=0:dim-2
    m=N/2-n;
    diagSm(n+1) = sqrt((N/2)*(N/2+1)-(m-1)*m);
end

diagSz = zeros(dim,1);
for n=1:dim
    m = N/2 +1 -n;
    diagSz(n) = m;
end

Sm = spdiags(diagSm,-1,dim,dim);
Sp = Sm';
Sz = spdiags(diagSz,0,dim,dim);
Sx = (Sm+Sp)/2;
Sy = -1i*(Sp-Sm)/2;
SxSq = (Sp^2 + Sm^2 + Sp*Sm + Sm*Sp)/4;
HamSqT = SxSq; %Hamiltonian for entanglement of A and B H = Sx^2

%%Now we compute the Hamiltonian terms for the evolution of subsystem A in
%%the basis |JA M_JA>

diagSmA = zeros(dimA-1,1);
for n=0:dimA-2
    m=N/4-n;
    diagSmA(n+1) = sqrt((N/4)*(N/4+1)-(m-1)*m);
end

diagSzA = zeros(dimA,1);
for n=1:dimA
    m = N/4 +1 -n;
    diagSzA(n) = m;
end

SmA = spdiags(diagSmA,-1,dimA,dimA);
SpA = SmA';
SzA = spdiags(diagSzA,0,dimA,dimA);
SxA = (SmA+SpA)/2;
SyA = -1i*(SpA-SmA)/2;
 
 HamA = SxA^2+SzA;  %Use this for H_A = H_{TF}
% HamA = SxA^2; %Use this for H_A = H_{OAT}

%Next we compute the diagonalizations of both Hamiltonians. For this, we
%assume that we have enough space in memory to diagonalize the Hamiltonians
%as dense matrices. If larger N is considered, the evolution can be
%computed using matrix exponentiation techniques for sparse matrices such
%as expv or expvm.

[VSqT,DSqT] = eig(full(HamSqT),'vector');

[DnT,ind] = sort(real(DSqT));
DSqT = DnT;
VSqT = VSqT(:,ind);

[VA,DA] = eig(full(HamA),'vector');

[DnTA,ind] = sort(real(DA));
DA = DnTA;
VA = VA(:,ind);


%Since OAT evolution is periodic and we take \Omega=1, it is enough for us
%to sample the evolution between 0 and pi/2

Nsamp = 300; %Number of time steps
tsim = linspace(0,pi/2,Nsamp); %Total simulation time
deltat = tsim(2)-tsim(1); %Delta t

datatime = zeros(Nsamp,2); 
ind1=1;

dimj2=0;
for j2=N/4
    for m2=-j2:j2
        dimj2=dimj2+1;
    end
end

psit=state0;

for t = tsim
    datatime(ind1,1)=t;
    datatime(ind1,2)=psit'*Sz*psit;
    rho_coupled = T*psit; %Here we change from the Dicke basis of the total subsystem to the Dicke basis of subsystem A.
    rho_coupled = rho_coupled*rho_coupled';

    redrho = sparse(dimj2,dimj2);
    tic
    for jj=1:dimj2
        stateTr = sparse(dimj2,1);
        stateTr(jj) = 1;
        redrho = redrho + (kron(speye(dimj2),stateTr))'*rho_coupled*(kron(speye(dimj2),stateTr));
    end

    datatime(ind1,3) = (2^(N/2))*(1-trace(redrho^2))/(2^(N/2)-1); %Linearized entanglement entropy S_AB
    Rho = psit*psit';
    
    %Next we compute the spin squeezing. To find the minimum variance
    %direction, we use the results presented in Ma, J., Wang, X., Sun, C. P., & Nori, F. (2011). Quantum spin squeezing. Physics Reports, 509(2-3), 89-165.

    %We express the average direction in spherical coordinates

    AvSpin = [real(trace(Rho*Sx)),real(trace(Rho*Sy)),real(trace(Rho*Sz))];
    AvSpin = AvSpin/norm(AvSpin);
    
    
    AvTet = acos(AvSpin(3));
    
    if AvSpin(1) > 0
       if AvSpin(2) > 0
           AvPhi = atan(AvSpin(2)/AvSpin(1));
       else
           AvPhi = 2*pi-atan(-AvSpin(2)/AvSpin(1));
       end
    else
       if AvSpin(2) > 0
           AvPhi = pi-atan(AvSpin(2)/-AvSpin(1));
       else
           AvPhi = atan(-AvSpin(2)/-AvSpin(1))+pi;
       end
    end
    
    n1 = [-sin(AvPhi),cos(AvPhi),0];
    n2 = [cos(AvPhi)*cos(AvTet),cos(AvTet)*sin(AvPhi),-sin(AvTet)];
    S1 = n1(1)*Sx + n1(2)*Sy + n1(3)*Sz;
    S2 = n2(1)*Sx + n2(2)*Sy + n2(3)*Sz;
    
    Acoeff = trace(Rho*(S1^2-S2^2));
    Bcoeff = trace(Rho*(S1*S2+S2*S1));
    
    if Bcoeff <= 0
       varphi = 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
    else
       varphi = pi - 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
    end
    
    SOpt = cos(varphi)*S1 + sin(varphi)*S2;
    VarOpt = trace(SOpt^2*Rho) - trace(SOpt*Rho)^2;
    
    AvSpin = [real(trace(Rho*Sx)),real(trace(Rho*Sy)),real(trace(Rho*Sz))];
    datatime(ind1,4) = 4*VarOpt/N; %This is \xi_AB (the Kitagawa Ueda squeezing parameter for AB)
    
    %Now we compute the evolution of subsystem A, since we want to find
    %the minimal spin squeezing in A we use optimization methods.

    InitialNum = 500; %Number of initial conditions per evolution time of AB
    
    VarEvol = zeros(InitialNum,2);
   
   for n = 1:InitialNum 
    incond = 5000*rand; %We choose the initial time, in this case a random number between 0 and 5000 \Omega t
    squeeze_func = @(t) SqueezingA(t,redrho,N,VA,DA,SxA,SyA,SzA); %We define the evolution under the squeezing Hamiltonian H_A here
    options = optimoptions('fmincon', 'Algorithm', 'sqp','Display','off');
    [t_min, min_squeezing] = fmincon(squeeze_func, incond,[],[],[],[],0,5000,[],options); %We force a bound in \Omega t, in this case to be positive and less than 5000
    VarEvol(n,1) =t_min; %Minimization time
    VarEvol(n,2) =min_squeezing; %Minimal spin squeezing for A (Kitagawa and Ueda)
    n;
   end
   [minvar,indminvar] = min(VarEvol(:,2)); %We find the minimum squeezing of A over all initial conditions.
   datatime(ind1,5) = VarEvol(indminvar,1); %Minimization time of the squeezing on A
   datatime(ind1,6) = minvar;%Minimal squeezing of A
   psit = UnitaryEvolution(psit,deltat,VSqT,DSqT); %Note that we update the state such that we always evolve only by deltat, this is not too relevant in the case of exact diagonalization, but its a good practice for exponentiation.
   ind1=ind1+1
end
%% This would generate the plot of min squeezing of A versus S_AB
scatter(datatime(1:end,6),datatime(1:end,3))




