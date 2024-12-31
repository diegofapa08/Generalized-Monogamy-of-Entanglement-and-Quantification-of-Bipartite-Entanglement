%First we define the system size, in this case both subsystems A and B have
%4 qubits each

N=256;

Nsamp = 150;%Number of time steps
tsim = linspace(0,pi/2,Nsamp); %Time range
Results = zeros(Nsamp, 13); %We will store, SAB, the total Spin Squeezing in the optimal direction, and the optimal perp direction vector (see Eq16)

%We define the spin operators in the product basis

Sz = 1/2*[1,0;0,-1];
Sx = 1/2*[0,1;1,0];
Sy = 1/2*[0,-1i;1i,0];

SxT = kron(Sx,eye(128)) + kron(kron(eye(2),Sx),eye(64)) + kron(kron(eye(4),Sx),eye(32)) +  kron(kron(eye(8),Sx),eye(16))  + kron(kron(eye(16),Sx),eye(8))+  kron(kron(eye(32),Sx),eye(4)) + kron(kron(eye(64),Sx),eye(2))+ kron(eye(128),Sx);
SyT = kron(Sy,eye(128)) + kron(kron(eye(2),Sy),eye(64)) + kron(kron(eye(4),Sy),eye(32)) +  kron(kron(eye(8),Sy),eye(16))  + kron(kron(eye(16),Sy),eye(8))+  kron(kron(eye(32),Sy),eye(4)) + kron(kron(eye(64),Sy),eye(2))+ kron(eye(128),Sy);
SzT = kron(Sz,eye(128)) + kron(kron(eye(2),Sz),eye(64)) + kron(kron(eye(4),Sz),eye(32)) +  kron(kron(eye(8),Sz),eye(16))  + kron(kron(eye(16),Sz),eye(8))+  kron(kron(eye(32),Sz),eye(4)) + kron(kron(eye(64),Sz),eye(2))+ kron(eye(128),Sz);

%We define the operators for subsystem A

SzA = kron(Sz,eye(8)) + kron(kron(eye(2),Sz),eye(4))+kron(kron(eye(4),Sz),eye(2))+kron(eye(8),Sz) ;
SxA =  kron(Sx,eye(8)) + kron(kron(eye(2),Sx),eye(4))+kron(kron(eye(4),Sx),eye(2))+kron(eye(8),Sx);
SyA =  kron(Sy,eye(8)) + kron(kron(eye(2),Sy),eye(4))+kron(kron(eye(4),Sy),eye(2))+kron(eye(8),Sy);

%We define the Hamiltonians we will use to evolve the total system and
%subsystem A

%HamSqA = SxA^2-SyA^2; %Hamiltonian for subsytem A
HamSqA = SxA^2+SzA; %TF for subsystem A
%HamSqA = SxA^2; %OAT for subsystem A
HamSqT = SxT^2; %Hamiltonian for entanglement of A and B

%We diagonalize the Hamiltonians to perform the unitary evolution

[VSqA,DSqA] = eig(HamSqA,'vector');
[VSqT,DSqT] = eig(HamSqT,'vector');

[DSqA,ind] = sort(real(DSqA));
DnA = DSqA;
VSqA = VSqA(:,ind);

[DSqT,ind] = sort(real(DSqT));
DnT = DSqT;
VSqT = VSqT(:,ind);


for n=1:Nsamp   
    %First we initialize the system in a product state with all spins down 
    spdown = [0;1];
    spup = [1;0];
    state0 = spdown;
    for jj=1:7
        state0 = kron(state0,spdown);
    end
    RhoD = state0*conj(transpose(state0));
    %Next we evolve the system using the Hamiltonian evolution, for that we
    %can simply change basis to the basis of Hamiltonian eigenstates
    vin0t = VSqT\RhoD*VSqT;
    tevol =tsim(n);
    vint = zeros(16,16);

    for jj=1:256
        for ll=1:256
            vint(jj,ll) = exp(-1i*(DnT(jj)-DnT(ll))*tevol)*vin0t(jj,ll);
        end
    end

    %We transform back to the original basis

    Rho = VSqT*vint/VSqT;
    

   %We get the reduced density matrix first
   RedRho = zeros(16,16);

   for jj = 1:16
       psi0 = zeros(1,16);
       psi0(1,jj) = 1;
       RedRho = RedRho + kron(eye(16),psi0)*Rho*kron(eye(16),ctranspose(psi0));
   end
   %Then the entanglement entropy is computed using the eigenvalues

   eigenred = eig(RedRho);
   SAB = 0;
   for jj=1:16
       if eigenred(jj) > 0
           SAB = SAB - eigenred(jj)*log(eigenred(jj));
       end
   end

   Results(n,1) = SAB;%Entanglement entropy
   Results(n,13) = 16/15*(1-trace(RedRho^2)); %Linear Entropy Eq 4

   %Now we compute the average spin direction and find the optimum
   %squeezing, using the details on Ref[29] about spin squeezing

   AvSpin = [real(trace(Rho*SxT)),real(trace(Rho*SyT)),real(trace(Rho*SzT))];
   AvSpin = AvSpin/norm(AvSpin);

   %Now we express the average direction in spherical coordinates

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
   S1 = n1(1)*SxT + n1(2)*SyT + n1(3)*SzT;
   S2 = n2(1)*SxT + n2(2)*SyT + n2(3)*SzT;

   Acoeff = trace(Rho*(S1^2-S2^2));
   Bcoeff = trace(Rho*(S1*S2+S2*S1));

   if Bcoeff <= 0
       varphi = 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
   else
       varphi = pi - 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
   end

   Results(n,2) = cos(varphi)*n1(1) + sin(varphi)*n2(1);
   Results(n,3) = cos(varphi)*n1(2) + sin(varphi)*n2(2);
   Results(n,4) = cos(varphi)*n1(3) + sin(varphi)*n2(3);   

   SOpt = cos(varphi)*S1 + sin(varphi)*S2;
   VarOpt = trace(SOpt^2*Rho) - trace(SOpt*Rho)^2;

   AvSpin = [real(trace(Rho*SxT)),real(trace(Rho*SyT)),real(trace(Rho*SzT))];
   Results(n,5) = 4*VarOpt/8;
   %Results(n,5) = VarOpt;
   Results(n,12) = tevol;

   %Let's also keep track of the squeezing in x and y directions
   Results(n,6) = 4*(trace(Rho*SxT*SxT)-trace(Rho*SxT)^2)/(trace(Rho*SzT)^2);
   Results(n,7) = 4*(trace(Rho*SyT*SyT)-trace(Rho*SyT)^2)/(trace(Rho*SzT)^2);

   %Now we want to see how much we can squeeze individually subsystem A,

   vin0tA = VSqA\RedRho*VSqA;
   vintA =zeros(16,16);
   VarEvol = zeros(100,4);
   indE = 1;
   for t=0:0.01:100 %We evolve each state obtained initially for this time and keep track of the minimal squeezing we can get
    for jj=1:16
        for ll=1:16
            vintA(jj,ll) = exp(-1i*(DnA(jj)-DnA(ll))*t)*vin0tA(jj,ll);
        end
    end
    rhoA = VSqA*vintA/VSqA;

   %Now we compute the average spin direction and find the optimum
   %squeezing, using the details on Ref[29] about spin squeezing

   AvSpin = [real(trace(rhoA*SxA)),real(trace(rhoA*SyA)),real(trace(rhoA*SzA))];
   AvSpin = AvSpin/norm(AvSpin);

   %Now we express the average direction in spherical coordinates

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
   S1 = n1(1)*SxA + n1(2)*SyA + n1(3)*SzA;
   S2 = n2(1)*SxA + n2(2)*SyA + n2(3)*SzA;

   Acoeff = trace(rhoA*(S1^2-S2^2));
   Bcoeff = trace(rhoA*(S1*S2+S2*S1));

   if Bcoeff <= 0
       varphi = 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
   else
       varphi = pi - 1/2*acos(-Acoeff/sqrt(Acoeff^2+Bcoeff^2));
   end

   SOpt = cos(varphi)*S1 + sin(varphi)*S2;
   VarOpt = trace(SOpt^2*rhoA) - trace(SOpt*rhoA)^2;
    
    VarEvol(indE,1) = 4*VarOpt/4;
    VarEvol(indE,2) =t;
    indE = indE +1;
   end

   Results(n,8) = min(VarEvol(:,1));
   n
end

%% We signal the data where the one-to-one relationship between linear entropy and minimum squeezing in A breaks
zoneS = zeros(1,5);
indS=1;
for jj=1:length(Results)
    if Results(jj,12) >= 0.8 %range of times where the one-to-one relationship breaks
        if Results(jj,12) <= pi-0.8
            zoneS(indS,1) = Results(jj,13);
            zoneS(indS,2) = Results(jj,5);
            zoneS(indS,3) = Results(jj,8);
            zoneS(indS,4) = Results(jj,12);
            zoneS(indS,5) = Results(jj,9);
            indS=indS+1;
        end
    end
end
%% Stores the generated data
%save('8qubitTF.txt','Results','-ascii')


%% Generates plot equivalent to Fig 7
plot(Results(:,12),Results(:,13),'LineWidth',2,'Color', [71 194 109]/255)
set(gca,'Fontsize',14)
hold on
plot(Results(:,12),Results(:,5),'LineWidth',2,'Color', [39 127 142]/255,'LineStyle',':')
plot(Results(:,12),Results(:,8),'LineWidth',2,'Color', [68 1 84]/255,'LineStyle','--')
scatter(zoneS(:,4),zoneS(:,3),'filled','MarkerFaceColor',[0 0 0])
scatter(zoneS(:,4),zoneS(:,1),'filled','MarkerFaceColor',[0 0 0])
xlabel('$\Omega t$','Interpreter','latex','FontSize',20)
legend('$S_{AB}$','$\xi^2_{KU,AB}$','$\min(\xi^2_{KU,A})$','Interpreter','latex','Fontsize',20,'location','southeast','box','off')
ylim([0,1.3])
hold off
%exportgraphics(gca,'8qubitTevol.jpg','Resolution',300)


