N=8; %Dimension of the Hilbert Space

Nsamp = 3000; %Number of states

Results=zeros(Nsamp,4);

for jj=1:Nsamp
    psi1 = [1,0,0,0,0,0,0,0];
    psi2 = [0,1,0,0,0,0,0,0];
    psi3 = [0,0,1,0,0,0,0,0];
    psi4 = [0,0,0,1,0,0,0,0];
    psi5 = [0,0,0,0,1,0,0,0];
    psi6 = [0,0,0,0,0,1,0,0];
    psi7 = [0,0,0,0,0,0,1,0];
    psi8 = [0,0,0,0,0,0,0,1];
    
    %The random numbers p_1 are used to have different ranks in the
    %generated density matrix

    p1 = rand;
    p2 = rand;
    p3 = rand;
    p4 = rand;
    p5 = rand;
    p6 = rand;
    p7 = rand;
    p8 = rand;
    
    if p1 < 0.5
        p1 = 0;
    else
        p1 = rand+1i*rand;
    end

    if p2 < 0.5
        p2 = 0;
    else
        p2 = rand+1i*rand;
    end

    if p3 < 0.5
        p3 = 0;
    else
        p3 = rand+1i*rand;
    end

    if p4 < 0.5
        p4 = 0;
    else
        p4 = rand+1i*rand;
    end

    if p5 < 0.5
        p5 = 0;
    else
        p5 = rand+1i*rand;
    end

    if p6 < 0.5
        p6 = 0;
    else
        p6 = rand+1i*rand;
    end

    if p7 < 0.5
        p7 = 0;
    else
        p7 = rand+1i*rand;
    end

    if p8 < 0.5
        p8 = 0;
    else
        p8 = rand+1i*rand;
    end

    psi0 = p1*psi1 +  p2*psi2 +  p3*psi3 +  p4*psi4 +  p5*psi5 +  p6*psi6 +  p7*psi7 +  p8*psi8; 
    psi0 = psi0/norm(psi0); 
    rho0 = ctranspose(psi0)*psi0; %Density Matrix

    %The reduced density matrix for subsystem A is defined below

    rhoA = kron(eye(4),[1,0])*rho0*kron(eye(4),[1;0]) + kron(eye(4),[0,1])*rho0*kron(eye(4),[0;1]);

    lambda1 = max(eig(rhoA)); %We define lambda_1 (Eq 8 in the main text)

    CAB2 = 4*(lambda1 - lambda1^2); %Eq 9 in the main text
    
    CBmax1 = 1/2*(1+sqrt(1-CAB2)); %Eq 10 in the main text
    CBmax2 = 1/2*(1-sqrt(1-CAB2));

    %Now we concurrence as defined in Eq 1

    S = kron([0,-1i;1i,0],[0,-1i;1i,0]);
    A = S*conj(rhoA)*S;
    CMat = rhoA*A;
    eigC = eig(CMat);
    eigC = sort(eigC,'descend');

    CB = max([0,sqrt(eigC(1))-sqrt(eigC(2))-sqrt(eigC(3))-sqrt(eigC(4))]);
    
    Results(jj,1) = CAB2;
    Results(jj,2) = CB;
    Results(jj,3) = CBmax1;
    Results(jj,4) = CBmax2;
end

ResultsOrd = sortrows(Results,1);


%% Generate a plot equivalent to Fig 2
scatter(ResultsOrd(:,1),ResultsOrd(:,2),'MarkerFaceColor', [74,194,109]/255, 'MarkerEdgeColor',[74,194,109]/255,'MarkerFaceAlpha',0.7,'SizeData',20)
hold on
plot(ResultsOrd(:,1),ResultsOrd(:,3),'LineWidth',3,'LineStyle','--','Color',[0 0 0])
set(gca,'Fontsize',14)
xlabel('$C_{AB}$','Interpreter','latex','FontSize',20)
ylabel('$C_{A_1A_2}$','Interpreter','latex','FontSize',20)
legend('Random states','$C_{A_1A_2}^{max}$ (analytic)','Interpreter','latex','box','off', 'fontsize',17)
ylim([0,1.1])
hold off
exportgraphics(gca,'ConcurrenceBound.pdf','contenttype','image','Resolution',300)

