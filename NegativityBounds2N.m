x1 = linspace(1,0.5,200); %We sample the possible values of lambda_1 for 2+1 systems
N12MaxN1 = sqrt(x1.^2 +(1-x1).^2)-(1-x1); %Generates the 2+1 result in Eq.13
NABN1 = 2*sqrt(x1.*(1-x1)); %NAB for 2+1 systems
dataN1 = [NABN1;N12MaxN1]; 

%%2+N systems:

Nsamp = 300; %Number of sample states for each rank 
dataN2 = zeros(11*Nsamp,6);
rankrange = linspace(0,0.99,10); %This is implemented to change the rank of the generated density matrix
ind1=1;
for q = rankrange
for jj=1:Nsamp
    p1 = rand; p2=rand;p3=rand;p4=rand;

    if p1>q
        p1 = 0;
    end
    if p2>q
        p2 = 0;
    end
    if p3>q
        p3 = 0;
    end
    if p4>q
        p4 = 0;
    end
    
    if p1==0 && p2==0 && p3==0 && p4==0
        p1=1;
    end
    par=[p1,p2,p3,p4]/(p1+p2+p3+p4); 
    par = sort(par,'descend'); %Eigenvalues lambda_i of the reduced density matrix in descendent order

    dataN2(ind1,1) = sqrt(par(1)*par(2))+sqrt(par(1)*par(3))+sqrt(par(2)*par(3))+sqrt(par(1)*par(4))+sqrt(par(2)*par(4))+sqrt(par(3)*par(4)); %Eq 15 in the main text
    dataN2(ind1,2) = max([0,sqrt((par(1)-par(3))^2+(par(2)-par(4))^2)-par(2)-par(4)]); %Eq 12 in the main text
    dataN2(ind1,3) = par(1);
    dataN2(ind1,4) = par(2);
    dataN2(ind1,5) = par(3);
    dataN2(ind1,6) = par(4);
    ind1=ind1+1
end
end

for jj=Nsamp*10+1:Nsamp*11 %We add this extra set of data on a range that is hard to sample in the loop above
    p1 = rand; p2=rand;p3=0.1*rand;p4=0.1*rand;
    %p4=0;
    par=[p1,p2,p3,p4]/(p1+p2+p3+p4);
    par = sort(par,'descend');

    dataN2(jj,1) = sqrt(par(1)*par(2))+sqrt(par(1)*par(3))+sqrt(par(2)*par(3))+sqrt(par(1)*par(4))+sqrt(par(2)*par(4))+sqrt(par(3)*par(4));
    dataN2(jj,2) = max([0,sqrt((par(1)-par(3))^2+(par(2)-par(4))^2)-par(2)-par(4)]);
    dataN2(jj,3) = par(1);
    dataN2(jj,4) = par(2);
    dataN2(jj,5) = par(3);
    dataN2(jj,6) = par(4);
end
%% Here we classify the generated states depending on whether one or more eigenvalues lambda_i are zero

dataN2pure = zeros(1,6);
dataN2p1p1 = zeros(1,6);
dataN2p1p2p3=zeros(1,6);
dataN2p1p2p3p4=zeros(1,6);
ind1=1;ind2=1;ind3=1;ind4=1;
for jj=1:3300
    if dataN2(jj,3)==1
        dataN2pure(ind1,:)= dataN2(jj,:);
        ind1=ind1+1;
    end

    if dataN2(jj,4) ~=0 && dataN2(jj,5)==0 && dataN2(jj,6)==0
        dataN2p1p2(ind2,:)= dataN2(jj,:);
        ind2=ind2+1;
    end

    if dataN2(jj,5) ~=0 && dataN2(jj,6)==0
        dataN2p1p2p3(ind3,:)= dataN2(jj,:);
        ind3=ind3+1;
    end

    if dataN2(jj,3)>0 && dataN2(jj,4)>0 && dataN2(jj,5)>0 && dataN2(jj,6)>0
        dataN2p1p2p3p4(ind4,:)= dataN2(jj,:);
        ind4=ind4+1;
    end
end
%% Here we generate a plot equivalent to Fig3 without the red boundaries, for generating those boundaries see the file BoundariesNeg.m
ydomain=linspace(0,1.1,100); 
xdomain = zeros(100,1)+1/3+1/sqrt(3);%Threshold value of N_{AB} found explicitly in the main text
hold on
scatter(dataN2p1p2p3p4(:,1)/1.5,dataN2p1p2p3p4(:,2),'MarkerFaceColor', [39 127 142]/255, 'MarkerEdgeColor',[39 127 142]/255,'MarkerFaceAlpha',0.7,'SizeData',30,'marker','s')
scatter(dataN2p1p2p3(:,1)/1.5,dataN2p1p2p3(:,2),'MarkerFaceColor',  [74 194 109]/255, 'MarkerEdgeColor', [74 194 109]/255,'MarkerFaceAlpha',0.7,'SizeData',30,'marker','d')
scatter(dataN2p1p2(1:100,1)/1.5,dataN2p1p2(1:100,2),'MarkerFaceColor', [70 50 127]/255, 'MarkerEdgeColor',[70 50 127]/255,'MarkerFaceAlpha',0.7,'SizeData',30,'marker','o')
plot(NABN1./3,N12MaxN1,'LineWidth',2,'LineStyle','--','Color',[1 0 0]) % Eq 13 result for 2+1 systems
scatter(0,1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],'marker','pentagram','SizeData',140,'LineWidth',2)
scatter(0,1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'marker','pentagram','SizeData',100)
scatter(0.3333,0.20,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],'marker','d','SizeData',100,'LineWidth',1.2)
scatter(0.666,0,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],'marker','s','SizeData',100,'LineWidth',1.2)
scatter(1,0,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],'marker','o','SizeData',100,'LineWidth',1.2)
plot(xdomain,ydomain,'LineStyle',':','Color',[0 0 0],'LineWidth',2)
text(0.93,0.5,'$N_{AB}^{TH}$','Interpreter','latex','FontSize',16)
text(0.03,0.5,'$-- 2+1$','Interpreter','latex','FontSize',16,'Color','red')
set(gca,'Fontsize',14)
ylim([-0.02,1.1])
xlim([-0.02,1.05])
xlabel('$N_{AB}$','Interpreter','latex','FontSize',20)
ylabel('$N_{A_1A_2}^{max}$','Interpreter','latex','FontSize',20)
legend('$\lambda_4 \neq 0$','$\lambda_3\neq 0, \lambda_4=0$','$\lambda_2\neq0, \lambda_3=\lambda_4=0$','','','','','Interpreter','latex','box','on','Fontsize',16,'location','northeast', 'numColumns',1)
hold off
%exportgraphics(gca,'NegativityBounds.pdf','contenttype','image','Resolution',300)





