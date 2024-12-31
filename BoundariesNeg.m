f = @(x) maxN(x(1), x(2), x(3), x(4)); %To maximimize the function a minus sign can be added -maxN(...)
datafinal=zeros(100,2);
ind1=1;
for N0=linspace(0,1,100) %Defines the range of N_{AB}* defined in the Appendix A
    g = @(x) Neg(x(1), x(2), x(3), x(4)) - N0; %constraint equivalent to g_2 in Appendix A
    
    eq_constraint = @(x) deal([f(x)], [g(x); sum(x) - 1]); %This additionally enforces g_1 in Appendix A
    
    lb = [0.25, 0, 0, 0]; %These two arrays define the minimum and maximum value of each eigenvalue lambda_i
    ub = [1, 1/2, 1/3, 1/4];
    
    %The variables A and b are used to define the inequalities h_1, h_2,
    %and h_3 in Appendix A

    A = [-1,  1,  0,  0;  % x2 - x1 <= 0
          0, -1,  1,  0;  % x3 - x2 <= 0
          0,  0, -1,  1]; % x4 - x3 <= 0
    b = [0; 0; 0];
    
    data=zeros(100,3)
    for nn=1:100
        %To optimize we consider multiple random initial states 
        x0=[rand(),rand(),rand(),rand()];
        x0=sort(x0,'descend');
        x0 = x0/(x0(1)+x0(2)+x0(3)+x0(4));
        %Below we define the parameters for the optimization
        options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','StepTolerance', 1e-6,'FunctionTolerance', 1e-4,'ConstraintTolerance', 1e-6);
        if maxN(x0(1),x0(2),x0(3),x0(4)) > 0 %We check that initially the maximum negativity within A is a positive value
            [x_max, y_max_neg] = fmincon(f, x0, A, b, [], [], lb, ub, eq_constraint, options);
            data(nn,1)=y_max_neg;
            data(nn,2)=maxN(x_max(1),x_max(2),x_max(3),x_max(4));
            data(nn,3)=Neg(x_max(1),x_max(2),x_max(3),x_max(4));
            if abs(data(nn,3)- N0) >= 1e-3 %Here we confirm that the optimization is yielding a value consistent with g_2
               data(nn,1) =100; 
            end
        else
            data(nn,1)=100; %We use this to reject initial values that yield a negative negativity (see Eq 12)
            data(nn,2)=100;
            data(nn,3)=100;
        end  
    end
    datafinal(ind1,1) =N0;
    datafinal(ind1,2)=min(data(:,1)); %To prevent getting trapped in local minima we take the minimum value over all initial states
    ind1=ind1+1;
end


% Here we process the data to eliminate all data with invalid initial conditions
newdata=zeros(1,2);
indn=1;
for jj=1:length(datafinal)
    if datafinal(jj,2)<100
        newdata(indn,1) = datafinal(jj,1);
        newdata(indn,2) = datafinal(jj,2);
        indn=indn+1;
    end
end

%Below we store the data generated
%save('BoundariesNegative.txt','newdata','-ascii')



% This function defines the maximum negativity within subsystem A defined in
% Eq 12 of the main text
function y = maxN(x1,x2,x3,x4)
y = sqrt((x1-x3)^2+(x2-x4)^2)-x2-x4;
end
%This function defines the negativity in Eq 15 of the main text
function y = Neg(x1,x2,x3,x4)
y = 2/3*(sqrt(x1*x2)+sqrt(x1*x3)+sqrt(x1*x4)+sqrt(x2*x3)+sqrt(x2*x4)+sqrt(x3*x4));
end