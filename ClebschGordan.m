function cg = ClebschGordan(j1,m1,j2,m2,J,M)
%%
%We check the validity of the inputs
if m1+m2 ~= M
    fprintf('m1 and m2 dont add up to M');
    cg=0;
    return
end

if j1 <0
    fprintf('j1 needs to be a non-negative integer / half-integer');
    cg=0;
    return
end

if j2<0
    fprintf('j2 needs to be a non-negative integer / half-integer');
    cg=0;
    return
end

if J<0
    fprintf('J needs to be a non-negative integer / half-integer');
    cg=0;
    return
end

if abs(j1-j2) > J 
    fprintf('J is out of bounds');
    cg=0;
    return
end

if abs(j1+j2) < J 
    fprintf('J is out of bounds');
    cg=0;
    return
end

if abs(m1) > j1
    fprintf('m1 is out of bounds');
    cg=0;
    return
end

if abs(m2) > j2
    fprintf('m1 is out of bounds');
    cg=0;
    return
end

if abs(M) > J
    fprintf('m1 is out of bounds');
    cg=0;
    return
end

if mod(abs(2*m1),1)~=0
    fprintf('m1 is not an integer / half-integer');
    cg=0;
    return
end

if mod(abs(2*m2),1)~=0
    fprintf('m2 is not an integer / half-integer');
    cg=0;
    return
end

if mod(abs(2*M),1)~=0
    fprintf('M is not an integer / half-integer');
    cg=0;
    return
end

if mod(abs(2*j1),1)~=0
    fprintf('j1 is not an integer / half-integer');
    cg=0;
    return
end

if mod(abs(2*j2),1)~=0
    fprintf('j2 is not an integer / half-integer');
    cg=0;
    return
end

if mod(abs(2*J),1)~=0
    fprintf('J is not an integer / half-integer');
    cg=0;
    return
end


%%
%If all the checks are completed we compute the CG coefficients, we use
%Gamma functions which are more stable than factorials for large numbers

term1 = sqrt((2*J+1)*gamma(J+j1-j2+1)*gamma(J-j1+j2+1)*gamma(j1+j2-J+1))/sqrt(gamma(j1+j2+J+2));
term2 = sqrt(gamma(J+M+1)*gamma(J-M+1)*gamma(j1-m1+1)*gamma(j1+m1+1)*gamma(j2-m2+1)*gamma(j2+m2+1));

lowL = max([0,j2-J-m1,j1-J+m2]);
highL = min([j1+j2-J,j1-m1,j2+m2]);

term3=0;

for k=lowL:highL
    term3=term3+((-1)^k)/(gamma(k+1)*gamma(j1+j2-J-k+1)*gamma(j1-m1-k+1)*gamma(j2+m2-k+1)*gamma(J-j2+m1+k+1)*gamma(J-j1-m2+k+1));
end

cg=term1*term2*term3;
end