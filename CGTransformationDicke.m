function [T, Prod_index]= CGTransformation(N)
%Constructs the transformation matrix <j1,m1,j2,m2 | J, M > for an ensemble
%of N spins, in this case we consider that Jmax = N/2 and j1max=j2max=N/4

num_Dicke = N + 1; %number of columns
num_Product = 0; %number of rows, will be updated

Prod_index = containers.Map;

index =1;

m1range = linspace(N/4,-N/4,2*N/4 + 1);
m2range = linspace(N/4,-N/4,2*N/4 +1);

%This first loop calculates the number of product states and assign an
%index to them linked to the values of j1,m1,j2,m2. Here we consider that
%we are in the Dicke manifold of N spins.

for j1 = N/4
    for j2 = N/4
        for m1=m1range
            for m2=m2range
                Prod_index(sprintf('%d_%d_%d_%d',j1,m1,j2,m2))=index;
                num_Product = num_Product +1;
                index = index +1;
            end
        end
    end
end

%Here we compute the CG coefficient for each entry

T = sparse(num_Product,num_Dicke);
MRange = linspace(N/2,-N/2,N+1);
for M = MRange
    dicke_idx = -M + N/2 + 1; %The first state of this basis is  |N/2,N/2>
    for j1 = N/4
        for j2 = N/4
            for m1 = m1range
                m2 = M - m1;
                if abs(m2) > j2
                    continue; % Skip invalid states
                end

                % Compute Clebsch-Gordan coefficient
                CG_coeff = ClebschGordan(j1,m1,j2,m2,N/2,M);

                if abs(CG_coeff) > 1e-10
                    j1
                    j2
                    coupled_idx = Prod_index(sprintf('%d_%d_%d_%d', j1, m1, j2, m2));
                    T(coupled_idx, dicke_idx) = CG_coeff;
                end
            end
        end
    end
end

end