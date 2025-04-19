function sq = SqueezingA(t,rho0,N,V,D,SxA,SyA,SzA)
%Computes the squeezing on subsystem A in an equivalent way to the way we
%do it for AB in the main code.

rhoA = UnitaryEvolutionRho(rho0,t,V,D);

 %Now we compute the average spin direction and find the optimum
   %squeezing, using the details on Ma paper about spin squeezing

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
    
   sq = real(4*VarOpt/(N/2));
end
