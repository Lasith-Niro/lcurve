function [icorner,rho_icorner,eta_icorner,jelim]=l_c_trian(rho,eta,p)
% This program computes the corner of the parametric L-curve given by
% of || A xreg - b ||, || L xreg ||
% Input
%  rho: norm of residuals || A xreg - b ||
%  eta: norm of regularized solutions || L xreg ||
%
% Output
%  icorner: index of the corner point
%  rho_icorner: norm of residual || A xreg - b || corresponding to the corner
%  point
%  eta_icorner: norm of regularized solution || L xreg || corresponding 
%  to the corner point
% 
%  Authors: J. Longina Castellanos and Valia Guerra
% 
%            Group of Numerical Methods
%            Institute of Cybernetics, Mathematics and Physics
%            Calle E No. 309 esq. 15, Vedado, Ciudad de La Habana
%            Cuba, C.P. 10400
%
%  Last Revision: December, 2004
%
%  Reference:
% 
%  Castellanos J.L., G�mez S., Guerra V., (2002), "The Triangle method for 
%  finding the corner of the L curve". Applied Numerical Mathematics 43 (2002),
%  359-373
% 
npoints= length(rho);
if (length(rho) ~= length(eta))
  error('rho and eta do not have the same dimension')
end
%  Initialize the corner (ICORNER) as the last point 
icorner = npoints;
% Convert to logarithms.
lrho = log(rho); leta = log(eta); 
% Check the monotonicity conditions on log(//x//) (increasing) and
% log(//Ax-b//) (decreasing) and eliminate the pair that does not accom-
% plish both conditions at the same time.
% 
i1=1;jelim=0;
lrhor(1)=lrho(1);
letar(1)=leta(1);
for i=2:npoints
    if lrho(i)-lrhor(i1) < sqrt(eps) & leta(i)-letar(i1) > sqrt(eps)
      i1=i1+1;
      lrhor(i1)=lrho(i);
      letar(i1)=leta(i);
    else
      jelim=jelim+1;
      elim(jelim)=i;
    end
end
npointsr=length(lrhor);
%  There must be at least 3 points. Otherwise the last point is
%  selected as corner point
if (length(lrhor) < 3)
  error('Too few data points for the Triangle method')
end
% Scale the points to the box [-1,1]x[-1,1]
nux = 1;
nlx = -1;
nuy = 1;
nly = -1;
t1x = nux - nlx;
t1y = nuy - nly;      
t2x = lrhor(1)*nlx - lrhor(npointsr)*nux;
t3x = max(2.22e-16,(lrhor(1)-lrhor(npointsr)));
t2y = letar(npointsr)*nly - letar(1)*nuy;
t3y = max(2.22e-16,(letar(npointsr)-letar(1)));
for i=1:npointsr
         rh(i) = (t1x*lrhor(i) + t2x)/t3x;
         et(i) = (t1y*letar(i) + t2y)/t3y;
end
% Make the choice of the corner 
icorner = npointsr;
cosant = -2;
cosmax = -2;
arneg = -eps^(1/3);
cte = cos(7*pi/8);
for i = 1:npointsr-2
    for j = i:npointsr-2
% 
%    B --> P(i)
%    A --> P(j+1)
%   C --> P(itlin)
%   
%    Triangle ABC
% 
%         C 
%          x
%    
%       A x
%                x B
% 
% Compute vector BA (B-A)
% 
     bmax = rh(i) - rh(j+1);
     bmay = et(i) - et(j+1);
%   Compute vector AC (A-C) 
     amcx = rh(j+1) - rh(npointsr);
     amcy = et(j+1) - et(npointsr);
%  Compute the scalar product BA(trans)AC
     atb = (bmax*(-amcx) + bmay*(-amcy));
%  Compute //AB// , //CA//
     nbma = sqrt(bmax*bmax + bmay*bmay);
     namc = sqrt(amcx*amcx + amcy*amcy);
% Compute cosine of angle (CAB) and set in COSFI
     cosfi = atb/(nbma*namc);
% Compute the area of the triangle as det(BA,AC)
     area = bmax*amcy - bmay*amcx;
     if (area < arneg) & (cosfi >= cte) & (cosfi < 0.0) & (cosfi >= cosmax)
% 
%  The cosine of current angle CAB is negative and greater than
%  the greatest negative cosine with negative area of the triangle 
%  obtained before the current triangle area is sufficiently negative
%  then.
%  Set point A=P(j+1) as the posible corner: ICORNER = J+1
   icorner = j+1;
% 
% Save COSFI in COSMAX, the greatest cosine value of all the angles
% whose corresponding triangle area is negative
% 
   cosmax = cosfi;
   end
  end
end
%  Restore the correct index of the corner in the original
%  points
cicor = icorner;
if jelim ~= 0
  for i=1:jelim
              if cicor >= elim(i) 
                    cicor = cicor+1;
                end
     end
end
%  Set ICORNER to the index of the original points 
icorner = cicor;
rho_icorner=rho(icorner);
eta_icorner=eta(icorner);
if (nargin > 2)
% Plot
plot(lrho,leta,'.:')
hold on
plot(lrho(icorner),leta(icorner),'*r')
title('Discrete L-curve and the corner point')
xlabel('log(residual)')
ylabel('log(solution norm)')
hold off
end

