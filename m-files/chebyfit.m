function chebyfit;
% generate noisy data, fit by Chebyshev polynomial, determine appropriate
% degree
clear all;
close all;
format short e;
solnrm = []; resnrm=[];
seed = input('random seed: ');
randn('state',seed);
% deetermine matrix
n = input('Number of Chebyshev points: ');
% generate zeros z(k) of Chebyshev polynomial of degree n
for k=1:n
  z(k)=cos(pi*(2*k-1)/(2*n));
end
% generate orthogonal matrix determined by a tabulation of all Chebyshev
% polynomials of degree less than n at the nodes z(k)
for k=1:n
  A(k,1)=1; A(k,2)=z(k);
  for j=2:n-1
    A(k,j+1)=2*z(k)*A(k,j)-A(k,j-1);
  end
end
% orthonormalize matrix
A=A*sqrt(2/n); A(:,1)=A(:,1)/sqrt(2);
% define noise-free function f
f=zeros(n,1);
for k=1:n
  f(k)=exp(z(k));
end
% exact Fourier coefficients
  cf=A'*f;
% define noisy function ftilde
eta = input('eta (noise level roughly 10^(-eta):  (eta=0 -> no noise) ');
if eta==0
  fdelta=f; err=zeros(n,1);
else
  w=randn(length(f),1);
  err=norm(f,2)*w/(10^eta*sqrt(n));
  noiselevel=norm(err)/norm(f,2)
  fdelta = f + err;
end
normnoise=norm(err)
% Fourier coefficients of noise
cerr=A'*err;
% Fourier coefficients of noisy function
cfdelta=A'*fdelta;
% unnormalize matrix and determine coefficients for unnormalized Chebyshev 
% polynomials
A=A/sqrt(2/n); A(:,1)=A(:,1)*sqrt(2);
cfdelta=cfdelta*sqrt(2/n); cfdelta(1)=cfdelta(1)/sqrt(2);
% compute least-squares polynomials: column j of matrices polyn contains
% a tabulation of the least-squares polynomial of degree j-1.
polyn=cfdelta(1)*A(:,1);
for k=2:n
  polyn(:,k)=polyn(:,k-1)+cfdelta(k)*A(:,k); 
end
% compute residual error for least-squares polynomials of degrees 0,1,...,p-1
p=3*n/4;
for k=1:p
   reserr(k)=norm(polyn(:,k)-fdelta);
end
% compute approximation error of least-squares polynomial on fine grid
% with m=10*n points
m=10*n;
% generate m Chebyshev points
for k=1:m
  zm(k)=cos(pi*(2*k-1)/(2*m));
end
% generate orthogonal matrix determined by a tabulation of all Chebyshev
% polynomials of degree less than n at the nodes zm(k)
for k=1:m
  Am(k,1)=1; Am(k,2)=zm(k);
  for j=2:n-1
    Am(k,j+1)=2*zm(k)*Am(k,j)-Am(k,j-1);
  end
end
% tabulate function at m nodes
fm=zeros(m,1);
for k=1:m
  fm(k)=exp(zm(k));
end
% evaluate least-squares polynomial on fine grid
polm=cfdelta(1)*Am(:,1);  
for k=2:n
  polm(:,k)=polm(:,k-1)+cfdelta(k)*Am(:,k); 
end
% compute least-squares error on fine grid. Note weighting factor.
for k=1:n
  err(k)=norm(polm(:,k)-fm)/sqrt(m);
end
[minerr,minidx]=min(err);
minidx
minerr
figure
semilogy([0:n-1],err,'k.:'); hold on
semilogy(minidx-1,err(minidx),'r*'); hold off
title('errors in computed approximate solutions');
%figure
%semilogy([0:n-1],err,'k.:'); hold on
%semilogy(minidx-1,err(minidx),'r*'); hold off
% residual L-curve
[icorner,rho_icorner,eta_icorner,jelim]=l_c_trian(reserr,[1:length(reserr)]);
corner_index_for_residual_L_curve_by_triangle_method=icorner
reserr_icorner=reserr(icorner);
figure;
semilogy([0:p-1],reserr,'k.:');
hold on;
semilogy(icorner-1,reserr(icorner),'*r');
title('Discrete residual L-curve and the corner point by triangle method')
xlabel('residual index');
ylabel('log(residual)');
hold off
%figure;
%semilogy([0:p-1],reserr,'k.:');
%hold on;
%semilogy(icorner-1,reserr(icorner),'*r');
%hold off
% pruning algorithm
[k_corner,info]=corner(reserr,[1:length(reserr)]);
corner_index_for_residual_L_curve_by_pruning_method=k_corner
figure;
semilogy([0:p-1],reserr,'k.:');
hold on;
semilogy(k_corner-1,reserr(k_corner),'*r');
title('Discrete residual L-curve and the corner point by pruning method')
xlabel('residual index');
ylabel('log(residual)');
hold off
%figure;
%semilogy([0:p-1],reserr,'k.:');
%hold on;
%semilogy(k_corner-1,reserr(k_corner),'*r');
%hold off
%
disp('solnmbr    error    reserror')
[[1:p]',err(1:p),reserr(1:p)']
%
solnmbr=input('plot computed solution nmbr (>0) = ');
while solnmbr>0
   figure;
   plot(zm,polm(:,solnmbr),'-k');
   hold on
   plot(zm,fm,'-.b');
   hold off
   legend('computed polynomial','error-free function');
   solnmbr=input('plot computed solution nmbr (>0) = ');
end
