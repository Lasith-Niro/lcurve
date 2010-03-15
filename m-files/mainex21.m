function [solnrm,resnrm,sol,exactsol,b] = mainex21;
% [solnrm,resnrm,sol,xx,b] = main
% test stopping criterion for tsvd
close all;
format short e;
solnrm = []; resnrm=[];
seed = input('random seed: ');
randn('state',seed);
% deetermine matrix
n = input('Enter order of the matrix: ');
zk=ones(n,1); A=zk;
z=exp(2*pi*i*[0:n-1]'/n);
for k=1:n-1
  zk=zk.*z;
  A=[A,zk]; 
end
a=input('approx of 1/(x-a): a=');
% determine exact rhs, tabulate function
  f=1./(a*ones(n,1)-z);
% exact Fourier coefficients
  cf=A'*f/n;
% determine noisy rhs
eta = input('eta (noise level roughly 10^(-eta):  (eta=0 -> no noise) ');
if eta==0
  ferr=f; err=zeros(n,1);
else
  w=randn(length(f),1);
  err=norm(f,2)*w/(10^eta*sqrt(n));
  noiselevel=norm(err)/norm(f,2)
  ferr = f + err;
end
normnoise=norm(err)/sqrt(n)
% Fourier coefficientr of noise
cerr=A'*err;
% Fourier coefficients of noisy function
cferr=A'*ferr/n;
% compute solution 
sol=[]; sol(:,1)=cferr(1)*A(:,1); 
for k=2:n
  sol(:,k)=sol(:,k-1)+cferr(k)*A(:,k); 
end
% compute residual error
p=3*n/4
for k=1:p
   reserr(k)=norm(sol(:,k)-ferr)/sqrt(n);
end
% compute approximations of on fine grid with m=10*n points
m=10*n;
zk=ones(m,1); Atall=zk;
z=exp(2*pi*i*[0:m-1]'/m); 
for k=1:n-1
  zk=zk.*z;
  Atall=[Atall,zk]; 
end
% tabulate function
ftall=1./(a*ones(m,1)-z);
%
soltall=[]; soltall(:,1)=cferr(1)*Atall(:,1); 
for k=2:n
  soltall(:,k)=soltall(:,k-1)+cferr(k)*Atall(:,k); 
end
% compute uniform norm approximation error on fine grid
for k=1:n
  err(k)=norm(soltall(:,k)-ftall)/sqrt(m);
end
[minerr,minidx]=min(err);
minidx
minerr
figure
semilogy([0:n-1],err,'k.:'); hold on
semilogy(minidx-1,err(minidx),'r*'); hold off
title('errors in computed approximate solutions');
figure
semilogy([0:n-1],err,'k.:'); hold on
semilogy(minidx-1,err(minidx),'r*'); hold off
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
figure;
semilogy([0:p-1],reserr,'k.:');
hold on;
semilogy(icorner-1,reserr(icorner),'*r');
hold off
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
figure;
semilogy([0:p-1],reserr,'k.:');
hold on;
semilogy(k_corner-1,reserr(k_corner),'*r');
hold off
%
disp('solnmbr    error    reserror')
[[1:p]',err(1:p),reserr(1:p)']
%
tval=[0:pi/(5*n):2*pi]';
a
ftall=[ftall;1/(a-1)];
sizetval=size(tval);
sizesoltall=size(soltall);
soltall=[soltall;soltall(1,:)];
sizesoltall=size(soltall);
solnmbr=input('plot computed solution nmbr (>0) = ');
while solnmbr>0
   figure;
   plot(tval,real(soltall(:,solnmbr)),'-k');
   hold on
   plot(tval,real(ftall),'-.b');
   hold off
   legend('real part of computed approx. solution','real part of exact solution');
   figure;
   plot(tval,real(soltall(:,solnmbr)),'-k');
   hold on
   plot(tval,real(ftall),'-.b');
   hold off
%
   figure;
   plot(tval,imag(soltall(:,solnmbr)),'-k');
   hold on
   plot(tval,imag(ftall),'-.b');
   hold off
   legend('imag part of computed approx. solution','imag part of exact solution');
   figure;
   plot(tval,imag(soltall(:,solnmbr)),'-k');
   hold on
   plot(tval,imag(ftall),'-.b');
   hold off
   solnmbr=input('plot computed solution nmbr (>0) = ');
end
