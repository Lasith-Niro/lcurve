function [solnrm,resnrm,sol,exactsol,b,U] = main;
% [solnrm,resnrm,sol,xx,b] = main
% test stopping criterion for tsvd
close all;
format short e;
solnrm = []; resnrm=[];
seed = input('random seed: ');
randn('state',seed);
n = input('Enter order of the matrix: ');
method=input('1->phillips,2->shaw,3->baart,4->deriv2,5->inv.lapl.trnsf. ');
if method==1
   [A,b,exactsol]=phillips(n);
   % A*exactsol-b is not tiny
   b=A*exactsol;
   disp('matrix symmetric');
   condnumb=cond(A)
elseif method==2
   [A,b,exactsol]=shaw(n); b=A*exactsol;
   disp('matrix symmetric');
   condnumb=cond(A)
elseif method==3
   [A,b,exactsol]=baart(n); b=A*exactsol;
   disp('matrix nonsymmetric');
   condnumb=cond(A)
elseif method==4
   ex=input('1->sol=t, 2->sol=exp(t), 3->sol p.w. linear ');
   [A,b,exactsol]=deriv2(n,ex); b=A*exactsol;
   % for ex=1 and ex=2 
   disp('matrix symmetric');
   condnumb=cond(A)
else
   % inverse Laplace transform
   ex=input('1->sol=exp(-t/2), 2->sol=1-exp(-t/2), 3->t^2*exp(-t/2). 4->sol p.w.constant ');
   [A,b,exactsol]=ilaplace(n,ex); b=A*exactsol;
   disp('matrix nonsymmetric');
   condnumb=cond(A)   
end
% determine rhs
errc = input('noise level: ');
err = randn(length(b),1);
err = errc*err*norm(b,2)/norm(err,2);
b = b + err;
normnoise=norm(err)
% TSVD
[U,S,V]=svd(A); 
singval=diag(S);
% extrapolation test
s=0;
gamma=[];
UTb=U'*b;
for k=1:n
    s=s+(singval(k)/UTb(k))^2;
    gamma(k)=1/sqrt(s);
end
figure;
semilogy(gamma,'k.:')
% end test
% matrix is n by n
p=input('number of tsvd solutions=');
%p=n-1;
[sol,resnrm,solnrm]=tsvd(U,singval,V,b,[1:p]);
err=[]; reserr=[];
for k=1:p
   err(k)=norm(sol(:,k)-exactsol);
   reserr(k)=norm(A*sol(:,k)-b);
end
for k=1:p-1
    deltasol(k)=norm(sol(:,k+1)-sol(:,k));
end
deltasol(p)=0;
figure
semilogy(err,'k'); hold; semilogy(err,'k+');
title('errors in computed approximate solutions');
figure
semilogy(deltasol,'k'); hold; semilogy(deltasol,'k+');
title('norm of differences of consecutive iterates');
% residual L-curve
[icorner,rho_icorner,eta_icorner,jelim]=l_c_trian(reserr,[1:p]);
corner_index_for_residual_L_curve_by_triangle_method=icorner
reserr_icorner=reserr(icorner);
figure;
semilogy(reserr,'k.:');
hold on;
semilogy(icorner,reserr(icorner),'*r');
title('Discrete residual L-curve and the corner point by triangle method')
xlabel('residual index');
ylabel('log(residual)');
hold off
% backtrack to min deltasol
if icorner==p,
  icorner=p-1;
end
deltasolmin=deltasol(icorner);
while deltasol(icorner-1)<=deltasolmin & icorner>1
    icorner=icorner-1;
    deltasolmin=deltasol(icorner);
end
corner_index_after_backtracking=icorner
% forward track to min deltasol
deltasolmin=deltasol(icorner);
while deltasol(icorner+1)<=deltasolmin & icorner+1<p
    icorner=icorner+1;
    deltasolmin=deltasol(icorner);
end
corner_index_after_forward_tracking=icorner
% pruning algorithm
[k_corner,info]=corner(reserr,[1:p]);
conrner_info=info
corner_index_for_residual_L_curve_by_pruning_method=k_corner
figure;
semilogy(reserr,'k.:');
hold on;
semilogy(k_corner,reserr(k_corner),'*r');
title('Discrete residual L-curve and the corner point by pruning method')
xlabel('residual index');
ylabel('log(residual)');
hold off
% backtrack to min deltasol
if k_corner==p
  k_corner=p-1;
end
deltasolmin=deltasol(k_corner);
while deltasol(k_corner-1)<=deltasolmin & k_corner>1
    k_corner=k_corner-1;
    deltasolmin=deltasol(k_corner);
end
corner_index_after_backtracking=k_corner
% forward track to min deltasol
deltasolmin=deltasol(k_corner);
while deltasol(k_corner+1)<=deltasolmin
    k_corner=k_corner+1;
    deltasolmin=deltasol(k_corner);
end
corner_index_after_forward_tracking=k_corner
% standard L-curve
logreserr=log(reserr'); logsolnrm=log(solnrm);
figure;
plot(logsolnrm,logreserr,'k'); hold
for k=1:p
    plot(logsolnrm(k),logreserr(k),'r+'); pause(0.5);
end
title('L-curve');
% determine corner by triangle method
[icorner,rho_icorner,eta_icorner,jelim]=l_c_trian(reserr,solnrm);
corner_index_for_L_curve_by_triangle_method=icorner
figure;
plot(logsolnrm,logreserr,'k.:');
hold on;
plot(logsolnrm(icorner),logreserr(icorner),'*r');
title('Discrete L-curve and the corner point by triangle method')
xlabel('log(solnrm)');
ylabel('log(residual)');
hold off
% backtrack to min deltasol
if icorner==p
  icorner=p-1;
end
deltasolmin=deltasol(icorner);
while deltasol(icorner-1)<=deltasolmin & icorner>1
    icorner=icorner-1;
    deltasolmin=deltasol(icorner);
end
corner_index_after_backtracking=icorner
% forward track to min deltasol
deltasolmin=deltasol(icorner);
while deltasol(icorner+1)<=deltasolmin & icorner+1<p
    icorner=icorner+1;
    deltasolmin=deltasol(icorner);
end
corner_index_after_forward_tracking=icorner
% pruning algorithm
[k_corner,info]=corner(reserr,solnrm);
conrner_info=info
corner_index_for_L_curve_by_pruning_method=k_corner
figure;
plot(logsolnrm,logreserr,'k.:');
hold on;
plot(logsolnrm(k_corner),logreserr(k_corner),'*r');
title('Discrete L-curve and the corner point by pruning method')
xlabel('log(solnrm)');
ylabel('log(residual)');
hold off
% backtrack to min deltasol
if k_corner==p,
  k_corner=p-1;
end
deltasolmin=deltasol(k_corner);
while deltasol(k_corner-1)<=deltasolmin & k_corner>1
    k_corner=k_corner-1;
    deltasolmin=deltasol(k_corner);
end
corner_index_after_backtracking=k_corner
% forward track to min deltasol
deltasolmin=deltasol(k_corner);
while deltasol(k_corner+1)<=deltasolmin & k_corner+1<p
    k_corner=k_corner+1;
    deltasolmin=deltasol(k_corner);
end
corner_index_after_forward_tracking=k_corner
%
disp('L-curve points: solnmbr  logsolnrm  logreserr');
[[1:p]',logsolnrm,logreserr]
disp('solnmbr    error    reserror   deltasol')
[[1:p]',err',reserr',deltasol']
%
solnmbr=input('plot computed solution nmbr (>0) = ');
while solnmbr>0
   figure;
   plot(sol(:,solnmbr),'-k');
   hold on
   plot(exactsol,'-.b');
   hold off
   legend('computed approx. solution','exact solution');
   solnmbr=input('plot computed solution nmbr (>0) = ');
end
