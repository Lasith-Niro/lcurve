% define globals for our inner functions
global b;
global A;
global S;
global bh;

% clear everything
clf;

% random seeds
seed = 1;
randn('state',seed);

% size of our matrix / solutions
n = 40;

% build A and b according to the following methods
%method=input('1->phillips,2->shaw,3->baart,4->deriv2,5->inv.lapl.trnsf. ');
method = 5;

% Phillips
if method==1
   [A,b,exactsol]=phillips(n);
   % A*exactsol-b is not tiny
   b=A*exactsol;
   disp('matrix symmetric');
   condnumb=cond(A)
% Shaw
elseif method==2
   [A,b,exactsol]=shaw(n); b=A*exactsol;
   disp('matrix symmetric');
   condnumb=cond(A)
% Baart
elseif method==3
   [A,b,exactsol]=baart(n); b=A*exactsol;
   disp('matrix nonsymmetric');
   condnumb=cond(A)
% Deriv2
elseif method==4
   ex=input('1->sol=t, 2->sol=exp(t), 3->sol p.w. linear ');
   [A,b,exactsol]=deriv2(n,ex); b=A*exactsol;
   % for ex=1 and ex=2 
   disp('matrix symmetric');
   condnumb=cond(A)
% Inverse Laplace Transform
else
   ex=input('1->sol=exp(-t/2), 2->sol=1-exp(-t/2), 3->t^2*exp(-t/2). 4->sol p.w.constant ');
   [A,b,exactsol]=ilaplace(n,ex); b=A*exactsol;
   disp('matrix nonsymmetric');
   condnumb=cond(A)   
end

% create noise vector err and add it to b
eta = input('eta (noise level roughly 10^(-eta): ');
err = rand(length(b),1);
err = err / norm(err);
err = err * norm(b) * 10^(-eta);

% remember the exact b signal
bexact    = b;
b         = b + err;
normnoise = norm(err);

% SVD
[U,S,V]=svd(A); 

% for finding p(x) from Dr. Reichel's method for min(Ax-b)
global bh;
bh = U'*b;

% p(x) function for efficiently finding solution norms min(Ax - b)
function [p] = p(x)
    global bh
    global S
    global A
    p = 0;
    for j = 1:length(A)
        p = p + (x /(S(j,j)^2 + x))^2 * bh(j)^2;
    end
end

% preal function for getting the actual value using matrices. depricated
function [preal] = preal(x)
    global A
    global b
    preal = norm(A*(A'*A + x*eye(length(A)))^-1 * A'*b - b)^2;
end

% derived p'(x) function for Newton's Method... Pretty sure it works
function [dp] = dp(x)
    global A
    global bh;
    global S;
    dp = 0;
    for j = 1:length(A)
        dp = dp + bh(j)^2 * ((2*x)/(S(j,j)^2 + x)^2 \
                         - (2*x^2)/(S(j,j)^2 +x)^3);
    end
end

hold on

% we can plot p(x) and see where it intersects normnoise
%for l = 0.01:.1:1
%   plot(l,p(l),'o')
%end

% plot normnoise
%plot([0,1],[normnoise^2,normnoise^2])

% Newton's method to find the zero of p(x) - normnoise ^2
x = .0001;
xprev = 0;
while (abs(p(x) - normnoise^2) > 1e-10) && (abs(xprev - x) > 1e-10)
    xprev = x;
    x = x - (p(x) - normnoise^2)/dp(x);
end

% using lambda, we can make our solution x using tikhonov
lambda = x;
xsol = (A'*A + lambda *eye(n))^-1 * A' *b;

% plot our final answers
hold on
plot(exactsol,'b')
plot(xsol,'r');
finalerror = norm(exactsol - xsol)
