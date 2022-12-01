clear all
close all
format compact

tol=10^-7; %tolerance on the residual norm - stoppin criteria

es=1;  % considered example
% NOTE : New examples can be added by modifing file "Coeff_Matrix.m"

n=7;  % amount of mesh poitns : 2^NN-1


w=0.75;  % weight of Jacobi
method=3; % method to compute the optimal weight for jacobi 
          %1 = numerically optimal one obtained through many tests (slow anddoes not always work)
          %0 = uses w in input
          %2 = Provides optimal weight with minimizes the residual after 1 iteration of TGM
          %3 = computed through the approach described in the paper:
          %"Multigrid for two-sided fractional differential equations discretized by finite volume elements on graded meshes" by Marco Donatelli, Rolf Krause, Mariarosa Mazza, Ken Trotti
mesh='composite' % 'composite' = mesh with power(n) points in the first interval (see below)
                 % otherwise standard graded mesh mapped by a function

ss=0.2; %length of the graded part of interval [0 1]
beta=0.9;  %fractional derivative order 2-beta, 0<beta<1
g=0.5;  %gamma in [0,1]
smoother='j'; % j= jacobi,  'gmres' = GMRES
MGlvl=inf;  % levels of MGM, 1 = two-grid, 2-three grid, ... , inf = as many levels as possible
%NOTE: lowest level where coefficient matrix has size < 3x3
band=0;  % band approximation of the coefficient matrix? 0 = no approx, 3,5,7,... means tridiagonal, pentadiagonal,...
length_connection=0.05; %length of the smooth connection between the graded part and the linear one NOTE: length_connection + ss < 1
force_power=0; %0 or 1 (on or off): force the use of a particular power q of the mapping x^q, without concerning about machine precision

maxit=100; %maximum number of iterations



%do not change
galerkin=0;  % 1 = galerkin approach (expensive and does not work for the moment), 0 = rediscretization in place of galerkin



if strcmp(mesh,'composite')
    if 1==2
        %if you want an amount of mesh points in the first interval which is a
        %function of N use this:
        power=@(n)log2(2^n-1);%sqrt(n); %amount of mesh point in the first interval [0,h] as a function of N
        
        % do not change
        grid=3; % 3= composite mesh with amount of points in the first interval given by "power(N)"
        N=2^n-1; %number of unknowns
    else
        %if you want fixed mesh points :  uniform mesh with N2 points, N1 in
        %the first interval
        n=[4,8];
        
        % do not change
        global N1 N2
        grid=[];power=[];
        N1=2^n(1);N2=2^n(2);N=N1+N2;
    end
else
    grid=1; % 1 = singular mesh, 2=non singular mesh (inefficient)
    power=[]; %power q of the mapping x^q. If empty is computed through the rule in the paper
    
    % do not change
    N=2^n-1; %number of unknowns
end




pre=smoother; % methods of the pre and post smoothers
post=pre;


if strcmp(pre,'j')
    nu1=1;  % 1 iteration of pre smoother
    nu2=1;  % 1 iteration of post smoother
else
    nu1=3;   % GMRES, more iterations are needed
    nu2=3;
end
    



s=[ss,0];


% generating coefficient matrix
[A1,HH,HHder]=Coeff_Matrix(beta,N,es,s,band,g,grid,power,length_connection,force_power);
stenX=[1 2 1];




if method~=0
    [A,P,R,w2] = V_SetupPhase_1D(A1,pre,post,galerkin,stenX,MGlvl,method);
else
    [A,P,R,w2] = V_SetupPhase_1D(A1,pre,post,galerkin,stenX,MGlvl,0);
    w2=w;
end

% MGM as solver
b=A1.b; uex=A1.U; A1.Bandx=0;
[x,it] = V(A,zeros(size(b)),b,w2,tol,P,R,nu1,nu2,maxit,pre,post);
norm2=norm(uex-x)/norm(uex);
norminf=norm(uex-x,inf);
fprintf('Iterations to tolerance of MGM: %.0f\rerror in 2-norm: %0.3e\rerror in inf-norm: %0.3e\r\r',it,norm2,norminf)
res=norm(b-A1*x)/norm(b)


% preconditioned GMRES with MGM
f=@(r)V(A,zeros(size(b)),r,w2,-1,P,R,nu1,nu2,1,pre,post);
[x,flag,RELRES,it,RESVEC] = SOLVER(@(x)A1*x,b,tol,maxit,zeros(size(b)),1,f);
norm2=norm(uex-x)/norm(uex);
norminf=norm(uex-x,inf);
fprintf('Iterations to tolerance of GMRES with MGM as preconditioner: %.0f\rerror in 2-norm: %0.3e\rerror in inf-norm: %0.3e\r',it,norm2,norminf)
res=norm(b-A1*x)/norm(b)


