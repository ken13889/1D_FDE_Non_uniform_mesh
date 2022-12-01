function [w] = W_JAC_gen(A,method,N,coarsening)
% ---------------- Generates the optimal weight for Jacobi ----------------
%
% METHOD 2 and 3 work for uniform grids only
%
% ISTRUCTIONS:
% In order for this function to work: 
% -> "A" must be a class whose name is BUILDER and its first 3 input must be:
% (1) a dummy variable, if it is empty than it is just useless, but if
% a matrix of the same class is in input than it copies the whole class
% except for the two following parameters in input
% (2) number of points [Nx,Ny] if the equation is 2D or just Nx if is 1D
% (3) the scale parameters [r,s] if 2D or [r] if 1D
%
%
% -> The class must have a built-in funciton named diag(obj), which
% returns the diagonal of the matrix whitout assembling it
%
% -> The class must have the following properties: 
% 1D case: N (grid size), dimension (=1 for 1D), wJ (weight of Jacobi)
% 2D case: Nx,Ny (grid size), dimension (=2 for 2D), wJ (weight of Jacobi)
% (pre and post-smoother have the same weight)
%
% There must be a function called "V_SetupPhase" whose input parameters are
% (1) the matrix class, (2) pre smoother 'j', (3) post smoother 'j', 
% (4) galerkin (1 for galerking, 0 for geometric),
% (5) stenX (stencil of projectors in direction x,
% (6) stenY (stencil of projectors in direction y,
% (7) MGMlvl (levels of the hierarchy),
% (8) method (method to compute the weight of Jacobi.. if is 'skip' than it 
%     does not compute them).


global W_GEN_WARNING1 W_GEN_WARNING2 W_GEN_WARNING3 c

% if isobject(A) && method==4
% %     disp('matrix A is class and method 4 is selected.. using method 2 which is the same but faster')
%     method=2;
% end
% if isobject(A)==0 && method==3
%     disp('matrix A is class and method 4 is selected.. using method 2 which is the same but faster')
%     method=4;
% end

if method~=round(method) || method==0
    w=method;
elseif method==1 %Find best weights through many tests  (SLOW) 
    if isobject(A)
        if A.dimension==1
            b=ones(A.Nx,1);
        elseif A.dimension==2
            b=ones(A.Nx*A.Ny,1);
        else
            disp('not implemented')
        end
    else
        b=ones(length(A),1);
    end
    if isempty(W_GEN_WARNING1)
        if isobject(A)==0
            warning("The selected method could be slow... (Jacobi's weight is being computed by running w-Jacobi many times)")
        else
            warning("Trying to estimate 'W' for Jacobi but GALERKIN method is selected")
        end
        W_GEN_WARNING1=1;
    end
    W=0.1:0.1:1;
    for w=1:length(W)
        x=zeros(size(b));
        for i=1:10
            x=x+W(w)*(b-A*x)./diag(A);
        end
        RES(w)=norm(b-A*x)/norm(b);
    end
    [val,pos]=min(RES);
    if val>1
        warning('Minimum residual is >1... none of the tested weights is working')
        w=NaN;
    else
        w=W(pos);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif method==2  % OPTIMAL WEIGHT CHOSEN BY RUNNING TGM (SLOW)
    n=4;  % If n>4 the estimate of the eigenvalues could be slow
    if isobject(A)==0 % matrice ma non classe
        error('"A" deve essere una classe.. altrimenti il metodo non funziona... si potrebbe modificare')
    end
    if isempty(W_GEN_WARNING3)
        warning('Using projectors with stencil [1 2 1] in the Jacobi weight generator (METHOD=2).')
        W_GEN_WARNING3=1;
    end
    stenX=[1 2 1]; stenY=stenX;
    if A.dimension==1
        % smaller linear system such that N<=2^(2n)
        if A.Nx<2^(2*n), N=A.Nx; else N=2^(2*n)-1; end
        try r=A.r; catch ERR, r=[]; end
        if isempty(r)
            A1=BUILDER(A,N);
        else
            A1=BUILDER(A,N,r);
        end
        [A,P,R,~] = V_SetupPhase_1D(A1,'j','j',0,stenX,1,0);
%         P{1}=spdiags(ones(L,1)*stenX,-1:1,L,L); P{1}=P{1}(:,2:2:end); P{1}=P{1}/norm(P{1},inf); R{1}=P{1}';
        % ALTERNATIVE WHICH DOES NOT USE V_SETUP FUNCTION..............
%         P{1}=spdiags(ones(L,1)*stenX,-1:1,L,L); P{1}=P{1}(1,2:2:end); P{1}=P{1}/norm(P{1},inf); R{1}=P{1};
%         A1{1}=A1; A1{2}=R{1}*Print(A1)*P{1};
    elseif A.dimension==2 % 2D equation
        % smaller linear system such that Nx*Ny<=2^(2n)
        if A.Nx<2^n, Nx=A.Nx; cc=round(2^n/Nx); else, Nx=2^n; cc=1; end
        if A.Ny<2^n
            Ny=A.Ny; cc=round(2^n/Ny);
            if A.Nx>2^n, Nx=min(2^n*cc,A.Nx); end
        else
            Ny=min(2^n*cc,A.Ny);
        end
        
        %make the size odd
        if mod(Nx,2)==0, Nx=Nx-1;end,  if mod(Ny,2)==0, Ny=Ny-1;end 
        if strcmp(coarsening,'x')
            stenX=[1 2 1];
            stenY=[]; 
        elseif strcmp(coarsening,'y')
            stenX=[];
            stenY=[1 2 1];
        elseif strcmp(coarsening,'xy')
            stenX=[1 2 1];
            stenY=[1 2 1];
        end
        A1=BUILDER(A,[Nx,Ny],[A.r,A.s]); AA{1,1}=A1;
        [P{1},R{1},nx,ny] = Proj_Gen2D(stenX,stenY,Nx,Ny);
        gal=0;
        if gal==0
            AA{2,2}=BUILDER(A1,[nx,ny]);
        else
            AA{2,2}=R{1}*Print(A1)*P{1};
        end
        A=AA;
    else
        error('not implemented')
    end
    W=0;
    W1=[1:-0.05:0.1, 0.08, 0.06, 0.05, 0.04, 0.02];
    maxit=20; res_old=inf; w_old=W1(1); res=[];
    % THIS ASSUMES THAT BY DECREASING THE WEIGHT FROM 1 to 0 THE ITERATIONS
    % SHOULD DECREASE AND INCREASE AGAIN... THEN IT FINDS THE MINIMUM
    it=[];
    for w=W1
        b=A1.b; x=zeros(size(A1.b)); tol=10^-7; nb=norm(b); iter=0; Res=1;
        while Res(end)>tol && iter<maxit %TWO GRID
            iter=iter+1;
            x=x+w*(b-A{1,1}*x)./diag(A{1,1});
            r=R{1,1}*(b-A{1,1}*x);
            e=A{2,2}\r;
            x=x+P{1,1}*e;
            RES=(b-A{1,1}*x);
            Res(end+1)=norm(RES)/nb;
            x=x+w*RES./diag(A{1,1});
        end
        it(end+1)=iter;
        res(end+1)=norm(A1*x-b)/norm(b);
%         if iter==maxit && res(end)>res_old && res(end)<1
%             w=w_old;
%             break
%         end
%         w_old=w;
%         res_old=res(end);
%         maxit=iter;
    end
    pos=find(it==min(it));
    w=W1(pos);
    pos=find(res(pos)==min(res(pos)));
    w=w(pos);
    if res(end)>1 && w==W1(end)
        warning('W_JAC_gen   did not succeed in generating an optimal weight... V-cycle never converges')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif method==3  % Computing optimal weight through the idea in first paper (need function f)
    % Eigenvalues of the iteration matrix of Jacobi are inside the region
    % delimited by function f
    if isempty(W_GEN_WARNING2)
%         warning("The selected method could malfunction... (function 'f' for the eigenvalues region has to be accurate)")
        W_GEN_WARNING2=1;
    end
    n=4;  % If n>4 the estimate of the eigenvalues could be slow
    
    if isobject(A)
        if A.dimension==1
            % smaller linear system such that N<=2^(2n)
            if A.Nx<2^(2*n), N=A.Nx; else, N=2^(2*n)-1; end
            try
            A1=BUILDER(A,N);
            catch ERR
                stop=0; N=A.Nx;
                while stop==0
                    if N>2^(2*n)
                        N=(N-mod(N,2))/2;
                    else
                        stop=1;
                    end
                end
                A1=BUILDER(A,N);
            end
            temp=Print(A1)./diag(A1);
        else % 2D equation
            % smaller linear system such that Nx*Ny<=2^(2n)
            if A.Nx<2^n, Nx=A.Nx; cc=round(2^n/Nx); else, Nx=2^n; cc=1; end
            if A.Ny<2^n
                Ny=A.Ny; cc=round(2^n/Ny);
                if A.Nx>2^n, Nx=min(2^n*cc,A.Nx); end
            else
                Ny=min(2^n*cc,A.Ny);
            end
            %make the size odd
            if mod(Nx,2)==0, Nx=Nx-1;end,  if mod(Ny,2)==0, Ny=Ny-1;end 
            A1=BUILDER(A,[Nx,Ny],[A.r,A.s]);
            temp=Print(A1)./diag(A1);
        end
    else
        if length(A)<=2^(2*n)
            temp=A./diag(A);
        else
            error('Not implemented')
        end
    end
    W=0; i=0; maxval=1;
    Eig=eig(full(temp));
    for w=[1:-0.05:0.1,0.09:-0.02:0.03]
        i=i+1;
        eigs=1-w*Eig;
        rho=max(abs(eigs));
        ReIm=[real(eigs),imag(eigs)];
        if exist('f','var')==0
            if W_GEN_WARNING2<2
                disp('Using pre-set function "f" to compute optimal "w" for Jacobi')
                W_GEN_WARNING2=2;
            end
            if exist('c','var')==0
                c=0.95; %c=0.8 must be >0... it decreases the curve
            else
                if isempty(c)
                    c=0.95; %c=0.8 must be >0... it decreases the curve
                end
            end
            ff=@(x,c)(x<=1).*c.*x/2-c/2;
            f=@(x,c)sqrt(1-x.^2)+ff(x,c);%-c/2*x+c/2;
%             fplot(@(x)f(x,c))
%             ylim([0,1])
%             xlim([-1,1])
%             hold on
%             plot(ReIm(:,1),ReIm(:,2))
%             'asd'
        end
        hh=min(f(ReIm(:,1),c)-abs(ReIm(:,2)));
        if  (hh<maxval && 0<=hh)  && rho<1
            W=w;
            break
        end
    end

    if W==0
        w=0.25;
        warning(['No weights available, setting w=',num2str(w)])
    else
        w=W;
    end
elseif method==4  % REDUCING THE SIZE OF THE MATRIX MANUALLY AND THEN USING TGM TO FIND OPTIMAL W
    n=4;
    if isobject(A)
        warning('obj input with this method is not implemented... converting to matrix through Print function')
        try N=[A.Nx,A.Ny]; catch ERR, N=[A.Nx]; end
        A=Print(A);
    end
    if exist('N','var')==0
        error('N must be an input when choosing method 4')
    end
    stenX=[1 2 1]; % STENCIL OF PROJECTORS
    if length(N)==1  % 1D
        len=N;
        if len>2^(2*n)
            t=0;
            while len>2^(2*n) 
                len=(len-mod(len,2))/2;
                t=t+1;
            end
            L=length(2^t:2^t:N);
            AA=zeros(L,L);
            for i=1:L
                AA(i,:)=A(2^t*i,(2^t*i-i+1):(L+2^t*i-i));
            end
        else
            AA=A;
            L=N;
        end
        P{1}=spdiags(ones(L,1)*stenX,-1:1,L,L); P{1}=P{1}(:,2:2:end); P{1}=P{1}/norm(P{1},inf); R{1}=P{1}';
        b=ones(L,1);
    elseif length(N)==2  % 2D
        % This works ONLY if the 2D matrix is of the form IoT+ToI, where o
        % is the kronecker product, T a toeplitz and I the identity matrix
        % (there could be non constand diagonals anywhere.. still works)
        %
        % THIS COULD BE IMPROVED BY USING SEMICORASENING INSTEAD OF FULL
        % COARSENING... IT COULD COMPUTE MORE SUITABLE WEIGHTS FOR SPECIFIC
        % CASES OF SEMICOARSENING
        len1=N(1); len2=N(2);
        t1=0; t2=0; turn=0;
        while len1*len2>2^(2*n) 
            if turn==0 %halves N(1) then N(2) then N(1),...
                len1=(len1-mod(len1,2))/2;
                t1=t1+1; turn=1;
            else
                len2=(len2-mod(len2,2))/2;
                t2=t2+1; turn=0;
            end
        end
        if t1>0 || t2>0
            L1=length(2^t1:2^t1:N(1)); L2=length(2^t2:2^t2:N(2));
            AA=zeros(L1*L2,L1*L2); 
            START_A=(2^t2-1)*N(1)+1:2^t2*N(1);
            START_AA=1:L1;
            for i=1:L2 %block index (row)
                index1A=(i*2^t2-1)*N(1)+1:i*2^t2*N(1);
                index1AA=(i-1)*L1+1:i*L1;
            for j=1:L2 %block index (colum)
                index2A=(START_A)+(j-1+(2^t2-1)*(i-1))*N(1);
                index2AA=(START_AA)+(j-1)*L1;
                block=A(index1A,index2A);
                temp=zeros(L1,L1);
                for k=1:L1
                    temp(k,:)=block(2^t1*k,(2^t1*k-k+1):(L1+2^t1*k-k));
                end
                AA(index1AA,index2AA)=temp;
            end
            end
        else
            L1=len1; L2=len2; AA=A;
        end
        if strcmp(coarsening,'x')
            stenX=[1 2 1];
            stenY=[]; 
        elseif strcmp(coarsening,'y')
            stenX=[];
            stenY=[1 2 1];
        elseif strcmp(coarsening,'xy')
            stenX=[1 2 1];
            stenY=[1 2 1];
        end
        
        nx=L1; ny=L2;
        [P1,R1,~,~] = Proj_Gen2D(stenX,stenY,nx,ny);
        P{1}=P1; R{1}=R1;
        b=ones(L1*L2,1);
    else
        error('not implemented')
    end
 
    A1{1}=AA; A1{2}=R{1}*AA*P{1};
    i=0; W1=[1:-0.05:0.1, 0.08, 0.06, 0.05, 0.04, 0.02];
    maxit=20;
    % THIS ASSUMES THAT BY DECREASING THE WEIGHT FROM 1 to 0 THE ITERATIONS
    % SHOULD DECREASE AND INCREASE AGAIN... THEN IT FINDS THE MINIMUM
     res=[]; it=[]; w_old=W1(1); res_old=inf;
    for w=W1
        i=i+1;
%         ww(1,1,1)=w; ww(2,1,1)=w;
%         [x,iter] = V(A1,zeros(size(b)),b,ww,10^-7,P,R,1,1,maxit);
        x=zeros(size(b)); tol=10^-7; nb=norm(b); iter=0; Res=1;
        while Res>tol && iter<maxit %TWO GRID
            iter=iter+1;
            x=x+w*(b-A1{1}*x)./diag(A1{1});
            r=R{1}*(b-A1{1}*x);
            e=A1{2}\r;
            x=x+P{1}*e;
            RES=(b-A1{1}*x);
            Res=norm(RES)/nb;
            x=x+w*RES./diag(A1{1});
        end
        res(end+1)=norm(A1{1}*x-b)/norm(b);
        it(end+1)=iter;
%         if iter==maxit && res(end)>res_old
%             w=w_old;
%             break
%         end
%         w_old=w;
%         res_old=res(end);
%         res(end)
%         maxit=iter;
    end
    pos=find(it==min(it));
    w=W1(pos);
    pos=find(res(pos)==min(res(pos)));
    w=w(pos);
elseif method==5  % REDUCING THE SIZE OF THE MATRIX MANUALLY AND THEN USING METHOD 3 TO FIND OPTIMAL W (CHEAPER THAN 4)
    n=4;
    if isobject(A)
        warning('obj input with this method is not implemented... converting to matrix through Print function')
        try N=[A.Nx,A.Ny]; catch ERR, N=[A.Nx]; end
        A=Print(A);
    end
    if exist('N','var')==0
        error('N must be an input when choosing method 4')
    end
    stenX=[1 2 1]; % STENCIL OF PROJECTORS
    if length(N)==1  % 1D
        len=N;
        if len>2^(2*n)
            t=0;
            while len>2^(2*n) 
                len=(len-mod(len,2))/2;
                t=t+1;
            end
            L=length(2^t:2^t:N);
            AA=zeros(L,L);
            for i=1:L
                AA(i,:)=A(2^t*i,(2^t*i-i+1):(L+2^t*i-i));
            end
        else
            AA=A;
        end
    elseif length(N)==2  % 2D
        % This works ONLY if the 2D matrix is of the form IoT+ToI, where o
        % is the kronecker product, T a toeplitz and I the identity matrix
        % (there could be non constand diagonals anywhere.. still works)
        %
        % THIS COULD BE IMPROVED BY USING SEMICORASENING INSTEAD OF FULL
        % COARSENING... IT COULD COMPUTE MORE SUITABLE WEIGHTS FOR SPECIFIC
        % CASES OF SEMICOARSENING
        len1=N(1); len2=N(2);
        t1=0; t2=0; turn=0;
        while len1*len2>2^(2*n) 
            if turn==0 %halves N(1) then N(2) then N(1),...
                len1=(len1-mod(len1,2))/2;
                t1=t1+1; turn=1;
            else
                len2=(len2-mod(len2,2))/2;
                t2=t2+1; turn=0;
            end
        end
        if t1>0 || t2>0
            L1=length(2^t1:2^t1:N(1)); L2=length(2^t2:2^t2:N(2));
            AA=zeros(L1*L2,L1*L2); 
            START_A=(2^t2-1)*N(1)+1:2^t2*N(1);
            START_AA=1:L1;
            for i=1:L2 %block index (row)
                index1A=(i*2^t2-1)*N(1)+1:i*2^t2*N(1);
                index1AA=(i-1)*L1+1:i*L1;
            for j=1:L2 %block index (colum)
                index2A=(START_A)+(j-1+(2^t2-1)*(i-1))*N(1);
                index2AA=(START_AA)+(j-1)*L1;
                block=A(index1A,index2A);
                temp=zeros(L1,L1);
                for k=1:L1
                    temp(k,:)=block(2^t1*k,(2^t1*k-k+1):(L1+2^t1*k-k));
                end
                AA(index1AA,index2AA)=temp;
            end
            end
        else
            L1=len1; L2=len2; AA=A;
        end
    else
        error('not implemented')
    end
    [w] = W_JAC_gen(AA,3,N);
else
    w=0;
end

end







function [P,R,nx,ny] = Proj_Gen2D(stenX,stenY,nx,ny)
    if fix(length(stenY)/2)==length(stenY)/2
        centerY=-length(stenY)/2:(length(stenY)-2)/2;
    else
        centerY=-(length(stenY)-1)/2:(length(stenY)-1)/2; 
    end
    if fix(length(stenX)/2)==length(stenX)/2
        centerX=-length(stenX)/2:(length(stenX)-2)/2;
    else
        centerX=-(length(stenX)-1)/2:(length(stenX)-1)/2; 
    end

    if isempty(centerX)
        Px=speye(nx);
    else
        Px=spdiags(ones(nx,1)*stenX,centerX,nx,nx);
        Px=Px(2:2:end,:); 
    end
    if isempty(centerY)
        Py=speye(ny);
    else
        Py=spdiags(ones(ny,1)*stenY,centerY,ny,ny);
        Py=Py(2:2:end,:);
    end
    ny=min(size(Py)); nx=min(size(Px));
    R=kron(Py,Px);
%     if gal==0
        P=R'/norm(R,1);
        R=R/norm(R,inf);
%     else
%         P=R';
%     end
end