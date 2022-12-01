classdef BUILDER
% Functions:
% Print(obj)  = computes whole matrix (used in multigrid on lowest level)
% mtimes = overloads the classic matrix vector operator '*'
%          product '*' works with the j-th main diag. block if obj.Blick=j
% subsindex(j) = gives the j-th main diagonal block. Use: M(j) 

properties
    Nx; dimension; Bandx; b; Kfun; K; Ix; uL; uR;
    HH; h; ffun; gamma; beta; x; uex;
    Center_Toep_r; Center_Toep_c;
    row_Toep_R; row_Toep_L; Nzx;
    col_Toep_R; col_Toep_L;
    FFT_ToepL;  FFT_ToepR;
    START; END; A; Diag;
    KL; KR; U; H; bandedA;
end
methods
function obj = BUILDER(A,N,HH,f,g,beta,u_ex,uL,uR,k,Ix,Bandx) %constructor
    if nargin==2
        if abs(A.Nx-N)<10^-10
            obj=A;
            return
        end
        % Build any exact copy of A but with size N
        obj.dimension=A.dimension; obj.Bandx=A.Bandx; %obj.Bandy=Bandy;
        obj.Kfun=A.Kfun;  
        obj.Ix=A.Ix;      Ix=A.Ix; 
        obj.uex=A.uex; 
        obj.uL=A.uL;   obj.uR=A.uR;
        obj.ffun=A.ffun;
        obj.gamma=A.gamma; obj.beta=A.beta;
        obj.x=A.x;    obj.h=A.h;
        obj.HH=A.HH;  HH=obj.HH;
        if length(HH)>1  % grid=3
%             p=floor(log2(A.Nx/N));
%             start=floor(2^(sqrt(p)));
% %             start=2^p;
%             if length(     [0,start:2^p:A.Nx-2^p+2,A.Nx+1])==N+2
%                 grid=obj.x([0,start:2^p:A.Nx-2^p+2,A.Nx+1]);
%             elseif length( [0,start+1:2^p:A.Nx-2^p+2,A.Nx+1])==N+2
%                 grid=obj.x([0,start+1:2^p:A.Nx-2^p+2,A.Nx+1]);
%             elseif length( [0,start+2:2^p:A.Nx-2^p+2,A.Nx+1])==N+2
%                 grid=obj.x([0,start+2:2^p:A.Nx-2^p+2,A.Nx+1]);
%             elseif length( [0,start+3:2^p:A.Nx-2^p+2,A.Nx+1])==N+2
%                 grid=obj.x([0,start+3:2^p:A.Nx-2^p+2,A.Nx+1]);
%             elseif length( [0,2:2^p:A.Nx-2^p+2,A.Nx+1])==N+2
%                 'BUILDER line 49..'
%                 grid=obj.x([0,2:2^p:A.Nx-2^p+2,A.Nx+1]);
%             else
%                 error('Errore di lunghezza della griglia')
%             end
            grid=obj.x([0:1:A.Nx+1]);
            while length(grid)~=N+2
%                 grid=grid([1,2:2:end-1,end]);
                grid=grid([1,2:2:end-1,end]);
                if length(grid)==N+2
                    break
                elseif length(grid)==N+3
                    grid=grid([1:end-2,end]);
                    break
                elseif length(grid)<N+2
%                     error('problems with the grid')
                    x2=linspace(0,1,N+2);
                    grid=interp1(linspace(0,1,A.Nx+2),obj.x([0:1:A.Nx+1]),x2,'linear');
%                     plot(linspace(0,1,A.Nx+2),obj.x([0:1:A.Nx+1]))
%                     hold on
%                     plot(x2,grid,'.')
                    break
                end
            end
%             if length([0,2:2^p:A.Nx-1,A.Nx+1])==N+2
%                 grid=obj.x([0,2:2^p:A.Nx-1,A.Nx+1]);
%             elseif length([0,3:2^p:(A.Nx-2^p),A.Nx+1])==N+2
%                 grid=obj.x([0,3:2^p:(A.Nx-2^p),A.Nx+1]);
%             elseif length( [0,2^p:2^p:A.Nx-1,A.Nx+1])==N+2
%                 grid=obj.x([0,2^p:2^p:A.Nx-1,A.Nx+1]);
%             elseif length( [0,2^p:2^p:(A.Nx-2^p),A.Nx+1])==N+2
%                 grid=obj.x([0,2^p:2^p:(A.Nx-2^p),A.Nx+1]);
%             elseif length( [0,2^p:2^p:A.Nx,A.Nx+1])==N+2
%                 grid=obj.x([0,2^p:2^p:A.Nx,A.Nx+1]);
%             else
%                 error('Errore di lunghezza della griglia')
%             end
%             global qq
%             if isempty(qq)
%                 qq=1;
%             else
%                 qq=qq+1;
%             end
%             plot(grid,ones(size(grid))*qq,'.')
%             xlim([0,1])
%             hold on
            skip=1;
        else
            skip=0;
        end
        % MEW MATRIX SIZE
        obj.Nx=N; 
        f=[];
        g=obj.gamma;
        uL=obj.uL;   uR=obj.uR;
        beta=obj.beta;
        k=obj.Kfun;
        u_ex=obj.uex;
        Bandx=A.Bandx;
        obj.Bandx=A.Bandx;
    else
        skip=0;
        if mod(Bandx,2)==0 && Bandx~=0  % odd number
            error('Bandx must be an odd number')
        end
        obj.dimension=1; obj.Bandx=Bandx; %obj.Bandy=Bandy;
        obj.Kfun=k;  obj.Ix=Ix;      obj.uex=u_ex; 
        obj.uL=uL;   obj.uR=uR;
        obj.Nx=N;    obj.ffun=f;
        obj.gamma=g; obj.beta=beta;
        obj.HH=HH;
    end
    if length(HH)==1
        h=1/(N+1);
        unif_grid=(Ix(1)+h):h:(Ix(2)-h);
        grid=[Ix(1),HH(unif_grid),Ix(2)]';
    else
        if skip==0
            stop=0;
            global N1 N2
            if isempty(N1)==0 && isempty(N2)==0
                m=N1;
                N=N1+N2;
                if N~=N1+N2
                    error('N must be the sum between N1 and N2')
                end
            else
                m=max(ceil(HH{1}(N)),1);
            end
%             while stop==0
                N=N-m+1;
                h=1/(N+1);
                unif_grid=(Ix(1)+h):h:(Ix(2)-h);
                grid=[];
                for i=1:m-1
                    grid(i)=h/2^(m-i);
                end
                grid=[Ix(1),grid,unif_grid,Ix(2)]';
                if min(grid(2:end))<eps
                    warning('minimum stepsize is <10^-16')
                end
%                 if min(grid(2:end))<10^-16
%                     m=m-1;
%                     disp('Decreasing the amount of grid point in the first part of the refined mesh, in order to have the smallest step > 10^-16')
%                 else
%                     stop=1;
%                 end
%             end
            N=N+m-1;
        end
    end
    %building variable spacing
    H=zeros(N,1);
    XX=zeros(2*N+1,1);
%     for i=2:N+2
%         H(i-1)=grid(i)-grid(i-1);  % x(1)=Ix(1)  
%     end
    H=grid(2:N+2)-grid(1:N+1);
    XX(2:2:2*N+1)=grid(2:end-1);
    XX(1:2:2*N+1)=grid(2:end)-H/2;


    x=@(i)XX(i*2.*(i>0 & i<N+1)+(i<=0 | i>=N+1)).*...
        [(i>0 & i<N+1)]'+[Ix(1)*(i==0)+Ix(2)*(i==N+1)]'; % X(0.5) gives the mid point between x(1) and x(0)
    obj.x=x;
    obj.U=u_ex(x(1:N));
    h=@(i)H(i.*(i>0 & i<N+1)+(i<=0 | i>=N+1)).*(i'>0 & i'<N+1)+(x(N+1)-x(N))*(i'==N+1)+999*(i'<=0 | i'>N+1);
    obj.H=spdiags((1./h(1:N)),0,N,N);
%     obj.H=spdiags(ones(N,1),0,N,N);
    obj.h=h;
    K=@(i)k(x(i));

    Dphi=zeros(N,2);
    for i=1:N % derivative of the hat basis _/\_
        Dphi(i,1)=1/h(i);  %_/
        Dphi(i,2)=-1/h(i+1);  %\_
    end
    temp=@(i,j)(g*D_caputo(obj,beta,Dphi(j,:),[x(j-1),x(j),x(j+1)],x(i),'L')...
        +(1-g)*D_caputo(obj,beta,Dphi(j,:),[x(j-1),x(j),x(j+1)],x(i),'R'));
    
    % Find the toeplitz structure inside matrix A
    counter=0; START=0; END=0; STOP=0;
%     for i=1:N
%         if abs(h(i+1)-h(i))<10^-14 && STOP==0 && (length(HH)==1 || (length(HH)==2 && i>1))%if h(i-1)==h(i)
%             START=i; STOP=1;
%         elseif abs(h(i+1)-h(i))<10^-14 && ( STOP==1 || STOP==2 )
%             END=i+1; STOP=2;
%         elseif abs(h(i+1)-h(i))>=10^-14 && STOP==2
%             break
%         end
%     end
    for i=N:-1:1
        if abs(h(i+1)-h(i))<10^-14 && STOP==0 && (length(HH)==1 || (length(HH)==2 && i>1))%if h(i-1)==h(i)
            END=i; STOP=1;
        elseif abs(h(i+1)-h(i))<10^-14 && ( STOP==1 || STOP==2 )
            START=i+1; STOP=2;
        elseif abs(h(i+1)-h(i))>=10^-14 && STOP==2
            break
        end
    end
    END=min(N,END-1); % needed to avoid problems with non uniform grid on
                      % the right side of the interval
    if END<0, END=0; end
    % Building first row and  column of the toeplitz
    if START~=0 && END~=0
        len=END-START+1;
        A=spalloc(N,N,(N-len)*(N+len)) ; KR=zeros(len,1); KL=KR;
        D=zeros(N,1);
        row_Toep_R=zeros(len,1); row_Toep_L=row_Toep_R;
        col_Toep_R=row_Toep_R; col_Toep_L=row_Toep_R;
        i=START; k=0;
        for j=START:END
            k=k+1;
            row_Toep_L(k)=temp(i-1/2,j); row_Toep_R(k)=temp(i+1/2,j);
        end
        j=START; k=0;
        for i=START:END
            k=k+1;
            col_Toep_L(k)=temp(i-1/2,j); col_Toep_R(k)=temp(i+1/2,j);
            KL(k)=K(i-1/2);  KR(k)=K(i+1/2);
        end
        % Build the remaining elements (non-uniform grid)
        for i=[1:START-1,END+1:N]
            for j=1:N
                A(i,j)=K(i-1/2)*temp(i-1/2,j)-K(i+1/2)*temp(i+1/2,j);
                if i==j
                    D(i)=K(i-1/2)*temp(i-1/2,i)-K(i+1/2)*temp(i+1/2,i);
                end
            end
        end
        for i=[START:END]
            for j=[1:START-1,END+1:N]
            	A(i,j)=K(i-1/2)*temp(i-1/2,j)-K(i+1/2)*temp(i+1/2,j);
            end
        end
        for i=1:N
            D(i)=K(i-1/2)*temp(i-1/2,i)-K(i+1/2)*temp(i+1/2,i);
        end
        obj.row_Toep_R=row_Toep_R; obj.row_Toep_L=row_Toep_L;
        obj.col_Toep_R=col_Toep_R; obj.col_Toep_L=col_Toep_L;
        obj.START=START;
        obj.END=END;
        [cL,~]=emb_toep_in_circ(col_Toep_L,row_Toep_L,obj);
        [cR,nz]=emb_toep_in_circ(col_Toep_R,row_Toep_R,obj);
        obj.Nzx=nz;
        obj.FFT_ToepL=fft(cL,nz);
        obj.FFT_ToepR=fft(cR,nz);
        obj.KR=KR; obj.KL=KL;
        obj.Diag=obj.H*D;
%         norm(obj.Diag-obj.H*diag(A))
%         'asd'
    else % there is no toeplitz structure inside matrix A
        obj.START=[];          obj.END=[];
        A=zeros(N,N);
        for i=1:N
            for j=1:N
                A(i,j)=K(i-1/2)*temp(i-1/2,j)-K(i+1/2)*temp(i+1/2,j);
            end
        end
        obj.Diag=diag(obj.H).*diag(A);
    end
    obj.A=A; 
    obj.bandedA=[];
    if Bandx>0
        obj.bandedA=Print(obj);
    end
    if isempty(f) && isempty(u_ex)==0
        obj.b=obj*u_ex(x(1:N));
    elseif isempty(f)==0
        temp1=@(i)K(i)*(g*D_caputo(obj,beta,[0,-1/h(1)],[Ix(1),Ix(1),x(1)],x(i),'L')...
            +(1-g)*D_caputo(obj,beta,[0,-1/h(1)],[Ix(1),Ix(1),x(1)],x(i),'R'));
        % when j=N+1
        temp2=@(i)K(i)*(g*D_caputo(obj,beta,[1/h(N+1),0],[x(N),Ix(2),Ix(2)],x(i),'L')...
            +(1-g)*D_caputo(obj,beta,[1/h(N+1),0],[x(N),Ix(2),Ix(2)],x(i),'R'));
        B=zeros(N,1);
%         warning('Usare Simpson per calcolare vettore b - BUILDER riga 164')
        for i=1:N
            B(i)=(x(i+1/2)-x(i-1/2))*f((x(i+1/2)+x(i-1/2))/2)...   %mid point rule (second order approx)
                +temp1(i+1/2)*uL+temp2(i+1/2)*uR...
                -temp1(i-1/2)*uL-temp2(i-1/2)*uR;
        end
        obj.b=obj.H*B;
    end
    [r1,c1]=size(obj.b);   if c1>r1, obj.b=obj.b'; end
end


function [cA,nz] = emb_toep_in_circ(cA,rA,obj)
%creates the first column of the circulant 2^t-dimenional in which the
%toeplitz is embedded
%NOTE: works with any dimension of the toeplitz
%__________________________________________________________________________
% OUTPUT
% cA = first column of the circuland
% nz = number of zeros needed to reach a power of 2: 2*length(cA)+z=2^t
%__________________________________________________________________________
% INPUT
% cA = first column of the toeplitz
% rA = first row of the toeplitz

[r1,c1]=size(cA);
len=c1;
if r1>c1
    cA=cA'; % cA is now row vector
    len=r1;
end
[r2,c2]=size(rA);
if r2>c2
    rA=rA'; % rA is now row vector
end
if r1~=r2 || c1~=c2
    error('Row and column dimension of the toeplitz does not match')
end
d=length(cA);
nz=-1; 
t=1.1;      %needed to get into the while cycle
while t~=fix(t)
    nz=nz+1;
    t=log2(2*d+nz);
end
cA(d+1:d+nz+1)=zeros(nz+1,1);
rA=rA(2:end);
cA(d+2+nz:2^t)=fliplr(rA); %fliplr flips row vectors
cA=cA';
nz=length(cA);
end



function [I] = D_caputo(obj,beta,phi,supp,x,LR)
    %Compute the symbolic left and right caputo derivatives  
    % 0_D_x^(-alpha) and x_D_1^(-alpha), 0<alpha<1
    %
    % works with constant phi only. Can easily be extended to polynomial phi

    % alpha = FD order
    % phi = constant [phi_1,phi_2]  hat basis /\
    % supp = interval containing the suppport of phi   (x_1,x_2,x_3)  where x_2 is the mid point of the hat base /\
    % x = variable where to evaluate the derivative
    % LR = 'L' or 'R' (left or right caputo FD)

    if beta<0 || beta>1
        error('This function requires 0<alpha<1')
    end


    if strcmp(LR,'L') % Ix(1)...x     %int_0^x  phi*(x-s)^alpha  ds
        if x<=supp(1)  % integrating when phi=0
            I=0;
        elseif x>supp(1) && x<=supp(2) %integrating in the middle of the support of phi (when phi is /)
            I=-phi(1)/gamma(beta+1)*(0-(x-supp(1))^beta);
        elseif x>supp(2) && x<=supp(3) %integrating in the middle of the support of phi (when phi is /\)
            I=-phi(1)/gamma(beta+1)*((x-supp(2))^beta-(x-supp(1))^beta)+...  %/
                -phi(2)/gamma(beta+1)*(0-(x-supp(2))^beta);   %\
        else %integrating after the support of phi (phi=0 and does not contribute to the area)
            I=-phi(1)/gamma(beta+1)*((x-supp(2))^beta-(x-supp(1))^beta)+...
                -phi(2)/gamma(beta+1)*( (x-supp(3))^beta-(x-supp(2))^beta );
        end
    elseif strcmp(LR,'R')  % x...Ix(2)   %int_x^1  phi*(s-x)^alpha  ds
        if supp(3)<=x  % integrating when phi=0
            I=0;
        elseif x>supp(2) && x<=supp(3) %integrating in the middle of the support of phi (when phi is \ )
            I=phi(2)/gamma(beta+1)*(supp(3)-x)^beta;
        elseif x>supp(1) && x<=supp(2) %integrating in the middle of the support of phi (when phi is \ )
            I=phi(2)/gamma(beta+1)*((supp(3)-x)^beta-(supp(2)-x)^beta)+... %phi is \ 
                phi(1)/gamma(beta+1)*((supp(2)-x)^beta-0);                      %phi is /
        else %integrating after the support of phi (phi=0 and does not contribute to the area)
            I=phi(2)/gamma(beta+1)*((supp(3)-x)^beta-(supp(2)-x)^beta)+... %phi is \ 
                phi(1)/gamma(beta+1)*((supp(2)-x)^beta-(supp(1)-x)^beta);  %phi is /
        end
    else
        error('wrong LR input')
    end
end



function matrix = Print(obj,r_new,s_new) %,r_new,s_new  
    band=obj.Bandx;
    if isempty(obj.bandedA) || obj.Bandx==0
        if strcmp(band,'L')
            % Scaled laplacian (check if it is correct)
            error('not implemented yet')
             %  REMEMBER THE obj.H*A
            warning('needs a check')
        else
            nx=obj.Nx;
            A2=sparse(nx,nx);
            if isempty(obj.START) && isempty(obj.END)
                matrix=obj.A;
            else
                matrix=obj.A;
                if obj.START*obj.END>sqrt(obj.Nx)
                    matrix=full(matrix);
                end
                r_Toep_R=bander(obj.row_Toep_R,band,obj);
                c_Toep_R=bander(obj.col_Toep_R,band,obj);
                r_Toep_L=bander(obj.row_Toep_L,band,obj);
                c_Toep_L=bander(obj.col_Toep_L,band,obj);
                matrix(obj.START:obj.END,obj.START:obj.END)=...
                    obj.KL.*toeplitz(c_Toep_L,r_Toep_L)-obj.KR.*toeplitz(c_Toep_R,r_Toep_R);
            end
            matrix=obj.H*matrix;
            if band>0
                band=(band-1)/2;
                for i=-band:band
                    temp=diag(matrix,i);
                    if i<0
                        A2=A2+spdiags([temp;zeros(nx-length(temp),1)],i,nx,nx);
                    else
                        A2=A2+spdiags([zeros(nx-length(temp),1);temp],i,nx,nx);
                    end
                end
                matrix=A2;
            end
        end
    else
        matrix=obj.bandedA;
    end
    if nargin>1
        warning('useless selected inputs..')
    end
end

function vec2=bander(vec,band,obj)
    [r,c]=size(vec);
    vec2=sparse(r,c);
    if strcmp(band,'L') %laplacian approx
        vec2(1)=2;
        if length(vec2)>1
            vec2(2)=-1;
        end
%         vec2=vec2/2; %in order to have GL+GR which is a correctly scaled Laplacian 
        warning('pls scale correctly the laplacian approx in bander funciton.. /4 at each level.. not depending on alpha')
    elseif band>=1
        band=(band+1)/2;
        vec2(1:min(end,band))=vec(1:min(end,band));
    elseif band==0
        vec2=vec;
    else
        error('band must be "L" (laplacian approximation) or a number >=0')
    end
end



function SOL = mtimes(obj,x)
%Performs the product matrix*(vector or matrix) considering the structure
nx=obj.Nx;
x_center=x(obj.START:obj.END);
[row,col]=size(x_center);

if numel(x)~=nx
    error(strcat('Dimension mismatch, size of x=[',num2str(size(x))...
        ,'], size of the matrix is [',num2str(nx),']'))
end

if obj.Bandx==0
    if nx>2^3% && obj.Bandx==0 && obj.Bandy==0 
        xhat=fft(x_center,obj.Nzx);
        SOL=zeros(numel(x),1);
        temp1=ifft(obj.FFT_ToepL.*xhat); temp1=temp1(1:row);
        temp2=ifft(obj.FFT_ToepR.*xhat); temp2=temp2(1:row);
        SOL(obj.START:obj.END)=obj.KL.*temp1-obj.KR.*temp2;
        SOL=obj.H*(SOL+obj.A*x);
    else
        SOL=Print(obj)*x;
    end
else
    SOL=obj.bandedA*x;
end
end


function x = mldivide(obj,b)
    x=Print(obj)\b;
end


function x = diag(obj)
    if strcmp(obj.Bandx,'L')
        error('da implementare')
    else
        x=obj.Diag;
    end
end



function [B,obj] = subsref(obj,S) %prints the j-th main diagonal block 
    switch S(1).type
    case '()'
%         B = obj(S.subs);  
        error('nothing here')
    case '.'
%     Toep_cx=obj.Toep_cx;    Toep_rx=obj.Toep_rx;
%     Toep_cy=obj.Toep_cy;    Toep_ry=obj.Toep_ry;
        if length(S)==1
            B = obj.(S.subs);   
        elseif length(S(2).subs)==1
            B = obj.(S(1).subs)(S(2).subs{1});
        elseif length(S(2).subs)==2
            B = obj.(S(1).subs)(S(2).subs{1},S(2).subs{2});
        else
            error('not implemented')
        end
    %case '{}'
    %  B = obj.(S.subs);
    end
end


end
end