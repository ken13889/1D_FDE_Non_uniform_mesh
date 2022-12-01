function [x,M] = Smoother(A,x,b,it,w,smoother,GMRES,Iter_matrix)
%possible choices are
% j = jacobi
% gsl = lower gauss seidel
% gsu = upper gauss seidel
% sorl = lower successive over relaxation
% soru = upper successive over relaxation
% b+number = band approx of A. Number must be odd (example b5)
% ilu = incomplete LU factorization
% ilup = incomplete LU factorization with pivoting of [-1 0 1] (beta=1 and gamma=1)
% 
%
% GMRES = 0/1 means Smoother as GMRES preconditioner
% Iter_matrix = 0/1 implies the computation of the iteration matrix and
% eigenvalues

n=length(A);
pre2=[];

if Iter_matrix==1 && n>2^10
    error('Iteration matrix is too big ( N > 2^10 )')
end



if strcmp(smoother,'j')  % JACOBI
    if GMRES==0
        W=diag(A);
        for i=1:it
            x=x+w*(b-A*x)./W;
        end
    else
        W=spdiags(1./diag(A),0,n,n);
        pre1=@(A,x,b)x+w*W*(b-A*x);
        [x]=gmres(A,b,[],[],it,pre1,[],x);
    end
    if Iter_matrix==1, M=speye(length(A))-w*spdiags(1./diag(A),0,n,n)*A; end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strcmp(smoother(1:2),'gs') % LOWER GAUSS SEIDEL 
    if abs(w-1)>10^-10
        warning('In the case of Lower Gauss-Seidel, "w" must be set to "1"')
    end
    if strcmp(smoother(3),'u')
        W=triu(A);
    else
        W=tril(A); % DEFAULT
    end
    if GMRES==0
        for i=1:it
            x=x+w*(W\(b-A*x));
        end
    else
        pre1=@(A,x,b)x+w*(W\(b-A*x));
        [x]=gmres(A,b,[],[],it,pre1,[],x);
    end
    if Iter_matrix==1, M=speye(n)-W\A; end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strcmp(smoother,'gmres')
    try
        [x,flag]=gmres(@(x)A*x,b,[],eps,it,[],[],x);
    catch err
        [x,flag]=gmres(A,b,[],eps,it,[],[],x);
    end
	if Iter_matrix==1, M='GMRES'; end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strcmp(smoother(1:3),'sor')
    if strcmp(smoother(4),'u')
        W=triu(A);
    else
        W=tril(A); % DEFAULT
    end
    if GMRES==0
        for i=1:it
            x=(1-w)*x+w*(W\(b-A*x));
        end
    else
        pre1=@(A,x,b)(1-w)*x+w*(W\(b-A*x));
        [x]=gmres(A,b,[],[],it,pre1,[],x);
    end
    if Iter_matrix==1, M=(1-w)*speye(n)-w*(W\A); end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strcmp(smoother(1),'b') || strcmp(smoother,'ilu')
    GMRES=1;
    if strcmp(smoother,'ilu')
        warning('ilu is applied to a band approximation... not to the full matrix')
        len=11;  m=(len-1)/2;
    else
        len=round(str2double(post(2:end)));  m=(len-1)/2;
        if mod(len,2)==0 % ven number
            error('band must be an odd number')
        end
    end
    band=zeros(n,len);
    for i=-m:m
        temp=diag(A,i);
        band(1:length(temp),i+m+1)=temp;
    end
    W=spdiags(band,-m:m,n,n);    
    if strcmp(smoother,'ilu')
        setup.type = 'crout';
        setup.milu = 'row';
        setup.droptol = 0.1;
        [L,U]=ilu(W,setup);
        [x]=gmres(A,b,[],[],it,L,U,x);
        if Iter_matrix==1, M=speye(n)-w*(W\A); end
    else
        [x]=gmres(A,b,[],[],it,W,[],x);
        if Iter_matrix==1, M=speye(n)-w*(W\A); end
    end
    
    
end




end

