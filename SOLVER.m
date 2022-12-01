function [x,flag,RELRES,it,RESVEC] = SOLVER(A,b,tol,maxit,x0,sol,f,Forced_Res_Tol)
% flag=1 -> qualcosa non va
%GRANTS REQUIRED RESIDUAL TOLERANCE !!!!!

% if isobject(A)==1, n=maxit; else, n=length(A); end
% n=maxit;
if nargin<8
    Forced_Res_Tol=1;
end

% RESTART=100; %to reduce ram requirements
% MAX_RAM=25*10^3; %MegaByte  %server has >100gb of ram
% try size = whos('x0'); size_var = size.bytes/1000^2; %in MegaByte
%     RESTART=min(round(MAX_RAM/size_var),RESTART);
% end


tol1=tol; maxitt=maxit;
tol=tol/10;
if isempty(f) || nargin<7
    f=[];
end
it=0; x=x0; res=1; tryagain=0; Forced_EXIT=1; flag=0;
while (res>tol1 && it<maxitt) && Forced_EXIT && flag~=1
    maxit=maxitt-it;
    if sol==1 %gmres
        [x,flag,RELRES,iter2,RESVEC]=gmres(A,b,[],tol,min(maxit,length(b)),f,[],x);
        if iter2(2)==0 && flag==0
            [x,~,RELRES,iter2,RESVEC]=gmres(A,b,[],eps,1,f,[],x);
        end
%         len=length(RESVEC)-1;
        len=iter2(2);%+RESTART*(iter2(1)-1);
    elseif sol==2 %bicgstab
        [x,flag,RELRES,iter,RESVEC]=bicgstab(A,b,tol,maxit,f,[],x);
        if iter==0 && flag==0
            [x,~,RELRES,iter,RESVEC]=bicgstab(A,b,eps,1,f,[],x);%gmres(A,b,[],eps,1,f,[],x);
        end
        len=ceil((length(RESVEC)-1)/2);
    elseif sol==3 %cgnr
        [x,flag,RELRES,iter,RESVEC]=pcg(A,b,tol,maxit,f,[],x);
        len=length(RESVEC)-1;
    else
        error('only 3 solver')
    end
    if tol<eps
        flag=1;
        it=maxitt;
        break
    else
        it=it+len;
    end
    
    if flag==3  && tryagain==0  %stagnate
        tryagain=1;
    end
    res=norm(A(x)-b)/norm(b);
    if Forced_Res_Tol==1
        if flag~=3 && flag~=4 %if it stagnates there is no need to reduce the tolerance 
            tol=tol/10;
%             tol=tol/(res/tol);
        end
    else
        if tryagain==1 || flag~=3
            Forced_EXIT=0;
        end
    end   
end



end

