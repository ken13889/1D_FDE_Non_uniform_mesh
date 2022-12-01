function [x,it,Res] = V(A,x,b,w,tol,P,R,nu1,nu2,maxit,pre,post)

if exist('pre','var')==0
    pre='j'; 
end
preGMRES=0;
if exist('post','var')==0
    post='j'; 
end
postGMRES=0;

lvl=-1;
[rr,cc]=size(A); temp=cell(length(A));
if rr==1 || cc==1
    for i=1:length(A)
        temp{i,i}=A{i};
    end
    A=temp;
end
for i=1:length(A)
    try
    if isempty(A{i,i})==0
        lvl=lvl+1;
    end
    end
end
for i=1:lvl
    n(i)=length(A{i});
end
if numel(w)==1
    W=ones(2,lvl,lvl);
    w=w*W;
elseif numel(w)==2
    W=ones(2,lvl,lvl);
    W(1,:,:)=w(1);
    W(2,:,:)=w(2);
    w=W;
end

if tol<=0 || tol>=1
    maxit=abs(tol);
    tol=-1;
end

Res=inf; e{1}=x; r{1}=b;
while Res(end)>tol*norm(b) && length(Res)<=maxit
    for i=1:lvl
        if i>1
            e{i}=zeros(size(r{i}));
        end
        % PRE 
%         if strcmp(pre,'j')
%             for k=1:nu1
%                 res=r{i}-A{i}*e{i};
%                 e{i}=e{i}+w(1,i,i)*spdiags(1./diag(A{i}),0,n(i),n(i))*res;
%             end
%         else % GS
%             for k=1:nu1
%                 res=r{i}-A{i}*e{i};
%                 e{i}=e{i}+w(1,i,i)*tril(A{i})\res;
%             end
%         end
        [e{i}] = Smoother(A{i,i},e{i},r{i},nu1,w(1,i,i),pre,preGMRES,0);
        % CGC
        r{i+1}=r{i}-A{i,i}*e{i};
        if i==1
            Res(end+1)=norm(r{i+1});
        end
        r{i+1}=R{i,i}*r{i+1};
    end
    e{lvl+1}=A{lvl+1,lvl+1}\r{lvl+1};
    for i=lvl:-1:1
        e{i}=e{i}+P{i,i}*e{i+1};
        % POST
%         if strcmp(post,'j')
%             for k=1:nu2
%                 res=r{i}-A{i}*e{i};
%                 e{i}=e{i}+w*spdiags(1./diag(A{i}),0,n(i),n(i))*res;
%             end
%         else % GS
%             for k=1:nu2
%                 res=r{i}-A{i}*e{i};
%                 e{i}=e{i}+w*tril(A{i})\res;
%             end
%         end
        [e{i}] = Smoother(A{i,i},e{i},r{i},nu2,w(2,i,i),post,postGMRES,0);
    end
%     plot(x)
%     drawnow
%     pause(0.1)
end


it=length(Res)-1;
Res=Res(2:end);

if it==maxit && tol>0
%     disp('maxit reached')
end

x=e{1};
if any(isnan(x))
    it=maxit;
    disp('Diverges!')
end



end