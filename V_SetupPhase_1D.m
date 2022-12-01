function [A,P,R,w] = V_SetupPhase_1D(A1,pre,post,galerkin,stenX,MGlvl,method,RECALL,A2)
% pre is 'j' 'gs' ,... or a cell array {'j','gs',stenX,MGlvl} which implies
% a pre smoother V-Cycle in direction x (the only possible since
% dimension=1) whose pre smoother is 'j' and post is 'gs'
% same for post

%RECALL,A2  is needed in order not to compute again the optimal twice 
%when Setup_Phase is called due to pre o post smoother V-Cycle

global V_SETUP_PHASE_REMINDER
if exist('RECALL','var')==0
    RECALL=0;
end

%Check
preV=0; postV=0;
if iscell(pre),    preV=1;  end
if iscell(post),   postV=1; end





min=3; %solution on the coarsest grid min * min

if isempty(V_SETUP_PHASE_REMINDER)
    V_SETUP_PHASE_REMINDER=1;
    fprintf(['Setup Phase of V-Cycle: Solution on the coarsest level where the matrix has minimum size ',...
        '%1.0f x %1.0f\r'],min,min)
end

if RECALL==0
    A{1,1}=A1; %on the diag of A there is the hierarchy of the main V-Cycle
    nx=A1.Nx;
else
    A{1,1}=A2;
    nx=RECALL;
end
%tic
i=0;
while nx>min && i<MGlvl
    i=i+1;
    if RECALL==1 && i==1
        method1=method;
        method=0;
    elseif RECALL==1 && i>1
        method=method1;
    end
    %pre-smoother
    if preV==1  % pre smoother is V-Cycle
        [tempA,tempP,tempR,tempW] = V_SetupPhase_1D(A1,pre{1},pre{2},galerkin,pre{3},pre{4},method,nx,A{i,i});
        for k=1:length(tempA{2:end}), A{i+k,i}=tempA{k+1};end
        for k=1:length(tempP), P{i+k,i}=tempP{k};end
        for k=1:length(tempR), R{i+k,i}=tempR{k};end
        w(1,i,i+1:length(tempW(1,2:end))+i)=tempW(1,2:end);
    else %its a standard smoother... add to A the correct weight to make it converge if Jacobi or SOR
        if isempty(pre) || all(pre==0)
            % THERE IS NO PRE SMOOTHER
        else
            if RECALL==0, tempA=A{i,i}; else, tempA=A{i}; end
            if strcmp(pre,'j')
                W = W_JAC_gen(tempA,method,nx);
            elseif strcmp(pre,'gmres')
                W = 0;
            else
                'da implementare'
                W=0;
            end
            if RECALL==0   
                w(1,i,i)=W;
            else    % CHIAMATA RICORSIVA DELLA SETUP PHASE
                w(1,i)=W;
            end
        end
    end
    %post-smoother
    if postV==1  % pre smoother is V-Cycle
        [tempA,tempP,tempR,tempW] = V_SetupPhase_1D(A1,post{1},post{2},galerkin,post{3},post{4},method,nx,A{i,i});
        for k=1:length(tempA(2:end)), A{i,i+k}=tempA{k+1};end
        for k=1:length(tempP), P{i,i+k}=tempP{k};end
        for k=1:length(tempR), R{i,i+k}=tempR{k};end
        w(2,i,i+1:(length(tempW(2,2:end))+i))=tempW(2,2:end);
    else %its a standard smoother... add to A the correct weight to make it converge if Jacobi or SOR
        if isempty(post) || all(post==0)
        % THERE IS NO PRE SMOOTHER
        else
            if strcmp(pre,post)
                if RECALL==0   
                    w(2,i,i)=w(1,i,i);
                else    % CHIAMATA RICORSIVA DELLA SETUP PHASE
                    w(2,i)=w(1,i);
                end
            else
                if RECALL==0, tempA=A{i,i}; else, tempA=A{i}; end
                if strcmp(post,'j')
                    W = W_JAC_gen(tempA,method,nx);
                elseif strcmp(post,'gmres')
                    W = 0;
                else
                    'da implementare'
                    W=0;
                end
            end
            if RECALL==0   
                w(2,i,i)=W;
            else    % CHIAMATA RICORSIVA DELLA SETUP PHASE
                w(2,i)=W;
            end
        end
    end
    % CENTER
%     [tempP,tempR,nx] = Proj_Gen1D(stenX,nx);
    
proj=2;
    if mod(A{i,i}.Nx,2)==1 % se dispari
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TT{i}=A{i,i};
        tt=1./diag(TT{i}.H); % h_1,...,h_n
        ll=length(tt(1:2:end))+(1-mod(TT{i}.Nx,2));
        temp1=spdiags([[tt(2:2:end);0],[0;tt(1:2:end-1)]],0:1,ll,ll);
        sums=full(sum(temp1,2));
    %     sums(1)=1;
        temp1=temp1./sums;
        temp1=temp1(:,2:end);
        temp1(end,:)=0;
        temp1(end,end)=0.5;
        temp2=kron(temp1,[1;0]);
        temp2=temp2(1:end-1,:);
        temp3=kron(speye(ll),[0;1]);
        temp3=temp3(1:end-1,1:end-1);

        tempP=temp2+temp3;
        if size(tempP,1)==length(tt)+1
            tempP=tempP(2:end,:);
        end
        tempR=1/2*tempP';
        nx=size(tempP,2);
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nx=(A{i,i}.Nx-mod(A{i,i}.Nx,2))/2;
%         tempA=BUILDER(A1,nx);
%         grid1=A{i,i}.x(1:A{i,i}.Nx);
%         grid2=tempA.x(1:nx);
        tempP=sparse(A{i,i}.Nx,nx);
        for j=1:A{i,i}.Nx
            if mod(j,2)==1 %dispari
                tempP(j,(j+1)/2)=1;
            else
%                 if mod(A{i,i}.Nx,2)==0 % se pari
                    h1=A{i,i}.h(j);
                    h2=A{i,i}.h(j+1);
%                 else
%                     h1=A{i,i}.h(j+1);
%                     h2=A{i,i}.h(j+2);
%                 end   
                tempP(j,j/2)=h2/(h1+h2);
                tempP(j,j/2+1)=h1/(h1+h2);
            end            
        end
        if abs(norm(tempP,inf)-1)>10^-15
            error('molto strano...')
        end
        if size(tempP,2)==nx+1
            tempP=tempP(:,1:nx);
        elseif size(tempP,2)>nx+1
            error('strano..')
        end
        tempR=tempP'/norm(tempP,1);
    end
%     size(tempP)
%     'asd'
%     if size(tempP,1)~=TT{i}.Nx || size(tempP,2)~=TT{i+1}.Nx
%         TT{i}.Nx
%         TT{i+1}.Nx
%         error('dimension do not match')
%     end
    %%%%%%%%%%%%%%%%%%%%
    
    
    if RECALL==0
        tempA=A{i,i};
    else
        tempA=A{i};
    end
    if galerkin==1
        if i==1 && (isobject(tempA)==1)
            tempA=tempR*Print(tempA)*tempP;
        else
            tempA=tempR*tempA*tempP;
        end
    else
        tempA=BUILDER(A1,nx);  % this yields better results in Paper 4 (caputo) with grid 3
%         tempA=BUILDER(A{i,i},nx);
    end
%     'asd'
    if i==1
        TT{1}=A1;
        TT{2}=tempA;
    else
        TT{i+1}=tempA;
    end
    
    
    if RECALL==0
        P{i,i}=tempP;
        R{i,i}=tempR;
        A{i+1,i+1}=tempA;
    else
        P{i}=tempP;
        R{i}=tempR;
        A{i+1}=tempA;
    end
end

% If band preconditioner is enabled than we assemble it here to speed up V
try
if all(A1.Bandx~=0) && galerkin==0
    [r,c]=size(A);
    for i=1:r
    for j=1:c
        try
            A{i,j}=Print(A{i,j});
        end
    end
    end
end
catch ERR
    warning('Couldnt assemble the band preconditioner because property .BAND of the class BUILDER is missing')
end

end





function [P,R,nx] = Proj_Gen1D(stenX,nx)
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
    nx=min(size(Px));
    R=Px;
%     if gal==0
        P=R'/norm(R,1);
        R=R/norm(R,inf);
%     else
%         P=R';
%     end
end

function v2=DIAG(v)
    len=min(size(v));
    v2=cell(len,1);
    for i=1:len
        v2{i}=v{i,i};
    end
end


