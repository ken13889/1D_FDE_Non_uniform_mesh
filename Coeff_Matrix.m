function [A,HH,Hder]=Coeff_Matrix(beta,N,ex,s,Band,g,grid,gg,length_connection,force_g)

% beta = Fractional derivative
% ex = example
% s = spacing [0.1 0.2] means 0.1 of non uniform grid left and 0.2 right
% global length_connection
% if isempty(smallesth)
    smallesth=10^(-16);
% end
if grid>2
    if isa(s, 'function_handle')
        temp=s; c=1;
    else
        temp=@(x)sqrt(x); c=1;
    end
else
	logg=@(x,base)log(x)/log(base);
    base=1/(N+1);
    c=logg(smallesth,base);
%     gg_min=min(c,gg);
    temp=(2-(1-beta))/(1-beta);
    temp=min(c,temp);
end
if nargin<=6 %no grid gg and length_connection
    grid=1;
elseif nargin<=7 %no gg and length_connection
    gg=temp;
    length_connection=0; force_g=0;
elseif nargin<=8 %no length_connection
    if isempty(gg)
        gg=temp;
    end
    length_connection=0; force_g=0;
elseif nargin<=9 %no force_g
    if isempty(gg)
        gg=temp;
    end
    if isempty(length_connection)
        length_connection=0;
    end
    force_g=0;
else
    if isempty(gg)
        gg=temp;
    end
    if isempty(length_connection)
        length_connection=0;
    end
    if isempty(force_g)
        force_g=0;
    end
end

if grid==1
if abs(gg-min(c,gg))>eps*5 && force_g==0
    disp(['changing grid from power=',num2str(gg),' to power=',num2str(c)])
    gg=c;
end
end



if ex==1
    if g<1
        f=@(x)(1-g).*(1-beta)./(gamma(beta).*x.*(1-x).^(1-beta));%+x.*(1-x);
    else
        f=@(x)x*0;%+x.*(1-x);
    end
    u_ex=@(x)x.^(1-beta);
    uL=0;
    uR=1;
    k=@(x)1+x*0;
    Ix=[0 1];  % interval
elseif ex==2
    if nargin<6
        g=0.5; %gamma
    end
    f=[];
    u_ex=@(x)(1-x).^(1-beta).*(x).^(1-beta);
    uL=0;
    uR=0;
    k=@(x)1+x*0;
    Ix=[0 1];  % interval
elseif ex==3
    'fixed gamma'
    g=0.5;
    d2au=@(x,a)-2*((x.^a+(1-x).^a)/(gamma(a+1))...
        -(6*(x.^(a+1)+(1-x).^(a+1)))/gamma(a+2)...
        +(12*(x.^(a+2)+(1-x).^(a+2)))/gamma(a+3));%d^(2-a)u1
    u1=@(x)x.^2.*(1-x).^2;
    u_ex=@(x)x.^2.*(1-x).^2;
    uL=0;
    uR=0;
    k=@(x)1+x*0;
    f=@(x,y)1/2*k(x).*d2au(x,beta);
    Ix=[0 1];  % interval
    
elseif ex==4
	eta=@(a)1;%1/(2*cos((1-a)*pi/2));
    
    d1au=@(x,a)-2*((x.^(a+1)-(1-x).^(a+1))/(gamma(a+2))...
        -(6*(x.^(a+2)-(1-x).^(a+2)))/gamma(a+3)...
        +(12*(x.^(a+3)-(1-x).^(a+3)))/gamma(a+4));%d^(2-a)u1
    d2au=@(x,a)-2*((x.^a+(1-x).^a)/(gamma(a+1))...
        -(6*(x.^(a+1)+(1-x).^(a+1)))/gamma(a+2)...
        +(12*(x.^(a+2)+(1-x).^(a+2)))/gamma(a+3));%d^(2-a)u1
    u1=@(x)x.^2.*(1-x).^2;
    
    
    
    Kx=@(x,y)dx+0*x;
    dxKx=@(x,y)0*x;%derivative of Kx with respect to x
  
    Ky=@(x,y)dy+0*x;
    dyKy=@(x,y)0*x;%derivative of Kx with respect to x
    
    v=@(x,y)eta(alpha)*Kx(x,y).*d2au(x,a);
%     v=NaN;
    u_exact=@(x,y)u1(x);
    %intervallo di definizione
    Ix=[0,1]; %intervallo spaziale
else
    error('Wrong example')
end



if grid==1 %singular grid
    [HH,Hder] = GridGen(gg,Ix,s,length_connection);
elseif grid==2
    bb=beta;%+0.05
    HH=@(x)bb*x.^3+(1-bb)*x;
    Hder=@(x)bb*3*x.^2+(1-bb);
else
    if iscell(s)
        s=s{1};
    elseif isa(gg, 'function_handle')
        s=gg;
    end
    HH={s,1}; % graded mesh with the first interval [0,h] subdivided in half, than half,...  a "grid" amount of times
    Hder=NaN;
end

A=BUILDER([],N,HH,f,g,beta,u_ex,uL,uR,k,Ix,Band);









end


