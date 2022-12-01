function [H,Hder] = GridGen(power,Ix,s,length_connection)
%generates grid with a lot of points near the edges of the domain Ix
% the bigger "power" is the more points we have near the edges of Ix
%
%Ix=[a,b] (domain), s=zone near the boundaries where to put more grid points
if power<1
    error('power must be >=1')
end

if nargin<=3
    lc=0.05;
else
    if isempty(length_connection)
        lc=0.05;
    else
        lc=length_connection;
    end
end
if length(s)==3
    lc=s(2);
    s=[s(1),s(3)];
end
x0=Ix(1);
x1=Ix(1)+s(1);
x2=Ix(2)-s(2);
x3=Ix(2);
if x1>x2
    error('s(1)+s(2)>I(2)-I(1)... It must be lower than the length of the interval')
end

% f1=@(x,a,b)(x-a)/(b-a); %from [a,b] to [0,1]
% f2=@(x,a,b)2*x-1;       %from [0,1] to [-1,1]
% f3=@(x,a,b)-x;          %from [-1,1] to [-1,1], mirror!!
% f4=@(x,a,b)(x+1)/2;     %from [-1,1] to [0,1]
% f5=@(x,a,b)x*(b-a)+a;   %from [0,1] to [a,b]
% temp=@(x)p1(f4(f3(f2(f1(x,x2,Ix(2))))));
% p3=@(x)f5(temp,x2,Ix(2));  % NOT WORKING

if s(1)==0 && s(2)==0 %uniform grid
	H=@(x)x.*(x>=Ix(1) & x<=Ix(2));
    Hder=@(x)1+0*x;
elseif s(1)==0 && s(2)~=0
    p1=@(x)x.^power;  
    mod1=@(x)x-x2+Ix(1); % x in [x2,Ix(2)]-> [Ix(1),x1]
    mod2=@(x)(Ix(1)+s(2))-x; %inversion of the input
    p3=@(x)Ix(2)-p1(mod2(mod1(x)));
    p2=@(x)(x-x1)*(p3(x2)-x1)/(x2-x1)+x1; % bind p1 to p3 with a line
    H=@(x)p2(x).*(x>=Ix(1) & x<x2)+p3(x).*(x>=x2 & x<=Ix(2));
    d3=matlabFunction(diff(sym(p3)));
    Hder=@(x)((p3(x2)-x1)/(x2-x1)).*(x>=Ix(1) & x<x2)+d3(x).*(x>=x2 & x<=Ix(2));
elseif s(1)~=0 && s(2)==0
    if x1>=x2 || lc==0  % THIS generates a non smooth function
        p1=@(x)(x-Ix(1)).^power+Ix(1);  
        if s(1)<Ix(2)
            p2=@(x)(x-x1)*(Ix(2)-p1(x1))/(x2-x1)+p1(x1); % bind p1 to p3 with a line
        else
            p2=@(x)x;
        end
        H=@(x)p1(x).*(x>=Ix(1) & x<=x1)+p2(x).*(x>x1 & x<=x2);
        d1=matlabFunction(diff(sym(p1)));
        Hder=@(x)d1(x).*(x>=Ix(1) & x<=x1)+feval(matlabFunction(diff(sym(p2)))).*(x>x1 & x<=x2);
    else
        % x11 is between x1 and x2, determines the length of the second degree polynomial which connect the first poly with the line
        x11=x1+lc;
        p1=@(x)(x-Ix(1)).^power+Ix(1);     dp1=@(x)power*(x-Ix(1)).^(power-1);  %derivative of p1
        V=[x1^2     x1      1    0    0;
           2*x1      1      0    0    0;
           x11^2     x11      1   -x11  -1;
           2*x11      1      0    -1   0;
             0       0      0    x3   1];
        b=[p1(x1) dp1(x1) 0 0 x3]';
        abc=V\b;
        p2=@(x)abc(1)*x.^2+abc(2)*x+abc(3);  dp2=@(x)2*abc(1)*x+abc(2);
        line=@(x)abc(4)*x+abc(5);           dline=@(x)abc(4);
        H=@(x)p1(x).*(x>=Ix(1) & x<=x1)+p2(x).*(x>x1 & x<=x11)+line(x).*(x>x11 & x<=x3);
%         fplot(H,[0,1])
        Hder=@(x)dp1(x).*(x>=Ix(1) & x<=x1)+dp2(x).*(x>x1 & x<=x11)+dline(x).*(x>x11 & x<=x3);
    end
else
    p1=@(x)(x-Ix(1)).^power+Ix(1);  
    p11=@(x)x.^power;  
    mod1=@(x)x-x2+Ix(1); % x in [x2,Ix(2)]-> [Ix(1),x1]
    mod2=@(x)Ix(1)+s(2)-x; %inversion of the input
    mod3=@(x)x;%(x-x2)/(s(2)/s(1))+x2;  %from [x2,Ix(2)] to [x2,Ix(2)+(s(1)-s(2))]
    p3=@(x)Ix(2)-p11(mod2(mod1(x)));
    p2=@(x)(x-x1)*(p3(x2)-p1(x1))/(x2-x1)+p1(x1); % bind p1 to p3 with a line
    H=@(x)p1(x).*(x>=Ix(1) & x<=x1)+p2(x).*(x>x1 & x<x2)+p3(x).*(x>=x2 & x<=Ix(2));
    d1=matlabFunction(diff(sym(p1)));
    d2=feval(matlabFunction(diff(sym(p2))));
    d3=matlabFunction(diff(sym(p3)));
    Hder=@(x)d1(x).*(x>=Ix(1) & x<=x1)+d2.*(x>x1 & x<x2)+d3(x).*(x>=x2 & x<=Ix(2));
end




% grid=Ix(1):0.01:Ix(2);
% % plot(x2:0.01:x3,p3(x2:0.01:x3))
% plot(grid,H(grid))
% 'asd'

end

