function pdex4
m = 0;
x = linspace(0,1,101);
t = linspace(0,240000,51);


sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

surf(x,t,u) 
title('Concentration of Dye in a Closed Pipe')
xlabel('Distance x')
ylabel('Time t')

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
D = 4.34*10^-6;
c = 1/0.1;
f = D*DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = pdex1ic(x)
if (0.45 < x) && (x < 0.55)
    u0 = 1;
else
    u0 = 0;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;