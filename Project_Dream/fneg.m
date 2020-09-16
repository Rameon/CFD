%% CALCULATE NEGATIVE FLUX VECTOR %%
function f = fneg(V,n)

rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

gamma = 1.4;
a     = speedsound(P,rho);
conV  = u*n(1) + v*n(2);        % Contravariant velocity
M     = conV/a;                 % Contravarient Mach number

f(1)  = -rho*a/4*(M-1)^2;
f(2) = f(1) * (u + n(1)*(-conV-2*a)/gamma);
f(3) = f(1) * (v + n(2)*(-conV-2*a)/gamma);
f(4) = f(1) * (h_0(V) - a^2*(M+1)^2/(gamma + 1));