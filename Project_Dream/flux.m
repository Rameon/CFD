%% CALCULATE FULL FLUX VECTOR %%
function f = flux(V,n)

rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Contravariant velocity
conV  = u*n(1) + v*n(2);

f(1) = rho*conV;
f(2) = rho*u*conV + P*n(1);
f(3) = rho*v*conV + P*n(2);
f(4) = rho*h_0(V)*conV;