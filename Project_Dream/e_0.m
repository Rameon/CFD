%% CALCULATE TOTAL ENERGY %%
function e = e_0(V)

rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Total energy
e = P/((1.4-1)*rho)+.5*(u^2+v^2);