%% CALCULATE TOTAL ENTHALPY %%
function h = h_0(V)

rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Total enthalpy
h = e_0(V)+P/rho;