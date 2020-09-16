%% CALCULATE SPEED OF SOUND %%
function c = speedsound(P,rho)
% speed of sound = sqrt(gamma*pressure/rho)
c = sqrt(1.4*P/rho);