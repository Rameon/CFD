clear;clc;

%% LOAD/GENERATE GRID %%
loadGrid = 1;                       %1 - load grid from file
tic
if loadGrid == 1
    fid = fopen('grid/p_0012.grd', 'r');

    npi = fscanf(fid, '%d', 1);
    npj = fscanf(fid, '%d', 1);
    npk = fscanf(fid, '%d', 1);
    
    x = fscanf(fid, '%f', [npi,npj]);
    y = fscanf(fid, '%f', [npi,npj]);
    z = fscanf(fid, '%f', [npi,npj]);
    disp('Grid read successfully');
    
    fclose(fid);
end
x(1,1)=2;

 x=rot90(x,1);
 y=rot90(y,1);

a= npi;
npi=npj;
npj=a;

mesh(x,y,zeros(npi,npj))
axis image;view(2);drawnow;

disp('Grid generated');
toc;disp(' ');

%% GRID METRICS %%
tic
z_width = [0,0,1];
nci = npi-1;                         % Number of cells
ncj = npj-1;

for i = 1:nci
    for j = 1:ncj
        e_xlen = x(i+1,j+1)-x(i+1,j);
        e_ylen = y(i+1,j+1)-y(i+1,j);
        
        n_xlen = x(i,j+1)-x(i+1,j+1);
        n_ylen = y(i,j+1)-y(i+1,j+1);
        
        w_xlen = x(i,j)-x(i,j+1);
        w_ylen = y(i,j)-y(i,j+1);
        
        s_xlen = x(i+1,j)-x(i,j);
        s_ylen = y(i+1,j)-y(i,j);
        
        % For Computational Plane
        xmid(i,j) = (x(i,j) + x(i+1,j))/2;
        ymid(i,j) = (y(i,j) + y(i,j+1))/2;
        
        volume(i,j) = abs(dot(z_width,cross(-1*[s_xlen,s_ylen,0],[e_xlen,e_ylen,0])));
        
        sE(i,j) = sqrt((e_xlen)^2 + (e_ylen)^2);
        sN(i,j) = sqrt((n_xlen)^2 + (n_ylen)^2);
        sW(i,j) = sqrt((w_xlen)^2 + (w_ylen)^2);
        sS(i,j) = sqrt((s_xlen)^2 + (s_ylen)^2);
        
        temp_nE = cross([e_xlen,e_ylen, 0]/sE(i,j), z_width);
        temp_nN = cross([n_xlen,n_ylen, 0]/sN(i,j), z_width);
        temp_nW = cross([w_xlen,w_ylen, 0]/sW(i,j), z_width);
        temp_nS = cross([s_xlen,s_ylen, 0]/sS(i,j), z_width);
        
        nE{i,j} = [temp_nE(1) temp_nE(2)];
        nN{i,j} = [temp_nN(1) temp_nN(2)];
        nW{i,j} = [temp_nW(1) temp_nW(2)];
        nS{i,j} = [temp_nS(1) temp_nS(2)];
        
        clear temp_nE temp_nN temp_nW temp_nS e_xlen n_xlen w_xlen s_xlen e_ylen n_ylen w_ylen s_ylen
    end
end
disp('Grid Metrics calculated');
toc;disp(' ');

%% SET FREESTREAM CONDITIONS %%
free_rho      = 1.2;                            % kg/m^3
free_T        = 300.00;                         % Kelvin
free_P        = 100000;                         % Pa
free_M        = 0.5;                            % Mach
free_aoa      = 3;                              % deg

free_a        = speedsound(free_P,free_rho);    % m/s
free_u        = free_M*free_a*cosd(free_aoa);   % m/s
free_v        = free_M*free_a*sind(free_aoa);   % m/s
free_vel      = [free_u free_v];

% Freestream Primitive State Vector (PSV)
V_free        = [free_rho free_u free_v free_P];

fprintf('Freestream Conditions:\n');
fprintf('Mach:         %10.2f\n',free_M);
fprintf('Flow AOA:     %10.2f deg\n',free_aoa);
fprintf('u Velocity:   %10.2f m/s\n',free_u);
fprintf('v Velocity:   %10.2f m/s\n',free_v);
fprintf('Pressure:     %10.2f Pa\n',free_P);
fprintf('Temperature:  %10.2f K\n',free_T);
fprintf('Density:      %10.2f kg/m^3\n\n',free_rho);

%% SET ITERATION VARIABLES %%
iterations   = 100;
gbl_timestep = 0;
timestep     = 1e-5;    
CFL          = 0.5;
m_stage      = 4;               % m-stage time stepping
freq         = 10;
plotcontours = 0;
% 0 - off
% 1 - Density
% 2 - U Velocity
% 3 - V Velocity
% 4 - Pressure

plots_on     = 1 ;  

fprintf('Iteration Variables:\n');
fprintf('Iterations:   %5d\n',iterations);
fprintf('CFL:          %5.2f\n',CFL);
fprintf('M-Stage:      %5d\n',m_stage);
if gbl_timestep == 0
    fprintf('Timestepping  %5s\n\n','local');
else
    fprintf('Timestepping %5s\n\n','global');
end

%% INITIALIZE SOLUTION %%
tic
resid_i      = 0;
resid_0      = 0;
start_iter   = 0;
end_iter     = 0;
residReduced = 0;
fm = 1;

normals = {nE nN nW nS};
areas   = {sE sN sW sS};

con_density = zeros([nci,ncj]);
con_uvel    = zeros([nci,ncj]);
con_vvel    = zeros([nci,ncj]);
con_pres    = zeros([nci,ncj]);

for i = 1:nci
    for j = 1:ncj
        V{i,j}     = V_free;
        U{i,j}     = convV_U(V{i,j});
        resid{i,j} = [0 0 0 0];
    end
end

if loadGrid == 0
    V{10,5}(4) = 1.1*free_P;f
    U{10,5}    = convV_U(V{10,5});
end

diary off
disp('Solution initialized');
toc

%% MAIN LOOP %%
start_iter = start_iter + 1;

% Create plot window
if plotcontours == 1
    hc = figure('name','Density Contour');
elseif plotcontours == 2
    hc = figure('name','U Velocity Contour');
elseif plotcontours == 3
    hc = figure('name','V Velocity Contour');
elseif plotcontours == 4
    hc = figure('name','Pressure Contour');
end
tic

for iter = start_iter:(end_iter + iterations)
    ti1 = cputime;
    resid_i = 0;
    U0 = U;
    
    % M-stage timestepping
    for m = 1:m_stage
        resid = calcResid(V, V_free, normals, areas, nci, ncj);

        for i = 1:nci
            for j = 1:ncj
                if gbl_timestep == 0
                    vel = [V{i,j}(2) V{i,j}(3)];
                    cell_a = speedsound(V{i,j}(4),V{i,j}(1));
                    dt(1) = CFL * sE(i,j)/(abs(vel(1)*nE{i,j}(1) + vel(2)*nE{i,j}(2))+cell_a);
                    dt(2) = CFL * sN(i,j)/(abs(vel(1)*nN{i,j}(1) + vel(2)*nN{i,j}(2))+cell_a);
                    dt(3) = CFL * sW(i,j)/(abs(vel(1)*nW{i,j}(1) + vel(2)*nW{i,j}(2))+cell_a);
                    dt(4) = CFL * sS(i,j)/(abs(vel(1)*nS{i,j}(1) + vel(2)*nS{i,j}(2))+cell_a);
                    timestep = min(dt);
                end

                U{i,j} = U0{i,j} - 1/(m_stage-(m-1))*timestep/volume(i,j)*resid{i,j};
                V{i,j} = convU_V(U{i,j});
                
                con_density(i,j) = V{i,j}(1);
                con_uvel(i,j)    = V{i,j}(2);
                con_vvel(i,j)    = V{i,j}(3);
                con_pres(i,j)    = V{i,j}(4);
                
                resid_i = resid_i + resid{i,j}.^2;
            end
        end
    end
    
    resid_i = (resid_i).^.5/(nci*ncj);
    
    if iter == 1
        resid_0 = resid_i;
        fprintf('\nIter  cont. resid  x-mom resid  y-mom resid  energy resid  Time left\n');
    end
    
    if isnan(resid_i/resid_0)
        break;
        disp('Solution corrupt.');
    end
    
    if ((resid_i(2)/resid_0(2)) >= (1e1))
        if residReduced == 0
            CFL = CFL/2;
            notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
            disp(notice);
            residReduced = residReduced + 1;
        end
    end

    ti2 = cputime-ti1;
    
    if mod(iter,freq) == 0
        fprintf('%4d  %11.2e  %11.2e  %11.2e   %11.2e  %7s\n',iter,...
            resid_i(1)/resid_0(1),...
            resid_i(2)/resid_0(2),...
            resid_i(3)/resid_0(3),...
            resid_i(4)/resid_0(4),...
            fixTime(ti2*(end_iter+iterations-iter)));
    end
end

start_iter = iter;
end_iter = iter;
% toc


%% PLOT PRIMITIVE STATE VECTOR CONTOURS %%
% Plot 1: Density (kg/m^3)
% Plot 2: U Velocity (m/s)
% Plot 3: V Velocity (m/s)
% Plot 4: Pressure (Pa)

if plots_on == 1
    figure('name','Primitive State Variables');
    
    subplot(221);
    contourf(xmid',ymid',con_density');
    axis image;
    colorbar('peer',gca,'SouthOutside');
    title('Density (kg/m^3)')
    
    subplot(222);
    contourf(xmid',ymid',con_uvel');
    axis image;
    colorbar('peer',gca,'SouthOutside');
    title('U Velocity (m/s)')
    
    subplot(223);
    contourf(xmid',ymid',con_vvel');
    axis image;
    colorbar('peer',gca,'SouthOutside');
    title('V Velocity (m/s)')
    
    subplot(224);
    contourf(xmid',ymid',con_pres');
    axis image;
    colorbar('peer',gca,'SouthOutside');
    title('Pressure (Pa)')
end








