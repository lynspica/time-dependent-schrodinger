% Numerical analysis of Time-Dependent Schrödinger equation
% Written by: Ahmet Burak Yıldırım, 2023

% i hbar d(psi(r,t))/dt = -(hbar^2)/(2m) nabla^2 psi(r,t) + V(r) psi(r,t)

% 1 dimensional -> r = (x):
% i hbar d(psi(x,t))/dt = -(hbar^2)/(2m) (d^2/dx^2) psi(x,t) + V(x) psi(x,t)

% 2 dimensional -> r = (x,y):
% i hbar d(psi(x,y,t))/dt = -(hbar^2)/(2m) (d^2/dx^2 + d^2/dy^2) psi(x,y,t) + V(x,y) psi(x,y,t)

%% Time-dependent Schrödinger equation
% Parabolic PDE

% OBSERVATIONS:
% EE: dt = 1e-4, tmax = 10, oscillations for finite-well, 
% 1 - Explicit Euler 
% 2 - Crank-Nicolson
%   2a - Jacobi
%   2b - Gauss-Seidel
%   2c - SOR
%   2d - ADI

clear; close all; clc;

% Initialization of parameters

dimension = 1; % 1 or 2

% Potential choices are applicable according to selected dimension
% List of available potentials:
% FOR 1D: harmonic / finite-well /  tunneling  / free-particle
% FOR 2D: harmonic / finite-well / double-slit / free-particle

potential1d = "tunneling"; % harmonic / finite-well / tunneling / free-particle
potential2d = "single-slit";   %double-slit / harmonic / finite-well / tunneling / free-particle / single-slit

m    = 1;         % mass of the particle (in amu)
hbar = 1;         % reduced Planck's constant (in eV*fs)

% Parameter table 
% 1D - ExplicitEuler       : dt: TRUE: 2e-5 / CASES: 8e-5, 2e-4, 8e-4, 2e-3
% 1D - Crank-Nicolson      : dt: TRUE: 2e-5 / CASES: 8e-5, 2e-4, 8e-4, 2e-3
% 2D - ExplicitEuler       : dt: TRUE: 2e-5 / CASES: 8e-5, 2e-4, 8e-4, 2e-3
% 2D - Crank-Nicolson (ADI): dt: TRUE: 2e-5 / CASES: 8e-5, 2e-4, 8e-4, 2e-3

dt      = 2e-4;          % time step size (in fs)
dt_true = 2e-5;          % time step size for "true" solution (in fs)

tmax = 1;                % maximum time (in fs)

t = 0:dt:tmax;           % time grid
t_true = 0:dt_true:tmax; % time grid

if dimension == 1
    dx = 0.1;      % spatial step size in the x direction (in nm)
    Lx = 10;       % size of the grid in the x direction (in nm)
    x = -Lx/2:dx:Lx/2; 
    X = reshape(x,[size(x,2) 1]); % 1D grid

    kx = 0.2;         % wave vector: for momentum 
    r0 = -2; sigma = 0.5;
    psi_0 = exp(-0.5*(((X-r0))/sigma).^2).*exp(1i*(kx*(X)));
    psi_0 = Normalizewavefunction1D(psi_0);

    psi = zeros(length(x),length(t));
    % set the initial wave function
    psi(:,1) = psi_0; 

    V = potential1D(potential1d,x);
    
    % TRUE SOLUTION (using smaller time step)

    psi_true = zeros(length(x),length(t_true));
    % To compare true solution psi_true with psi, creating a new arraw
    % and filling same time steps
    psi_true(:,1) = psi_0; 
 
    % Solve the time-dependent Schrödinger equation using the finite difference method
    disp("Solving TDSE using the finite difference method")
    tic
    for k = 2:length(t)

        % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
        % v = (i hbar)/(2m) -> diffusion coefficient
        % f = psi(r,t)
        % df/dt = v (d^2f/dx^2) + V(x)/(i hbar) f

        v = (1i*hbar)/(2*m); % diffusion coefficient
        psi_t = psi(:,k-1);

        pot = V/(1i*hbar); 
        % Available solvers: ExplicitEuler1D / CrankNicolson1D
        psi(:,k) = CrankNicolson1D(psi_t, hbar, m, dx, dt, pot); % replace function
        psi(:,k) = Normalizewavefunction1D(psi(:,k));

    end
    toc
    disp("Solving for the true solution")
    tic
    % Solve for the true solution
    for k = 2:length(t_true)

        % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
        % v = (i hbar)/(2m) -> diffusion coefficient
        % f = psi(r,t)
        % df/dt = v (d^2f/dx^2) + V(x)/(i hbar) f

        v = (1i*hbar)/(2*m); % diffusion coefficient
        psi_t = psi_true(:,k-1);

        pot = V/(1i*hbar); 
        psi_true(:,k) = ExplicitEuler1D(psi_t, hbar, m, dx, dt_true, pot);
        psi_true(:,k) = Normalizewavefunction1D(psi_true(:,k));

    end
    toc
    psi_true_c = psi_true(:,1:length(t_true)/length(t):end);
end

if dimension == 2
    dx = 0.1;      % spatial step size in the x direction (in nm)
    dy = 0.1;      % spatial step size in the y direction (in nm)
    Lx = 10;       % size of the grid in the x direction (in nm)
    Ly = 10;       % size of the grid in the y direction (in nm)
    x = -Lx/2:dx:Lx/2; 
    y = -Ly/2:dy:Ly/2;
    [X,Y] = meshgrid(x,y); % 2D grid

    k = [10, 0];   % wave vector: for momentum 
    r0 = [-3, 0]; sigma = 0.5;
    psi_0 = exp(-((X-r0(1)).^2+(Y-r0(2)).^2)/(2*sigma)).*exp(1i*(k(1)*X+k(2)*Y));
    psi_0 = Normalizewavefunction2D(psi_0);
    psi = zeros(length(x),length(y),length(t));
    % set the initial wave function
    psi(:,:,1) = psi_0; 

    V = potential2D(potential2d,x,y);

    % TRUE SOLUTION (using smaller time step)

    psi_true = zeros(length(x),length(y),length(t_true));
    % To compare true solution psi_true with psi, creating a new arraw
    % and filling same time steps
    psi_true(:,:,1) = psi_0; 

    % Solve the time-dependent Schrödinger equation using the finite difference method
    disp("Solving TDSE using the finite difference method")
    tic
    for k = 2:length(t)

        % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
        % v = (i hbar)/(2m) -> diffusion coefficient
        % f = psi(r,t)
        % df/dt = v (d^2f/dx^2 + d^2f/dy^2) + V(x,y)/(i hbar) f

        v = (1i*hbar)/(2*m); % diffusion coefficient
        psi_t = psi(:,:,k-1);
        
        pot = V/(1i*hbar);
        % Available solvers: ExplicitEuler2D / CrankNicolson2D (ADI)
        psi(:,:,k) = ExplicitEuler2D(psi_t, hbar, m, dx, dy, dt, pot);
        psi(:,:,k) = Normalizewavefunction2D(psi(:,:,k));

    end
    toc
    disp("Solving for the true solution")
    tic
    % Solve for the true solution
    for k = 2:length(t_true)

        % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
        % v = (i hbar)/(2m) -> diffusion coefficient
        % f = psi(r,t)
        % df/dt = v (d^2f/dx^2 + d^2f/dy^2) + V(x,y)/(i hbar) f

        v = (1i*hbar)/(2*m); % diffusion coefficient
        psi_t = psi_true(:,:,k-1);

        pot = V/(1i*hbar);
        psi_true(:,:,k) = ExplicitEuler2D(psi_t, hbar, m, dx, dy, dt_true, pot);
        psi_true(:,:,k) = Normalizewavefunction2D(psi_true(:,:,k));

    end
    toc
    psi_true_c = psi_true(:,:,1:length(t_true)/length(t):end);
end

disp("Stability condition: alpha = " + dt*abs(v)/dx^2)

%% Save data: Dimension = 1, Method = Explicit Euler
% Run after each simulation

if dt == 8e-5
    X_1DEULER_dt8e5 = X;
    V_1DEULER_dt8e5 = V;
    t_1DEULER_dt8e5 = t;
    psi_1DEULER_dt8e5 = psi;
    psi_true_c_1DEULER_dt8e5 = psi_true_c;
    save 1DEULER_dt8e5.mat V_1DEULER_dt8e5 X_1DEULER_dt8e5 t_1DEULER_dt8e5 ...
        psi_1DEULER_dt8e5 psi_true_c_1DEULER_dt8e5
    clear
elseif dt == 2e-4
    X_1DEULER_dt2e4 = X;
    V_1DEULER_dt2e4 = V;
    t_1DEULER_dt2e4 = t;
    psi_1DEULER_dt2e4 = psi;
    psi_true_c_1DEULER_dt2e4 = psi_true_c;
    save 1DEULER_dt2e4.mat V_1DEULER_dt2e4 X_1DEULER_dt2e4 t_1DEULER_dt2e4 ...
        psi_1DEULER_dt2e4 psi_true_c_1DEULER_dt2e4
    clear
elseif dt == 8e-4
    X_1DEULER_dt8e4 = X;
    V_1DEULER_dt8e4 = V;
    t_1DEULER_dt8e4 = t;
    psi_1DEULER_dt8e4 = psi;
    psi_true_c_1DEULER_dt8e4 = psi_true_c;
    save 1DEULER_dt8e4.mat V_1DEULER_dt8e4 X_1DEULER_dt8e4 t_1DEULER_dt8e4 ...
        psi_1DEULER_dt8e4 psi_true_c_1DEULER_dt8e4
    clear
elseif dt == 2e-3
    X_1DEULER_dt2e3 = X;
    V_1DEULER_dt2e3 = V;
    t_1DEULER_dt2e3 = t;
    psi_1DEULER_dt2e3 = psi;
    psi_true_c_1DEULER_dt2e3 = psi_true_c;
    save 1DEULER_dt2e3.mat V_1DEULER_dt2e3 X_1DEULER_dt2e3 t_1DEULER_dt2e3 ...
        psi_1DEULER_dt2e3 psi_true_c_1DEULER_dt2e3
    clear
end

%% Save data: Dimension = 1, Method = Crank-Nicolson
% Run after each simulation

if dt == 8e-5
    X_1DCN_dt8e5 = X;
    V_1DCN_dt8e5 = V;
    t_1DCN_dt8e5 = t;
    psi_1DCN_dt8e5 = psi;
    psi_true_c_1DCN_dt8e5 = psi_true_c;
    save 1DCN_dt8e5.mat V_1DCN_dt8e5 X_1DCN_dt8e5 t_1DCN_dt8e5 psi_1DCN_dt8e5 psi_true_c_1DCN_dt8e5
    clear
elseif dt == 2e-4
    X_1DCN_dt2e4 = X;
    V_1DCN_dt2e4 = V;
    t_1DCN_dt2e4 = t;
    psi_1DCN_dt2e4 = psi;
    psi_true_c_1DCN_dt2e4 = psi_true_c;
    save 1DCN_dt2e4.mat V_1DCN_dt2e4 X_1DCN_dt2e4 t_1DCN_dt2e4 psi_1DCN_dt2e4 psi_true_c_1DCN_dt2e4
    clear
elseif dt == 8e-4
    X_1DCN_dt8e4 = X;
    V_1DCN_dt8e4 = V;
    t_1DCN_dt8e4 = t;
    psi_1DCN_dt8e4 = psi;
    psi_true_c_1DCN_dt8e4 = psi_true_c;
    save 1DCN_dt8e4.mat V_1DCN_dt8e4 X_1DCN_dt8e4 t_1DCN_dt8e4 psi_1DCN_dt8e4 psi_true_c_1DCN_dt8e4
    clear
elseif dt == 2e-3
    X_1DCN_dt2e3 = X;
    V_1DCN_dt2e3 = V;
    t_1DCN_dt2e3 = t;
    psi_1DCN_dt2e3 = psi;
    psi_true_c_1DCN_dt2e3 = psi_true_c;
    save 1DCN_dt2e3.mat V_1DCN_dt2e3 X_1DCN_dt2e3 t_1DCN_dt2e3 psi_1DCN_dt2e3 psi_true_c_1DCN_dt2e3
    clear
end

%% Figures for the Report: 1D TDSE: Euler method: ERROR

load('1DEULER_dt8e5.mat')
load('1DEULER_dt2e4.mat')
load('1DEULER_dt8e4.mat')
load('1DEULER_dt2e3.mat')

error_1DEULER_dt8e5 = abs(trapz(psi_1DEULER_dt8e5-psi_true_c_1DEULER_dt8e5,1));
error_1DEULER_dt2e4 = abs(trapz(psi_1DEULER_dt2e4-psi_true_c_1DEULER_dt2e4,1));
error_1DEULER_dt8e4 = abs(trapz(psi_1DEULER_dt8e4-psi_true_c_1DEULER_dt8e4,1));
error_1DEULER_dt2e3 = abs(trapz(psi_1DEULER_dt2e3-psi_true_c_1DEULER_dt2e3,1));

set(0, 'DefaultLineLineWidth', 2); close all;
tl = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
t1 = nexttile(1); hold on; grid on

plot(t_1DEULER_dt8e5,error_1DEULER_dt8e5)
plot(t_1DEULER_dt2e4,error_1DEULER_dt2e4)
plot(t_1DEULER_dt8e4,error_1DEULER_dt8e4)
plot(t_1DEULER_dt2e3,error_1DEULER_dt2e3)

legend("$\Delta t = 8\times10^{-5}$","$\Delta t = 2\times10^{-4}$",...
    "$\Delta t = 8\times10^{-4}$","$\Delta t = 2\times10^{-3}$")

t2 = nexttile(2); hold on; grid on 

plot(t_1DEULER_dt8e5,error_1DEULER_dt8e5)
plot(t_1DEULER_dt2e4,error_1DEULER_dt2e4)
plot(t_1DEULER_dt8e4,error_1DEULER_dt8e4)
plot(t_1DEULER_dt2e3,error_1DEULER_dt2e3)

xlim([0 0.4])
ylim([0 0.2*1e-2])

title(t1,"(a) Error for explicit Euler method in 1D")
title(t2,"(b) Close-up to initial behavior of error")
title(tl,"POTENTIAL: Free-particle", 'Interpreter', 'latex')

xlabel(tl,'$\textrm{Time step}$','Interpreter','latex');
ylabel(tl,'$\textrm{Error}= |\int \psi - \psi_{true} \; dx|$','Interpreter','latex');

picturewidth = 3*34.4/3; hw_ratio = 0.2; fs = 13; 
set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
pos = get(gcf, 'Position');
set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'northeast')
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

%% Figures for the Report: 1D TDSE: Crank-Nicolson method: ERROR

load('1DCN_dt8e5.mat')
load('1DCN_dt2e4.mat')
load('1DCN_dt8e4.mat')
load('1DCN_dt2e3.mat')

error_1DCN_dt8e5 = abs(trapz(psi_1DCN_dt8e5-psi_true_c_1DCN_dt8e5,1));
error_1DCN_dt2e4 = abs(trapz(psi_1DCN_dt2e4-psi_true_c_1DCN_dt2e4,1));
error_1DCN_dt8e4 = abs(trapz(psi_1DCN_dt8e4-psi_true_c_1DCN_dt8e4,1));
error_1DCN_dt2e3 = abs(trapz(psi_1DCN_dt2e3-psi_true_c_1DCN_dt2e3,1));

set(0, 'DefaultLineLineWidth', 2); close all;
tl = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
t1 = nexttile(1); hold on; grid on

plot(t_1DCN_dt8e5,error_1DCN_dt8e5)
plot(t_1DCN_dt2e4,error_1DCN_dt2e4)
plot(t_1DCN_dt8e4,error_1DCN_dt8e4)
plot(t_1DCN_dt2e3,error_1DCN_dt2e3)

legend("$\Delta t = 8\times10^{-5}$","$\Delta t = 2\times10^{-4}$",...
    "$\Delta t = 8\times10^{-4}$","$\Delta t = 2\times10^{-3}$")

t2 = nexttile(2); hold on; grid on 

plot(t_1DCN_dt8e5,error_1DCN_dt8e5)
plot(t_1DCN_dt2e4,error_1DCN_dt2e4)
plot(t_1DCN_dt8e4,error_1DCN_dt8e4)
plot(t_1DCN_dt2e3,error_1DCN_dt2e3)

xlim([0 0.4])
ylim([0 0.2*1e-2])

title(t1,"(a) Error for explicit Crank-Nicolson method in 1D")
title(t2,"(b) Close-up to initial behavior of error")
title(tl,"POTENTIAL: Free-particle", 'Interpreter', 'latex')

xlabel(tl,'$\textrm{Time step}$','Interpreter','latex');
ylabel(tl,'$\textrm{Error}= |\int \psi - \psi_{true} \; dx|$','Interpreter','latex');

picturewidth = 3*34.4/3; hw_ratio = 0.2; fs = 13; 
set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
pos = get(gcf, 'Position');
set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'northwest')
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

%% Figures for the Report: 1D TDSE: Euler method: SNAPSHOTS

load('1DEULER_dt8e5.mat')
load('1DEULER_dt2e4.mat')
load('1DEULER_dt8e4.mat')
load('1DEULER_dt2e3.mat')

set(0, 'DefaultLineLineWidth', 2); lw = 2;

close all; f = figure;
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_1DEULER_dt8e5; V = V_1DEULER_dt8e5;
V(1,1) = 1; V(1,end) = 1; % for visualization purposes

t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
potplot1 = plot(X, V,'-r','LineWidth',lw);
psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e5(:,1)),'-k','LineWidth',lw);
realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e5(:,1)),'-b','LineWidth',lw*0.5);
imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e5(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
potplot2 = plot(X, V,'-r','LineWidth',lw);
psiplot2 = plot(X,abs(psi_1DEULER_dt2e4(:,1)),'-k','LineWidth',lw);
realplot2 = plot(X,real(psi_1DEULER_dt2e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot2 = plot(X,imag(psi_1DEULER_dt2e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
potplot3 = plot(X, V,'-r','LineWidth',lw);
psiplot3 = plot(X,abs(psi_1DEULER_dt8e4(:,1)),'-k','LineWidth',lw);
realplot3 = plot(X,real(psi_1DEULER_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot3 = plot(X,imag(psi_1DEULER_dt8e4(:,1)),'-g','LineWidth',lw*0.5);

ylim([-0.5 0.5])
legend("V(x)","|$\psi$|","Re($\psi$)","Im($\psi$)")
picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4; fs = 13; 

no_snap = 5;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(realplot1)
    delete(imagplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(realplot2)
    delete(imagplot2)
    delete(potplot3)
    delete(psiplot3)
    delete(realplot3)
    delete(imagplot3)

    ts_1 = (length(psi_true_c_1DEULER_dt8e5)-1) / (no_snap-1);
    ts_2 = (length(psi_1DEULER_dt2e4)-1) / (no_snap-1);
    ts_3 = (length(psi_1DEULER_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
    potplot1 = plot(X, V,'-r','LineWidth',lw);
    psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e5(:,(snap-1)*ts_1+1)),'-k','LineWidth',lw);
    realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e5(:,(snap-1)*ts_1+1)),'-b','LineWidth',lw*0.5);
    imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e5(:,(snap-1)*ts_1+1)),'-g','LineWidth',lw*0.5);
    
    t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
    potplot2 = plot(X, V,'-r','LineWidth',lw);
    psiplot2 = plot(X,abs(psi_1DEULER_dt2e4(:,(snap-1)*ts_2+1)),'-k','LineWidth',lw);
    realplot2 = plot(X,real(psi_1DEULER_dt2e4(:,(snap-1)*ts_2+1)),'-b','LineWidth',lw*0.5);
    imagplot2 = plot(X,imag(psi_1DEULER_dt2e4(:,(snap-1)*ts_2+1)),'-g','LineWidth',lw*0.5);
    
    t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
    potplot3 = plot(X, V,'-r','LineWidth',lw);
    psiplot3 = plot(X,abs(psi_1DEULER_dt8e4(:,(snap-1)*ts_3+1)),'-k','LineWidth',lw);
    realplot3 = plot(X,real(psi_1DEULER_dt8e4(:,(snap-1)*ts_3+1)),'-b','LineWidth',lw*0.5);
    imagplot3 = plot(X,imag(psi_1DEULER_dt8e4(:,(snap-1)*ts_3+1)),'-g','LineWidth',lw*0.5);

    xlabel('x (nm)');
    tx = ['Time step = ', num2str(4/(no_snap-1)*snap-1)];
    
    legend("V(x)","$|\psi|$","Re($\psi$)","Im($\psi$)",'NumColumns',2)
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "1DEULER_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% Figures for the Report: 1D TDSE: Crank-Nicolson method: SNAPSHOTS

load('1DCN_dt8e5.mat')
load('1DCN_dt2e4.mat')
load('1DCN_dt8e4.mat')
load('1DCN_dt2e3.mat')

set(0, 'DefaultLineLineWidth', 2); lw = 2;

close all; f = figure;
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_1DCN_dt8e5; V = V_1DCN_dt8e5;
V(1,1) = 1; V(1,end) = 1; % for visualization purposes

t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
potplot1 = plot(X, V,'-r','LineWidth',lw);
psiplot1 = plot(X,abs(psi_true_c_1DCN_dt8e5(:,1)),'-k','LineWidth',lw);
realplot1 = plot(X,real(psi_true_c_1DCN_dt8e5(:,1)),'-b','LineWidth',lw*0.5);
imagplot1 = plot(X,imag(psi_true_c_1DCN_dt8e5(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
potplot2 = plot(X, V,'-r','LineWidth',lw);
psiplot2 = plot(X,abs(psi_1DCN_dt2e4(:,1)),'-k','LineWidth',lw);
realplot2 = plot(X,real(psi_1DCN_dt2e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot2 = plot(X,imag(psi_1DCN_dt2e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
potplot3 = plot(X, V,'-r','LineWidth',lw);
psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,1)),'-k','LineWidth',lw);
realplot3 = plot(X,real(psi_1DCN_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,1)),'-g','LineWidth',lw*0.5);

ylim([-0.5 0.5])
legend("V(x)","|$\psi$|","Re($\psi$)","Im($\psi$)")
picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4; fs = 13; 

no_snap = 5;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(realplot1)
    delete(imagplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(realplot2)
    delete(imagplot2)
    delete(potplot3)
    delete(psiplot3)
    delete(realplot3)
    delete(imagplot3)

    ts_1 = (length(psi_true_c_1DCN_dt8e5)-1) / (no_snap-1);
    ts_2 = (length(psi_1DCN_dt2e4)-1) / (no_snap-1);
    ts_3 = (length(psi_1DCN_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
    potplot1 = plot(X, V,'-r','LineWidth',lw);
    psiplot1 = plot(X,abs(psi_true_c_1DCN_dt8e5(:,(snap-1)*ts_1+1)),'-k','LineWidth',lw);
    realplot1 = plot(X,real(psi_true_c_1DCN_dt8e5(:,(snap-1)*ts_1+1)),'-b','LineWidth',lw*0.5);
    imagplot1 = plot(X,imag(psi_true_c_1DCN_dt8e5(:,(snap-1)*ts_1+1)),'-g','LineWidth',lw*0.5);
    
    t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
    potplot2 = plot(X, V,'-r','LineWidth',lw);
    psiplot2 = plot(X,abs(psi_1DCN_dt2e4(:,(snap-1)*ts_2+1)),'-k','LineWidth',lw);
    realplot2 = plot(X,real(psi_1DCN_dt2e4(:,(snap-1)*ts_2+1)),'-b','LineWidth',lw*0.5);
    imagplot2 = plot(X,imag(psi_1DCN_dt2e4(:,(snap-1)*ts_2+1)),'-g','LineWidth',lw*0.5);
    
    t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
    potplot3 = plot(X, V,'-r','LineWidth',lw);
    psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-k','LineWidth',lw);
    realplot3 = plot(X,real(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-b','LineWidth',lw*0.5);
    imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-g','LineWidth',lw*0.5);

    xlabel('x (nm)');
    tx = ['Time step = ', num2str(4/(no_snap-1)*snap-1)];
    
    legend("V(x)","$|\psi|$","Re($\psi$)","Im($\psi$)",'NumColumns',2)
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "1DCN_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% Figures for the Report: 1D TDSE: Tunneling: SNAPSHOTS

load('1DEULER_dt8e4.mat')
load('1DCN_dt8e4.mat')

set(0, 'DefaultLineLineWidth', 2); lw = 2;

close all; f = figure;
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_1DEULER_dt8e4; V = V_1DEULER_dt8e4;
V(1,1) = 1; V(1,end) = 1; % for visualization purposes

t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
potplot1 = plot(X, V,'-r','LineWidth',lw);
psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e4(:,1)),'-k','LineWidth',lw);
realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
potplot2 = plot(X, V,'-r','LineWidth',lw);
psiplot2 = plot(X,abs(psi_1DEULER_dt8e4(:,1)),'-k','LineWidth',lw);
realplot2 = plot(X,real(psi_1DEULER_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot2 = plot(X,imag(psi_1DEULER_dt8e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
potplot3 = plot(X, V,'-r','LineWidth',lw);
psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,1)),'-k','LineWidth',lw);
realplot3 = plot(X,real(psi_1DCN_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,1)),'-g','LineWidth',lw*0.5);

ylim([-0.5 0.5])
legend("V(x)","|$\psi$|","Re($\psi$)","Im($\psi$)")
picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4; fs = 13; 

no_snap = 5;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(realplot1)
    delete(imagplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(realplot2)
    delete(imagplot2)
    delete(potplot3)
    delete(psiplot3)
    delete(realplot3)
    delete(imagplot3)

    ts_1 = (length(psi_true_c_1DEULER_dt8e4)-1) / (no_snap-1);
    ts_2 = (length(psi_1DEULER_dt8e4)-1) / (no_snap-1);
    ts_3 = (length(psi_1DCN_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
    potplot1 = plot(X, V,'-r','LineWidth',lw);
    psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-k','LineWidth',lw);
    realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-b','LineWidth',lw*0.5);
    imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-g','LineWidth',lw*0.5);
    
    t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
    potplot2 = plot(X, V,'-r','LineWidth',lw);
    psiplot2 = plot(X,abs(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-k','LineWidth',lw);
    realplot2 = plot(X,real(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-b','LineWidth',lw*0.5);
    imagplot2 = plot(X,imag(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-g','LineWidth',lw*0.5);
    
    t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
    potplot3 = plot(X, V,'-r','LineWidth',lw);
    psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-k','LineWidth',lw);
    realplot3 = plot(X,real(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-b','LineWidth',lw*0.5);
    imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-g','LineWidth',lw*0.5);

    xlabel('x (nm)');
    tx = ['Time step = ', num2str(4/(no_snap-1)*snap-1)];
    
    legend("V(x)","$|\psi|$","Re($\psi$)","Im($\psi$)",'NumColumns',2)
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "1Dcomparison_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% GIF: 1D

load('1DEULER_dt8e4.mat')
load('1DCN_dt8e4.mat')

set(0, 'DefaultLineLineWidth', 2); lw = 2;

close all; f = figure; filename = "1d_tunneling.gif"; set(gcf,'color','w');
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_1DEULER_dt8e4; V = V_1DEULER_dt8e4;
V(1,1) = 1; V(1,end) = 1; % for visualization purposes

t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
potplot1 = plot(X, V,'-r','LineWidth',lw);
psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e4(:,1)),'-k','LineWidth',lw);
realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
potplot2 = plot(X, V,'-r','LineWidth',lw);
psiplot2 = plot(X,abs(psi_1DEULER_dt8e4(:,1)),'-k','LineWidth',lw);
realplot2 = plot(X,real(psi_1DEULER_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot2 = plot(X,imag(psi_1DEULER_dt8e4(:,1)),'-g','LineWidth',lw*0.5);
ylim([-0.5 0.5])

t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
potplot3 = plot(X, V,'-r','LineWidth',lw);
psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,1)),'-k','LineWidth',lw);
realplot3 = plot(X,real(psi_1DCN_dt8e4(:,1)),'-b','LineWidth',lw*0.5);
imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,1)),'-g','LineWidth',lw*0.5);

ylim([-0.5 0.5])
legend("V(x)","|$\psi$|","Re($\psi$)","Im($\psi$)")
picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4; fs = 13; 

no_snap = 101;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(realplot1)
    delete(imagplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(realplot2)
    delete(imagplot2)
    delete(potplot3)
    delete(psiplot3)
    delete(realplot3)
    delete(imagplot3)

    ts_1 = (length(psi_true_c_1DEULER_dt8e4)-1) / (no_snap-1);
    ts_2 = (length(psi_1DEULER_dt8e4)-1) / (no_snap-1);
    ts_3 = (length(psi_1DCN_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; grid on; set(gcf,'color','w');
    potplot1 = plot(X, V,'-r','LineWidth',lw);
    psiplot1 = plot(X,abs(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-k','LineWidth',lw);
    realplot1 = plot(X,real(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-b','LineWidth',lw*0.5);
    imagplot1 = plot(X,imag(psi_true_c_1DEULER_dt8e4(:,(snap-1)*ts_1+1)),'-g','LineWidth',lw*0.5);
    title("True solution, (Euler, $\Delta t = 2 \times 10^{-5}$)")

    t2 = nexttile(2); hold on; grid on; set(gcf,'color','w');
    potplot2 = plot(X, V,'-r','LineWidth',lw);
    psiplot2 = plot(X,abs(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-k','LineWidth',lw);
    realplot2 = plot(X,real(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-b','LineWidth',lw*0.5);
    imagplot2 = plot(X,imag(psi_1DEULER_dt8e4(:,(snap-1)*ts_2+1)),'-g','LineWidth',lw*0.5);
    title("Euler method, $\Delta t = 8\times 10^{-4}$")

    t3 = nexttile(3); hold on; grid on; set(gcf,'color','w');
    potplot3 = plot(X, V,'-r','LineWidth',lw);
    psiplot3 = plot(X,abs(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-k','LineWidth',lw);
    realplot3 = plot(X,real(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-b','LineWidth',lw*0.5);
    imagplot3 = plot(X,imag(psi_1DCN_dt8e4(:,(snap-1)*ts_3+1)),'-g','LineWidth',lw*0.5);
    title("Crank-Nicolson method, $\Delta t = 8\times 10^{-4}$")

    xlabel('x (nm)');
    tx = ['Time step = ', num2str(8e-4*(snap-1)*ts_2)];
    
    legend("V(x)","$|\psi|$","Re($\psi$)","Im($\psi$)",'NumColumns',2)
    title(tl, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    drawnow

   if snap == 1
       frame = getframe(f);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf);
   else
       frame = getframe(f);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append');
   end

end

%% Save data: Dimension = 2, Method = Explicit Euler
% Run after each simulation

if dt == 8e-5
    X_2DEULER_dt8e5 = X;
    Y_2DEULER_dt8e5 = Y;
    V_2DEULER_dt8e5 = V;
    t_2DEULER_dt8e5 = t;
    psi_2DEULER_dt8e5 = psi;
    psi_true_c_2DEULER_dt8e5 = psi_true_c;
    save 2DEULER_dt8e5.mat V_2DEULER_dt8e5 X_2DEULER_dt8e5 Y_2DEULER_dt8e5 t_2DEULER_dt8e5 ...
        psi_2DEULER_dt8e5 psi_true_c_2DEULER_dt8e5
    clear
elseif dt == 2e-4
    X_2DEULER_dt2e4 = X;
    Y_2DEULER_dt2e4 = Y;
    V_2DEULER_dt2e4 = V;
    t_2DEULER_dt2e4 = t;
    psi_2DEULER_dt2e4 = psi;
    psi_true_c_2DEULER_dt2e4 = psi_true_c;
    save 2DEULER_dt2e4.mat V_2DEULER_dt2e4 X_2DEULER_dt2e4 Y_2DEULER_dt2e4 t_2DEULER_dt2e4 ...
        psi_2DEULER_dt2e4 psi_true_c_2DEULER_dt2e4
    clear
elseif dt == 8e-4
    X_2DEULER_dt8e4 = X;
    Y_2DEULER_dt8e4 = Y;
    V_2DEULER_dt8e4 = V;
    t_2DEULER_dt8e4 = t;
    psi_2DEULER_dt8e4 = psi;
    psi_true_c_2DEULER_dt8e4 = psi_true_c;
    save 2DEULER_dt8e4.mat V_2DEULER_dt8e4 X_2DEULER_dt8e4 Y_2DEULER_dt8e4 t_2DEULER_dt8e4 ...
        psi_2DEULER_dt8e4 psi_true_c_2DEULER_dt8e4
    clear
elseif dt == 2e-3
    X_2DEULER_dt2e3 = X;
    Y_2DEULER_dt2e3 = Y;
    V_2DEULER_dt2e3 = V;
    t_2DEULER_dt2e3 = t;
    psi_2DEULER_dt2e3 = psi;
    psi_true_c_2DEULER_dt2e3 = psi_true_c;
    save 2DEULER_dt2e3.mat V_2DEULER_dt2e3 X_2DEULER_dt2e3 Y_2DEULER_dt2e3 t_2DEULER_dt2e3 ...
        psi_2DEULER_dt2e3 psi_true_c_2DEULER_dt2e3
    clear
end

%% Save data: Dimension = 2, Method = Crank-Nicolson
% Run after each simulation

if dt == 8e-5
    X_2DCN_dt8e5 = X;
    Y_2DCN_dt8e5 = Y;
    V_2DCN_dt8e5 = V;
    t_2DCN_dt8e5 = t;
    psi_2DCN_dt8e5 = psi;
    psi_true_c_2DCN_dt8e5 = psi_true_c;
    save 2DCN_dt8e5.mat V_2DCN_dt8e5 X_2DCN_dt8e5 Y_2DCN_dt8e5 t_2DCN_dt8e5 psi_2DCN_dt8e5 psi_true_c_2DCN_dt8e5
    clear
elseif dt == 2e-4
    X_2DCN_dt2e4 = X;
    Y_2DCN_dt2e4 = Y;
    V_2DCN_dt2e4 = V;
    t_2DCN_dt2e4 = t;
    psi_2DCN_dt2e4 = psi;
    psi_true_c_2DCN_dt2e4 = psi_true_c;
    save 2DCN_dt2e4.mat V_2DCN_dt2e4 X_2DCN_dt2e4 Y_2DCN_dt2e4 t_2DCN_dt2e4 psi_2DCN_dt2e4 psi_true_c_2DCN_dt2e4
    clear
elseif dt == 8e-4
    X_2DCN_dt8e4 = X;
    Y_2DCN_dt8e4 = Y;
    V_2DCN_dt8e4 = V;
    t_2DCN_dt8e4 = t;
    psi_2DCN_dt8e4 = psi;
    psi_true_c_2DCN_dt8e4 = psi_true_c;
    save 2DCN_dt8e4.mat V_2DCN_dt8e4 X_2DCN_dt8e4 Y_2DCN_dt8e4 t_2DCN_dt8e4 psi_2DCN_dt8e4 psi_true_c_2DCN_dt8e4
    clear
elseif dt == 2e-3
    X_2DCN_dt2e3 = X;
    Y_2DCN_dt2e3 = Y;
    V_2DCN_dt2e3 = V;
    t_2DCN_dt2e3 = t;
    psi_2DCN_dt2e3 = psi;
    psi_true_c_2DCN_dt2e3 = psi_true_c;
    save 2DCN_dt2e3.mat V_2DCN_dt2e3 X_2DCN_dt2e3 Y_2DCN_dt2e3 t_2DCN_dt2e3 psi_2DCN_dt2e3 psi_true_c_2DCN_dt2e3
    clear
end

%% Figures for the Report: 2D TDSE: Euler method: ERROR

load('2DEULER_dt8e5.mat')
load('2DEULER_dt2e4.mat')
load('2DEULER_dt8e4.mat')
load('2DEULER_dt2e3.mat')

error_2DEULER_dt8e5 = abs(trapz(trapz(psi_2DEULER_dt8e5-psi_true_c_2DEULER_dt8e5,1)));
error_2DEULER_dt2e4 = abs(trapz(trapz(psi_2DEULER_dt2e4-psi_true_c_2DEULER_dt2e4,1)));
error_2DEULER_dt8e4 = abs(trapz(trapz(psi_2DEULER_dt8e4-psi_true_c_2DEULER_dt8e4,1)));
error_2DEULER_dt2e3 = abs(trapz(trapz(psi_2DEULER_dt2e3-psi_true_c_2DEULER_dt2e3,1)));

error_2DEULER_dt8e5 = reshape(error_2DEULER_dt8e5,[length(error_2DEULER_dt8e5) 1]);
error_2DEULER_dt2e4 = reshape(error_2DEULER_dt2e4,[length(error_2DEULER_dt2e4) 1]);
error_2DEULER_dt8e4 = reshape(error_2DEULER_dt8e4,[length(error_2DEULER_dt8e4) 1]);
error_2DEULER_dt2e3 = reshape(error_2DEULER_dt2e3,[length(error_2DEULER_dt2e3) 1]);

set(0, 'DefaultLineLineWidth', 2); close all;
tl = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
t1 = nexttile(1); hold on; grid on

plot(t_2DEULER_dt8e5,error_2DEULER_dt8e5)
plot(t_2DEULER_dt2e4,error_2DEULER_dt2e4)
plot(t_2DEULER_dt8e4,error_2DEULER_dt8e4)
plot(t_2DEULER_dt2e3,error_2DEULER_dt2e3)

legend("$\Delta t = 8\times10^{-5}$","$\Delta t = 2\times10^{-4}$",...
    "$\Delta t = 8\times10^{-4}$","$\Delta t = 2\times10^{-3}$")

t2 = nexttile(2); hold on; grid on 

plot(t_2DEULER_dt8e5,error_2DEULER_dt8e5)
plot(t_2DEULER_dt2e4,error_2DEULER_dt2e4)
plot(t_2DEULER_dt8e4,error_2DEULER_dt8e4)
plot(t_2DEULER_dt2e3,error_2DEULER_dt2e3)

xlim([0 0.1])
%ylim([0 0.2*1e-2])

title(t1,"(a) Error for explicit Euler method in 2D")
title(t2,"(b) Close-up to initial behavior of error")
title(tl,"POTENTIAL: Free-particle", 'Interpreter', 'latex')

xlabel(tl,'$\textrm{Time step}$','Interpreter','latex');
ylabel(tl,'$\textrm{Error}= |\int \int \psi - \psi_{true} \; dx dy|$','Interpreter','latex');

picturewidth = 3*34.4/3; hw_ratio = 0.2; fs = 13; 
set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
pos = get(gcf, 'Position');
set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'northeast')
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

%% Figures for the Report: 2D TDSE: Crank-Nicolson method: ERROR

load('2DCN_dt8e5.mat')
load('2DCN_dt2e4.mat')
load('2DCN_dt8e4.mat')
load('2DCN_dt2e3.mat')

error_2DCN_dt8e5 = abs(trapz(trapz(psi_2DCN_dt8e5-psi_true_c_2DCN_dt8e5,1)));
error_2DCN_dt2e4 = abs(trapz(trapz(psi_2DCN_dt2e4-psi_true_c_2DCN_dt2e4,1)));
error_2DCN_dt8e4 = abs(trapz(trapz(psi_2DCN_dt8e4-psi_true_c_2DCN_dt8e4,1)));
error_2DCN_dt2e3 = abs(trapz(trapz(psi_2DCN_dt2e3-psi_true_c_2DCN_dt2e3,1)));

error_2DCN_dt8e5 = reshape(error_2DCN_dt8e5,[length(error_2DCN_dt8e5) 1]);
error_2DCN_dt2e4 = reshape(error_2DCN_dt2e4,[length(error_2DCN_dt2e4) 1]);
error_2DCN_dt8e4 = reshape(error_2DCN_dt8e4,[length(error_2DCN_dt8e4) 1]);
error_2DCN_dt2e3 = reshape(error_2DCN_dt2e3,[length(error_2DCN_dt2e3) 1]);

set(0, 'DefaultLineLineWidth', 2); close all;
tl = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
t1 = nexttile(1); hold on; grid on

plot(t_2DCN_dt8e5,error_2DCN_dt8e5)
plot(t_2DCN_dt2e4,error_2DCN_dt2e4)
plot(t_2DCN_dt8e4,error_2DCN_dt8e4)
plot(t_2DCN_dt2e3,error_2DCN_dt2e3)

legend("$\Delta t = 8\times10^{-5}$","$\Delta t = 2\times10^{-4}$",...
    "$\Delta t = 8\times10^{-4}$","$\Delta t = 2\times10^{-3}$")

t2 = nexttile(2); hold on; grid on 

plot(t_2DCN_dt8e5,error_2DCN_dt8e5)
plot(t_2DCN_dt2e4,error_2DCN_dt2e4)
plot(t_2DCN_dt8e4,error_2DCN_dt8e4)
plot(t_2DCN_dt2e3,error_2DCN_dt2e3)

xlim([0 0.1])
%ylim([0 0.2*1e-2])

title(t1,"(a) Error for Crank-Nicolson method in 2D")
title(t2,"(b) Close-up to initial behavior of error")
title(tl,"POTENTIAL: Free-particle", 'Interpreter', 'latex')

xlabel(tl,'$\textrm{Time step}$','Interpreter','latex');
ylabel(tl,'$\textrm{Error}= |\int \int \psi - \psi_{true} \; dx dy|$','Interpreter','latex');

picturewidth = 3*34.4/3; hw_ratio = 0.2; fs = 13; 
set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
pos = get(gcf, 'Position');
set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'northwest')
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

%% Figures for the Report: 2D TDSE: Euler method: SNAPSHOTS

load('2DEULER_dt8e5.mat')
load('2DEULER_dt2e4.mat')
load('2DEULER_dt8e4.mat')
load('2DEULER_dt2e3.mat')

x = X_2DEULER_dt2e3(1,:); y = Y_2DEULER_dt2e3(:,1);

set(0, 'DefaultLineLineWidth', 2); lw = 2; fs = 13;

close all; f = figure;
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_2DEULER_dt8e5; Y = Y_2DEULER_dt8e5; V = V_2DEULER_dt8e5;
V(:,1) = 1; V(:,end) = 1; V(1,:) = 1; V(end,:) = 1; % for visualization purposes

t1 = nexttile(1); hold on; set(gcf,'color','w');
potplot1 = imagesc(x,y,V);
psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt2e3(:,:,1)));

t2 = nexttile(2); hold on; set(gcf,'color','w');
potplot2 = imagesc(x,y,V);
psiplot2 = imagesc(x,y,abs(psi_2DEULER_dt2e4(:,:,1)));

t3 = nexttile(3); hold on; set(gcf,'color','w');
potplot3 = imagesc(x,y,V);
psiplot3 = imagesc(x,y,abs(psi_2DEULER_dt8e4(:,:,1)));
c = colorbar; clim([0 25]); c.Location = "southoutside";

picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4;

no_snap = 6;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(potplot3)
    delete(psiplot3)

    ts_1 = (length(psi_true_c_2DEULER_dt2e3)-1) / (no_snap-1);
    ts_2 = (length(psi_2DEULER_dt2e4)-1) / (no_snap-1);
    ts_3 = (length(psi_2DEULER_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; set(gcf,'color','w');
	potplot1 = imagesc(x,y,V);
	psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt2e3(:,:,(snap-1)*ts_1+1)));

    t2 = nexttile(2); hold on; set(gcf,'color','w');
	potplot2 = imagesc(x,y,V);
	psiplot2 = imagesc(x,y,abs(psi_2DEULER_dt2e4(:,:,(snap-1)*ts_2+1)));
    ylabel('y (nm)');

    t3 = nexttile(3); hold on; set(gcf,'color','w');
	potplot3 = imagesc(x,y,V);
	psiplot3 = imagesc(x,y,abs(psi_2DEULER_dt8e4(:,:,(snap-1)*ts_3+1)));
    %c = colorbar; clim([0 25]); c.Location = "southoutside";

    xlabel('x (nm)'); 
    tx = ['Time step = ', num2str(1/(no_snap-1)*(snap-1))];
    
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "2DEULER_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% Figures for the Report: 2D TDSE: Crank-Nicolson method: SNAPSHOTS

load('2DEULER_dt2e3.mat')
load('2DCN_dt8e5.mat')
load('2DCN_dt2e4.mat')
load('2DCN_dt8e4.mat')
load('2DCN_dt2e3.mat')

x = X_2DCN_dt2e3(1,:); y = Y_2DCN_dt2e3(:,1);

set(0, 'DefaultLineLineWidth', 2); lw = 2; fs = 13;

close all; f = figure;
tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

X = X_2DCN_dt8e5; Y = Y_2DCN_dt8e5; V = V_2DCN_dt8e5;
V(:,1) = 1; V(:,end) = 1; V(1,:) = 1; V(end,:) = 1; % for visualization purposes

t1 = nexttile(1); hold on; set(gcf,'color','w');
potplot1 = imagesc(x,y,V);
psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt2e3(:,:,1)));

t2 = nexttile(2); hold on; set(gcf,'color','w');
potplot2 = imagesc(x,y,V);
psiplot2 = imagesc(x,y,abs(psi_2DCN_dt2e4(:,:,1)));

t3 = nexttile(3); hold on; set(gcf,'color','w');
potplot3 = imagesc(x,y,V);
psiplot3 = imagesc(x,y,abs(psi_2DCN_dt8e4(:,:,1)));
c = colorbar; clim([0 25]); c.Location = "southoutside";

picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4;

no_snap = 6;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(potplot2)
    delete(psiplot2)
    delete(potplot3)
    delete(psiplot3)

    ts_1 = (length(psi_true_c_2DEULER_dt2e3)-1) / (no_snap-1);
    ts_2 = (length(psi_2DCN_dt2e4)-1) / (no_snap-1);
    ts_3 = (length(psi_2DCN_dt8e4)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; set(gcf,'color','w');
	potplot1 = imagesc(x,y,V);
	psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt2e3(:,:,(snap-1)*ts_1+1)));

    t2 = nexttile(2); hold on; set(gcf,'color','w');
	potplot2 = imagesc(x,y,V);
	psiplot2 = imagesc(x,y,abs(psi_2DCN_dt2e4(:,:,(snap-1)*ts_2+1)));
    ylabel('y (nm)');

    t3 = nexttile(3); hold on; set(gcf,'color','w');
	potplot3 = imagesc(x,y,V);
	psiplot3 = imagesc(x,y,abs(psi_2DCN_dt8e4(:,:,(snap-1)*ts_3+1)));
    %c = colorbar; clim([0 25]); c.Location = "southoutside";

    xlabel('x (nm)'); 
    tx = ['Time step = ', num2str(1/(no_snap-1)*(snap-1))];
    
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "2DCN_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% Figures for the Report: 2D TDSE: Double/single-slit: SNAPSHOTS

load('2DEULER_dt8e5_doubleslit.mat')

X_2DEULER_dt8e5_doubleslit = X_2DEULER_dt8e5;
Y_2DEULER_dt8e5_doubleslit = Y_2DEULER_dt8e5;
V_2DEULER_dt8e5_doubleslit = V_2DEULER_dt8e5;
psi_true_c_2DEULER_dt8e5_doubleslit = psi_true_c_2DEULER_dt8e5;
psi_2DEULER_dt8e5_doubleslit = psi_2DEULER_dt8e5;

clear X_2DEULER_dt8e5 Y_2DEULER_dt8e5 V_2DEULER_dt8e5 psi_true_c_2DEULER_dt8e5 psi_2DEULER_dt8e5

load('2DEULER_dt8e5_singleslit.mat')

X_2DEULER_dt8e5_singleslit = X_2DEULER_dt8e5;
Y_2DEULER_dt8e5_singleslit = Y_2DEULER_dt8e5;
V_2DEULER_dt8e5_singleslit = V_2DEULER_dt8e5;
psi_true_c_2DEULER_dt8e5_singleslit = psi_true_c_2DEULER_dt8e5;
psi_2DEULER_dt8e5_singleslit = psi_2DEULER_dt8e5;

clear X_2DEULER_dt8e5 Y_2DEULER_dt8e5 V_2DEULER_dt8e5 psi_true_c_2DEULER_dt8e5 psi_2DEULER_dt8e5

x = X_2DEULER_dt8e5_doubleslit(1,:); y = Y_2DEULER_dt8e5_doubleslit(:,1);

set(0, 'DefaultLineLineWidth', 2); lw = 2; fs = 13;

close all; f = figure;
tl = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

%V_2DEULER_dt8e5_doubleslit(:,1) = 1; V_2DEULER_dt8e5_doubleslit(:,end) = 1; 
%V_2DEULER_dt8e5_doubleslit(1,:) = 1; V_2DEULER_dt8e5_doubleslit(end,:) = 1; % for visualization purposes
%V_2DEULER_dt8e5_singleslit(:,1) = 1; V_2DEULER_dt8e5_singleslit(:,end) = 1; 
%V_2DEULER_dt8e5_singleslit(1,:) = 1; V_2DEULER_dt8e5_singleslit(end,:) = 1; % for visualization purposes

t1 = nexttile(1); hold on; set(gcf,'color','w');
potplot1 = imagesc(x,y,V_2DEULER_dt8e5_singleslit,'AlphaData', 0.5);
psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,1)),'AlphaData', 0.5);
c = colorbar; c.Location = "southoutside";
%clim([0 25]);

t2 = nexttile(2); hold on; set(gcf,'color','w');
potplot2 = imagesc(x,y,V_2DEULER_dt8e5_doubleslit,'AlphaData', 0.5);
psiplot2 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,1)),'AlphaData', 0.5);

c = colorbar; c.Location = "southoutside";
% clim([0 25]);

picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4;

no_snap = 6;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(potplot2)
    delete(psiplot2)

    ts_1 = (length(psi_2DEULER_dt8e5_singleslit)-1) / (no_snap-1);
    ts_2 = (length(psi_2DEULER_dt8e5_doubleslit)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; set(gcf,'color','w');
	potplot1 = imagesc(x,y,V_2DEULER_dt8e5_singleslit,'AlphaData', 0.5);
	psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,(snap-1)*ts_1+1)),'AlphaData', 0.5);
    caxis([0 max(max(abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,(snap-1)*ts_1+1))))])

    t2 = nexttile(2); hold on; set(gcf,'color','w');
	potplot2 = imagesc(x,y,V_2DEULER_dt8e5_doubleslit,'AlphaData', 0.5);
	psiplot2 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,(snap-1)*ts_2+1)),'AlphaData', 0.5);
    caxis([0 max(max(abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,(snap-1)*ts_1+1))))])

    ylabel(tl, 'y (nm)', 'Interpreter', 'latex');

    xlabel('x (nm)'); 
    tx = ['Time step = ', num2str(1/(no_snap-1)*(snap-1))];
    
    title(t1, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    ftx = "2DSLIT_" + snap;
    drawnow
    print(f,ftx,'-dpdf','-vector', '-fillpage')

end

%% GIF: 2D

load('2DEULER_dt8e5_doubleslit.mat')

X_2DEULER_dt8e5_doubleslit = X_2DEULER_dt8e5;
Y_2DEULER_dt8e5_doubleslit = Y_2DEULER_dt8e5;
V_2DEULER_dt8e5_doubleslit = V_2DEULER_dt8e5;
psi_true_c_2DEULER_dt8e5_doubleslit = psi_true_c_2DEULER_dt8e5;
psi_2DEULER_dt8e5_doubleslit = psi_2DEULER_dt8e5;

clear X_2DEULER_dt8e5 Y_2DEULER_dt8e5 V_2DEULER_dt8e5 psi_true_c_2DEULER_dt8e5 psi_2DEULER_dt8e5

load('2DEULER_dt8e5_singleslit.mat')

X_2DEULER_dt8e5_singleslit = X_2DEULER_dt8e5;
Y_2DEULER_dt8e5_singleslit = Y_2DEULER_dt8e5;
V_2DEULER_dt8e5_singleslit = V_2DEULER_dt8e5;
psi_true_c_2DEULER_dt8e5_singleslit = psi_true_c_2DEULER_dt8e5;
psi_2DEULER_dt8e5_singleslit = psi_2DEULER_dt8e5;

clear X_2DEULER_dt8e5 Y_2DEULER_dt8e5 V_2DEULER_dt8e5 psi_true_c_2DEULER_dt8e5 psi_2DEULER_dt8e5

x = X_2DEULER_dt8e5_doubleslit(1,:); y = Y_2DEULER_dt8e5_doubleslit(:,1);

set(0, 'DefaultLineLineWidth', 2); lw = 2; fs = 13;

close all; f = figure; filename = "2d_slit.gif"; set(gcf,'color','w');
tl = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

%V_2DEULER_dt8e5_doubleslit(:,1) = 1; V_2DEULER_dt8e5_doubleslit(:,end) = 1; 
%V_2DEULER_dt8e5_doubleslit(1,:) = 1; V_2DEULER_dt8e5_doubleslit(end,:) = 1; % for visualization purposes
%V_2DEULER_dt8e5_singleslit(:,1) = 1; V_2DEULER_dt8e5_singleslit(:,end) = 1; 
%V_2DEULER_dt8e5_singleslit(1,:) = 1; V_2DEULER_dt8e5_singleslit(end,:) = 1; % for visualization purposes

t1 = nexttile(1); hold on; set(gcf,'color','w');
potplot1 = imagesc(x,y,V_2DEULER_dt8e5_singleslit,'AlphaData', 0.5);
psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,1)),'AlphaData', 0.5);
c = colorbar; c.Location = "southoutside";
%clim([0 25]);

t2 = nexttile(2); hold on; set(gcf,'color','w');
potplot2 = imagesc(x,y,V_2DEULER_dt8e5_doubleslit,'AlphaData', 0.5);
psiplot2 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,1)),'AlphaData', 0.5);

c = colorbar; c.Location = "southoutside";
% clim([0 25]);

picturewidth = 3*34.4/3/4; hw_ratio = 0.2*3*4;

no_snap = 101;
for snap = 1:no_snap
    
    delete(potplot1)
    delete(psiplot1)
    delete(potplot2)
    delete(psiplot2)

    ts_1 = (length(psi_2DEULER_dt8e5_singleslit)-1) / (no_snap-1);
    ts_2 = (length(psi_2DEULER_dt8e5_doubleslit)-1) / (no_snap-1);

    t1 = nexttile(1); hold on; set(gcf,'color','w'); xlim([-5 5]); ylim([-5 5]);
	potplot1 = imagesc(x,y,V_2DEULER_dt8e5_singleslit,'AlphaData', 0.5);
	psiplot1 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,(snap-1)*ts_1+1)),'AlphaData', 0.5);
    caxis([0 max(max(abs(psi_true_c_2DEULER_dt8e5_singleslit(:,:,(snap-1)*ts_1+1))))])
    title("Single-slit, (Euler, $\Delta t = 8 \times 10^{-5}$)")

    t2 = nexttile(2); hold on; set(gcf,'color','w'); xlim([-5 5]); ylim([-5 5]);
	potplot2 = imagesc(x,y,V_2DEULER_dt8e5_doubleslit,'AlphaData', 0.5);
	psiplot2 = imagesc(x,y,abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,(snap-1)*ts_2+1)),'AlphaData', 0.5);
    caxis([0 max(max(abs(psi_true_c_2DEULER_dt8e5_doubleslit(:,:,(snap-1)*ts_1+1))))])
    title("Double-slit, (Euler, $\Delta t = 8 \times 10^{-5}$)")

    ylabel(tl, 'y (nm)', 'Interpreter', 'latex');

    xlabel('x (nm)'); 
    tx = ['Time step = ', num2str(8e-5*(snap-1)*ts_2)];
    
    title(tl, tx, 'Interpreter', 'latex')
    set(findall(gcf,'-property', 'FontSize'), 'FontSize', fs) % never change fontsize anymore!
    set(findall(gcf,'-property', 'Box'), 'Box', 'on') % optional
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(gcf, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(gcf, 'Units', 'centimeters', 'Position', [2 1 picturewidth hw_ratio*picturewidth])
    pos = get(gcf, 'Position');
    set(findobj(gcf, 'Type', 'legend'), 'FontSize', fs-3, 'Location', 'southoutside')
    set(gca,'LooseInset',get(gca,'TightInset'))
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
    drawnow

   if snap == 1
       frame = getframe(f);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf);
   else
       frame = getframe(f);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append');
   end

end

%% Functions

function f_next = ExplicitEuler1D(f, hbar, m, dx, dt, V)

    % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
    % v = (i hbar)/(2m) -> diffusion coefficient
    % f = psi(r,t)
    % df/dt = v (d^2f/dx^2) + V(x)/(i hbar) f
        
    f_next = 0*f;

    for x_ = 2:size(f,1)-1
        laplacian_x = 1/dx^2 * (f(x_-1) - 2*f(x_) + f(x_+1));
        f_next(x_) = f(x_) + (-hbar^2/(2*m))/(1i*hbar)*dt*(laplacian_x) + dt*V(x_)*f(x_)/(1i*hbar);
    end
end

function f_next = RungeKutta41D(f, hbar, m, dx, dt, V)

    % d(psi(r,t))/dt = (i hbar)/(2m) nabla^2 psi(r,t) + 1/(i hbar) V(r) psi(r,t)
    % v = (i hbar)/(2m) -> diffusion coefficient
    % f = psi(r,t)
    % df/dt = v (d^2f/dx^2) + V(x)/(i hbar) f
        
    f_next = 0*f;

    for x_ = 2:size(f,1)-1
        laplacian_x = 1/dx^2 * (f(x_-1) - 2*f(x_) + f(x_+1));
        f_next(x_) = f(x_) + (-hbar^2/(2*m))/(1i*hbar)*dt*(laplacian_x) + dt*V(x_)*f(x_)/(1i*hbar);
    end
end

function f_new = CrankNicolson1D(f, hbar, m, dx, dt, V)
    SoR = 1.8;
    
    N = length(f);
    f_temp = 0*f;
    f_new = f_temp;
    
    a1 = -(1i * hbar * dt)/(4*m*dx^2);
    b1 = (2*1i*hbar*dt)/(4*m*dx^2) + 1 - (V*dt)/(2*1i*hbar);
    c1 = a1;
    
    a2 = (1i * hbar * dt)/(4*m*dx^2);
    b2 = 1 - (2*1i*hbar*dt)/(4*m*dx^2) + (V*dt)/(2*1i*hbar);
    c2 = a2;
    
    A = zeros(N-2,N-2);
    B = zeros(N-2,N-2);

    for x_ = 1:N-2
        for y_ = 1:N-2
            if x_-y_ == -1
                A(x_,y_) = a1;
                B(x_,y_) = a2;
            elseif x_-y_ == 0
                A(x_,y_) = b1(x_);
                B(x_,y_) = b2(x_);
            elseif x_-y_ == +1
                A(x_,y_) = c1;
                B(x_,y_) = c2;
            end
        end
    end
    
    tol = 1e-6;
    max_iter = 1;
    
    for i = 1:max_iter
        RHS = B*f(2:N-1);
        f_new(2:N-1) = A\RHS;
        if norm(f_new - f, 2) < tol
            break
        end
    end
end

function f_new = CrankNicolson2D(f, hbar, m, dx, dy, dt, V)
    
    N = size(f,1);
    f_new = 0*f;

    a1_hs = -(1i * hbar * dt)/(4*m*dx^2);
    b1_hs = 1 + (2*1i*hbar*dt)/(4*m*dx^2) - (V*dt)/(4*1i*hbar);
    c1_hs = a1_hs;
    a2_hs = (1i * hbar * dt)/(4*m*dy^2);
    b2_hs = 1 - (2*1i*hbar*dt)/(4*m*dy^2) + (V*dt)/(4*1i*hbar);
    c2_hs = a2_hs; 

    a1_fs = -(1i * hbar * dt)/(4*m*dy^2);
    b1_fs = 1 + (2*1i*hbar*dt)/(4*m*dy^2) - (V*dt)/(4*1i*hbar);
    c1_fs = a1_fs;
    a2_fs = (1i * hbar * dt)/(4*m*dx^2);
    b2_fs = 1 - (2*1i*hbar*dt)/(4*m*dx^2) + (V*dt)/(4*1i*hbar);
    c2_fs = a2_fs; 
    
    A_hs = zeros(N-2,N-2);
    B_hs = zeros(N-2,N-2);
    A_fs = zeros(N-2,N-2);
    B_fs = zeros(N-2,N-2);
    
    for x_ = 1:N-2
        for y_ = 1:N-2
            if x_-y_ == -1
                A_hs(x_,y_) = A_hs(x_,y_) + a1_hs;
                B_hs(x_,y_) = B_hs(x_,y_) + a2_hs;
                
                A_fs(x_,y_) = A_fs(x_,y_) + c1_fs;
                B_fs(x_,y_) = B_fs(x_,y_) + c2_fs;
            elseif x_-y_ == 0
                A_hs(x_,y_) = A_hs(x_,y_) + b1_hs(x_,y_);
                B_hs(x_,y_) = B_hs(x_,y_) + b2_hs(x_,y_);
                
                A_fs(x_,y_) = A_fs(x_,y_) + b1_fs(x_,y_);
                B_fs(x_,y_) = B_fs(x_,y_) + b2_fs(x_,y_);
            elseif x_-y_ == +1
                A_hs(x_,y_) = A_hs(x_,y_) + c1_hs;
                B_hs(x_,y_) = B_hs(x_,y_) + c2_hs;
                
                A_fs(x_,y_) = A_fs(x_,y_) + a1_fs;
                B_fs(x_,y_) = B_fs(x_,y_) + a2_fs;
            end
        end
    end
    
    tol = 1e-6;
    max_iter = 1; %1000
    
    
    % Not yet iterative
    for i = 1:max_iter
        % First x-direction
        RHS_hs = B_hs*f(2:N-1,2:N-1)';
        % Then y-direction
        f_hs = (A_hs\RHS_hs);
        %f_hs = f(2:N-1,2:N-1);
        RHS_fs = B_fs*f_hs';
        f_new(2:N-1,2:N-1) = (A_fs\RHS_fs);
        f_new = f_new;
%         if norm(f_new - f, 2) < tol
%             break
%         end
        f = f_new;
    end
    
end

function f_next = ExplicitEuler2D(f, hbar, m, dx, dy, dt, V)
    f_next = 0*f;
    for x_ = 2:size(f,1)-1
        for y_ = 2:size(f,2)-1
            laplacian_x = 1/dx^2 * (f(x_-1,y_) - 2*f(x_,y_) + f(x_+1,y_));
            laplacian_y = 1/dy^2 * (f(x_,y_-1) - 2*f(x_,y_) + f(x_,y_+1));
            f_next(x_,y_) = f(x_,y_) + (-hbar^2/(2*m))/(1i*hbar)*dt*(laplacian_x + laplacian_y) + ...
                dt*V(x_,y_)*f(x_,y_)/(1i*hbar);
        end
    end
end

function psi = Normalizewavefunction1D(psi)
    A = sqrt(trapz(psi.^2));
    psi = psi/A;
end

function psi = Normalizewavefunction2D(psi)
    A = sqrt(trapz(trapz(psi.^2)));
    psi = psi/A;
end

function V = potential1D(type, x)
    switch type
        case "harmonic"
            % Define the potential energy function
            omega = 0.1; % frequency of the harmonic oscillator (in eV)
            m = 1;
            X = reshape(x,[size(x,2) 1]); % 1D grid
            V = m*omega^2*(X.^2)/2; % harmonic oscillator potential
        case "finite-well"
            a = 8;
            V = zeros(size(x));
            V(abs(x) <= a/2) = 0;
            V(abs(x) > a/2) = 1; %inf;
        case "tunneling"
            a = 1;
            width = 0.5;
            V = zeros(size(x));
            V(x >= a) = 0.5;
            V(x > a+width) = 0;
        case "free-particle"
            V = zeros(size(x));
    end
end

function V = potential2D(type, x, y)
    switch type
        case "harmonic"
            % Define the potential energy function
            omega = 0.1; % frequency of the harmonic oscillator (in eV)
            m = 1;
            [X,Y] = meshgrid(x,y); % 2D grid
            V = m*omega^2*(X.^2+Y.^2)/2; % harmonic oscillator potential
        case "finite-well"
            a = 8;
            V = zeros(length(x),length(y));
            V(abs(x) <= a/2, abs(y) <= a/2) = 0;
            V(abs(x) > a/2, abs(y) > a/2) = 1; %inf;
        case "tunneling"
            a = 1;
            width_x = 0.5;
            V = zeros(length(x), length(y));
            V(x >= a, y >= a) = 5;
            V(x > a+width_x, y > a+width_x) = 0;
        case "single-slit"
            a = 1;
            dy = y(2) - y(1);
            width_x = 0.5;
            width_y = 1;
            V = zeros(length(x), length(y));
            V(x >= a, :) = 100;
            V(x >= a, round((end+1)/2-width_y/dy/2):round((end+1)/2+width_y/dy/2)) = 0;
            V(x > a+width_x, :) = 0;
            V = transpose(V);
        case "double-slit"
            a = 1;
            %dx = x(2) - x(1);
            dy = y(2) - y(1);
            width_x = 0.5;
            width_y = 1;
            V = zeros(length(x), length(y));
            V(x >= a, :) = 100;
            V(x >= a, round((end+1)/2-width_y/dy):round((end+1)/2-width_y/dy/2)) = 0;
            V(x >= a, round((end+1)/2+width_y/dy/2):round((end+1)/2+width_y/dy)) = 0;
            V(x > a+width_x, :) = 0;
            V = transpose(V);
        case "free-particle"
            V = zeros(length(x), length(y));
    end
end