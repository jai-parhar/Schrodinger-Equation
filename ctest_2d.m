% Convergence test 2D

idtype = 0;
vtype = 0;
idpar = [2, 3];
tmax = 0.05;
lambda = 0.05;

% Level 6
[x_6, y_6, t_6, psi_6, psire_6, psiim_6, psimod_6, v_6] = ...
sch_2d_adi(tmax, 6, lambda, idtype, idpar, vtype, []);

% Level 7
[x_7, y_7, t_7, psi_7, psire_7, psiim_7, psimod_7, v_7] = ...
sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, []);

% Level 8
[x_8, y_8, t_8, psi_8, psire_8, psiim_8, psimod_8, v_8] = ...
sch_2d_adi(tmax, 8, lambda, idtype, idpar, vtype, []);

% Level 9
[x_9, y_9, t_9, psi_9, psire_9, psiim_9, psimod_9, v_9] = ...
sch_2d_adi(tmax, 9, lambda, idtype, idpar, vtype, []);


% Difference between 2 discretization levels
dpsi_6 = psi_6 - psi_7(1:2:end, 1:2:end, 1:2:end);
dpsi_7 = psi_7 - psi_8(1:2:end, 1:2:end, 1:2:end);
dpsi_8 = psi_8 - psi_9(1:2:end, 1:2:end, 1:2:end);

% RMS values
dpsi_6_rms = rms(dpsi_6, [2, 3]);
dpsi_7_rms = rms(dpsi_7, [2, 3]);
dpsi_8_rms = rms(dpsi_8, [2, 3]);

% Calculate exact solution
[t_grid, x_grid, y_grid] = meshgrid(t_9, x_9, y_9);
E = exp(-1i*(idpar(1)^2 + idpar(2)^2)*(pi^2)*t_grid) .* sin(idpar(1)*pi*x_grid) .* sin(idpar(2)*pi*y_grid);

% Difference between exact solution and discretization level
dE_6 = E(1:8:end, 1:8:end, 1:8:end) - psi_6;
dE_7 = E(1:4:end, 1:4:end, 1:4:end) - psi_7;
dE_8 = E(1:2:end, 1:2:end, 1:2:end) - psi_8;
dE_9 = E - psi_9;

% RMS values
dE_6_rms = rms(dE_6, [2, 3]);
dE_7_rms = rms(dE_7, [2, 3]);
dE_8_rms = rms(dE_8, [2, 3]);
dE_9_rms = rms(dE_9, [2, 3]);


% Plot convergence test of levels
hold on;
plot(t_6, dpsi_6_rms);
plot(t_7, 4 * dpsi_7_rms);
plot(t_8, (4^2) * dpsi_8_rms);
title('Convergence Test of 2D ADI Solver, Levels l = 6, 7, 8, 9')
xlabel('t')
ylabel('||d\psi^l||_2')
legend('||d\psi^{6}||_2', '4||d\psi^{7}||_2', '4^2||d\psi^{8}||_2', "Location","northwest")
drawnow;

%{
% Plot convergence test of exact solution
hold on;
plot(t_6, dE_6_rms);
plot(t_7, dE_7_rms);
plot(t_8, dE_8_rms);
plot(t_9, dE_9_rms);
title('Convergence Test of 2D ADI Solver, Exact \psi - Numerical Approximation \psi^l')
xlabel('t')
ylabel('||E(\psi^l)||_2')
legend('||E(\psi^{6})||_2', '4||E(\psi^{7})||_2', '4^2||E(\psi^{8})||_2', '4^3||E(\psi^{9})||_2', "Location","northwest");
drawnow;
%}

