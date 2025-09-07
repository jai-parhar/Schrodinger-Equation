
% Convergence test 1
idtype = 0;
vtype = 0;

idpar = [3];
tmax = 0.25;
lambda = 0.1;

% Level 6
[x_6, t_6, psi_6, psire_6, psiim_6, psimod_6, prob_6, v_6] = ...
sch_1d_cn(tmax, 6, lambda, idtype, idpar, vtype, []);

% Level 7
[x_7, t_7, psi_7, psire_7, psiim_7, psimod_7, prob_7, v_7] = ...
sch_1d_cn(tmax, 7, lambda, idtype, idpar, vtype, []);

% Level 8
[x_8, t_8, psi_8, psire_8, psiim_8, psimod_8, prob_8, v_8] = ...
sch_1d_cn(tmax, 8, lambda, idtype, idpar, vtype, []);

% Level 9
[x_9, t_9, psi_9, psire_9, psiim_9, psimod_9, prob_9, v_9] = ...
sch_1d_cn(tmax, 9, lambda, idtype, idpar, vtype, []);

% Find deltas, skipping over every 2nd element in second array
dpsi_6 = psi_6 - psi_7(1:2:end, 1:2:end);
dpsi_7 = psi_7 - psi_8(1:2:end, 1:2:end);
dpsi_8 = psi_8 - psi_9(1:2:end, 1:2:end);

% Uncomment to make first plot
%{
% Plot convergence test
hold on;
plot(t_6, rms(dpsi_6, 2));
plot(t_7, 4*rms(dpsi_7, 2));
plot(t_8, (4^2)*rms(dpsi_8, 2));
title('Convergence Test 1 of Levels l = 6, 7, 8, 9')
xlabel('t')
ylabel('||d\psi^l||_2')
legend('||d\psi^{6}||_2', '4||d\psi^{7}||_2', '4^2||d\psi^{8}||_2', "Location","northwest")
drawnow;
%}

% Make grid to compute exact values.
[x_grid,t_grid] = meshgrid(x_9,t_9);

% Find exact solution
E = exp(-1i * (idpar(1)^2) * (pi^2) * t_grid) .* sin(idpar(1)*pi*x_grid);

% Uncomment to make second plot
%{
hold on;
plot(t_6, rms(E(1:8:end, 1:8:end) - psi_6, 2));
plot(t_7, 4 * rms(E(1:4:end, 1:4:end) - psi_7, 2));
plot(t_8, (4^2) * rms(E(1:2:end, 1:2:end) - psi_8, 2));
plot(t_9, (4^3) * rms(E - psi_9, 2));
title('Convergence Test 1 of Exact \psi - Numerical Approximation \psi^l')
xlabel('t')
ylabel('||E(\psi^l)||_2')
legend('||E(\psi^{6})||_2', '4||E(\psi^{7})||_2', '4^2||E(\psi^{8})||_2', '4^3||E(\psi^{9})||_2', "Location","northwest");
drawnow;
%}



% Convergence test 2
idtype = 1;
vtype = 0;

idpar = [0.50, 0.075, 0.0];
tmax = 0.01;
lambda = 0.01;

% Level 6
[x_6, t_6, psi_6, psire_6, psiim_6, psimod_6, prob_6, v_6] = ...
sch_1d_cn(tmax, 6, lambda, idtype, idpar, vtype, []);

% Level 7
[x_7, t_7, psi_7, psire_7, psiim_7, psimod_7, prob_7, v_7] = ...
sch_1d_cn(tmax, 7, lambda, idtype, idpar, vtype, []);

% Level 8
[x_8, t_8, psi_8, psire_8, psiim_8, psimod_8, prob_8, v_8] = ...
sch_1d_cn(tmax, 8, lambda, idtype, idpar, vtype, []);

% Level 9
[x_9, t_9, psi_9, psire_9, psiim_9, psimod_9, prob_9, v_9] = ...
sch_1d_cn(tmax, 9, lambda, idtype, idpar, vtype, []);

% Find deltas, skipping over every 2nd element in second array
dpsi_6 = psi_6 - psi_7(1:2:end, 1:2:end);
dpsi_7 = psi_7 - psi_8(1:2:end, 1:2:end);
dpsi_8 = psi_8 - psi_9(1:2:end, 1:2:end);

% Uncomment to make plot
%{
% Plot convergence test
clf;
hold on;
plot(t_6, rms(dpsi_6, 2));
plot(t_7, 4*rms(dpsi_7, 2));
plot(t_8, (4^2)*rms(dpsi_8, 2));
title('Convergence Test 2 of Levels l = 6, 7, 8, 9')
xlabel('t')
ylabel('||d\psi^l||_2')
legend('||d\psi^{6}||_2', '4||d\psi^{7}||_2', '4^2||d\psi^{8}||_2', "Location", "northwest")
drawnow;
%}
