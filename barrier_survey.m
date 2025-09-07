% Barrier Survey

% Parameter
tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1; % Boosted Gaussian
idpar = [0.40, 0.075, 20.0];
vtype = 1;

% Bounds we are checking between
x_1 = 0.8;
x_2 = 1.0;

% Should be 251 in final use case, put to lower for testing
n_lnV_0 = 251;

% Set up array of ln(V_0) and F_e (will take ln of that later, need to compute first)
lnV_0 = linspace(-2, 5, n_lnV_0);
F_e = zeros(n_lnV_0,1);

for lnV_0_curr = 1:length(lnV_0)
    
    % Set up parameters for the potential
    % Working with ln(V_0) so we need to take exponential to convert back
    % to potential
    vpar = [0.6, 0.8, exp(lnV_0(lnV_0_curr))];

    % Compute solution
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

    % Find index of closest element to x_1
    [ d, ix ] = min( abs( x - x_1 ) );

    % Calculate probabilities
    
    % Mean over all time of the whole probability dist. Should be 1 if
    % normalized. Used to normalize.
    P_nx = mean(prob(:,end));

    % Mean over all time of all probability up to x_1, normalized
    P_x1 = mean(prob(:,ix)) / P_nx;
    
    % Mean over all time of all probability up to x_2, which is the end, normalized
    % Formula uses P_x2 = mean(prob(:,end))/P_nx; but this will just be 1 if 
    % properly normalized, and we just normalized it so it must be
    F_e_curr = (1 - P_x1)/(x(end)-x(ix));
    F_e(lnV_0_curr) = F_e_curr;

end

lnFe = log(F_e);

hold on;
plot(lnV_0, lnFe);
title(['Log-log Plot of Excess Fractional Probability between x=0.8 and x=1 vs Potential Barrier Strength'])
xlabel('ln(V_0)')
ylabel('ln(F_e(0.8,1.0))')
drawnow;



