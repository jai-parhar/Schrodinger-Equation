function [x t psi psire psiim psimod prob v] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
    % Solves the 1D Schrodinger equation using CN method
    % Inputs
    %
    % tmax: Maximum integration time
    % level: Discretization level
    % lambda: dt/dx
    % idtype: Selects initial data type
    % idpar: Vector of initial data parameters
    % vtype: Selects potential type
    % vpar: Vector of potential parameters
    %
    % Outputs
    %
    % x: Vector of x coordinates [nx]
    % t: Vector of t coordinates [nt]
    % psi: Array of computed psi values [nt x nx]
    % psire Array of computed psi_re values [nt x nx]
    % psiim Array of computed psi_im values [nt x nx]
    % psimod Array of computed sqrt(psi psi*) values [nt x nx]
    % prob Array of computed running integral values [nt x nx]
    % v Array of potential values [nx]

    nx = 2^(level) + 1;
    dx = 2^(-level);

    dt = lambda * dx;
    nt = round(tmax/dt) + 1;

    x = linspace(0, 1, nx);
    t = linspace(0, tmax, nt);
    psi0 = getPsi0_1D(x, idtype, idpar);
    v = getV_1D(x, nx, vtype, vpar);

    % Set up tridiagonal matrices
    A_l_diag = (0.5 / dx^2) * ones(nx,1);
    A_m_diag = (1i/dt - 1/dx^2) * ones(nx,1) - (0.5*v);
    A_u_diag = A_l_diag;

    % Fix boundary conditions on tridiag matrix
    A_l_diag(end-1) = 0;
    A_m_diag(1) = 1;
    A_m_diag(end) = 1;
    A_u_diag(2) = 0;
    
    % Generate sparse matrix
    A = spdiags([A_l_diag, A_m_diag, A_u_diag], -1:1, nx, nx);

    % Generates a matrix that will calculate the integral
    % Assumes fixed dx (which is the case here)
    % Each element follows form of running sum, using this on the values of
    % |psi|^2 will give running integral
    trap_sum_mat = tril(ones(nx)) + tril(ones(nx), -1);
    trap_sum_mat(:, 1) = 1;
    
    % Matrix which turns some function f [nx array] into integral given
    % fixed dx
    prob_integral_mat = (dx/2) * trap_sum_mat;

    % Initialize output
    psi = zeros(nt, nx);
    psire = zeros(nt, nx);
    psiim = zeros(nt, nx);
    psimod = zeros(nt, nx);
    prob = zeros(nt, nx);

    psi(1, :) = psi0;
    psire(1, :) = real(psi0);
    psiim(1, :) = imag(psi0);
    psimod(1, :) = abs(psi0);
    prob(1, :) = prob_integral_mat * transpose(psimod(1, :).^2);

    for n = 1:nt-1
        % Find next psi
        
        % Term on right hand side, B * psi_n
        term  = zeros(nx,1);
        term(2:nx-1) = - 0.5*(psi(n, 1:nx-2) + psi(n, 3:nx))/dx^2 ...
        + (0.5*transpose(v(2:nx-1)) + (1/dx^2) + (1i/dt)) .* psi(n, 2:nx-1);
        
        % Boundary conditions for this term
        term(1) = 0;
        term(end) = 0;

        % Find the next psi by solving the system
        psi_next = A \ term;
        
        % Boundary conditions for psi
        psi_next(1) = 0;
        psi_next(nx) = 0;

        % Store values
        psi(n+1, :) = psi_next;
        psire(n+1, :) = real(psi_next);
        psiim(n+1, :) = imag(psi_next);
        psimod(n+1, :) = abs(psi_next);
        
        prob(n+1, :) = prob_integral_mat * transpose(psimod(n+1, :).^2);

    end
end




function [psi0] = getPsi0_1D(x, idtype, idpar)
    % Gives initial data for 1D Schrodinger equation based on project
    % requirements
    % Inputs
    % 
    % x: Vector of x coordinates [nx]
    % idtype: Selects initial data type
    % idpar: Vector of initial data parameters
    %
    % Outputs
    % 
    % psi0: Array of computed psi values [nx]

    if idtype == 0
        m = idpar(1);

        psi0 = sin(m*pi*x);
    elseif idtype == 1
        x0 = idpar(1);
        delta = idpar(2);
        p = idpar(3);
        
        psi0 = exp(1i*p*x) .* exp( -((x-x0)./delta).^2 );
    else
        error("Unknown intial data type %d", idtype);
    end

    % Boundary conditions
    psi0(1) = 0;
    psi0(end) = 0;

end


function [v] = getV_1D(x, nx, vtype, vpar)
    % Gives potential for 1D Schrodinger equation based on project
    % requirements
    % Inputs
    % 
    % x: Vector of x coordinates [nx]
    % nx: Number of x coordinates
    % vtype: Selects potential type
    % vpar: Vector of potential parameters
    %
    % Outputs
    % 
    % v Array of potential values [nx]
    if vtype == 0
        % No potential
        v = zeros(nx,1);
    elseif vtype == 1
        v = zeros(nx,1);
        % 1D Barrier/Well
        xmin = vpar(1);
        xmax = vpar(2);
        Vc = vpar(3);
        
        v((x >= xmin) & (x <= xmax)) = Vc;
    else
        error("Unknown potential type %d", idtype);
    end
end