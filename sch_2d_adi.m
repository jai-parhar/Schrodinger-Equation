
function [x, y, t, psi, psire, psiim, psimod, v] = ...
         sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
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
    % y: Vector of y coordinates [ny]
    % t: Vector of t coordinates [nt]
    % psi: Array of computed psi values [nt x nx x ny]
    % psire Array of computed psi_re values [nt x nx x ny]
    % psiim Array of computed psi_im values [nt x nx x ny]
    % psimod Array of computed sqrt(psi psi*) values [nt x nx x ny]
    % v Array of potential values [nx x ny]
    
    nx = 2^(level) + 1;
    ny = nx;
    dx = 2^(-level);
    dy = dx;

    dt = lambda * dx;
    nt = round(tmax/dt) + 1;

    x = linspace(0, 1, nx);
    y = x;
    t = linspace(0, tmax, nt);
    
    % Initialize psi
    psi = zeros(nt, nx, ny);
    psi(1, :, :) = getPsi0_2D(x, y, idtype, idpar);
    
    % Initialize potential
    v = getV_2D(x, y, nx, ny, vtype, vpar);
    
    % Tridiagonal matrices Set-up
    % Find terms for upper (u), main (m), and lower (l) diagonals for the
    % matrix A
    A_l_diag = -(1i*dt/(2*dx^2)) .* ones(nx, 1);
    A_m_diag = (1 + (1i*dt/(dx^2))) .* ones(nx, 1);
    A_u_diag = A_l_diag;

    % Find terms for upper (u) and lower (l) diagonals for the
    % matrix B, we're gonna find the main diagonal later
    B_l_diag = A_l_diag;
    B_u_diag = B_l_diag;

    % Fix the boundary conditions for upper, main, lower diagonals for A
    A_l_diag(nx-1) = 0.0;
    A_m_diag(1) = 1; 
    A_m_diag(nx) = 1;
    A_u_diag(2) = 0.0;
    
    % Fix the boundary conditions for upper and lower diag for B
    B_u_diag(2) = 0.0;
    B_l_diag(ny-1) = 0.0;

    % Sparse matrices
    A = spdiags([A_l_diag, A_m_diag, A_u_diag], -1:1, nx, nx);
    
    % Start at n = 1, go until second last n. Finding value for next n as
    % opposed to going from n=2 to nt because this is easier to
    % understand. Please don't change this I know the darkness in your 
    % heart and soul and your need to change stuff to make it "better" 
    % just please whatever you do just do leave it like this and Also 
    % change this comment later.
    for n = 1:nt-1
        % Stage 1: Finding n+1/2

        % The first half of the right hand side term
        term2 = zeros(nx,ny);
        term2(2:nx-1,2:ny-1) = (1i*dt/(2*dx^2)) .* squeeze((psi(n,2:nx-1,3:ny))) + ... 
                                (1 - (1i*dt/(dx^2)) - (1i*dt/(2*dx^2)) .* dx^2 .* v(2:nx-1,2:ny-1)) .* squeeze((psi(n,2:nx-1,2:ny-1))) + ... 
                                (1i*dt/(2*dx^2)) .* squeeze((psi(n,2:nx-1,1:ny-2)));
        
        % Full right hand side term 
        term = zeros(nx,ny);
        term(2:nx-1, 2:ny-1) = (1i*dt/(2*dx^2)) * term2(3:nx,2:ny-1) + ... 
                            (1-(1i*dt/(dx^2))) * term2(2:nx-1,2:ny-1) + ... 
                            (1i*dt/(2*dx^2)) * term2(1:nx-2,2:ny-1);

        % Actually solve for psi_n+1/2
        psi_half_step = zeros(nx,ny);
        for curr_y = 1:ny
            term_col = term(:,curr_y);

            % Boundary conds
            term_col(1) = 0.0;
            term_col(1) = 0.0;
            
            % Find that half step term
            psi_half_step(:,curr_y) = A \ term_col;
        end
    
        % Stage 2: Finding n+1
        for curr_x = 1:nx
            % Finally making that middle diagonal. This is a big day for
            % me.
            B_m_diag = (1i*dt/(dx^2)) + 1 + (1i*dt/(2*dx^2)).*dx^2 .* (v(curr_x,:)).'; 
            B_m_diag(1) = 1.0; B_m_diag(ny) = 1.0;

            % Making matrix B
            B = spdiags([B_l_diag, B_m_diag, B_u_diag], -1:1, nx, ny);

            % Finding the full right hand term for n+1 
            term = zeros(nx,1);
            term(2:ny-1) = (psi_half_step(curr_x,2:ny-1));

            % Boundary conds
            term(1) = 0.0;
            term(ny) = 0.0;

            % Store next value
            psi(n+1,curr_x,:) = B \ term;
        end
    end

    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);
end




function [psi0] = getPsi0_2D(x, y, idtype, idpar)
    % Gives initial data for 2D Schrodinger equation based on project
    % requirements
    % Inputs
    % 
    % x: Vector of x coordinates [nx]
    % y: Vector of y coordinates [ny]
    % idtype: Selects initial data type
    % idpar: Vector of initial data parameters
    %
    % Outputs
    % 
    % psi0: Array of computed psi values [nx x ny]
    
    [x_grid, y_grid] = meshgrid(x,y);

    if idtype == 0
        mx = idpar(1);
        my = idpar(2);

        psi0 = (sin(my*pi*y_grid)) .* sin(mx*pi*x_grid);
    elseif idtype == 1
        x0 = idpar(1); 
        y0 = idpar(2);
        deltax = idpar(3); 
        deltay = idpar(4);
        px = idpar(5); 
        py = idpar(6);

        
        psi0 = (exp(1i*py*y_grid) .* (exp( -((y_grid-y0)./deltay).^2 ))) ...
                .* (exp(1i*px*x_grid) .* exp( -((x_grid-x0)./deltax).^2 ));
    else
        error("Unknown intial data type %d", idtype);
    end
    
    % Boundary Conditions
    psi0(1, :) = 0;
    psi0(:, 1) = 0;
    psi0(:, end) = 0;
    psi0(end, :) = 0;

    psi0 = squeeze(psi0);

end

function [v] = getV_2D(x, y, nx, ny, vtype, vpar)
    % Gives potential for 2D Schrodinger equation based on project
    % requirements
    % Inputs
    % 
    % x: Vector of x coordinates [nx]
    % y: Vector of y coordinates [ny]
    % vtype: Selects potential type
    % vpar: Vector of potential parameters
    %
    % Outputs
    % 
    % v Array of potential values [nx x ny]

    [x_grid, y_grid] = meshgrid(x,y);

    if vtype == 0
        v = zeros(nx, ny);
    elseif vtype == 1
        % 2D Barrier/Well
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        Vc = vpar(5);
        
        v = Vc * ((x_grid >= xmin) & (x_grid <= xmax) ...
               & (y_grid >= ymin) & (y_grid <= ymax));
    elseif vtype == 2
        % 2D Double slit
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        Vc = vpar(5);

        j_prime = round((size(y, 2) - 1) / 4) + 1;

        v = zeros(nx, ny);

        v(:, j_prime) = Vc;
        v(:, j_prime+1) = Vc;
        v(x >= x1 & x <= x2 | x >= x3 & x <= x4, j_prime) = 0;
        v(x >= x1 & x <= x2 | x >= x3 & x <= x4, j_prime+1) = 0;

    else
        error("Unknown potential type %d", idtype);
    end
end
