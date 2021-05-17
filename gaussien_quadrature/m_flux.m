function e_norm = m_flux(u, dx, nx, q, qL, qR, X, f)
    e_norm = 0;
    % create u of right length
    uu = [0; u; 0];
    XX = [0, X, 1];
    
    for cell = 2:nx+1
        b1 = 1 / dx * (q(cell) - q(cell - 1));
        fi = f(XX(cell));
        e_norm = e_norm + dx * (fi - b1)^2;
    end
    
    %left_border
    b1 = 2 / dx * (q(1) - qL);
    fi = f(XX(1));
    e_norm = e_norm + dx/2 * (fi - b1)^2;
    
    %right border
    b1 = 2 / dx * (qR - q(end));
    fi = f(XX(nx + 2));
    e_norm = e_norm + dx * (fi - b1)^2;

    e_norm = sqrt(e_norm);
end