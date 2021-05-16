function e_norm = m_energy_norm(u, q, qL, qR, k, dx, E, DX_u, DX_q)
% ||k^(-1/2) (r + k grad(v))||

e_norm = 0.;

% flux interpolation (q)
a_0 = q(1:end - 1);
a_1 = DX_q * q;

% potential interpolation (u)
b_1 = DX_u * u;

for cell = 1:length(E) - 1
    e_norm = e_norm + integral(@(x) 1 ./ k(x) .* (a_0(cell) + a_1(cell) .* (x - E(cell)) + k(x) .* b_1(cell)).^2, E(cell), E(cell + 1));
end

a_0 = qL;
a_1 = (q(1) - qL) * 2 / dx;
b_1 = u(1) / dx;
e_norm = e_norm + integral(@(x) 1 ./ k(x) .* (a_0 + a_1 * x + k(x) * b_1).^2, 0, E(1));

a_0 = q(end);
a_1 = (qR - q(end)) * 2 / dx;
b_1 = -u(end) / dx;
e_norm = e_norm + integral(@(x) 1 ./ k(x) .* (a_0 + a_1 * (x - E(end)) + k(x) * b_1).^2, E(end), 1);

e_norm = sqrt(e_norm);

end