function [err] = e_k_flux(q_fabricated, q, qL, qR, k, E, dx, DX_q)
% || k^(-1/2)(q - r)||
%   Detailed explanation goes here


a_0 = q(1:end - 1);
a_1 = DX_q * q;

err = 0.;
for cell = 1:length(E) - 1
    err = err + integral(@(x) 1 ./ k(x) .* (q_fabricated(x) - a_0(cell) - a_1(cell) .* (x - E(cell))).^2, E(cell), E(cell + 1));
end

a_0 = qL;
a_1 = (q(1) - qL) * 2 / dx;

err = err + integral(@(x) 1 ./ k(x) .* (q_fabricated(x) - a_0 - a_1 .* x).^2, 0, E(1));

a_0 = q(end);
a_1 = (qR - q(end)) * 2 / dx;
err = err + integral(@(x) 1 ./ k(x) .* (q_fabricated(x) - a_0 - a_1 .* (x - E(end))).^2, E(end), 1);

err = sqrt(err);

end

