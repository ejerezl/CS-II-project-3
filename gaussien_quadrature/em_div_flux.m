function [err] = em_div_flux(f, q, qL, qR, E, dx, DX_q)
% || f - div r||

a_1 = DX_q * q;

err = 0.;
for cell = 1:length(E) - 1
    err = err + integral(@(x) (f(x) - a_1(cell)).^2, E(cell), E(cell + 1));
end

a_1 = (q(1) - qL) * 2 / dx;
err = err + integral(@(x) (f(x) - a_1).^2, 0, E(1));

a_1 = (qR - q(end)) * 2 / dx;
err = err + integral(@(x) (f(x) - a_1).^2, E(end), 1);


err = sqrt(err);

end

