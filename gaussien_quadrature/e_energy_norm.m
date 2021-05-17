function [err] = e_energy_norm(grad_u, u, k, E, dx, DX_u)
% computation of the true error in the energy norm
%   given by || k^1/2 * grad(u - v)||

% The grad operator on the averages is the central difference on the
% cell centers

% N = length(X);

a_1 = DX_u * u;

err = 0.;
for cell = 1:length(u)
    err = err + integral(@(x) k(x) .* (grad_u(x) - a_1(cell)).^2, E(cell), E(cell + 1));
end

a_l = u(1) / dx;
err = err + integral(@(x) k(x) .* (grad_u(x) - a_l).^2, 0, E(1));

a_r = -u(end) / dx;
err = err + integral(@(x) k(x) .* (grad_u(x) - a_r).^2, E(end), 1);

err = sqrt(err);

end
