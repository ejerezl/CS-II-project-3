
% Poisson equation: -Div(k(x,y)Grad(u)) = f

% ex = 3;
ex = 8;

if ex == 1   % Log case
    k1 = @(x) ones(1, length(x)); k_const = true;
    u_fabricated = @(x) (1 - x) .* x .* log(x + 1);
    q_fabricated = @(x) x .* log(x + 1) - (1 - x) .* log(x + 1) - ((1 - x) .* x)./(x + 1);
    f = @(x) 1 ./ (x + 1).^2 .* ((2 * x.^2 + 4 * x + 2) .* log(x + 1) + 3 * x.^2 + 3 * x - 2);
elseif ex == 2  % easy case
    k1 = @(x) ones(1, length(x)); k_const = true;
    u_fabricated = @(x) -0.5*x.^2 + 0.5*x;
    q_fabricated = @(x) x - 0.5;
    f = @(x) ones(1, length(x));
elseif ex == 3  % sin(x^2)*(1 - x) k konst case
    k1 = @(x) ones(1, length(x)); k_const = true;
    u_fabricated = @(x) sin(x) .* (1 - x);
    q_fabricated = @(x) (x - 1) .* cos(x) + sin(x);
    f = @(x) 2 * cos(x) - (x - 1) .* sin(x);
elseif ex == 4  % exp case // hard case
    k1 = @(x) 1/2 * x + 1; k_const = false;
    u_fabricated = @(x) (1 - exp(x)) .* (x - 1);
    grad_u_fab = @(x) (1 - exp(x)) - (x - 1) .* exp(x);
    q_fabricated = @(x) -1/2 * x - 1 + 1/2 * x.^2 .* exp(x) + x .* exp(x);
    f = @(x) 1/2 * ( ( x.^2 + 4 * x + 2) .* exp(x) - 1);
elseif ex == 5  % easy case with nonkonstant k
    k1 = @(x) 1/2 * x + 1; k_const = false;
    u_fabricated = @(x) - 1/2 * x.^2 + 1/2 * x;
    q_fabricated = @(x) 1/2 * x.^2 + 3/4 * x + 1/2;
    f = @(x) x + 3/4;
elseif ex == 6  % uncontinous k
    k1 = @(x) arrayfun(@incontin_k, x); k_konst = false;
    f = @(x) ones(1, length(x));
elseif ex == 7  % uncontinous k / known sol
    k1 = @(x) arrayfun(@incontin_k, x); k_konst = false;
    f = @(x) arrayfun(@incontin_f, x);
    u_fabricated = @(x) sin(pi * x);
    q_fabricated = @(x) arrayfun(@incontin_q, x);
elseif ex == 8
    k1 = @(x) 2 - x; k_const = false;
    f = @(x) 4*pi*pi*sin(2*pi*x).*(2-x) + 2*pi*cos(2*pi*x);
    u_fabricated = @(x) sin(2*pi*x);
    grad_u_fab = @(x) 2 * pi * cos(2*pi*x);
    q_fabricated = @(x) -2*pi*cos(2*pi*x).*(2-x);
end


% create n logspaced grid widths
n = 7;
nx_start = 1e1;
nx_end = 1e1;
nx_vec = round(logspace(log10(nx_start), log10(nx_end), n));

dx_vec = [];
star_vec = [];
maj_vec = [];
energy_err = [];
kflux_err = [];
divflux_err = [];

L2_err = [];

% ACHTUNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!11!1!elf1!!111!!!!!!!!1
ex = -1;

for nx = nx_vec

	% Build matrices
	[A, b, G, D, K, dx, X, E, DX_u, DX_q] = assemble1d(nx,f,k1);
	dx_vec = [dx_vec, dx];

	% Solve linear system and compute flux
	u = A\b;
    if ex == 1
        %[q, qL, qR] = flux1d_sym1(u, K, G, f, dx);
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    elseif ex == 2
        %[q, qL, qR] = flux1d_sym2(u, K, G, f, dx);
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    elseif ex == 3
        %[q, qL, qR] = flux1d_sym3(u, K, G, f, dx);
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    elseif ex == 4
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    elseif ex == 5
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    else
        [q, qL, qR] = flux1d(u, K, G, f, dx);
    end
    
    % Poincare constant
    if k_const
        C_omega = 1/(sqrt(k1(0.5)) * pi);
        %C_omega = 1 / sqrt(eigs(A, 1, 'smallestabs'));
    else
        % approximate as the smallest non-one eigenvalue of A
        % candidates = eigs(A,10,'smallestabs');
        % idx = find(candidates ~= 1)
        % C_omega = candidates(idx(1));
        C_omega = 1 / sqrt(eigs(A, 1, 'smallestabs'));
    end
    
    %Plot difference of solutions
if nx == nx_start 
    figure(1)
    u_fabricatedvec = u_fabricated(X).'; 
    subplot(1,2,1)
    plot(X,abs(u-u_fabricatedvec)./abs(u))
    xlabel('x', 'FontSize', 16)
    ylabel('relative difference','FontSize',16)
    title({'relative difference between analytical','and numerical potential'},'FontSize',14)

%     X_plot = zeros(nx+3,1); %new grid
%     X_plot(1) = 0;
%     X_plot(end) = 1;
%     for i=2:nx
%         X_plot(i) = 0.5*dx + dx*(i-2);
%     end
    q_fabricatedvec = q_fabricated(E).';
%     w = size(E);
%     q_vec = zeros(w(2),1).';
%     q_vec(1) = qL;
%     q_vec(end) = qR;
%     q_vec(2:end-1) = q;
    subplot(1,2,2)
    plot(E,abs(q_fabricatedvec-q)./abs(q))
    xlabel('x', 'FontSize', 16)
    ylabel('relative difference','FontSize',16)
    title({'relative difference between analytical','and numerical flux'},'FontSize',14)
end    
    if ex == -1
        energy_ev = e_energy_norm(grad_u_fab, u, k1, E, dx, DX_u);
        k_er = e_k_flux(q_fabricated, q, qL, qR, k1, E, dx, DX_q);
        div_er = em_div_flux(f, q, qL, qR, E, dx, DX_q);
        
        %maj_energy = m_energy_norm_rechner_new(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm(u, q, qL, qR, k1, dx, E, DX_u, DX_q);
	% error in energy norm
    elseif ex == 1
        energy_ev = e_energy_norm_sym1(u_fabricated, u, k1, nx, dx, X);
        k_er = e_k_flux(q_fabricated, q, qL, qR, E, k1, X, dx); % expression to long for integral calculator
        div_er = e_div_flux(q_fabricated, q, qL, qR, E, dx, D);
        
        %maj_energy = m_energy_norm_rechner_new(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_k_const(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 2
        energy_ev = e_energy_norm_sym2(u, k1, nx, dx, X);
        k_er = e_k_flux_sym2(q, qL, qR, k1, X, dx, nx);
        div_er = em_div_flux_sym2(q, qL, qR, dx, nx, X);
        
        %maj_energy = m_energy_norm_rechner_new(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_k_const(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 3
        energy_ev = e_energy_norm_sym3(u_fabricated, u, k1, nx, dx, X);
        k_er = e_k_flux_sym3(q, qL, qR, k1, X, dx, nx);
        div_er = em_div_flux_sym3(q, qL, qR, dx, nx, X);
        
        %maj_energy = m_energy_norm_rechner_new(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_k_const(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 4
        energy_ev = e_energy_norm_sym4(u, nx, dx, X);
        k_er = e_k_flux_sym4(q, qL, qR, k1, X, dx, nx);
        div_er = em_div_flux_sym4(q, qL, qR, dx, nx, X);
        
        %maj_energy = m_energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_sym4(u, dx, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 5
        energy_ev = e_energy_norm_sym5(u, nx, dx, X);
        k_er = e_k_flux_sym5(q, qL, qR, k1, X, dx, nx);
        div_er = em_div_flux_sym5(q, qL, qR, dx, nx, X);
        
        %maj_energy = m_energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_sym4(u, dx, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 6
        energy_ev = 0;
        k_er = 0;
        div_er = 0;
        
        %maj_energy = m_energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_k_const(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    elseif ex == 7
        energy_ev = 0;
        k_er = 0;
        div_er = 0;
        
        %maj_energy = m_energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
        maj_energy = m_energy_norm_rechner_k_const(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    else
        energy_ev = e_energy_norm(u_fabricated, u, k1, X, dx, DX_u);
        k_er = e_k_flux(q_fabricated, q, qL, qR, E, k1, X, dx);
        div_er = e_div_flux(q_fabricated, q, qL, qR, E, dx, D);
        
        maj_energy = m_energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
        %[maj_energy2, maj_energy, maj_energy1] = m_energy_norm_rechner_kaja(u_fabricated, u, q, qL, qR, k1, X, dx, nx, DX_u, E);
        %maj_energy = m_energy(u, q, k1, qL, qR, dx, X, DX_u, DX_q);
    end
    energy_err = [energy_err, energy_ev];
    
    kflux_err = [kflux_err, k_er];
    divflux_err = [divflux_err, div_er];

	star_err = energy_ev + k_er + C_omega * div_er;
	star_vec = [star_vec, star_err];


	% compute energy norm
	% e_norm = energy_norm_rechner(u, dx, k1, nx, q, qL, qR, X);
	% maj_energy - e_norm
	maj_flux = m_flux(u, dx, nx, q, qL, qR, X, f);

	% majorant = e_norm + C_omega * e_norm2
    if ex == 4 || ex == 5 || ex == 6
        majorant = maj_energy + C_omega * maj_flux;
    else
        majorant = maj_energy + C_omega * div_er;
    end
	maj_vec = [maj_vec, majorant];
    
    if ex == 6
        l2_error = 0;
    else
        l2_error = sqrt(dx) * norm(u - u_fabricated(transpose(X)));
    end
    L2_err = [L2_err, l2_error];
end

eff_index = 3 * maj_vec ./ star_vec;
eff_index2 = maj_vec ./ energy_err;

% plot
% blue = [0, 0.4470, 0.7410];
% orange = [0.8500, 0.3250, 0.0980];
% yellow = [0.9290, 0.6940, 0.1250];
% violet = [0.4940, 0.1840, 0.5560];
% green = [0.4660, 0.6740, 0.1880];
% cyan = [0.3010, 0.7450, 0.9330];
% winered = [0.6350, 0.0780, 0.1840];

violet = [0.4940, 0.1840, 0.5560];
winered = [0.6350, 0.0780, 0.1840];

blue = [51, 64, 182] / 255;
orange = [253, 143, 0] / 255;
green = [15, 192, 9] / 255;
red = [255, 0, 0] / 255;
black = [0,0,0];


figure(2)
yyaxis left
title({'Error Comparison'}, 'FontSize', 15);
loglog(dx_vec, 3*maj_vec, '-.^')
hold on
loglog(dx_vec, maj_vec, ':*')
loglog(dx_vec, star_vec, '--o')
% loglog(dx_vec, energy_err, '-sg')
loglog(dx_vec, energy_err, '-s')
%loglog(dx_vec, kflux_err, '-pc')
% loglog(dx_vec, divflux_err, '-hm')
% loglog(dx_vec, res_2, '-o')
%loglog(dx_vec, dx_vec, '--', 'color', violet)
% loglog(dx_vec, L2_err, '-hm')
% loglog(dx_vec, dx_vec.^2, '--c')

ax = gca;
ax.FontSize = 13;


%loglog(dx_vec, dx_vec.^2, '--c')
%loglog(dx_vec, dx_vec.^(1/2), '--c')

ylabel('Error', 'FontSize', 15)

yyaxis right
plot(dx_vec, eff_index, '-.*')
plot(dx_vec, eff_index2, '--o')
% yline(3.5)
yline(3, 'color', black)
yline(1, 'color', black)
yline(0.5)

ylabel('Efficiency Index', 'FontSize', 15)

xlabel('$\Delta x$', 'Interpreter', 'latex', 'FontSize', 15)

% legend({'e star norm', 'majorant', 'e energy norm', 'dx', 'dx^2'}, 'location', 'northwest')

% set(findall(gcf,'-property','FontSize'),'FontSize',15)
legend({'$3\mathcal{M}(v,r; f)$','$\mathcal{M}(v,r; f)$', '$\|(e_v, e_r)\|_*$', '$\|k^{1/2} \nabla (u - v)\|$', '$I_{v,r}^{eff}$', '$I^{eff}_v$'}, 'location', 'southeast', 'Interpreter', 'latex', 'FontSize', 12)
% xticks(dx_vec)
% xticks(dx_vec(end:-1:1))
% XTickLabels = cellstr(num2str(round(log10(dx_vec(end:-1:1))), '10^%d'));
% xticklabels(XTickLabels)

