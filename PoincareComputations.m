%Get permeability from separate m-file
%k = @(x,y) kSquare(x, y, 0.6, 0.3, 1e-6, 1);
k = @(x,y) kCross(x, y, 0.6, 0.3, 1e-6, 1);
%k = @(x,y) 1;

% Source term which is only needed to assemble (coule be any lambda)
f = @(x,y) -2*(y*y-y)-2*(x*x-x);


%Determine number of cells in each direction
nx = 257;
ny = 257;
dx = 1/(nx-1);
dy = 1/(ny-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);


%Computing smallest eigenvalues
[V,L]= eigs(A, 15, "smallestabs");
%Finding "correct" sallest eigenvalue
entry = 1;
for i=1:10
    if L(i,i)< (1-1e-6)
        entry = i;
        break;
    end
    if L(i,i) > (1+1e-6)
        entry = i;
        break;
    end
end

%Computing Poincare constant (C and C_o should be equal)
C = 1/sqrt(norm(L(:,entry),2))
v_min = V(:,entry);
dv_min = (sqrt(K).*G)*v_min;
C_o = norm(v_min,2)/norm(dv_min,2)

% Plotting Permeability and smallest egenfunction
[X,Y] = meshgrid(0:dy:1,0:dx:1);
U = reshape(v_min,nx,ny);
k_vect = zeros(nx*ny,1);
for i = 1:ny*nx
    k_vect(i) = k(cells(i,1),cells(i,2));
end
k_plot = reshape(k_vect,nx,ny);
figure(1)
set(gcf,'Position',[100 100 1200 500])
subplot(1,2,1);
surf(X,Y,U,'EdgeColor','none')
title('smallest eigenfunction')
colorbar()
axis([0 1 0 1])
subplot(1,2,2);
surf(X,Y,k_plot,'EdgeColor','none')
title('Permeability')
axis([0 1 0 1])
colorbar()