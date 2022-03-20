

function [Ex, Ey, eFlowx, eFlowy, Vmap] = getG (Box, num_boxes, V0, sigma_out, sigma_in, l, w)
% Dimesions
nx = l;
ny = w; 

cMap = (1/sigma_out) * ones(nx, ny);

% Boundary conditions
boundary_left = 0;
boundary_right = V0;

% distance between points
d = 1;
d2 = d^2;

for n = 1:num_boxes
cMap(Box{n}.x(1):Box{n}.x(2), Box{n}.y(1):Box{n}.y(2)) = (sigma_in); %Box positions of cMap assigned the desired internal sigma
end

% Matrices
G = sparse(nx*ny);
B = zeros(1,nx*ny);

%--------------------------------------------------------------------------
% Iterative Solution
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            B(n) = boundary_left;

        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
            B(n) = boundary_right;

        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/d2;
            rxp = (cMap(i,j) + cMap(i+1,j))/d2;
            rym = (cMap(i,j) + cMap(i,j-1))/d2;
            ryp = (cMap(i,j) + cMap(i,j+1))/d2;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

%--------------------------------------------------------------------------
%Potential Solution

V = G\B';
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        Vmap(i, j) = V(n);
    end
end


% Electric Field Solutions
for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

% Current Solutions

eFlowx = cMap .* Ex;
eFlowy = cMap .* Ey;
% 
% current = (eFlowx.^2+eFlowy.^2).^0.5;

end