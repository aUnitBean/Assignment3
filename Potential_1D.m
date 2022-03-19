%Sarah Dolan, ELEC 4700, February 2022
%% Question 1 ai)
% RUN THIS FROM PART_1_LAPLACE
% The purpose of this code is solve for the electro static potential in a
% rectangular region using the finite difference method. The y boundaries
% are not fixed.


function V = Potential_1D(nX, nY, V0)

% rectangular region
% number of points
nx = nX;
ny = nY; 
n_iterations = 10000;

% boundary conditions
BC_left = 0;
BC_right = V0;

% matrix formation
V = zeros(nx, ny);

for k = 1 : n_iterations
        for i = 1 : nx           
           %checking boundary conditions
 
           switch i
                case 1
                    V_left= BC_left;
                    V_right = V(i + 1, :);
                case nx
                    V_left = V(i - 1, :);
                    V_right = BC_right;
                otherwise 
                    V_left = V(i - 1, :);
                    V_right = V(i + 1, :);
            end
                     
            V(i, :) = (V_left + V_right) / 2;
        end
end


end
