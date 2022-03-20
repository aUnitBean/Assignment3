
function aveCurrent = boxes(passageWidth)
%Duration of simulation
num_steps = 400;

%Potential Maximum
V0 = 1;

%Silicon Temperature
T = 300;

%Constants
C.m0 = 9.11 *10 ^ (-31);
C.mn = 0.26 * C.m0;
C.k = 1.381 * 10 ^ (-23);
C.q = 1.60217662 * 10 ^ (-19);


%Number of Particles
num_part = 30000;
part_plot = 20;

% Sigma Map
sigma_out = 1;
sigma_in = 10^-2;

% Dimesions
%passageWidth = passageWidth;
passageLength = 40;
l = 200; % length
w = 100; % width
silicon = zeros(w, l);

%Boxes!
num_boxes = 2;
Box = {};
Box{1}.y =[0 1/2*(w-passageWidth) ];
Box{1}.x =[1/2*(l-passageLength) 1/2*(l+passageLength)];
Box{2}.y =[1/2*(w+passageWidth) w];
Box{2}.x =[1/2*(l-passageLength) 1/2*(l+passageLength)];

%Get fields
[Ex, Ey, eFlowx, eFlowy, Vmap] = getG(Box, num_boxes, V0, sigma_out, sigma_in, l, w);
%Get trajectories
[all_x_positions, all_y_positions, part, aveCurrent] = ElectronTrajectories(Ex, Box, num_boxes, num_steps, w, l, num_part, C);

ElectronDensity (num_part, part, silicon, l, w)
end
