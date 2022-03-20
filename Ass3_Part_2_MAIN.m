%Duration of simulation
num_steps = 300;

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
num_part = 50000;
part_plot = 10;

% Sigma Map
sigma_out = 1;
sigma_in = 10^-2;

% Dimesions
passageWidth = 20;
passageLength = 40;
l = 200; % length
w = 100; % width

%Boxes!
num_boxes = 2;
Box = {};
Box{1}.y =[1 1/2*(w-passageWidth) ];
Box{1}.x =[1/2*(l-passageLength) 1/2*(l+passageLength)];
Box{2}.y =[1/2*(w+passageWidth) w];
Box{2}.x =[1/2*(l-passageLength) 1/2*(l+passageLength)];
% 
% Box{3}.x =[140 160];
% Box{3}.y =[45 55];
% 
% Box{4}.x =[40 60];
% Box{4}.y =[45 55];

%Get fields
[Ex, Ey, eFlowx, eFlowy, Vmap] = getG(Box, num_boxes, V0, sigma_out, sigma_in, l, w);
%Get trajectories
[all_x_positions, all_y_positions] = ElectronTrajectories(Ex, Box, num_boxes, num_steps, w, l, num_part, C);

ElectronDensity (num_part, part, silicon, length_silicon, width_silicon)

 %Every particle gets its own colour
    colours = jet(part_plot );
    myTitle = sprintf('Electron Trajectories in Silicon, Electric Field Present, E= %d V/nm', mean(Ex, "all"));
 
  %  This plots the linear trajectories of all the particles
    figure
    for m =1:part_plot 
    scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
    hold on
    end
    rectangle('Position',[Box{1}.x(1) Box{1}.y(1) (Box{1}.x(2)-Box{1}.x(1)) (Box{1}.y(2)-Box{1}.y(1)) ],'EdgeColor','b');
    hold on
    rectangle('Position',[Box{2}.x(1) Box{2}.y(1) (Box{2}.x(2)-Box{2}.x(1)) (Box{2}.y(2)-Box{2}.y(1)) ],'EdgeColor','b');
    hold off
    title(myTitle)
    axis([0 l 0 w])
    ylabel('y, (nm)')
    xlabel('x, (nm)')

%Potential and Current Plots

figure
% tiledlayout(2,1)
% ax1 = nexttile;
surf(Vmap');
shading interp;
colormap(cool(20));
title('Potential Map','FontSize', 12);
xlabel('x','FontSize', 12) 
ylabel('y','FontSize', 12) 
zlabel('Potential (V)','FontSize', 12) 

figure
% nexttile;
quiver(Ex', Ey');
hold on
rectangle('Position',[Box{1}.x(1) Box{1}.y(1) (Box{1}.x(2)-Box{1}.x(1)) (Box{1}.y(2)-Box{1}.y(1)) ],'EdgeColor','b');
hold on
rectangle('Position',[Box{2}.x(1) Box{2}.y(1) (Box{2}.x(2)-Box{2}.x(1)) (Box{2}.y(2)-Box{2}.y(1)) ],'EdgeColor','b');
title('Electric Field Vectors','FontSize', 12);
xlabel('x','FontSize', 12) 
ylabel('y','FontSize', 12) 
axis([0 l 0 w]);


