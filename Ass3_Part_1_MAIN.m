%Code for trajectories of electrons in silicon
%Sarah Dolan, 2022, ELEC 4700
set(0, 'DefaultFigureWindowStyle', 'docked')
close all

%Duration of simulation
num_steps = 500;

%Silicon Temperature
T = 300;

%Constants
C.m0 = 9.11 *10 ^ (-31);
C.mn = 0.26 * C.m0;
C.k = 1.381 * 10 ^ (-23);
C.q = 1.60217662 * 10 ^ (-19);

%Thermal velocity, mean time between collisions,
%Mean free path
v_Th = sqrt(2*C.k*T/C.mn)/10^-9; %Converted to nm/s
tau = 0.2 * 10 ^(-12);
lambda = 3.74 * 10 ^(-8);

%Time Step
delta_t = tau/150;

%Number of Particles
num_part = 30000;
part_plot = 10; % number of particles to be plotted
part_boundary = zeros(num_steps, 1); % number of particles that cross boundary

%Current 
current = zeros(num_steps, 1);

%Temperatures
temperatures = zeros(num_steps, 1);
temperature_map = T*ones(width_silicon, length_silicon);

%Silicon Dimensions
length_silicon = 200;
width_silicon = 100;
silicon = zeros(width_silicon, length_silicon);

%Potential, 1D
V0 = 0.1 / 10^-9; %Converted to V/nm
V_ai = Potential_1D(length_silicon, width_silicon,V0);

%Electric Field
[Ex, Ey] = gradient(V_ai);
%Force
F = Ey * C.q;


%particles are assigned random x and y positions
part.position = zeros(num_part, 2);
part.position(:,1) = length_silicon * rand (num_part, 1);
part.position(:,2) = width_silicon* rand (num_part, 1);

%every particle has a random angle 
part.phi = 2*pi* rand(num_part, 1); 

%particles are assigned x and y velocities
part.velocity = zeros(num_part, 2);
part.velocity(:,1) = cos(part.phi) * v_Th;
part.velocity(:,2) = sin(part.phi) * v_Th;


%particles are assigned x and y accelerations
acceleration = mean(F, "all")/(C.mn*(10^-9))

all_x_positions = zeros(num_part, num_steps);
all_y_positions = zeros(num_part, num_steps);
all_x_positions(:,1) = part.position(:,1);
all_y_positions(:,1) = part.position(:,2);

myTitle = sprintf('Electron Trajectories in Silicon, Electric Field Present, E= %d V/nm', mean(Ey, "all")*10^-9);

%------------------------------------------------------------------------------------------------------------------------------------------------
for i = 1:num_steps

%   This is the live plot
%     scatter(part.position(1:part_plot,1),part.position(1:part_plot,2), 'b');
%     title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)
  
    %Velocity Updates
    part.velocity(:,1) = part.velocity(:,1) + acceleration * delta_t;
    %Position Updates
    part.position = part.position + part.velocity * delta_t;
    
    %Checking Boundary Conditions
    for n = 1:num_part 
        if  part.position(n, 1) > length_silicon || part.position(n, 1) < 0
            if  part.position(n, 1) > length_silicon 
            part.position (n, 1) =  part.position (n, 1) - length_silicon; 
            part_boundary(i,1) = part_boundary(i,1) + 1; %Particle crosses positive boundary, increases current
            else 
            part.position(n, 1) = part.position(n, 1) + length_silicon; 
            part_boundary(i,1) = part_boundary(i,1) - 1; %Particle crosses negative boundary, decreases current
            end
        end
        if  part.position(n, 2) > width_silicon || part.position(n, 2) < 0 
            part.velocity(n,2) = -part.velocity(n,2);
            part.position (n,:) = part.position(n,:) + part.velocity(n,:) * delta_t;
        end
    end
       
    %Currents
    current(i) = part_boundary(i) *C.q / delta_t;

    % I have arrays for storing all of the positions so I can plot them all
    % at the end
    all_x_positions(:,i) = part.position(:,1);
    all_y_positions(:,i) = part.position(:,2);    

end
    %Every particle gets its own colour
    colours = jet(part_plot );
 
  %  This plots the linear trajectories of all the particles
    figure
    for m =1:part_plot 
    scatter(all_x_positions(m,:),all_y_positions(m,:),'.', 'color', colours(m))
    hold on
    end
    title(myTitle)
    axis([0 length_silicon 0 width_silicon])
    ylabel('y, (nm)')
    xlabel('x, (nm)')
% 
ElectronDensity (num_part, part, silicon, length_silicon, width_silicon)
TemperatureMap (num_part, part, temperature_map, C);

figure
surf(V_ai'*(10^-9));
title('Potential, Finite Difference Matrix Solution, 1D')
shading interp;
xlabel('x (nm)') 
ylabel('y (nm)') 
zlabel('Potential, (V)')

% 
% figure
% myTitle_Efield = sprintf('Electric Field Finite Difference Matrix Solution, 1D, E = %d V/nm', mean(Ey, 'all')*10^-9);
% surf(Ey, 'FaceColor','r');
% %surf(1 :width_silicon,1 :length_silicon, Ey','FaceColor','r');
% % zlim([9^6 10^7 ])
% xlabel('x (nm)') 
% ylabel('y (nm)') 
% zlabel('E-Field, (V/m)')
% title(myTitle_Efield)
% 
% figure
% myTitle_Force= sprintf('Electric Force, Finite Difference Matrix Solution, 1D, F = %d N', mean(F, 'all'));
% surf(1 :width_silicon,1 :length_silicon, F,'FaceColor','b');
% %zlim([ 0 10^-12])
% xlabel('x (nm)') 
% ylabel('y (nm)') 
% zlabel('E-Force, (N)')
% title(myTitle_Force)

figure
plot(delta_t * (1:num_steps), current)
title("Current in Silicon At a Given Time")
ylabel('Currnet, (A)')
xlabel('t, (s)')