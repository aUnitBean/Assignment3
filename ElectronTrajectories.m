
function [all_x_positions, all_y_positions, part, aveCurrent] = ElectronTrajectories(E, Box, num_boxes, num_steps, width_silicon, length_silicon, num_part, C)

part_plot = 10;
%Silicon Temperature
T = 300;

%Thermal velocity, mean time between collisions,
%Mean free path
v_Th = sqrt(2*C.k*T/C.mn)/10^-9; %Converted to nm/s
tau = 0.2 * 10 ^(-12);

%Time Step
delta_t = tau/100;
delta_time_step = 0;

%Number of Particles
% part_plot = 10; % number of particles to be plotted
part_boundary = zeros(num_steps, 1); % number of particles that cross boundary

%Current 
current = zeros(num_steps, 1);

%Temperatures
% temperatures = zeros(num_steps, 1);
% temperature_map = T*ones(width_silicon, length_silicon);

%silicon = zeros(width_silicon, length_silicon);

%Force
F = E * C.q/10^-9;



%particles are assigned random x and y positions
part.position = zeros(num_part, 2);
part.position(:,1) = length_silicon * rand (num_part, 1);
part.position(:,2) = width_silicon* rand (num_part, 1);

% Remove particles from box region
for b = 1: num_boxes
    q = 1;
    while q <num_part
        q = q+1;
        if part.position(q,1) > Box{b}.x (1) && part.position(q,1) < Box{b}.x (2) && part.position(q,2) > Box{b}.y (1) && part.position(q,2) < Box{b}.y (2)
               part.position(q,1) = length_silicon * rand;
               part.position(q,2) = width_silicon * rand;
               q = 1;
        end
    end
end


%every particle has a random angle 
part.phi = 2*pi* rand(num_part, 1); 

%particles are assigned x and y velocities
part.velocity = zeros(num_part, 2);
part.velocity(:,1) = cos(part.phi) * v_Th;
part.velocity(:,2) = sin(part.phi) * v_Th;

%particles are assigned x accelerations
acceleration = -F/(C.mn*(10^-9));
part.aceleration = zeros(num_part,1);
part.acceleration = getAcceleration (part, acceleration, num_part, length_silicon, width_silicon);

all_x_positions = zeros(num_part, num_steps);
all_y_positions = zeros(num_part, num_steps);
all_x_positions(:,1) = part.position(:,1);
all_y_positions(:,1) = part.position(:,2);

%For collisions
num_collisions = 0;
part.collisions = zeros(num_part,1);
duration = num_steps * delta_t;

for i = 1:num_steps
%   This is the live plot
%     scatter(part.position(1:part_plot,1),part.position(1:part_plot,2), 'b');
%     hold on
%     rectangle('Position',[Box{1}.x(1) Box{1}.y(1) (Box{1}.x(2)-Box{1}.x(1)) (Box{1}.y(2)-Box{1}.y(1)) ],'EdgeColor','b');
%     hold on
%     rectangle('Position',[Box{2}.x(1) Box{2}.y(1) (Box{2}.x(2)-Box{2}.x(1)) (Box{2}.y(2)-Box{2}.y(1)) ],'EdgeColor','b');
%     hold off
%  %   title(myTitle)
%     ylabel('y, (nm)')
%     xlabel('x, (nm)')
%     axis([0 length_silicon 0 width_silicon])
%     pause(0.01)
  
    %Velocity Updates
    part.velocity(:,1) = part.velocity(:,1) + part.acceleration * delta_t;
    %Position Updates
    part.position = part.position + part.velocity * delta_t;
    
    
    for n = 1:num_part 
        %Sattering
        P_scat = 1-exp(-(delta_time_step*delta_t + delta_t)/tau);
        delta_time_step = delta_time_step + 1;
      
        if (P_scat > 10*rand())
            part.collisions (n) = part.collisions (n)+1;
            delta_time_step = 0; %reset time step
            %velocity reassigned
            v = sqrt(part.velocity(n,1)^2+part.velocity(n,2)^2);
            new_random_velocity =  v * randn(1,1) + v; 
            new_random_phi = 2*pi* rand(); 
            part.velocity(n,1) = cos(new_random_phi) * new_random_velocity;
            part.velocity(n,2) = sin(new_random_phi) * new_random_velocity;
        end
        %Checking Boundary Conditions
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

       for b=1:num_boxes
           while part.position(n,1) > Box{b}.x (1) && part.position(n,1) < Box{b}.x (2) && part.position(n,2) > Box{b}.y (1) && part.position(n,2) < Box{b}.y (2)
              part.position (n,:) = part.position(n,:) - part.velocity(n,:) * delta_t;
      
                   if part.position (n,1) < Box{b}.x (1) || part.position(n,1) > Box{b}.x (2)
                      part.velocity(n,1) =-part.velocity(n,1) ;
                   end
                   if part.position (n,2) < Box{b}.y (1) || part.position(n,2) > Box{b}.y (2)
                      part.velocity(n,2) =-part.velocity(n,2) ;
                   end  
           end
       end

    end
    % I have arrays for storing all of the positions so I can plot them all
    % at the end

    %Currents
    current(i) = part_boundary(i) *C.q / delta_t;

    all_x_positions(:,i) = part.position(:,1);
    all_y_positions(:,i) = part.position(:,2);  
    part.aceleration = getAcceleration (part, acceleration, num_part, length_silicon, width_silicon);
end

aveCurrent = mean(current);
      
end
