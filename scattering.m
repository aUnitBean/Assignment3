function [velocity_x, velocity_y] = scattering (part, delta_time_step, delta_t, tau, v_Th)

P_scat = 1-exp(-(delta_time_step*delta_t+delta_t)/tau);
delta_time_step = delta_time_step + 1;
duration = num_steps * delta_t;
  
    if (P_scat > rand())
        part.collisions (n) = part.collisions (n)+1;
        delta_time_step = 0; %reset time step
        %velocity reassigned
        new_random_velocity = sqrt(part.velocity(n,1)^2 +part.velocity(n,2)^2* randn(1,1)) + v_Th; 
        new_random_phi = 2*pi* rand(); 
        velocity_x = cos(new_random_phi) * new_random_velocity;
        velocity_y = sin(new_random_phi) * new_random_velocity;
    end
  
end