 
function ElectronDensity (num_part, part, silicon, length_silicon, width_silicon)
%electron density
    for m =1:num_part      
        my_position_y = round(part.position(m,2));
        if my_position_y <= 0
            my_position_y = 1;
        end
        if my_position_y >width_silicon
            my_position_y = width_silicon;
        end
        my_position_x = round(part.position(m,1));
        if my_position_x <= 0
           my_position_x = 1;
        end
        if my_position_x >= length_silicon
           my_position_x = length_silicon;
        end
        silicon(my_position_y, my_position_x)  = silicon(my_position_y, my_position_x) + 1;
    end

    figure
    electron_density = silicon;
     surf(electron_density);
     title('Electron Density')
    ylabel('y (nm)')
    xlabel('x(nm)')
    zlabel('Electrons per Nano Metre')