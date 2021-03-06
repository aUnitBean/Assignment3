function TemperatureMap (num_part, part, temperature_map, C)

  for n = 1:num_part
        my_position_y = round(part.position(n,2));
        if my_position_y <= 0
            my_position_y = 1;
        end
        if my_position_y >100
            my_position_y = 100;
        end
        my_position_x = round(part.position(n,1));
        if my_position_x <= 0
            my_position_x = 1;
        end
        if my_position_x > 200
            my_position_x = 200;
        end
    v_tot_square = (part.velocity(n,1).^2 + part.velocity(n,2).^2);
    v_tot_square = v_tot_square *((10^-9)^2); % convert velocity to m/s
    KE_part = (1/2) * C.mn * v_tot_square; 
    temperature_map (my_position_y, my_position_x) = temperature_map (my_position_y, my_position_x) + KE_part / C.k;

  end

    figure
    surf(temperature_map);
    title('Temperature Map')
    ylabel('y (nm)')
    xlabel('x(nm)')
    zlabel('Temperature (K)')
end