  

function part_acceleration =  getAcceleration (part, acceleration, num_part, length_silicon, width_silicon)

part_acceleration = zeros (num_part, 1);
for n =1:num_part
            my_position_y = round(part.position(n,2));
            if my_position_y <= 0
                my_position_y = 1;
            end
            if my_position_y >width_silicon
                my_position_y = width_silicon;
            end
            my_position_x = round(part.position(n,1));
            if my_position_x <= 0
               my_position_x = 1;
            end
            if my_position_x >= length_silicon
               my_position_x = length_silicon;
            end
          part_acceleration(n) = acceleration(my_position_x, my_position_y);
end
end
