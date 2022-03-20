
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
passageWidth = zeros(5,1);

passageLength = 40;
l = 200; % length
w = 100; % width
silicon = zeros(w, l);

ave_current = zeros(5,1);

for i = 1:5
    passageWidth(i) = i*10;
    ave_current(i) = boxes(passageWidth(i));
end


plot(passageWidth,ave_current);
title("Current in Silicon vs Passage Width")
ylabel('Currnet, (A)')
xlabel('width, (nm)')







