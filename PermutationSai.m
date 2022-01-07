clear all
close all

% event position and time
nev = 3;
nsen= 3;
z = rand(nev,2);
% sensors
s = [0; 0.5; 1];

%vel
v = 0.5;

% arrival
y = abs(z(:,1)-s)/v+z(:,2)+ randn(3)/1000;
disp(z)
disp(y)

% permute
