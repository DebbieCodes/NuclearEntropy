clear all;
close all;
% domain of the problem is x=[0,1], with sensors at positions x=0,1 and one
% random position in between. The event occurs within x=[0,1].

% The goal is to determine the location of the event.

% We know: 
% - the time the 3 sensors receive a signal.
% - the signal created propagates in both directions (x and -x)
% We do not know:
% -  the propagation velocity of the signal created
% 

% Method: 
% - create an ensemble of events randomly scatterred around the domain (random initial guess).
% - calculate what time these events should reach the sensors.
% - compare measured arrival time (ym) with predicted arrival time (y_pred) for each
% ensemble member
% - calculate kalman gain based on ensemble.
% - update [time, pos, vel] based on gain for each member
% - repeat untill the error: ym - y_pred is minimized.

%% Set parameters
makePlot = true;
%Simulation parameters
nev = 1; nsen = 3; % Please don't try this with >1 event, the code needs to be modified.
ndim = 2;
%Random event
zt= rand(nev,ndim); %true location and time of event, between 0 an 1
%sensor placement
s = [0 (rand-0.5)/4+0.5 1]; % sensor on each edge and one somewhere in between.
%velocity model
v = 0.5; % universal velocity for all events (only 1 now)
% random noise -- 1\%
sig = 0.01;
%arrival
y=abs(zt(:,1)-s)/v+zt(:,2)+randn(nev,nsen)*sig; % arrival time at each sensor
disp(y)

xvec = linspace(0,1,100);

% plot locations sensors and true event
figure;
hold on
plot(s,zeros(3),'kx')
plot(zt(1),zt(2),'r*')
plot(xvec,abs(xvec-zt(1))/v+zt(2),':')
xline(s,'--')
xlim([-0.1 1.1])
ylim([0 1.1])
title('sensor locations and true event')
xlabel('x')
ylabel('t')

%% Ensemble method to find event
% Estimation 
% estimate number of ensemble members to use
nens = nsen*nev*ndim/sqrt(sig); % This is a kludge

z = [rand(nev,ndim,nens)]; % create ensemble of events randomly scattered around the domain in time and position.
vens=rand(1,nens); % also give them all their own random velocity, as we assume the velocity is unknown.
esprev = 0; %  initialize variance over the ensemble
des = inf; % initialize error, delta_e sample
jj = 1;

if(makePlot==true)
    figure;
    hold on
    plot(s,zeros(3),'kx')
    plot(zt(1),zt(2),'r*')
    xline(s,'--')
    xlim([-0.1 1.1])
    ylim([0 1.1])
    title('true event and ensemble')
    xlabel('x')
    ylabel('t')
    aa = axis;
end
    

while(jj<100 && max(abs(des))>1e-4) % need better convergence criterion. Now: loop for maximum 100, and while error is higher than 1e-4
   if(makePlot==true)
         z_pos = reshape(z(1,1,:),1,[]);
         z_time = reshape(z(1,2,:),1,[]);
         plot(z_pos, z_time,'.b')
         axis(aa)
         pause(1)
         drawnow
    end
    Z = [reshape(z,[nev*ndim nens]);vens]; % 3x60 matrix, for each ensemble member have: [pos, time, velocity]
    dZ = Z - mean(Z,2); % subtract ensemble mean
    
    Yh = ttModel(z,vens,s); % Arrival time at each sensor for each ensemble member using a model for propagation ( so size is 3xnens)
    dYh = Yh - mean(Yh,2); % subtract ensemble mean
    
    Czy = (dZ*dYh')/(nens -1 ); % Cross variance dZ and dY (3x3 matrix)
    Cyy  = cov(Yh');  %(3x3 matrix)
    
    G = Czy*pinv(Cyy + sig^2*eye(nsen)); %Kalman Gain, inefficient approach. (3x3 matrix)
    errv = (y' - Yh); % y wasn't perturbed. Difference between measured arrival time and predicted arrival time of each ensemble member
    zup = G*errv; % update z. 3x60 matrix. delta step to update [pos, time, velocity] for each ensemble member
    
    % use continued variance reduction as criterion
    % If all the ensemble members converge to one location (the answer we
    % are looking for), then the variance over the ensemble members will decrease
    es = var(zup,[],2);  %(1x3) matrix, variance for [pos, time, velocity]
    des = es - esprev; % delta variance:
    esprev = es; % store the variance of this step so we can compare with the next.
    
    %update
    Z = Z + zup; % update [pos, time, velocity] for each ensemble member
    % store [pos, time] and % vel separately
    z = reshape(Z(1:2,:),nev,ndim,nens); % get first to entries from Z [pos, time] and store as z
    vens = Z(3,:); % write velocity to separate vector
    %iterate
    jj = jj+1;
end
%answer
[mean(Z,2) [zt';v]]
disp(jj)



function y = ttModel(z,vens,s) % gives a (3xnens) output
for i = 1:size(z,3)
    y(:,i) = abs(z(:,1,i)-s)./vens(i)+z(:,2,i);
    % y = abs(x - s)/v + t
    %% for the other approach
    % [sgn(x-s)/v        0        0; ...
    %  0          -abs(x-s)/v^2   0;...
    %  0                 0        1]
end
end