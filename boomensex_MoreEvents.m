clear all;
close all;

% domain of the problem is x=[0,1], with sensors at positions x=0,1 and one
% random position in between. 2 (or more?) events occurs within x=[0,1].

% The goal is to determine the location of the events.

% We know: 
% - the time the 3 sensors receive a signal (so also the number of events)
% - the signals created propagate in both directions (x and -x)
% 
% We do not know:
% -  the propagation velocity of the signal created


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
nev = 2; nsen = 3; % Please don't try this with >1 event, the code needs to be modified.
ndim = 2;
%Random event
zt= rand(nev,ndim); %true location and time of event, between 0 an 1
%sensor placement
s = [0 (rand-0.5)/4+0.5 1]; % sensor on each edge and one somewhere in between.
%velocity model
v = 0.5; % universal velocity for all events 
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
for(j=1:nev)
    plot(zt(j,1),zt(j,2),'r*')
    plot(xvec,abs(xvec-zt(j,1))/v+zt(j,2),':')
end
xline(s,'--')
xlim([-0.1 1.1])
ylim([0 1.1])
title('sensor locations and true event')
xlabel('x')
ylabel('t')

%% Ensemble method to find events
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
    for(j=1:nev)
        plot(zt(j,1),zt(j,2),'r*')
        plot(xvec,abs(xvec-zt(j,1))/v+zt(j,2),':')
    end
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
    Z = [reshape(z,[nev*ndim nens]);vens]; % (3*n_events +1) x n_ensmembers matrix, for each ensemble member have: [pos1, time1,pos2, time2,pos3, time3,vel]
    dZ = Z - mean(Z,2); % subtract ensemble mean
    
    Yh1 = ttModel(z,vens,s); % Arrival time at each sensor for each ensemble member using a model for propagation, for each event:
    %Yh = [n_sensors,n_events,n_ensemble], 
    % reshape to: [n_sensors*nevents,n_ensmble], so we can calculate the
    % variance and do transpose over 2d matrix.
    Yh = reshape(Yh1,[nev*nsen nens]);
    dYh = Yh - mean(Yh,2); % subtract ensemble mean
    
    Czy = (dZ*dYh')/(nens -1 ); % Cross variance dZ and dY 
    Cyy  = cov(Yh');  %
    
    G = Czy*pinv(Cyy + sig^2*eye(nsen*nev)); %Kalman Gain, inefficient approach. 
    y_resh = reshape(y,[1,nev*nsen]);
    errv = (y_resh' - Yh); % y wasn't perturbed. Difference between measured arrival time and predicted arrival time of each ensemble member
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
    z = reshape(Z(1:(nev*ndim),:),nev,ndim,nens); % get first to entries from Z [pos, time] and store as z
    vens = Z(nev*ndim+1,:); % write velocity to separate vector
    %iterate
    jj = jj+1;
end

%% Display results
for j = 1:nev
    [[mean(Z(2*j:2*j+1,:),2);mean(Z(end,:),2)] [zt(j,:)';v]]
end
%answer event 2:
disp(jj)



function y = ttModel(z,vens,s) % gives a n_sens x n_events x  n_ens output, but put into same shape as Z:
    n_events = size(z,1);
    n_ens = size(z,3);
    n_sens = size(s,2);
    y = zeros(n_sens,n_events,n_ens);
    for i = 1:n_events % loop over each ensemble member
        for ii = 1:n_ens
            y(:,i,ii) = abs(z(i,1,ii)-s)./vens(ii)+z(i,2,ii);
            % y = abs(x - s)/v + t
            %% for the other approach
            % [sgn(x-s)/v        0        0; ...
            %  0          -abs(x-s)/v^2   0;...
            %  0                 0        1]
        end
    end
end