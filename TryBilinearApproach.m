% Generate Truth
n = 3; % number of sensors
s = [0; 0.5; 1];
x_t = [rand(2, 1); 1]; % x = [x t 1]
h_t = [rand(1); 1]; % h = [w 1]
nu = rand(n, 1) / 1000; % noise

X_t = [abs(s-x_t(1)) repmat(x_t(2), size(s))]

H_t = [repmat([-h_t(1) 1], size(s)) s*h_t(1)];
y_t = H_t * x_t + nu;

% Solve
h_i = [rand(1); 1]; % initial guess

err = 1;
while err > 1e-3
    H_i = [repmat([-h_i(1) 1], size(s)) s*h_i(1)];
    x_i = pinv(H_i'*H_i)*H_i'*y_t;
    
    X_i = [s - x_i(1) repmat(x_i(2), size(s))];
    h_i = pinv(X_i'*X_i)*X_i'*y_t;
    
    y = [repmat([-h_i(1) 1], size(s)) s*h_i(1)]*x_i;
    err = sqrt(mean((y_t - y).^2));
end

