%% Problem 2: Maximum likelihood/Least squares
M = 1e4;  % Number of observations
b = 4;    % True value of the Rayleigh parameter
x = raylrnd(b, M, 1);  % Generate random Rayleigh-distributed data

% Plot histogram of the data using the updated hist_density function
hist_density(x, 40)  % Calls your function to plot a normalized histogram
hold on

% Maximum Likelihood Estimate (MLE) for b
my_est_ml = sqrt(sum(x.^2) / (2 * M));  

% Least Squares Estimate (LSE) for b
my_est_ls = (2 / (M * pi)) * sum(x);  

% Plot the ML and LS estimates as stars and the true b as a circle
plot(my_est_ml, 0, 'r*')  % Red star for MLE
plot(my_est_ls, 0, 'g*')  % Green star for LSE
plot(b, 0, 'ro')          % Red circle for the true value of b

plot(0:0.1:6, raylpdf(0:0.1:16, my_est_ml), 'r')
hold off
