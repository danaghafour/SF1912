%% Problem 1 - Simulation of confidence intervals

% Parameters:
n = 25;               % Number of measurements (sample size for each interval)
mu = 2;               % Expected value (mean) of the normal distribution
sigma = 1;            % Standard deviation of the normal distribution
alpha = 0.05;         % Significance level (used for 95% confidence interval)

% Simulation of n * 100 observations (n observations for each of 100 intervals)
x = normrnd(mu, sigma, n, 100);   % n x 100 matrix of observations, generated from N(mu, sigma)

% Estimation of mu by the mean of each sample
xbar = mean(x);       % Vector containing the means of 100 samples (1 mean for each sample of 25 observations)

% Computation of upper and lower limits for the confidence intervals
lowerl = xbar - norminv(1 - alpha/2) * sigma / sqrt(n);  % Lower limit of the confidence intervals
upperl = xbar + norminv(1 - alpha/2) * sigma / sqrt(n);  % Upper limit of the confidence intervals

% Plot all the intervals, marking those that do not cover the true value in red
figure(1)             % Create a new figure window
hold on               % Hold the plot to overlay multiple intervals
for k = 1:100
    if upperl(k) < mu          % If the upper limit is below the true value of mu
        plot([lowerl(k) upperl(k)], [k k], 'r')   % Plot the interval in red
    elseif lowerl(k) > mu      % If the lower limit is above the true value of mu
        plot([lowerl(k) upperl(k)], [k k], 'r')   % Plot the interval in red
    else
        plot([lowerl(k) upperl(k)], [k k], 'b')   % Otherwise, plot the interval in blue
    end
end

% Adjust the axis to improve the visual appearance of the plot
b1 = min(xbar - norminv(1 - alpha/2) * sigma / sqrt(n));  % Minimum value for the x-axis
b2 = max(xbar + norminv(1 - alpha/2) * sigma / sqrt(n));  % Maximum value for the x-axis
axis([b1 b2 0 101])   % Set the axis limits to minimize unused space in the figure

% Plot the true value of mu as a vertical green line
plot([mu mu], [0 101], 'g')

hold off              % Release the hold on the plot
