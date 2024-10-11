%% Problem 3: Confidence interval for Rayleigh distribution
load wave_data.mat
% Plot the signal
subplot(211), plot(y(1:end))
title('Wave Data Signal')

% Plot histogram of the signal
subplot(212), hist_density(y, 40)
title('Histogram of Wave Data')


% Number of data points in the wave_data
n = length(y);
% Maximum Likelihood Estimate (MLE) a method of estimating the parameters of an assumed probability distribution, given some observed data
my_est = sqrt(sum(y.^2) / (2 * n));
% Least Squares Estimate (LSE) 
my_est_ls = (2 / (n * pi)) * sum(y);

% Confidence interval calculation (95% confidence level)
alpha = 0.05;  % Confidence level (5% significance level)
z_alpha_2 = norminv(1 - alpha/2);  % z-value for 95% confidence interval

% Standard error of the estimate
std_error = std(y) / sqrt(n); 
upper_bound = my_est - norminv(1 - alpha/2) * std_error;  % Lower limit of the confidence intervals
lower_bound = my_est + norminv(1 - alpha/2) * std_error;  % Upper limit of the confidence intervals

% Plot the confidence interval and estimates
hold on
plot(lower_bound, 0, 'g*')  % Green star for lower bound
plot(upper_bound, 0, 'g*')  % Green star for upper bound
plot(my_est, 0, 'r*')       % Red star for MLE estimate
plot(0:0.1:6, raylpdf(0:0.1:6, my_est), 'r')
hold off

title('MLE and Confidence Interval for Rayleigh Parameter')