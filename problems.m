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

%% Problem 4 - Comparison of distributions in different populations
% loading data about 747 women who gave birth to their first 
% child in Malmö during 1991-1993.

load birth.dat

childWeight = birth(:, 3);
motherAge = birth(:, 4);
motherWeight = birth(:, 15);
motherHeight = birth(:, 16);

figure;

% Histogram for birth weight
subplot(2, 2, 1);
hist_density(childWeight, 40);
title('Birth weight distributions');
xlabel('Weight (g)');
ylabel('Density');

% Histogram for mother's age
subplot(2, 2, 2);
hist_density(motherAge, 40);
title('Mother Age Distribution');
xlabel('Age (years)');
ylabel('Density');

% Histogram for mother's weight
subplot(2, 2, 3);
hist_density(motherWeight, 40);
title('Mother Weight Distribution');
xlabel('Age (years)');
ylabel('Density');

% Histogram for mother's height
subplot(2, 2, 4);
hist_density(motherHeight, 40);
title('Mother Height Distribution');
xlabel('Age (years)');
ylabel('Density');

%% Problem 4: Distributions of given data
load birth.dat
% Icke-rökande eller slutat
x = birth(birth(:, 20) < 3, 3);
% Rökande
y = birth(birth(:, 20) == 3, 3);

subplot(2,2,1), boxplot(x),
axis([0 2 500 5000])
title('not smoking during pregnancy');
subplot(2,2,2), boxplot(y),
title('smoking during pregnancy');
axis([0 2 500 5000])

subplot(2,2,3:4), ksdensity(x),
hold on
[fy, ty] = ksdensity(y);
title('probability density function (PDF)');
plot(ty, fy, 'r')
hold off


%% Problem 4.2: Distributions of given data 
load birth.dat
% Svensk
x = birth(birth(:, 7) ==1 , 3);
% Icke-svensk
y = birth(birth(:, 7) == 0, 3);

subplot(2,2,1), boxplot(x),
axis([0 2 500 5000])
title('Icke-svensk');
subplot(2,2,2), boxplot(y),
title('Svensk');
axis([0 2 500 5000])

subplot(2,2,3:4), ksdensity(x),
hold on
[fy, ty] = ksdensity(y);
title('probability density function (PDF)');
plot(ty, fy, 'r')
hold off

%% Problem 5 - Testing for normality

% Load the data
load birth.dat

% Extract variables
childWeight = birth(:, 3);
motherAge = birth(:, 4);
motherWeight = birth(:, 15);
motherHeight = birth(:, 16);

% Create a figure with subplots
figure;

% Normal Probability Plot for Birth Weight
subplot(2, 2, 1);
normplot(childWeight);
title('Normal Probability Plot - Birth Weight');

% Normal Probability Plot for Mother\'s Age
subplot(2, 2, 2);
normplot(motherAge);
title('Normal Probability Plot - Mothers Age');

% Normal Probability Plot for Mother\'s Weight
subplot(2, 2, 3);
normplot(motherWeight);
title('Normal Probability Plot - Mothers Weight');

% Normal Probability Plot for Mother\'s Height
subplot(2, 2, 4);
normplot(motherHeight);
title('Normal Probability Plot - Mothers Height');


% JB-test for Birth Weight
motherAgeTest= jbtest(motherAge,0.05);

% JB-test for Birth Weight
motherHeightTest= jbtest(motherHeight,0.05);

% JB-test for Birth Weight
motherWeightTest= jbtest(motherWeight,0.05);

% JB-test for Birth Weight
childWeightTest= jbtest(childWeight,0.05);


% Store the variable names and JB test results in cell arrays
variables = {'Mothers Age:', 'Mothers Height:', 'Mothers Weight:', 'Childs Weight:'};
jb_results = [motherAgeTest, motherHeightTest, motherWeightTest, childWeightTest];

% Loop through each variable and print the results
for i = 1:length(variables)
    if jb_results(i) == 1
        fprintf('\n%s Given the significance level 0.05 and jbtest outputs 1, \nThe test has rejected the null hypothesis. Therefore not normally distributed., \n', variables{i});
    else
        fprintf('\n%s Given the significance level 0.05 and jbtest outputs 0, \nThe test has accepted the null hypothesis. Therefore normally distributed.\n', variables{i});
    end
end

%% Problem 6 - Simple linear regression
load moore.dat
x= moore(:,1); % the year
y=moore(:,2); % number of trasistors
w = log(y); % the logarithm of the number of transistors/unit 
% area should increase linearly over time.

% Matrix with ones in the columne on left side, and data values on the right side
X = [ones(size(x)) x];
[beta_hat,bint,r,rint,stats] = regress(w,X); 
% the logarithm of the numbers of transistors increase linearly

% Call the estimate w_circumflex 
w_circumflex = X*beta_hat ;

figure
plot(x,w)
hold on
plot(x, w_circumflex)

% Comparing w to w_circumflex to the data y 
res= w_circumflex-w;
disp(res)
subplot(2,1,1), normplot(res)
subplot(2,1,2), hist(res)

%R² value from regress function
r2 = stats(1) 

% prediction of the number of transistors for year 2025 β = (β0, β1)T 
numOfTrans2025 = exp(beta_hat(1)+beta_hat(2)*2025);
fprintf("number of transistor in 2025: %d \n", numOfTrans2025);
fprintf("R²: %d \n", r2);
 