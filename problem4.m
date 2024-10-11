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
