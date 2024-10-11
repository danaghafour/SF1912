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