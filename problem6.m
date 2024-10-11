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
 