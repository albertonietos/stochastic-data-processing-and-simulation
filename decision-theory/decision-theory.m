%% Lab 2, created by Alberto Nieto Sandino

%% Assignment 1.1 - Assumptions about the data
clear all; close all; clc

% Read the data with stock values
load('stockdata.tsv');

% Calculate the log-returns
X = zeros(size(stockdata,1)-1,size(stockdata,2)-1);
for i = 2:size(stockdata,1) 
    X(i,:) = log(stockdata(i,2:8)) - log(stockdata(i-1,2:8));
end

% Create cell array for stock names
stocks = {'AstraZeneca','Electrolux','Ericsson','Gambio','Nokia',...
          'Swedish Match','Svenska Handelsbanken'};

% Plot log of stock values over time
figure();
for j = 1:7
    subplot(4,2,j)
    plot(stockdata(:,1), X(:,j));
    xlim([min(stockdata(:,1)) max(stockdata(:,1))])
    title(stocks(j));
    xlabel('time (days)');
    ylabel('X(t)');
end

% Plot the histograms and qq-plots for every stock value
figure(); % Histograms
for k = 1:7
    subplot(4,2,k)
    histogram(X(:,k));
    title(stocks(k));
end

figure()
for k = 1:7
    subplot(4,2,k)
    qqplot(X(:,k));
    title(stocks(k));
end
clear i j k

% Formal goodness of fit for each stock
h = zeros(7,1);
p = zeros(7,1);
for i = 1:7
    [h(i), p(i), stats(i)] = chi2gof(X(:,i),'nbins',20);
end
   
% Normal assumptions for the log-returns do not seem possible since the
% goodness of fit rejects the null hypothesis.

% Estimate the autocorrelation function of log-returns and its absolute
% values
figure()
for j = 1:7
    subplot(4,2,j)
    autocorr(X(:,j));
    title(stocks(j));
%     ylabel('r_h');
end

figure()
for j = 1:7
    subplot(4,2,j)
    autocorr(abs(X(:,j)));
    title(stocks(j));
    ylabel('r_h');
end
% The sample ACF has significant autocorrelation for some of the lags in
% all the stocks.

% Estimate the mean and standard deviation of the log-returns
mu = zeros(7,1); 
sd = zeros(7,1);
for k = 1:7
    mu(k) = mean (X(:,k));
    sd(k) = std (X(:,k));
end

% Estimate the correlation between the log-returns
clear i j k
r_sign = 0.8; % Minimum value of the correlation coeficient to be 
% considered significant
R = corrcoef(X);
for i = 1:7
    for j = 1:7
        if i~=j
            r_temp = R(i,j);
            if r_temp > r_sign
                disp (['Correlation between stocks ' num2str(i)...
                    ' and ' num2str(j) ' is likely,'...
                    'correlations coefficients are significant']);
            else
                disp (['Correlation between stocks is ' num2str(i)...
                    ' and ' num2str(j) ' unlikely,'...
                    'correlations coefficients are not significant']);
            end
        end
    end
end

% All the correlation coeficients are not significant, thus we can conclude 
% it is unlikely that there is a linear dependency between the stock
% values.
% It's expected that they're not correlated since they change based on
% different basis depending on how each company is doing.


%% Assignment 1.2 - Utility and expected utility

close all;

% Explore the utility function for different k's
k = 0.00001; 
figure()
for i = 1:6
    u = @(x)(1 - exp(-k*x));
    subplot(2,3,i);
    fplot(u,'linewidth',1.5);
    title(['k = ' num2str(k)]);
    k = k*10;
end

% The smaller the k, the more neutral the behaviour. The bigger the k, the
% more risk averse it is.

clc;
w = [1; zeros(6,1)]; % Considering 1st stock at a time
mu = mean (X).'; % Mean vector of log-returns
cov_mat = cov(X); % Covariance matrix of log-returns

z = mu.'*w; % Mean of 1st stock
sigma = sqrt(w.'*cov_mat*w); % Standard deviation for 1st stock

clear x
% Expected utility for stock nbr 1
% x = X(:,1); % utility
k = 1;
u = @(x)(1 - exp(-k*x));
pi_x = @(x)(1/(sqrt(2*pi)*sigma)).*exp(-((x-z).^2/(2*sigma.^2)));
y = @(x)((1 - exp(-k*x)).*(1/(sqrt(2*pi)*sigma)).*exp(-((x-z).^2/(2*sigma.^2))));

figure();
subplot(2,2,1:2)
    fplot(pi_x,[-1 1],'g');
    title('\pi_x(y)');
subplot(2,2,3)
    fplot(u,[-1 1],'b');
    title('u(y)');
subplot(2,2,4)
    fplot(y,[-1 1], 'r');
    title('\pi_x(y)*u(y)');

U = zeros(7,5);
for i = 1:7 % Loop thru stocks
    k = 0.0001;
    w = zeros(7,1);
    w(i) = 1; % Select each stock
    z = mu.'*w; % Mean of each stock
    sigma = sqrt(w.'*cov_mat*w); % Standard deviation for each stock
    for j = 1:5 % Loop thru k's
        y = @(x) u_AND_pi (x, k, z, sigma);
        U(i,j) = quad(y, min(X(:)), max(X(:)));
        k = k*10;
    end
end

figure();
mesh(U);
xlabel('k values');
ylabel('Stocks');
zlabel('Expected utility, U(w)');

figure()
k = 0.0001;
for i = 1:5
    subplot(3,2,i)
    stem(U(:,i))
    title(['k = ' num2str(k)]);
    k = k*10;
    xlabel('Stocks');ylabel('U(w)');
end

% For the most neutral behaviour, the best choice would be to invest in
% stock #4. As the behaviour gets more averse to risk #4 keeps being the
% best option for k = {0.0001, 0.01, 0.1}, except for k = 0.001 were #3 is
% better. For the most risk averse case (k=1), the most recommended would
% be #6 with #4 following closely.

% k = 0.01 gives the highest results considering all stocks.

%% Assignment 3 - Optimization of two stocks

clear all; close all; clc

% Read the data with stock values
load('stockdata.tsv');

% Calculate the log-returns
X = zeros(size(stockdata,1)-1,size(stockdata,2)-1);
for i = 2:size(stockdata,1) 
    X(i,:) = log(stockdata(i,2:8)) - log(stockdata(i-1,2:8));
end

% Take out the stock values of Ericsson (#3) and Gambio (#4)
X = X(:,3:4);

% Estimate the mean vector and the covariance matrix
mu = mean (X).'; % Mean vector of log-returns
cov_mat = cov(X); % Covariance matrix of log-returns

% Calculate the expected utility
k = 1;
w_1 = linspace(0, 1);
w_2 = 1 - w_1;
w = [w_1; w_2]; % Matrix of every pair of weights
U = zeros(length(w), 1);

for i = 1:length(w)
    z = mu.'*w(:,i); % Mean of i-th row
    sigma = sqrt(w(:,i).'*cov_mat*w(:,i)); % Standard deviation for i-th row
    y = @(x) u_AND_pi (x, k, z, sigma);
    U(i) = quad(y, min(X(:)), max(X(:)));
end

% Plot U(w) against w_1
figure()
% plot(w_1, U)
max_U = find(U == max(U));
plot(w_1, U,'-p','MarkerIndices',max_U, 'MarkerFaceColor','red',...
    'MarkerSize',15)
xlabel('Stock weights: left = 100% Gambio, right = 100% Ericsson');
ylabel('Expected utility, U(w)');


% Find the optimum point more accurately with fmincon
k = 1; % we keep k = 1, for now
clear w

% x0 = [0.15; 0.85]; 
x0 = [0.2; 0.8];% Starting point in approx. solution
A = [];
b = [];
Aeq = [1,1];
beq = 1;
lb = [0, 0];
ub = [1, 1];

% k = 0.00001;
% The solver fmincon will try to minimize the function, thus we try to
% minimize the opposite of the Eq. 2.2 which we want to maximize


fun = @(w) maximizefun(w, k, mu, cov_mat);
options.StepTolerance = 1e-20;
optimum_w = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);

% Evaluate the expected utility with the optimal weights
z = mu.'*optimum_w; % Mean of i-th row
sigma = sqrt(optimum_w.'*cov_mat*optimum_w); % Standard deviation for i-th row
y = @(x) u_AND_pi (x, k, z, sigma);

% When using the whole function it seems the mean is too powerful and seems
% to be where most of the influence lies.

% Using fmincon on function fun tries to maximize the expected return while
% minimizing the risk (standard deviation) of the investment. The risk can
% be decreased with only a small sacrifice of the return and that is what
% the algorithm does. It gives prevalence to the stock with higher expected
% return while diversificating, thus reducing the risk.

% The bigger the k, the more importance is given to the minimization of the
% risk (standard deviation), thus the investment should be more diverse
% with more equal weights.


% Find the optimal weights for several different k's
clear i k
k = linspace (0.00001, 1);
optima_w = zeros(2,100);
for i = 1:length(k)
    fun = @(w) maximizefun(w, k(i), mu, cov_mat);
    optima_w(:,i) = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
end

close all
figure()
yyaxis left
% title('Evolution of weights for different values of k');
plot(k, optima_w(1,:))
ylabel('Values for w_1')

yyaxis right
plot(k, optima_w(2,:))
xlabel('Values for k');
ylabel('Values for w_2')

figure()
risk = zeros(100,1);
for i = 1: size(optima_w, 2)
    risk (i) = optima_w(:,i)' * cov_mat * optima_w(:,i);
end
plot(k,risk)
title('Risk over k values')
xlabel('Values for k');
ylabel('Risk');

% Plotting the risk over the k values, we can see that as k increases, the
% risk aversion of the investment increases and the risk decreases as a
% consequence of our decision.

%% Assignment 4 - Optimization with 7 stocks

clear all; close all; clc

% Read the data with stock values
load('stockdata.tsv');

% Calculate the log-returns
X = zeros(size(stockdata,1)-1,size(stockdata,2)-1);
for i = 2:size(stockdata,1) 
    X(i,:) = log(stockdata(i,2:8)) - log(stockdata(i-1,2:8));
end

% Estimate the mean vector and the covariance matrix
mu = mean (X).'; % Mean vector of log-returns
cov_mat = cov(X); % Covariance matrix of log-returns

% Find the optimum point more accurately with fmincon

x0 = rand(7,1);
x0 = x0/sum(x0);
A = [];
b = [];
Aeq = ones(1,7);
beq = 1;
lb = zeros(1,7);
ub = ones(1,7);

% Define the function to maximize
fun = @(w) maximizefun(w, k, mu, cov_mat);

% Find the optimal weights for several different k's
clear i k
k = linspace (0.00001, 1);
optima_w = zeros(7,100);
fval = zeros(1,100);
for i = 1:length(k)
    fun = @(w) maximizefun(w, k(i), mu, cov_mat);
    [optima_w(:,i),fval(i)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
end

% Plot the final results of the weihgts
figure()
title('Evolution of weights for different values of k');
plot(k, optima_w)
ylabel('Values for the weights')
xlabel('Values for k');
legend('AstraZeneca','Electrolux','Ericsson','Gambio','Nokia','Swedish Match','Svenska Handelsbanken');

% Calculate the expected utility
for i = 1:length(optima_w)
    z = mu.'*optima_w(:,i); % Mean of i-th row
    sigma = sqrt(optima_w(:,i).'*cov_mat*optima_w(:,i)); % Standard deviation for i-th row
    y = @(x) u_AND_pi (x, k(i), z, sigma);
    U(i) = quad(y, min(X(:)), max(X(:)));
end

% Calculate the expected utility for a naive case
naive_w = (1/7)*ones(7,100);
for i = 1:length(naive_w)
    z = mu.'*naive_w(:,i); % Mean of i-th row
    sigma = sqrt(naive_w(:,i).'*cov_mat*naive_w(:,i)); % Standard deviation for i-th row
    y = @(x) u_AND_pi (x, k(i), z, sigma);
%     U(i) = quad(y, min(X(:)), max(X(:)));
    U_naive(i) = integral(y, -inf, inf);
end

% Calculate the expected utility for a all-in Gambio
allin_w = zeros(7,100);
allin_w(4,:) = 1;
for i = 1:length(allin_w)
    z = mu.'*allin_w(:,i); % Mean of i-th row
    sigma = sqrt(allin_w(:,i).'*cov_mat*allin_w(:,i)); % Standard deviation for i-th row
    y = @(x) u_AND_pi (x, k(i), z, sigma);
    U(i) = quad(y, min(X(:)), max(X(:)));
end

figure()
plot(k, U)
hold on
plot(k,U_naive)
plot(k,U_allin)
legend('Optimal','Naive','All Gambio')
xlabel('k values')
ylabel('U(w)')