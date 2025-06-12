% Verify Down-and-Out Call Option Pricing

% Parameters
S = 50;         % underlying price
K = 50;         % strike price
r = 0.035;      % risk-free rate
T = 1;          % time to maturity
sigma = 0.3;    % volatility
B = 45;         % barrier price

% Compute the Black-Scholes price of a vanilla call option
Cv_bs = blsprice(S, K, r, T, sigma);

% Compute the price of a down-and-out call option using the analytical formula
Cdo_theory = Cv_bs - (S/B)^(-2*r/sigma^2) * blsprice(B^2/S, K, r, T, sigma);

% Compute the price of a down-and-out call option using the barrier option pricing function
Settle = datetime(2015,1,1);
Maturity = datetime(2016,1,1);
RateSpec = intenvset('ValuationDate', Settle, 'StartDates', Settle, 'EndDates', Maturity, ...
                     'Rates', r, 'Compounding', -1, 'Basis', 1);
StockSpec = stockspec(sigma, S);

Cdo = barrierbybls(RateSpec, StockSpec, 'call', K, Settle, Maturity, 'DO', B);

% Display the results
fprintf('Price of vanilla call option               : %.4f\n', Cv_bs);
fprintf('Price of down-and-out call option (formula): %.4f\n', Cdo_theory);
fprintf('Price of down-and-out call option (function): %.4f\n', Cdo);
