% Down-and-Out Call Option Pricing Surface

% Parameters
K = 40;         % Strike price
B = 60;         % Barrier level
r = 0.05;       % Risk-free interest rate
sigma = 0.2;    % Volatility
T = 1;          % Time to maturity
Smin = 0;       % Minimum stock price
Smax = 200;     % Maximum stock price (for plotting)
tmin = 0;       % Start time
tmax = T;       % End time

% Define the grid
N = 100;                            % Number of grid points
S = linspace(Smin, Smax, N);        % Stock prices
t = linspace(tmin, tmax, N);        % Time points
[S, t] = meshgrid(S, t);            % Create 2D grid

% Initialize the option price matrix
Cdo = zeros(size(S));

% Compute down-and-out call option prices on the grid
for i = 1:N
    for j = 1:N
        tau = T - t(i, j);  % Time to maturity at each grid point
        if S(i,j) < B
            Cdo(i,j) = 0;   % Option is knocked out
        else
            Cv_bs = blsprice(S(i,j), K, r, tau, sigma);  % Vanilla call
            Cdo(i,j) = Cv_bs - (S(i,j)/B)^(-2*r/sigma^2) * ...
                       blsprice(B^2/S(i,j), K, r, tau, sigma);  % Down-and-out
        end
    end
end

% Plot the surface of the option price
surf(S, t, Cdo)
title('Down-and-Out Call Option Price')
xlabel('Stock Price')
ylabel('Time to Maturity')
zlabel('Option Price')
