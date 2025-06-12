function [V, ster, CPUt, varsc, eb] = ouscratch_bb(N, M_, seed, delta, varargin)
% Improved Monte Carlo simulation with Brownian Bridge correction for barrier crossing

if seed ~= -1
    rng(seed);
end

pinta = false;
if nargin == 5
    pinta = true;
end

% Model and option parameters
[S0, T, K, B0, sigma, kappa, theta, r] = deal(14, 2, 14, 13.5, 0.5, 2, 14, 0);

% Generate Brownian increments
Wfino = randn(N, 2 * M_);

for k = [1, 2]
    tic;

    if k == 1
        dW = (Wfino(:, 1:2:end-1) + Wfino(:, 2:2:end)) / sqrt(2);
    else
        dW = Wfino;
    end

    M = k * M_;
    h = T / M;
    S = NaN(N, M + 1);
    S(:, 1) = S0;
    B = B0;

    % Initialize HIT status
    HIT = ones(N, 1);

    % Simulate paths and detect barrier hits using Brownian Bridge
    for j = 1:M
        mu = S(:, j) + kappa * (theta - S(:, j)) * h;
        sigma_sqrt_h = sigma * sqrt(h);

        % Simulate next step
        S(:, j + 1) = mu + sigma_sqrt_h .* dW(:, j);

        % Brownian bridge adjustment
        bridge = brownianBridge(N, h);
        bridge_min = min(S(:, j), S(:, j + 1)) + bridge;
        HIT(bridge_min <= B) = 0;
    end

    % Compute payoffs
    score = exp(-r * T) * max(0, S(:, end) - K);
    score(HIT == 0) = 0;

    % Collect stats
    V(k) = mean(score);
    varsc(k) = var(score);
    ster(k) = 3 * sqrt(varsc(k) / N);
    CPUt(k) = toc;
end

% Richardson extrapolation for bias estimation
eb = (V(2) - V(1)) / (1 - 2^delta);

if pinta
    for k = [1, 2]
        fprintf('N=%d, h=%.5f: V=%.5f ± %.5f, CPUt=%.3f s (%.2f%% hit barrier)\n', ...
            N, T/(k*M_), V(k), ster(k), CPUt(k), 100*(1 - sum(HIT)/N));
    end

    fprintf('Estimated bias = %.5g\n', eb);
    fprintf('abs(bias) = %.5g\n', abs(eb));
    fprintf('2 * max(σ error) = %.5g\n', 2 * max([abs(ster(1)), abs(ster(2))]));

    if abs(eb) < 2 * max([abs(ster(1)), abs(ster(2))])
        fprintf('WARNING: Statistical error NOT negligible compared to bias. Unreliable estimate!\n');
        eb = NaN;
    end
end
end
