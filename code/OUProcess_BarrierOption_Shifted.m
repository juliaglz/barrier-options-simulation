function [V, ster, CPUt, varsc, eb] = ouscratch_shifted(N, M_, semilla, delta, varargin)
% Monte Carlo simulation of a down-and-out option under OU process
% using barrier shifting (weak convergence improvement)

if semilla ~= -1
    rng(semilla);
end

pinta = false;
if nargin == 5
    pinta = true;
end

% Parameters
[S0, T, K, B0, sigma, kappa, theta, r] = deal(14, 2, 14, 13.5, 0.5, 2, 14, 0);
Wfino = randn(N, 2 * M_);

for k = [1, 2]
    tic;

    % Generate Brownian increments
    if k == 1
        dW = (Wfino(:, 1:2:end-1) + Wfino(:, 2:2:end)) / sqrt(2);
    else
        dW = Wfino;
    end

    M = k * M_;
    h = T / M;
    S = NaN(N, M + 1);
    S(:, 1) = S0;

    % Shifted barrier (for weak convergence improvement)
    B = B0 + 0.5826 * sigma * sqrt(h);

    % Simulate paths
    for j = 1:M
        S(:, j + 1) = S(:, j) + kappa * (theta - S(:, j)) * h + sigma * sqrt(h) .* dW(:, j);
    end

    % Detect barrier hits
    HIT = ones(N, 1);
    for j = 1:N
        if min(S(j, :)) <= B
            HIT(j) = 0;
        end
    end

    % Compute payoffs
    score = exp(-r * T) * max(0, S(:, end) - K);
    score(HIT == 0) = 0;

    % Statistics
    V(k) = mean(score);
    varsc(k) = var(score);
    ster(k) = 3 * sqrt(varsc(k) / N);
    CPUt(k) = toc;
end

% Richardson extrapolation for bias
eb = (V(2) - V(1)) / (1 - 2^delta);

if pinta
    for k = [1, 2]
        fprintf('N=%d, h=%.5f: V=%.5f ± %.5f, CPUt=%.3f s (%.2f%% hit barrier)\n', ...
            N, T/(k*M_), V(k), ster(k), CPUt(k), 100 * (1 - sum(HIT)/N));
    end

    fprintf('Estimated bias = %.5g\n', eb);
    fprintf('abs(eb) = %.5g\n', abs(eb));
    fprintf('2 * max(σ error) = %.5g\n', 2 * max([abs(ster(1)), abs(ster(2))]));

    if abs(eb) < 2 * max([abs(ster(1)), abs(ster(2))])
        fprintf('WARNING: Statistical error NOT negligible compared to bias. Unreliable estimate!\n');
        eb = NaN;
    end
end
end
