function [V, ster, CPUt, varsc, eb] = ouscratch(N, M_, semilla, delta, varargin)
% [V, ster, CPUt, varsc, eb] = ouscratch(1e5, 128, -1, 0.5, true);

if semilla ~= -1
    rng(semilla);
end

pinta = false;
if nargin == 5
    pinta = true;
end

% Model and option parameters
[S0, T, K, B0, sigma, kappa, theta, r] = deal(14, 2, 14, 13.5, 0.5, 2, 14, 0);

% Generate Brownian increments
Wfino = randn(N, 2*M_);

for k = [1, 2]
    tic;

    if k == 1
        dW = (Wfino(:,1:2:end-1) + Wfino(:,2:2:end)) / sqrt(2);  % Brownian coarsening
    else
        dW = Wfino;  % Full resolution
    end

    M = k * M_;
    h = T / M;
    S = NaN(N, M+1);
    S(:,1) = S0;

    B = B0; % Use fixed barrier (no shifting)

    % Simulate all trajectories
    for j = 1:M
        S(:, j+1) = S(:, j) + kappa * (theta - S(:, j)) * h + sigma * sqrt(h) * dW(:, j);
    end

    % Detect paths that hit the barrier
    HIT = ones(N, 1);
    for j = 1:N
        if min(S(j, :)) <= B
            HIT(j) = 0;
        end
    end

    % Compute payoffs
    score = exp(-r * T) * max(0, S(:, end) - K);
    score(HIT == 0) = 0;  % Knocked-out paths get zero payoff

    % Store results
    V(k) = mean(score);                    % Option price
    varsc(k) = var(score);                 % Variance
    ster(k) = 3 * sqrt(varsc(k) / N);      % 3-sigma confidence interval
    CPUt(k) = toc;
end

% Richardson extrapolation bias estimate
eb = (V(2) - V(1)) / (1 - 2^delta);

if pinta
    for k = [1, 2]
        fprintf('N=%d, h=%.5f: V=%.5f ± %.5f, CPUt=%.3f s (%.2f%% hit barrier)\n', ...
            N, T/(k*M_), V(k), ster(k), CPUt(k), 100*(1 - sum(HIT)/N));
    end

    fprintf('Estimated bias = %.5g\n', eb);

    if abs(eb) < 2 * max([abs(ster(1)), abs(ster(2))])
        fprintf('WARNING: Statistical error NOT negligible relative to bias. Unreliable estimate!\n');
        eb = NaN;
    end
end
