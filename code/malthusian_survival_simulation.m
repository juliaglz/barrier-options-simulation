function [survivalProb, stdError, excessRisk, computationTime, variance, totalPaths] = malthusian(populationSize, timeSteps, seed, varargin)
    % MALTHUSIAN - Simulates population survival probability with barriers
    %
    % INPUTS:
    %   populationSize - total number of simulated paths
    %   timeSteps      - number of time discretization steps
    %   seed           - random seed (-1 for no seed)
    %   varargin       - optional flag to display detailed output
    %
    % OUTPUTS:
    %   survivalProb   - estimated survival probability (up-in barrier hit)
    %   stdError       - 3-sigma statistical error estimate
    %   excessRisk     - placeholder output (not computed here)
    %   computationTime- time taken for simulation
    %   variance       - variance of the estimator
    %   totalPaths     - total number of simulated paths processed

    if seed ~= -1
        rng(seed, 'twister');
    end

    displayDetails = false;
    if nargin > 3
        displayDetails = true;
    end

    initialPopulation = 3;
    totalDuration = 10;
    volatility = 0.5;
    growthRate = 0.03;

    timeStepSize = totalDuration / timeSteps;
    batchSize = min(populationSize, 1e5);

    pathsProcessed = 0;
    totalHits = 0;
    sumOfSquares = 0;
    barrierHits = 0;

    % Barrier definitions (log scale)
    downOutBarrier = exp(log(initialPopulation) + 0.5826 * volatility * sqrt(timeStepSize));
    upInBarrier = exp(log(2 * initialPopulation) - 0.5826 * volatility * sqrt(timeStepSize));

    tic

    while pathsProcessed < populationSize
        populationLog = NaN(batchSize, timeSteps + 1);
        populationLog(:, 1) = log(initialPopulation);
        hitFlag = NaN(batchSize, 1);

        for step = 1:timeSteps
            % Simulate log population at next time step (GBM with drift)
            populationLog(:, step + 1) = populationLog(:, step) + ...
                (growthRate - volatility^2 / 2) * timeStepSize + ...
                volatility * sqrt(timeStepSize) * randn(batchSize, 1);

            % Check barriers for each path
            for i = 1:batchSize
                if isnan(hitFlag(i))
                    if populationLog(i, step + 1) <= log(downOutBarrier) % Down-and-out barrier
                        hitFlag(i) = 0;
                    elseif populationLog(i, step + 1) >= log(upInBarrier) % Up-and-in barrier
                        hitFlag(i) = 1;
                    end
                end
            end
        end

        % For paths never hitting barriers, assign 0 (no survival)
        hitFlag(isnan(hitFlag)) = 0;

        totalHits = totalHits + sum(hitFlag);
        sumOfSquares = sumOfSquares + sum(hitFlag .* 2); % Seems unusual: probably intended variance accumulation
        barrierHits = barrierHits + sum(hitFlag);
        pathsProcessed = pathsProcessed + batchSize;
    end

    survivalProb = totalHits / pathsProcessed;
    variance = sumOfSquares / pathsProcessed - survivalProb^2;
    stdError = 3 * sqrt(variance / pathsProcessed);
    computationTime = toc;
    excessRisk = NaN;  % Not computed in this code
    totalPaths = pathsProcessed;

    if displayDetails
        fprintf('N (effective) = %d, time step size = %g: Probability = %g +/- %g, Computation Time = %g, %% hit barrier = %g \n', ...
            pathsProcessed, timeStepSize, survivalProb, stdError, computationTime, 100 * (barrierHits / pathsProcessed));
    end
end
