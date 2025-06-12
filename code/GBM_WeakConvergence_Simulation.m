% Problem parameters and exact solution
r = 0.1;
sig = 0.3;
T = 0.2;
S0 = 100;
K = 100;
B = 85;
Ve = 6.3076;

% Monte Carlo simulation parameters
M = 1e7;
M2 = 1e4;

for p = 1:7
    N = 2^p;
    h = T / N;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;

    for m = 1:M2:M
        m2 = min(M2, M - m + 1);
        S = S0 * ones(1, m2);
        S2 = S0 * ones(1, m2);

        for n = 1:N/2
            dW1 = sqrt(h) * randn(1, m2);
            S = S .* (1 + r*h + sig * dW1);

            dW2 = sqrt(h) * randn(1, m2);
            S = S .* (1 + r*h + sig * dW2);

            S2 = S2 .* (1 + r*2*h + sig * (dW1 + dW2));
        end

        P = exp(-r * T) * max(S - K, 0);
        P2 = exp(-r * T) * max(S2 - K, 0);

        sum1 = sum1 + sum(P);
        sum2 = sum2 + sum(P.^2);
        sum3 = sum3 + sum(P - P2);
        sum4 = sum4 + sum((P - P2).^2);
    end

    hh(p) = h;
    err1(p) = sum1 / M - Ve;
    err2(p) = 3 * sqrt((sum2 / M - (sum1 / M)^2) / (M - 1));
    err3(p) = sum3 / M;
    err4(p) = 3 * sqrt((sum4 / M - (sum3 / M)^2) / (M - 1));
end

% Plot weak convergence errors comparing to exact solution
figure;
pos = get(gcf, 'pos');
pos(3:4) = pos(3:4) * [0.8 0.8];
set(gcf, 'pos', pos);
loglog(hh, abs(err1), 'b-*', hh, err2, 'r-*');
title('Weak convergence -- comparison to exact solution');
xlabel('h');
ylabel('Error');
legend('Weak error', 'MC error', 'location', 'NorthWest');

% Plot weak convergence errors comparing difference from 2h approximation
figure;
pos = get(gcf, 'pos');
pos(3:4) = pos(3:4) * [0.8 0.8];
set(gcf, 'pos', pos);
loglog(hh, abs(err3), 'b-*', hh, err4, 'r-*');
title('Weak convergence -- difference from 2h approximation');
xlabel('h');
ylabel('Error');
legend('Weak error', 'MC error', 'location', 'NorthWest');
