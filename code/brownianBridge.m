function bridge = brownianBridge(N, T)
% Generates a Brownian bridge vector of size N over interval T
bridge = zeros(N, 1);
t = linspace(0, T, N + 2)';
t = t(2:end-1);  % Exclude 0 and T

bridge(1) = sqrt(T) * randn;
for i = 2:N
    dt = t(i) - t(i - 1);
    bridge(i) = bridge(i - 1) + sqrt(dt) * randn;
end
end
