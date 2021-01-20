close all
load('data.mat')

%% Specify common parameters and create containers
num_particles = 100;
num_steps = 625;
num_measurements = 78;
X_BestParticles = zeros(size(XTRUE));
X_AverageParticles = zeros(size(XTRUE));
X_WorstParticles = zeros(size(XTRUE));
predicted_X = repmat([0; 0; 0], 1, num_particles);
corrected_X = predicted_X;
max_index = randi(num_particles);
avg_index = randi(num_particles);
min_index = randi(num_particles);

%% Iterate through dataset to filter
for i = 1:num_steps
    % Sample motion model
    V = normrnd(VG(1, i), 0.5, [1, num_particles]);
    G = normrnd(VG(2, i), 3*pi/180, [1, num_particles]);
    predicted_X(1, :) = corrected_X(1, :) + 0.025*V.*cos(G + corrected_X(3, :));
    predicted_X(2, :) = corrected_X(2, :) + 0.025*V.*sin(G + corrected_X(3, :));
    predicted_X(3, :) = corrected_X(3, :) + 0.025*V.*sin(G)/4;
    predicted_X(3, :) = atan2(sin(predicted_X(3, :)), cos(predicted_X(3, :)));
    
    if (~mod(i, 8))
        % Calculate weights from measurement model
        sigma = [0.2 2*pi/180].^2;
        w = ones(num_particles, 1);
        for k = 1:size(Z, 2)
            j = Z(3, k, i);
            if (~isnan(j))
                mu = [sqrt((lm(1, j) - predicted_X(1, :)).^2 + (lm(2, j) - predicted_X(2, :)).^2);
                      atan2(lm(2, j) - predicted_X(2, :), lm(1, j) - predicted_X(1, :)) - predicted_X(3, :)]';
                mu(:, 2) = atan2(sin(mu(:, 2)), cos(mu(:, 2)));
                w_tmp = mvnpdf(Z(1:2, k, i)', mu, sigma);
                if (sum(w_tmp) > 0)
                    w = w.*w_tmp;
                end
                w = w/sum(w);
            end
        end
        max_index = find(w == max(w), 1);
        avg_index = find(abs(w - mean(w)) == min(abs(w - mean(w))), 1);
        min_index = find(w == min(w), 1);
        
        % Perform low-variance resampling
        J = num_particles;
        r = normrnd(0, 1/J);
        c = w(1);
        k = 1;
        for j = 1:J
            U = r + (j - 1)/J;
            while (U > c)
                k = k + 1;
                if (k > J)
                    k = k - J;
                end
                c = c + w(k);
            end
            corrected_X(:, j) = predicted_X(:, k);
        end
    else
        corrected_X = predicted_X;
    end
    
    % Save data for plotting
    X_BestParticles(:, i) = predicted_X(:, max_index);
    X_AverageParticles(:, i) = predicted_X(:, avg_index);
    X_WorstParticles(:, i) = predicted_X(:, min_index);
end

%% Calculate RMS of each trajectory with respect to the true path
error_odo = rms(XODO - XTRUE, 2);
error_max = rms(X_BestParticles - XTRUE, 2);
error_avg = rms(X_AverageParticles - XTRUE, 2);
error_min = rms(X_WorstParticles - XTRUE, 2);
disp('             error_x     error_y     error_theta')
disp('-------------------------------------------------')
disp(['error_odo', '    ', num2str(error_odo', 5)])
disp(['error_max', '    ', num2str(error_max', 5)])
disp(['error_avg', '    ', num2str(error_avg', 5)])
disp(['error_min', '    ', num2str(error_min', 5)])

%% Plot the trajectories
% Plot the trajectory of the best particle
figure('Position', get(groot,'ScreenSize'))
plot(XTRUE(1, :), XTRUE(2, :), '-.k')
hold on
plot(XODO(1, :), XODO(2, :), ':b')
hold on
plot(X_BestParticles(1, :), X_BestParticles(2, :), '-r')
legend({'True Path', 'Unfiltered Path', 'Filtered Path'}, 'Location', 'southwest')
title('Trajectory of Best Particles')

% Plot the trajectory of the average particle
figure('Position', get(groot,'ScreenSize'))
plot(XTRUE(1, :), XTRUE(2, :), '-.k')
hold on
plot(XODO(1, :), XODO(2, :), ':b')
hold on
plot(X_AverageParticles(1, :), X_AverageParticles(2, :), '-r')
legend({'True Path', 'Unfiltered Path', 'Filtered Path'}, 'Location', 'southwest')
title('Trajectory of Average Particles')

% Plot the trajectory of the worst particle
figure('Position', get(groot,'ScreenSize'))
plot(XTRUE(1, :), XTRUE(2, :), '-.k')
hold on
plot(XODO(1, :), XODO(2, :), ':b')
hold on
plot(X_WorstParticles(1, :), X_WorstParticles(2, :), '-r')
legend({'True Path', 'Unfiltered Path', 'Filtered Path'}, 'Location', 'southwest')
title('Trajectory of Worst Particles')