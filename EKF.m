close all
load('data.mat')

%% Specify common parameters and create containers
num_steps = 625;
R = diag([1 1 1])^2;
Q_beam = diag([0.2 2*pi/180])^2;
predicted_mu = zeros(3, 1);
predicted_sigma = eye(3);
corrected_mu = predicted_mu;
corrected_sigma = predicted_sigma;
X_EKF = zeros(size(XTRUE));

%% Iterate through dataset to filter
for i = 1:num_steps
    % Perform prediction
    predicted_mu(1) = corrected_mu(1) + 0.025*VG(1, i)*cos(corrected_mu(3) + VG(2, i));
    predicted_mu(2) = corrected_mu(2) + 0.025*VG(1, i)*sin(corrected_mu(3) + VG(2, i));
    predicted_mu(3) = corrected_mu(3) + 0.025*VG(1, i)*sin(VG(2, i))/4;
    predicted_mu(3) = atan2(sin(predicted_mu(3)), cos(predicted_mu(3)));
    G = eye(3) + 0.025*VG(1, i)*[zeros(3, 2) [-sin(corrected_mu(3) + VG(2, i));
                                               cos(corrected_mu(3) + VG(2, i));
                                               0]];
    predicted_sigma = G*corrected_sigma*G' + R;
    
    % Perform correction
    if (~mod(i, 8))
        num_beams = sum(Z(3, :, i) > 0);
        H = zeros(2*num_beams, 3);
        Q = diag(diag(repmat(Q_beam, num_beams, num_beams)));
        z = zeros(2*num_beams, 1);
        z_ = zeros(2*num_beams, 1);
        n = 0;
        for k = 1:size(Z, 2)
            j = Z(3, k, i);
            if (~isnan(j))
                dx = lm(1, j) - predicted_mu(1);
                dy = lm(2, j) - predicted_mu(2);
                q = dx^2 + dy^2;
                n = n + 1;
                H(n:n+1, :) = 1/q*[-sqrt(q)*dx  -sqrt(q)*dy     0; 
                                            dy          -dx    -q];
                z(n:n+1) = Z(1:2, k, i);
                z_(n:n+1) = [sqrt(q);
                             atan2(dy, dx) - predicted_mu(3)];
                z_(n+1) = atan2(sin(z_(n+1)), cos(z_(n+1)));
            end
        end
        K = predicted_sigma*H'/(H*predicted_sigma*H' + Q);
        corrected_mu = predicted_mu + K*(z - z_);
        corrected_sigma = (eye(3) - K*H)*predicted_sigma;
    else
        corrected_mu = predicted_mu;
        corrected_sigma = predicted_sigma;
    end
    
    % Save data for plotting
    X_EKF(:, i) = corrected_mu;
end

%% Calculate RMS of each trajectory with respect to the true path
XTRUE_norm = normc(XTRUE)/2;
error_odo = rms(XODO - XTRUE, 2);
error_ekf = rms(X_EKF - XTRUE, 2);
disp('             error_x     error_y     error_theta')
disp('-------------------------------------------------')
disp(['error_odo', '    ', num2str(error_odo', 5)])
disp(['error_ekf', '    ', num2str(error_ekf', 5)])

%% Plot the trajectories
figure('Position', get(groot,'ScreenSize'))
plot(XTRUE(1, :), XTRUE(2, :), '-.k')
hold on
plot(XODO(1, :), XODO(2, :), ':b')
hold on
plot(X_EKF(1, :), X_EKF(2, :), '-r')
legend({'True Path', 'Unfiltered Path', 'Filtered Path'}, 'Location', 'southwest')
title('Trajectory filtered by EKF')