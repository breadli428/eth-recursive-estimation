function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 4000; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
   
    radius = estConst.d * sqrt(rand(1, N_particles));
    angle = 2 * pi * rand(1, N_particles);
    center_pos = [estConst.pA; estConst.pB];
    init_center = (rand(1, N_particles) < 0.5) + 1;
    postParticles.x_r = center_pos(init_center, 1)' + radius .* cos(angle); % 1xN_particles matrix
    postParticles.y_r = center_pos(init_center, 2)' + radius .* sin(angle); % 1xN_particles matrix
    postParticles.phi = (rand(1, N_particles) * 2 - 1) * estConst.phi_0; % 1xN_particles matrix
    postParticles.kappa = (rand(1, N_particles) * 2 - 1) * estConst.l; % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!


% Prior Update:
v_f = (rand(1, N_particles) - 0.5) * estConst.sigma_f;
v_phi = (rand(1, N_particles) - 0.5) * estConst.sigma_phi;
priorParticles.x_r = prevPostParticles.x_r + (act(1) + v_f) .* cos(prevPostParticles.phi);
priorParticles.y_r = prevPostParticles.y_r + (act(1) + v_f) .* sin(prevPostParticles.phi);
priorParticles.phi = prevPostParticles.phi + act(2) + v_phi;
priorParticles.kappa = prevPostParticles.kappa;

% Posterior Update:
distance = zeros(1, N_particles);
density = zeros(1, N_particles);
for k = 1:N_particles
    distance(k) = compute_distance(priorParticles.x_r(k), priorParticles.y_r(k), priorParticles.phi(k), priorParticles.kappa(k), estConst.contour);
    density(k) = measurement_noise_density(sens - distance(k), estConst.epsilon);
end
while sum(density) == 0
    disp("Zero measurement likelihood encountered, generating new priors...")
    [priorParticles, density] = resample_prior(N_particles, priorParticles, sens, estConst);
end
beta = density / sum(density);

bin = resampling(beta);
for k = 1:N_particles
    postParticles.x_r(k) = priorParticles.x_r(bin(k));
    postParticles.y_r(k) = priorParticles.y_r(bin(k));
    postParticles.phi(k) = priorParticles.phi(bin(k));
    postParticles.kappa(k) = priorParticles.kappa(bin(k));
end
K = 0.01;
postParticles = roughening(postParticles, K);


    function density = measurement_noise_density(w, epsilon)
        w = abs(w);
        if w <= 2 * epsilon
            density =   -w / (5 * epsilon^2) + 2 / (5 * epsilon);
        elseif w <= 2.5 * epsilon
            density =  2 * w / (5 * epsilon^2) - 4 / (5 * epsilon);
        elseif w <= 3 * epsilon
            density = -2 * w/ (5 * epsilon^2) + 6 / (5 * epsilon);
        else
            density = 0;
        end
    end


    function distance = compute_distance(x_r, y_r, phi, kappa, contour)
        contour(8, 1) = kappa;
        contour(9, 1) = kappa;
        contour = [contour; contour(1,:)];
        ratio = zeros(1,10);
        dist = zeros(1,10);
        for i = 1:size(contour, 1) - 1
            x_1 = contour(i, 1);
            y_1 = contour(i, 2);
            x_2 = contour(i + 1, 1);
            y_2 = contour(i + 1, 2);
            % area of triangle r12 / dist
            d_area_r_1_2_d_dist = sin(phi) * (x_2 - x_1) + cos(phi) * (y_1 - y_2);
            % area of triangle rc2 / dist
            d_area_r_c_2_d_dist = sin(phi) * (x_2 - x_r) + cos(phi) * (y_r - y_2);
            % ratio = v_c2 / v_12 = area of triangle rc2 / area of triangle r12
            ratio(i) = d_area_r_c_2_d_dist / d_area_r_1_2_d_dist;
            % area of triangle r12
            area_r_1_2 = (y_2 - y_1) * (x_r - x_2) + (x_1 - x_2) * (y_r - y_2);
            dist(i)  = area_r_1_2 / d_area_r_1_2_d_dist;
        end
        detection_flag = ratio >= 0 & ratio <= 1 & dist >= 0;
        detected_dists = dist .* detection_flag;
        distance = min(detected_dists(detected_dists~=0));
        if isempty(distance)
            distance = inf;
        end
    end

    function [priorParticles, density_new] = resample_prior(N, priorParticles_0, sens, estConst)
        N_mass = 20 * N;
        priorParticles.x_r = mean(priorParticles_0.x_r) + (rand(1, N_mass) * 2 - 1) * 0.5;
        priorParticles.y_r = mean(priorParticles_0.y_r) + (rand(1, N_mass) * 2 - 1) * 0.5;
        priorParticles.phi = mean(priorParticles_0.phi) + (rand(1, N_mass) * 2 - 1) * pi / 3;
        priorParticles.kappa = mean(priorParticles_0.kappa) + (rand(1, N_mass) * 2 - 1) * 0.02;
        distance_new = zeros(1, N_mass);
        density_new = zeros(1, N_mass);
        for j = 1:N_mass
            distance_new(j) = compute_distance(priorParticles.x_r(j), priorParticles.y_r(j), priorParticles.phi(j), priorParticles.kappa(j), estConst.contour);
            density_new(j) = measurement_noise_density(sens - distance_new(j), estConst.epsilon);
        end
        [~, index] = maxk(density_new, N);
        priorParticles.x_r = priorParticles.x_r(index);
        priorParticles.y_r = priorParticles.y_r(index);
        priorParticles.phi = priorParticles.phi(index);
        priorParticles.kappa = priorParticles.kappa(index);
        density_new = density_new(index);
    end


    function bin = resampling(beta)
        gamma = cumsum([0 beta]);
        [~, ~, bin] = histcounts(rand(1, length(beta)), gamma);
    end

    function postParticles = roughening(postParticles, K)
        N = length(postParticles.x_r);
        d = 4;
        E_x_r = max(postParticles.x_r) - min(postParticles.x_r);
        E_y_r = max(postParticles.y_r) - min(postParticles.y_r);
        E_phi = max(postParticles.phi) - min(postParticles.phi);
        E_kappa = max(postParticles.kappa) - min(postParticles.kappa);
        sigma_x_r = K * E_x_r * N^(-1 / d);
        sigma_y_r = K * E_y_r * N^(-1 / d);
        sigma_phi = K * E_phi * N^(-1 / d);
        sigma_kappa = K * E_kappa * N^(-1 / d);
        postParticles.x_r = postParticles.x_r + sigma_x_r * randn(1, N);
        postParticles.y_r = postParticles.y_r + sigma_y_r * randn(1, N);
        postParticles.phi = postParticles.phi + sigma_phi * randn(1, N);
        postParticles.kappa = postParticles.kappa + sigma_kappa * randn(1, N);
    end

end % end estimator