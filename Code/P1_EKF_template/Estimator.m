function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst = zeros(1, 2); % 1x2 matrix
    linVelEst = zeros(1, 2); % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix
    
    % initial state variance
    posVar = estConst.StartRadiusBound^2 / 4 * ones(1, 2); % 1x2 matrix
    linVelVar = zeros(1, 2); % 1x2 matrix
    oriVar = estConst.RotationStartBound^2 / 3; % 1x1 matrix
    windVar = estConst.WindAngleStartBound^2 / 3; % 1x1 matrix
    driftVar = estConst.GyroDriftStartBound^2 / 3; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar, linVelVar, oriVar, windVar, driftVar]);
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst]';
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.
% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% prior update
tspan = [tm - dt, tm];
num_state = 7;
prior_init = [estState.xm; reshape(estState.Pm, [], 1)];

[~, prior] = ode45(@(t, prior) prior_update(prior, actuate, estConst), tspan, prior_init);

    function sys_col = prior_update(prior, actuate, estConst)
        x_hat = prior(1:7);
        s_x = x_hat(3);
        s_y = x_hat(4);
        phi = x_hat(5);
        rho = x_hat(6);
        
        u_t = actuate(1);
        u_r = actuate(2);
        
        v_d = 0;
        v_r = 0;
        v_rho = 0;
        v_b = 0;
        
        q = [s_x;
            s_y;
            cos(phi) * (tanh(u_t) - estConst.dragCoefficientHydr * (s_x^2 + s_y^2) * (1 + v_d)) - estConst.dragCoefficientAir * (s_x - estConst.windVel * cos(rho)) * sqrt((s_x - estConst.windVel * cos(rho))^2 + (s_y - estConst.windVel * sin(rho))^2);
            sin(phi) * (tanh(u_t) - estConst.dragCoefficientHydr * (s_x^2 + s_y^2) * (1 + v_d)) - estConst.dragCoefficientAir * (s_y - estConst.windVel * sin(rho)) * sqrt((s_x - estConst.windVel * cos(rho))^2 + (s_y - estConst.windVel * sin(rho))^2);
            estConst.rudderCoefficient * u_r * (1 + v_r);
            v_rho;
            v_b];
        
        A = [0, 0,                                                                                                                                                                                                                              1,                                                                                                                                                                                                                              0,                                                       0,                                                                                                                                                                                                                                                                  0, 0;
            0, 0,                                                                                                                                                                                                                              0,                                                                                                                                                                                                                              1,                                                       0,                                                                                                                                                                                                                                                                  0, 0;
            0, 0, - (3*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2))/50 - (s_x*cos(phi)*(v_d + 1))/5 - ((2*s_x - (3*cos(rho))/2)*((3*s_x)/50 - (9*cos(rho))/200))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)),                                                                      - (s_y*cos(phi)*(v_d + 1))/5 - (((3*s_x)/50 - (9*cos(rho))/200)*(2*s_y - (3*sin(rho))/2))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)), -sin(phi)*(tanh(u_t) - (s_x^2/10 + s_y^2/10)*(v_d + 1)), - (9*sin(rho)*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2))/200 - (((3*s_x)/50 - (9*cos(rho))/200)*((3*sin(rho)*(s_x - (3*cos(rho))/4))/2 - (3*cos(rho)*(s_y - (3*sin(rho))/4))/2))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)), 0;
            0, 0,                                                                      - (s_x*sin(phi)*(v_d + 1))/5 - ((2*s_x - (3*cos(rho))/2)*((3*s_y)/50 - (9*sin(rho))/200))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)), - (3*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2))/50 - (s_y*sin(phi)*(v_d + 1))/5 - ((2*s_y - (3*sin(rho))/2)*((3*s_y)/50 - (9*sin(rho))/200))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)),  cos(phi)*(tanh(u_t) - (s_x^2/10 + s_y^2/10)*(v_d + 1)),   (9*cos(rho)*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2))/200 - (((3*s_y)/50 - (9*sin(rho))/200)*((3*sin(rho)*(s_x - (3*cos(rho))/4))/2 - (3*cos(rho)*(s_y - (3*sin(rho))/4))/2))/(2*((s_x - (3*cos(rho))/4)^2 + (s_y - (3*sin(rho))/4)^2)^(1/2)), 0;
            0, 0,                                                                                                                                                                                                                              0,                                                                                                                                                                                                                              0,                                                       0,                                                                                                                                                                                                                                                                  0, 0;
            0, 0,                                                                                                                                                                                                                              0,                                                                                                                                                                                                                              0,                                                       0,                                                                                                                                                                                                                                                                  0, 0;
            0, 0,                                                                                                                                                                                                                              0,                                                                                                                                                                                                                              0,                                                       0,                                                                                                                                                                                                                                                                  0, 0];
        
        L = [                          0,     0, 0, 0;
            0,     0, 0, 0;
            -cos(phi)*(s_x^2/10 + s_y^2/10),     0, 0, 0;
            -sin(phi)*(s_x^2/10 + s_y^2/10),     0, 0, 0;
            0, 2*u_r, 0, 0;
            0,     0, 1, 0;
            0,     0, 0, 1];
        
        
        P = reshape(prior(8:end), [7 7]);
        Q_c = diag([estConst.DragNoise, estConst.RudderNoise, estConst.WindAngleNoise, estConst.GyroDriftNoise]);
        P_dot = (A * P) + (P * A') + (L * Q_c * L');
        
        sys_col = [q; reshape(P_dot, [], 1)];
    end

x_hat_p = prior(end, 1:num_state)';
P_p = reshape(prior(end, num_state+1:end)', [num_state num_state]);

p_x = x_hat_p(1);
p_y = x_hat_p(2);
phi_p = x_hat_p(5);
b_p = x_hat_p(7);

w_a = 0;
w_b = 0;
w_c = 0;
w_g = 0;
w_n = 0;

h = [w_a + ((p_x + 1000)^2 + (p_y - 1000)^2)^(1/2);
    w_b + ((p_x - 2000)^2 + p_y^2)^(1/2);
    w_c + ((p_y - 2000)^2 + p_x^2)^(1/2);
    b_p + phi_p + w_g;
    phi_p + w_n];

H = [(2*p_x + 2000)/(2*((p_x + 1000)^2 + (p_y - 1000)^2)^(1/2)), (2*p_y - 2000)/(2*((p_x + 1000)^2 + (p_y - 1000)^2)^(1/2)), 0, 0, 0, 0, 0
              (2*p_x - 4000)/(2*((p_x - 2000)^2 + p_y^2)^(1/2)),                         p_y/((p_x - 2000)^2 + p_y^2)^(1/2), 0, 0, 0, 0, 0
                             p_x/((p_y - 2000)^2 + p_x^2)^(1/2),          (2*p_y - 4000)/(2*((p_y - 2000)^2 + p_x^2)^(1/2)), 0, 0, 0, 0, 0
                                                              0,                                                          0, 0, 0, 1, 0, 1
                                                              0,                                                          0, 0, 0, 1, 0, 0];

M = eye(5);                                                          
                                                          
R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);

if ~isfinite(sense(3))
    h(3) = [];
    sense(3) = [];
    H(3, :) = [];
    M = eye(4);
    R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.GyroNoise, estConst.CompassNoise]);
end

K = (P_p * H') / ((H * P_p * H') + (M * R * M'));
x_hat_m = x_hat_p + K * (sense' - h);
P_m = (eye(num_state) - (K * H)) * P_p;

% Get resulting estimates and variances
estState.xm = x_hat_m;
estState.Pm = P_m;

% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
windEst = estState.xm(6);
driftEst = estState.xm(7);

posVar = [estState.Pm(1, 1), estState.Pm(2, 2)];
linVelVar = [estState.Pm(3, 3), estState.Pm(4, 4)];
oriVar = estState.Pm(5, 5);
windVar = estState.Pm(6, 6);
driftVar = estState.Pm(7, 7);

end