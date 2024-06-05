function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    %defining the variables C,K & R before proceeding for the update
    C = [eye(6) zeros(6,9)];
    R = eye(6)*0.000021;
    K = covarEst*transpose(C)*inv(C*covarEst*transpose(C) + R);

    %updating the uCurr and covar_curr using K&C
    uCurr = uEst + K*(z_t - C*uEst);
    covar_curr = covarEst - K*C*covarEst;


end

