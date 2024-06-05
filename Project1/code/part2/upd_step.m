function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

I = [zeros(3,3),zeros(3,3),eye(3,3),zeros(3,3),zeros(3,3)];

%mesurement - z
z = I * uEst + zeros(3,1);

C = I ;

% Noise tuning
R = eye(3,3)*0.7;

%computing the Kalman Gain
K = covarEst * C' /(C * covarEst* C' + R );

%computing the mean and covariance
uCurr = uEst + K * (z_t - z );
covar_curr = covarEst- K* C * covarEst;
end