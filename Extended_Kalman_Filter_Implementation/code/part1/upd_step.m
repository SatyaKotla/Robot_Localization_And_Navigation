function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

I = [eye(3),zeros(3),zeros(3),zeros(3),zeros(3);
    zeros(3),eye(3),zeros(3),zeros(3),zeros(3)];

z = I * uEst + zeros(6,1);


C = I ;                                                                    

% Noise tuning - R
R = eye(6,6) * 0.7;

%Computing the kalman gain
K = covarEst * C' *inv(C * covarEst* C' + R );

%updating the mean and covariance with measurement inputs
uCurr = uEst + K * (z_t - z );
covar_curr = covarEst- K* C * covarEst;

end