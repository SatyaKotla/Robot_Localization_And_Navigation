clear; % Clear variables
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
Z = sampledVicon(1:6,:);%all the measurements that you need for the update
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %J ust for saving state his.
prevTime = 0; %last time step in real time
%write your code here calling the pred_step.m and upd_step.m functions
for i = 1:length(sampledTime)

   % sensor.omg for angular velocity
   angVel= sampledData(i).omg;
   %sensor.acc to get the accelaration
   acc = sampledData(i).acc;

   % to get the delta_t - difference between current time i.e. sampledTime
   % at ith step - Previous Time
   dt  = sampledTime(i) - prevTime;
   %PrevTime to be replaced by current time 
   prevTime = sampledTime(i);
   % getting the measurement at time
   z_t = Z(:,i);

   %prediction step
   [covarEst,uEst] =  pred_step(uPrev,covarPrev,angVel,acc,dt);

   %update step
   [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst);

   %updating the previous mean,covariance with currenct mean and covariance
   uPrev = uCurr;
   covarPrev = covar_curr;
   savedStates(:,i) = uCurr;
end
plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);