clear; % Clear variables
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
Z = sampledVicon(7:9,:);%all the measurements that you need for the update
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %J ust for saving state his.
prevTime = 0; %last time step in real time
%write your code here calling the pred_step.m and upd_step.m functions
for i = 1:length(sampledTime)
   %angular velocity from the sensor struct data - sensor.omg
   angVel= sampledData(i).omg;

   %accelaration from the sensor struct data - sensor.omg
   acc = sampledData(i).acc;

   %Time difference
   dt  = sampledTime(i) - prevTime;

   %updating the previous time with current time i.e sampledTime(i)
   prevTime = sampledTime(i);
   
   %measurement inputs at time t
   z_t = Z(:,i);

   %prediction step
   [covarEst,uEst] =  pred_step(uPrev,covarPrev,angVel,acc,dt);

   %update step
   [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst);

   %updating the previous mean and covariance variables
   uPrev = uCurr;
   covarPrev = covar_curr;
   savedStates(:,i) = uCurr;

end
plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);