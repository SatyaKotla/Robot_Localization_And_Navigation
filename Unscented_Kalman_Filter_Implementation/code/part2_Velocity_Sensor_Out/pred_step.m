 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
 %% Defining the tuning parameters to define the sigma points spread and noise
    alpha=0.001;
    k=1;
    beta = 2;
    n=27;
    Q_t=eye(12)*0.5;
    lambda=alpha^2*(n+k)-n; 

%%  mean and Covariance - augmented 

    mean_augmented=[uPrev;zeros(12,1)];
    covariance_augmented=[covarPrev zeros(15,12); zeros(12,15) Q_t];
    covariance_augmented=chol(covariance_augmented);

%% Step 1 : Compute sigma points with the augmented

    x_aug_zero = mean_augmented;
    x_aug_add=zeros(27,n);
    x_aug_sub=zeros(27,n);

    for i = 1:n
        x_aug_add(:,i) = mean_augmented+sqrt(n+lambda)*(covariance_augmented(:,i));
        x_aug_sub(:,i) = mean_augmented-sqrt(n+lambda)*(covariance_augmented(:,i));
    end
        
    sigma_points = [x_aug_zero x_aug_add x_aug_sub];

%%  Step2 - Propagating the sigma points throught the function - f i.e f_xun

process_model = pro_mod(sigma_points,angVel,acc,dt,n) ;

%% Step 3 - Computing  the predicted mean and covariance matrix
    
    Wc_0 = (lambda/(lambda+n))+(1-alpha^2+beta);
    Wc_i = 1/(2*(lambda+n));
    Wm_0 = lambda/(lambda+n);
    Wm_i = 1/(2*(lambda+n));
    uEst=[]; %initializing the uEst
    covarEst=[];%initializing the covarEst

    for j= 1:(2*n+1)
        if j == 1
            uEst = Wm_0*process_model(:,1);

        else 
            uEst = uEst + (Wm_i*process_model(:,j));
        end
    end

    for l = 1:(2*n+1)
        if l == 1
            covarEst = Wc_0*(process_model(:,1)-uEst)*(process_model(:,1)-uEst)';
        else
            covarEst=covarEst+(Wc_i*(process_model(:,l)-uEst)*(process_model(:,l)-uEst)');
        end
    end 
end

