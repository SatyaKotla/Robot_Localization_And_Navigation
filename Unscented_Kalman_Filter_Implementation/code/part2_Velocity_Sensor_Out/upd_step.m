function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    k = 1;
    alpha_bar = 0.001;
    beta_bar = 2;
    %n_bar = n + nv;
    n_bar = 15;
    lambda_bar = alpha_bar^2 * (n_bar + k) - n_bar;

    %% Step 1 : Compute sigma points with the augmented

    covariance_augmented = chol(covarEst,"lower");
    x_aug_zero = uEst;
    x_aug_add=zeros(15,n_bar);
    x_aug_sub=zeros(15,n_bar);
    for i = 1:n_bar
        x_aug_add(:,i) = uEst+sqrt(n_bar+lambda_bar)*(covariance_augmented(:,i));
        x_aug_sub(:,i) = uEst-sqrt(n_bar+lambda_bar)*(covariance_augmented(:,i));     
    end
        
    sigma_points = [x_aug_zero x_aug_add x_aug_sub];

    %% Step 2: Propagate Sigma Points through the nonlinear function 𝑔
    Rcb=[sqrt(2)/2, -sqrt(2)/2, 0; -sqrt(2)/2, -sqrt(2)/2, 0; 0, 0, -1];
    Rbc=transpose(Rcb);
    tbc=[-0.04; 0;-0.03];

    zt_i=zeros(3,length(2*(n_bar)+1));
     vt = normrnd(0,(0.00093),[3,1]);
     for i = 1:(2*n_bar+1)
         Rbw = eul2rotm(uEst(4:6)');
         Rwb = inv(Rbw);
         g_i =Rcb*Rwb'*[sigma_points(7,i);sigma_points(8,i);sigma_points(9,i)] - Rbc'*(cross(tbc',(Rbc'*z_t(4:6))'))'+vt;
         zt_i(:,i) = g_i + vt;       
     end

    %% Step 3: Compute the predicted mean, predicted covariance of the measurement, and predicted cross-covariance

    Wc_0_bar = lambda_bar/((lambda_bar+n_bar))+(1-alpha_bar^2+beta_bar);
    Wc_i_bar = 1/(2*(lambda_bar+n_bar));
    Wm_0_bar = lambda_bar/(lambda_bar+n_bar);
    Wm_i_bar = 1/(2*(lambda_bar+n_bar));
    
    %Computing the z_mu_t using vectorized approach
    weight_length = size(zt_i,2);
    weights = repmat(Wm_i_bar, 1, weight_length);
    weights(1) = Wm_0_bar;
    z_mu_t = zt_i * weights';

    %Calculating C and S
    C = zeros(15,3);
    S = zeros(3,3);
    R = eye(3)*0.011;
       
        for n=1:(2*n_bar+1)
            if n==1
                C=Wc_0_bar*(sigma_points(:,n)-(uEst))*((zt_i(:,n))-z_mu_t)';
                S = Wc_0_bar*((zt_i(:,n))-z_mu_t)*((zt_i(:,n))-z_mu_t)';
            else
                C=C+Wc_i_bar*(sigma_points(:,n)-(uEst))*((zt_i(:,n))-z_mu_t)';
                S = S+Wc_i_bar*((zt_i(:,n))-z_mu_t)*((zt_i(:,n))-z_mu_t)';
            end
        end
        
        S = S + R;
    
    % Computing the filter gain and the filtered state mean and covariance, conditional to the measurement

    K = C/(S);
    
    uCurr = uEst + K * (z_t(1:3) - z_mu_t);
    covar_curr = covarEst - (K * S * transpose(K));
end

