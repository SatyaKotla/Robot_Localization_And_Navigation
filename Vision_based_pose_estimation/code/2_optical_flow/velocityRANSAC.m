function [Vel] = velocityRANSAC(optV,optPos,Z,R_c2w,e)
%% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
    %% Input Parameter Description
    % optV = The optical Flow
    % optPos = Position of the features in the camera frame 
    % Z = Height of the drone
    % R_c2w = Rotation defining camera to world frame
    % e = RANSAC hyper parameter    
    
    %% Output Parameter Description
    % Vel = Linear velocity and angualr velocity vector
    p_success = 0.99; 
    M = 3;
    k = log(1 - p_success) / log(1 - (e^M));  
    
    max_inliers = 0; %instantiating the maximum number of inliers
    
    Vel = zeros(6, 1);  % Initialize Vel

    for current_inliers = 1:k
        % Randomly select 3 positions for RANSAC indices
        random_index = randperm(length(optPos), 3);
        index = optPos(random_index, :);

        % Creating the H matrix
        
        H = zeros(6, 6);
        for l = 1:3
            h_l = H_matrix(index(l, 1), index(l, 2), Z);
            H((2 * l - 1):(2 * l), :) = h_l;
        end
        
        % Extracting the corresponding optical flow vectors
        v_opt = optV([2 * random_index(1) - 1, 2 * random_index(1)...
                      2 * random_index(2) - 1, 2 * random_index(2)...
                      2 * random_index(3) - 1, 2 * random_index(3)]);

        % Solving for the velocity 
        vel_and_omega = pinv(H) * v_opt;

        % Count inliers
        current_inliers = 0;
        for n = 1:length(optPos)
            H_n = H_matrix(optPos(n, 1), optPos(n, 2), Z);
            P_n = optV([2 * n - 1, 2 * n]);
            diff = norm(H_n * vel_and_omega - P_n)^2;
            if (diff <= 0.0001) %setting the threshold
                current_inliers = current_inliers + 1;
            end
        end
        
        % Update the max_inliers and Vel variables if the current iteration
        % has more inliers than the previous max inliers
        if (current_inliers >= max_inliers)
            max_inliers = current_inliers;
            Vel = vel_and_omega;
        end
    end 
end