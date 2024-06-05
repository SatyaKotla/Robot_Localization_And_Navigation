%% PROJECT 2 VELOCITY ESTIMATION
close all;
clear all;
clc;
addpath('../data')

%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 1;
ransac = false;  % Flag for RANSAC

[sampledData, sampledVicon, sampledTime] = init(datasetNum);

%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION
K = [311.0520, 0, 201.8724;
     0, 311.3885, 113.6210;
     0, 0, 1];

t = [sampledData.t]; %extracting the time stamp
t = movmean(t, 10); % smoothing the time stamp using movmean

R_b2c = [sqrt(2)/2, -sqrt(2)/2, 0; -sqrt(2)/2, -sqrt(2)/2, 0; 0, 0, -1]; % rotation matrix body to camera
T_b2c = [0.04; 0; -0.03]; % translation vector body to camera
H_b2c = [R_b2c, T_b2c; 0 0 0 1]; % homogenous or transformation matrix body to camera

estimatedV = zeros(6, length(sampledData)); % initializing estimated velocity variable


for n = 2:length(sampledData)
    %% Initalize Loop load images
    prev_img = sampledData(n-1).img;
    cur_img = sampledData(n).img;
    d_t = t(n) - t(n-1);

    %% Detect good points
    %Using SURF Features algorithm to detect and match points
    prev_pts = detectSURFFeatures(prev_img);
    cur_pts = detectSURFFeatures(cur_img);
    [prev_features, prev_validPts] = extractFeatures(prev_img, prev_pts);
    [cur_features, cur_validPts] = extractFeatures(cur_img, cur_pts);
    index_pairs = matchFeatures(prev_features, cur_features);
    matched_prev_pts = prev_validPts(index_pairs(:, 1));
    matched_cur_pts = cur_validPts(index_pairs(:, 2));

    % Converting points to homogeneous coordinates
    matched_prev_pts = [matched_prev_pts.Location, ones(size(matched_prev_pts, 1), 1)]';
    matched_cur_pts = [matched_cur_pts.Location, ones(size(matched_cur_pts, 1), 1)]';
   
    %% Initalize the tracker to the last frame.
    pt_track = vision.PointTracker('MaxBidirectionalError', 1);

    %% Find the location of the next points;
    initialize(pt_track, matched_prev_pts(1:2, :)', prev_img);
    [pts, validity] = pt_track(cur_img);

    %% Calculate velocity
    % Use a for loop
    
    %% Calculate Height
    [position, orientation, R_c2w] = estimatePose(sampledData, n);

    R_w2b = eul2rotm(orientation);
    H_w2b = [R_w2b, position; 0, 0, 0, 1];
    H_b2w = inv(H_w2b);
    H_c2w = H_b2c\H_b2w; %inv(H_b2c)*H_b2w
    H_w2c = inv(H_c2w); %Homogenous matrix - world to camera
    R_w2c = H_w2c(1:3, 1:3);
    
    mpn = [];
    pn = [];
    for m = 1:length(pts)
        temp_pn = [pts(m, 1); pts(m, 2); 1];
        temp_mpn = [matched_prev_pts(1:2, m); 1];
        temp_pn = K\temp_pn; %inv(K)
        temp_mpn = K\temp_mpn;
        pn = [pn; transpose(temp_pn)];
        mpn = [mpn; transpose(temp_mpn)];
    end
    vel = [];
    for i = 1:length(pts)
        v_i = [(pn(i, 1) - mpn(i, 1)) / d_t; (pn(i, 2) - mpn(i, 2)) / d_t];
        vel = [vel; v_i];
    end
    %% RANSAC    
    % Write your own RANSAC implementation in the file velocityRANSAC

    if ransac == false
        H = [];
        for l = 1:length(pts)
            x = pn(l, 1);
            y = pn(l, 2);
            Z = position(3) / (dot([x; y; 1], -1 * R_w2c(:, 3)));
            H_i = H_matrix(x, y, Z);
            H = [H; H_i];
        end

        current_vel = pinv(H) * vel;
        Hb2c_skew = [0, -T_b2c(3, 1), T_b2c(2, 1); 
                     T_b2c(3, 1), 0, -T_b2c(1, 1); 
                     -T_b2c(2, 1), T_b2c(1, 1), 0];
        T_matrix = [R_w2b, zeros(3); zeros(3), R_w2b] * [R_b2c, -R_b2c * Hb2c_skew;zeros(3), R_b2c];
        Vel = T_matrix * current_vel;
    else
        for l = 1:length(pts) 
        x = pn(l,1);   
        y = pn(l,2);
        Z = position(3)/(dot([x;y;1],-1*R_w2c(:,3)));
        
        end
        e = 0.7;
        current_vel = velocityRANSAC(vel, pn, Z, R_c2w, e);
        Hb2c_skew = [0, -T_b2c(3, 1), T_b2c(2, 1); 
                     T_b2c(3, 1), 0, -T_b2c(1, 1); 
                     -T_b2c(2, 1), T_b2c(1, 1), 0];
        T_matrix = [R_w2b, zeros(3);zeros(3), R_w2b] * [R_b2c, -R_b2c * Hb2c_skew;zeros(3), R_b2c];
        Vel = T_matrix * current_vel;
    end
    %% Thereshold outputs into a range.
    % Not necessary
    
    %% Fix the linear velocity
    % Change the frame of the computed velocity to world frame
    
    %% ADD SOME LOW PASS FILTER CODE
    % Not neceessary but recommended 
    %estimatedV(:,n) = Vel;
    
    %% STORE THE COMPUTED VELOCITY IN THE VARIABLE estimatedV AS BELOW
    estimatedV(:,n) = Vel; % Feel free to change the variable Vel to anything that you used.
    % Structure of the Vector Vel should be as follows:
    % Vel(1) = Linear Velocity in X
    % Vel(2) = Linear Velocity in Y
    % Vel(3) = Linear Velocity in Z
    % Vel(4) = Angular Velocity in X
    % Vel(5) = Angular Velocity in Y
    % Vel(6) = Angular Velocity in Z
end

%Smoothing the output using movmean function with window size of 10
for i = 1:size(estimatedV, 1)
    estimatedV(i, :) = movmean(double(estimatedV(i, :)), 10);
end
plotData(estimatedV, sampledData, sampledVicon, sampledTime, datasetNum)
