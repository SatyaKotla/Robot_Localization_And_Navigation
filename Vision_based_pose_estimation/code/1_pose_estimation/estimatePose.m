function [position, orientation] = estimatePose(data, t)
%% CHANGE THE NAME OF THE FUNCTION TO estimatePose
% Please note that the coordinates for each corner of each AprilTag are
% defined in the world frame, as per the information provided in the
% handout. Ideally a call to the function getCorner with ids of all the
% detected AprilTags should be made. This function should return the X and
% Y coordinate of each corner, or each corner and the centre, of all the
% detected AprilTags in the image. You can implement that anyway you want
% as long as the correct output is received. A call to that function
% should made from this function.
    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    
    %% Output Parameter Defination
    
    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order ZYX
    
    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order ZYX
  
    % getting the tag ids from the data struct

    tag_id = data(t).id;

    %Calling the getCorner function to get the calculated corner
    %coordinates
    res = getCorner(tag_id);
   
    P0 = res(1:2,:);
    P1 = res(3:4,:);
    P2 = res(5:6,:);
    P3 = res(7:8,:);
    P4 = res(9:10,:);

    %getting the image coordinates from the given data
    P0_bar = data(t).p0;
    P1_bar = data(t).p1;
    P2_bar = data(t).p2;
    P3_bar = data(t).p3;
    P4_bar = data(t).p4;

    %step1 finding H from A
    %Initiating the A matrix

    A = [];

    for i = 1:length(tag_id)
        a = zeros(10,9);
        
        a(1,:)=[P0(1,i) P0(2,i) 1 0 0 0 -P0_bar(1,i)*P0(1,i) -P0_bar(1,i)*P0(2,i) -P0_bar(1,i)];
        a(2,:)=[0 0 0 P0(1,i) P0(2,i) 1 -P0_bar(2,i).*P0(1,i) -P0_bar(2,i)*P0(2,i) -P0_bar(2,i)];
        a(3,:)=[P1(1,i) P1(2,i) 1 0 0 0 -P1_bar(1,i)*P1(1,i) -P1_bar(1,i)*P1(2,i) -P1_bar(1,i)];
        a(4,:)=[0 0 0 P1(1,i) P1(2,i) 1 -P1_bar(2,i)*P1(1,i) -P1_bar(2,i)*P1(2,i) -P1_bar(2,i)];
        a(5,:)=[P2(1,i) P2(2,i) 1 0 0 0 -P2_bar(1,i)*P2(1,i) -P2_bar(1,i)*P2(2,i) -P2_bar(1,i)];
        a(6,:)=[0 0 0 P2(1,i) P2(2,i) 1 -P2_bar(2,i)*P2(1,i) -P2_bar(2,i)*P2(2,i) -P2_bar(2,i)];
        a(7,:)=[P3(1,i) P3(2,i) 1 0 0 0 -P3_bar(1,i)*P3(1,i) -P3_bar(1,i)*P3(2,i) -P3_bar(1,i)];
        a(8,:)=[0 0 0 P3(1,i) P3(2,i) 1 -P3_bar(2,i)*P3(1,i) -P3_bar(2,i)*P3(2,i) -P3_bar(2,i)];
        a(9,:)=[P4(1,i) P4(2,i) 1 0 0 0 -P4_bar(1,i)*P4(1,i) -P4_bar(1,i)*P4(2,i) -P4_bar(1,i)];
        a(10,:)=[0 0 0 P4(1,i) P4(2,i) 1 -P4_bar(2,i)*P4(1,i) -P4_bar(2,i)*P4(2,i) -P4_bar(2,i)];
        
        A=[A;a];
    end
    
    %A*h = 0

    [~, ~, V] = svd(A,"econ");

    % Forming the matrix HL from the svd
    h = [V(1,9) V(2,9) V(3,9);
             V(4,9) V(5,9) V(6,9);
             V(7,9) V(8,9) V(9,9)];
    
    %Normalization of matrix HL
    singularvalues = svd(h);

    H = h/singularvalues(2,1);
    H = H *sign(V(9,9));

    %Calculating the pose and translation after getting H
    K = [311.0520 0 201.8724; 0 311.3885 113.6210; 0 0 1];
    R_mat = K\H;   %K_inv * H
    R1_hat = R_mat(:,1);
    R2_hat = R_mat(:,2);
    T_hat = R_mat(:,3);
    
    USV = [R1_hat R2_hat cross(R1_hat,R2_hat)];
    
    [U,~,V] = svd(USV);
    
    R_c2w = U*[1,0,0;0,1,0;0,0,det(U*V')]*V';
    T_c2w = T_hat/norm(R1_hat);

    H_c2w = [R_c2w T_c2w; 0 0 0 1]; %world to camera
     
    R_b2c = [sqrt(2)/2, -sqrt(2)/2, 0; -sqrt(2)/2, -sqrt(2)/2, 0; 0, 0, -1];
    T_b2c = [0.04; 0; -0.03];
    H_b2c = [R_b2c T_b2c; 0 0 0 1];            % transformation matrix body to camera
    H_w2b = inv(H_c2w)*inv(H_b2c);
    
    position = H_w2b(1:3,4); %position vector
    orientation = rotm2eul(H_w2b(1:3,1:3),"ZYX"); %orientation vector converted to euler angles
end