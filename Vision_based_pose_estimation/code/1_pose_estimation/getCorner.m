function res = getCorner(id)
%% CHANGE THE NAME OF THE FUNCTION TO getCorner
    %% Input Parameter Description
    % id = List of all the AprilTag ids detected in the current image(data)
    
    %% Output Parameter Description
    % res = List of the coordinates of the 4 corners (or 4 corners and the
    % centre) of each detected AprilTag in the image in a systematic method

    count_id = length(id); % count of number of tag_ids available

    %instantiating the corners
    P0 = zeros(2, count_id); 
    P1 = zeros(2, count_id);
    P2 = zeros(2, count_id);
    P3 = zeros(2, count_id);
    P4 = zeros(2, count_id);
    for i = 1:count_id
        [p0, p1, p2, p3, p4] = getcorner_individual(id(i)); %calling the getcorner_individual function to calculate for each step
        P0(:,i) = p0;
        P1(:,i) = p1;
        P2(:,i) = p2;
        P3(:,i) = p3;
        P4(:,i) = p4;
    end 
    
    res = [P0; P1 ;P2;P3;P4];
end