function [R,G] = rot_G_matrix(angles)
    %computing the rotation matrix from the given euler angles

    
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);

    Rx = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    R = Rz*Ry*Rx;
    G = [-sin(theta) 0 1; cos(theta)*sin(phi) cos(phi) 0;cos(theta)*cos(phi) -sin(phi) 0];
    
end