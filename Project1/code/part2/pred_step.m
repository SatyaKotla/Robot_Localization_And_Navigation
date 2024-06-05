 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time
wx = angVel(1,1);
wy = angVel(2,1);
wz = angVel(3,1);

ax = acc(1,1);
ay = acc(2,1);
az = acc(3,1);

vx = uPrev(7,1);
vy = uPrev(8,1);
vz = uPrev(9,1);

roll = uPrev(4,1);
pitch = uPrev(5,1);
yaw = uPrev(6,1);

bgx = uPrev(10,1);
bgy = uPrev(11,1);
bgz = uPrev(12,1);

bax = uPrev(13,1);
bay = uPrev(14,1);
baz = uPrev(15,1);

% Rotation Matrix
% Rx = [1 0 0 ; 0 cos(roll) -sin(roll) ; 0 sin(roll) cos(roll)]
% Ry = [cos(pitch) 0 sin(pitch) ; 0 1 0 ; -sin(pitch) 0 cos(pitch)]
% Rz = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ;0 0 1]

% G(x2)^-1 =[(cos(yaw)*sin(pitch))/cos(pitch) (sin(pitch)*sin(yaw))/cos(pitch) 1 ; 
%             -sin(yaw) cos(yaw) 0;
%             cos(yaw)/cos(pitch) sin(yaw)/cos(pitch) 0];                                                             

% f = [x3 ; G(x2)^-1*Rot_zyx*(wm - x4 - ng) ; g + (Rot_zyx)*(am - x5 - na); nbg ; nba];

f = [vx; vy; vz; ((sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch))*(bgz - wz) - ((sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch) - (cos(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch))*(bgy - wy) - (bgx - wx)*(cos(yaw)^2 + sin(yaw)^2); (bgz - wz)*(cos(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + sin(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))) - (bgy - wy)*(cos(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) + sin(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll))); - (bgz - wz)*(cos(pitch)*cos(roll) + (cos(yaw)*sin(pitch)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch) - (sin(pitch)*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch)) - (bgy - wy)*(cos(pitch)*sin(roll) - (cos(yaw)*sin(pitch)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) + (sin(pitch)*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch)) - (bgx - wx)*(sin(pitch)*cos(yaw)^2 + sin(pitch)*sin(yaw)^2 - sin(pitch)); (az - baz)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) - (ay - bay)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)) + cos(pitch)*cos(yaw)*(ax - bax); (ay - bay)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - (az - baz)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + cos(pitch)*sin(yaw)*(ax - bax); cos(pitch)*cos(roll)*(az - baz) - sin(pitch)*(ax - bax) + cos(pitch)*sin(roll)*(ay - bay) - (981/100); (0); (0); (0); (0); (0); (0)];

uEst  = uPrev + dt* f;
    

%At = jacobian(f,x)
At = [(0), (0), (0), (0), (0), (0), (1), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (1), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (1), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), ((sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch))*(bgy - wy) + ((sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch) - (cos(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch))*(bgz - wz), - (bgz - wz)*(cos(roll)*cos(yaw)^2 + cos(roll)*sin(yaw)^2 + (cos(yaw)*sin(pitch)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch)^2 - (sin(pitch)*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch)^2) - (bgy - wy)*(cos(yaw)^2*sin(roll) + sin(roll)*sin(yaw)^2 - (cos(yaw)*sin(pitch)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch)^2 + (sin(pitch)*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch)^2), (0), (0), (0), (0), - cos(yaw)^2 - sin(yaw)^2, (cos(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) - (sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch), ((sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch)), (0), (0), (0); (0), (0), (0), (bgy - wy)*(cos(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + sin(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))) + (bgz - wz)*(cos(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) + sin(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll))), (0), (0), (0), (0), (0), (0), - cos(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - sin(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)), (cos(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + sin(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))), (0), (0), (0);
    (0), (0), (0), (bgz - wz)*(cos(pitch)*sin(roll) - (cos(yaw)*sin(pitch)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) + (sin(pitch)*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch)) - (bgy - wy)*(cos(pitch)*cos(roll) + (cos(yaw)*sin(pitch)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch) - (sin(pitch)*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch)), - (bgz - wz)*(cos(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) - cos(roll)*sin(pitch) - sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + cos(roll)*cos(yaw)^2*sin(pitch) + cos(roll)*sin(pitch)*sin(yaw)^2 + (cos(yaw)*sin(pitch)^2*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch)^2 - (sin(pitch)^2*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch)^2) - (bgy - wy)*(sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - cos(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)) - sin(pitch)*sin(roll) + cos(yaw)^2*sin(pitch)*sin(roll) + sin(pitch)*sin(roll)*sin(yaw)^2 - (cos(yaw)*sin(pitch)^2*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch)^2 + (sin(pitch)^2*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch)^2) - (bgx - wx)*(cos(pitch)*cos(yaw)^2 + cos(pitch)*sin(yaw)^2 - cos(pitch)), (0), (0), (0), (0), - sin(pitch)*cos(yaw)^2 - sin(pitch)*sin(yaw)^2 + sin(pitch), (cos(yaw)*sin(pitch)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) - cos(pitch)*sin(roll) - (sin(pitch)*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch), (sin(pitch)*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*sin(pitch)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch) - cos(pitch)*cos(roll), (0), (0), (0);
    (0), (0), (0), (ay - bay)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) + (az - baz)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)), cos(pitch)*cos(roll)*cos(yaw)*(az - baz) - cos(yaw)*sin(pitch)*(ax - bax) + cos(pitch)*cos(yaw)*sin(roll)*(ay - bay), (az - baz)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) - (ay - bay)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - cos(pitch)*sin(yaw)*(ax - bax), (0), (0), (0), (0), (0), (0), -cos(pitch)*cos(yaw), (cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)), - sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
    (0), (0), (0), - (ay - bay)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) - (az - baz)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)), cos(pitch)*cos(roll)*sin(yaw)*(az - baz) - sin(pitch)*sin(yaw)*(ax - bax) + cos(pitch)*sin(roll)*sin(yaw)*(ay - bay), (az - baz)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) - (ay - bay)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)) + cos(pitch)*cos(yaw)*(ax - bax), (0), (0), (0), (0), (0), (0), -cos(pitch)*sin(yaw), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), (cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw));
    (0), (0), (0), cos(pitch)*cos(roll)*(ay - bay) - cos(pitch)*sin(roll)*(az - baz), - cos(pitch)*(ax - bax) - cos(roll)*sin(pitch)*(az - baz) - sin(pitch)*sin(roll)*(ay - bay), (0), (0), (0), (0), (0), (0), (0), sin(pitch), -cos(pitch)*sin(roll), -cos(pitch)*cos(roll);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0)];

% F = I + dt*At
Ft = eye(15,15) + dt* At;

%Ut = jacobian(f,n)
Ut = [(0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    - cos(yaw)^2 - sin(yaw)^2, (cos(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) - (sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch), (sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch), (0), (0), (0), (0), (0), (0), (0), (0), (0); (0), - cos(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - sin(yaw)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)), cos(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) + sin(yaw)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    - sin(pitch)*cos(yaw)^2 - sin(pitch)*sin(yaw)^2 + sin(pitch), (cos(yaw)*sin(pitch)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/cos(pitch) - cos(pitch)*sin(roll) - (sin(pitch)*sin(yaw)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/cos(pitch), (sin(pitch)*sin(yaw)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/cos(pitch) - (cos(yaw)*sin(pitch)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/cos(pitch) - cos(pitch)*cos(roll), (0), (0), (0), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), -cos(pitch)*cos(yaw), (cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)), - sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), -cos(pitch)*sin(yaw), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), (cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), sin(pitch), -cos(pitch)*sin(roll), -cos(pitch)*cos(roll), (0), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (1), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (1), (0), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (1), (0), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (1), (0), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (1), (0);
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (1)];

Vt  = Ut;

Q = eye(12,12) * dt;

%covariance update
covarEst = Ft* covarPrev * (Ft)' + Vt* Q * Vt';

end

