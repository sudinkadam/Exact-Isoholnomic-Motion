function  states = KinematicsPlanarPurcell(Init,Control)
global steps k L stepLength

%Init.stateG = [1;0;0;0;0];
%Init.stateB = [1;0;1;0];
Init = [Init.stateB; Init.stateG];
prev_alpha1=Init(1);
prev_alpha2=Init(2);
theta=Init(3);
R1 = [cos(prev_alpha1),      -sin(prev_alpha1);
    sin(prev_alpha1),       cos(prev_alpha1)];
R2 = [cos(prev_alpha2),      -sin(prev_alpha2);
    sin(prev_alpha2),       cos(prev_alpha2)];
G = [cos(theta),    -sin(theta),    Init(4);
    sin(theta),     cos(theta),     Init(5);
    0,              0,        1];

% R1 = [Init.stateB(1),      -Init.stateB(2);
%       Init.stateB(2),       Init.stateB(1)];
% 
% R2 = [Init.stateB(3),      -Init.stateB(4);
%       Init.stateB(4),       Init.stateB(3)];
% 
% G = [Init.stateG(1),    -Init.stateG(2),     Init.stateG(3);
%      Init.stateG(2),     Init.stateG(1),     Init.stateG(4);
%                   0,                  0,     1]; 

%% Initialization
q_next=[R1,              zeros(2,2),     zeros(2,3);
        zeros(2,2),              R2,     zeros(2,3);
        zeros(3,2),       zeros(3,2),            G];
% states(:,1)=[q_next(1,1); q_next(2,1); q_next(3,3); q_next(4,3); q_next(5,5); q_next(6,5); q_next(5,7);  q_next(6,7)];
states(:,1)=[prev_alpha1; prev_alpha2; theta; q_next(5,7);  q_next(6,7)];
%% Control sequence for pure X
control_sequence = Control;
%% Simulation
for row=1:steps
    u = [control_sequence(1,row); control_sequence(2,row)];
    [CF,q_next]=state_rhs01(q_next,u);
    alpha1=atan2(q_next(2,1),q_next(1,1));
    alpha2=atan2(q_next(4,3),q_next(3,3));
    theta=atan2(q_next(6,5),q_next(5,5));
    states(:,row+1)=[alpha1; alpha2; theta; q_next(5,7);  q_next(6,7)];
end

