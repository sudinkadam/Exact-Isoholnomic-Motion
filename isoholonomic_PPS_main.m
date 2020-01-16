%% Solve Optimal Control Problem using CasAdi

clear; clc;
tic
import casadi.*

k = 1; % viscous drag coefficient
l = 0.1;
a = 0.25/100;
mu=0.95;
k = 4*pi*mu/(log(l/a));
l = 1; % Length of all the links
T = 4.0;          % Time horizon
% global dT  
dT = 0.02;
N = ceil(T/dT); % number of control intervals

alpha_dot_max = 80*pi/180;     %Max limb rotation speed in rad/s

StDim = 5;
ConDim = 2;
ControlUpBnd = [alpha_dot_max;  alpha_dot_max];
ControlLwBnd = -[alpha_dot_max;  alpha_dot_max];
StateUpBnd = [inf;inf;inf;inf;inf];
StateLwBnd = -[inf;inf;inf;inf;inf];

% Initial and final states 
% IntSt = [0; 0; 0; 0; 0];
% FnlSt = [0.001; -0.001; 0; 0.1; 0];
init_limb_angle=55*pi/180;
IntSt = [init_limb_angle; -init_limb_angle; 0; 0; 0];
FnlSt = [init_limb_angle; -init_limb_angle; 0.02686; 0; 0];
%% Symbolic equations

% Declare model variables
x = SX.sym('x',[StDim, 1]);
u = SX.sym('u',[ConDim, 1]);

% Model equations
CF = return_connection_Gutman(x,k,l);
xdot = [eye(2); CF]*u;

% Return the the power matrix term based on the expression in 
% Wiezel O and Or Y 2016 Optimization and small-amplitude analysis of Purcell’s three-link microswimmer model Proc. R. Soc. A 47220160425
% https://royalsocietypublishing.org/doi/10.1098/rspa.2016.0425
W = return_power_matrix(x,k,l);
%W=eye(2);
% %% Gutman power
gutman_control = [-2*55*pi/180;0];
gutman_x=[55*pi/180;0;0;0;0];
gutman_dt=0.01;
gutman_count=1;
for gutman_t=0:gutman_dt:1
    gutman_power=gutman_control'*return_power_matrix(gutman_x,k,l)*gutman_control;
    gutman_power_array(1,gutman_count)=gutman_power;
    gutman_x(1) = gutman_x(1)+gutman_control(1)*gutman_dt
    gutman_count=gutman_count+1;
end
% % plot(tgrid, gutman_power_array)
% % grid on;
% % ylim([0 100])
% % ylabel('gutman_power(W)')
% % xlabel('Time (sec)')
% % legend('Power')

%%
% Objective term 
% L = u(1)^2 + u(2)^2;
% L=0;
L = u'*W*u;

% Continuous time dynamics
% Formulate discrete time dynamics
% Fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = T/N/M;
f = Function('f', {x, u}, {xdot, L});
X0 = MX.sym('X0', StDim);
U = MX.sym('U', ConDim);
X = X0;
Q = 0;
%    T=0;
for j=1:M
    [k1, k1_q] = f(X, U);
    [k2, k2_q] = f(X + DT/2 * k1, U);
    [k3, k3_q] = f(X + DT/2 * k2, U);
    [k4, k4_q] = f(X + DT * k3, U);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});

% Evaluate at a test point
Fk = F('x0', IntSt, 'p', zeros(ConDim, 1));
disp(Fk.xf)
disp(Fk.qf)

% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g = {};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0', StDim);
w = {w{:}, Xk};
lbw = [lbw; IntSt];
ubw = [ubw; IntSt];
w0 = [w0; IntSt];
% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], ConDim);
    w = {w{:}, Uk};
    lbw = [lbw; ControlLwBnd];
    ubw = [ubw; ControlUpBnd];
    w0 = [w0;  zeros(ConDim,1)];

    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'p', Uk);
    Xk_end = Fk.xf;
    J=J+Fk.qf;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], StDim);
    w = [w, {Xk}];
    if k<N-1
%         lbw = [lbw; [-inf*ones(2,1);-ones(2,1);-inf*ones(4,1)]];
%         ubw = [ubw; [+inf*ones(2,1);+ones(2,1);+inf*ones(4,1)]];
%         w0 = [w0; [zeros(2,1);1;0;zeros(4,1)]];
        lbw = [lbw; StateLwBnd];
        ubw = [ubw; StateUpBnd];
        w0 = [w0; zeros(StDim,1)];
    else
        % Final Condition on States
        lbw = [lbw; FnlSt];
        ubw = [ubw;  FnlSt];
        w0 = [w0; FnlSt];
    end

    % Add equality constraint
    g = [g, {Xk_end-Xk}];
    lbg = [lbg; zeros(StDim,1)];
    ubg = [ubg; zeros(StDim,1)];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
w_opt = [w_opt; nan*ones(ConDim, 1)];
OptStCon = reshape(w_opt, StDim+ConDim, N+1);

%% Save optimal state-action trajectory to .csv file

% TrajIndex = (1:N)';
% Traj = OptStCon(:,TrajIndex)';
% OptTraj = [TrajIndex,Traj];
% % TrajHeader={'tInd','x','y','Ctheta','Stheta','alpha','alphaDot','Phi_Dot_l','Phi_Dot_r','ur','ul'};
% % csvwrite_with_headers('OptTraj.csv',OptTraj,TrajHeader)
% 
% Save data to MATLAB file

 Opt_St = OptStCon(1:StDim, :);
 Opt_Con = OptStCon(StDim+1:StDim+ConDim, :); 
 tgrid = linspace(0, T, N+1);
% save('Man.mat', 'tgrid','Opt_St','Opt_Con','IntSt','FnlSt')

%% Plot the solution
alpha1 = Opt_St(1,:);
alpha3 = Opt_St(2,:);
x = Opt_St(3,:);
y = Opt_St(4,:);
theta = Opt_St(5,:);

%% Plot the solution
% %colorSt=['k','b','c','m','r','g','y','k'];
% theta = Opt_St(3,:);
% thrust = 0.1*Opt_Con(1,:);
% % quiverX = 0.02*thrust.*sin(theta);
% % quiverY = 0.02*thrust.*cos(theta);
% quiverX = 0.02*sin(theta);
% quiverY = 0.02*cos(theta);
% quiverX1 = 0.01*Opt_St(4,:);
% quiverY1 = 0.01*Opt_St(5,:);
% figure
% plot(Opt_St(1,:),Opt_St(2,:),'*'); hold on
% quiver(Opt_St(1,:),Opt_St(2,:),quiverX,quiverY,'ShowArrowHead','off')
% quiver(Opt_St(1,:),Opt_St(2,:),quiverX1,quiverY1,'ShowArrowHead','off')
% xlabel('x'); ylabel('z')
% grid on
% axis equal

Init_Graphics_Style

subplot(2,2,1)
plot(tgrid, Opt_St(1,:)*180/pi,tgrid, Opt_St(2,:)*180/pi)
grid on
legend('alpha1','alpha2')
title('Limb angles (deg)')

subplot(2,2,2)
plot(tgrid, Opt_Con(1,:)*180/pi,tgrid, Opt_Con(2,:)*180/pi)
grid on
legend('u1','u2')
title('Limb rates/controls (deg/s)')

subplot(2,2,3)
plot(tgrid, Opt_St(3,:),tgrid, Opt_St(4,:))
grid on
legend('X','Y')
title('Translational position (m)')

subplot(2,2,4)
plot(tgrid, Opt_St(5,:)*180/pi)
grid on
legend('theta')
title('Angular position (deg)')

%% Individual plots
figure1=figure(1);
plot(tgrid, Opt_St(1,:)*180/pi,tgrid, Opt_St(2,:)*180/pi)
grid on
legend('\alpha_1','\alpha_2')
ylabel('Limb angles (deg)')
xlabel('Time (sec)')
% title('Limb angles (deg)')
saveas(figure1,'Limb_angles.jpg')

figure2=figure(2);
plot(tgrid, Opt_Con(1,:)*180/pi,tgrid, Opt_Con(2,:)*180/pi)
grid on
legend('u_1','u_2')
ylabel('Limb velocities (deg/s)')
xlabel('Time (sec)')
% title('Limb rates/controls (deg/s)')
saveas(figure2,'Limb_rates.jpg')

figure3=figure(3);
plot(tgrid, Opt_St(3,:),tgrid, Opt_St(4,:))
grid on
ylabel('Translational position (m)')
xlabel('Time (sec)')
axis([0 4 -0.3 0.4])
legend('x position','y position')
% title('Translational position (m)')
saveas(figure3,'translational_position.jpg')

figure4=figure(4);
plot(tgrid, Opt_St(5,:)*180/pi)
grid on
ylabel('Rotational position (deg)')
xlabel('Time (sec)')
legend('\theta')
% title('Angular position (deg)')
saveas(figure4,'angular_position.jpg')

%% Plot power
for count=1:length(Opt_St)-1
    power=Opt_Con(:,count)'*return_power_matrix(Opt_St(:,count),k,l)*Opt_Con(:,count);
    power_array(1,count)=power;
    count=count+1;
end
plot(tgrid, power_array)
grid on;
ylim([0 100])
ylabel('Power(W)')
xlabel('Time (sec)')
legend('Power')
saveas(figure4,'power.jpg')
%% Data for Tikz
%data_to_save=[tgrid',Opt_St(1,:)'*180/pi,Opt_St(2,:)'*180/pi,...
%              Opt_St(3,:)',Opt_Con(1,:),Opt_Con(2,:)'*180/pi, power_array'];

%Data = table(tgrid',Opt_St(1,:)'*180/pi,Opt_St(2,:)'*180/pi, ...
%            Opt_St(3,:), Opt_St(4,:), Opt_St(5,:)'*180/pi,Opt_Con(1,:)'*180/pi,Opt_Con(2,:)'*180/pi, power_array');
Data = table(tgrid',Opt_St(1,:)'*180/pi,Opt_St(2,:)'*180/pi, ...
            Opt_St(3,:), Opt_St(4,:), Opt_St(5,:)'*180/pi,Opt_Con(1,:)'*180/pi,Opt_Con(2,:)'*180/pi, power_array');

writetable(Data,'./Data.txt')

toc