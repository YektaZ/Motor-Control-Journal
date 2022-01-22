% Two link arm model 
% with optimal torque

% :: X = [theta1, theta2, dtheta1, dtheta2, T1, T2]
% :: u = [dT1; dT2]

% cost = u'*R*u

% T + J'*F= I(theta)*dd_theta + C(theta,d_theta)*theta 


function [X,u,x_pos,y_pos, vel_pos] = TwoLinkArm_min_dT(Parameters) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters

n = 6; % number of states
rr = 2; % number of inputs 

a1 = Parameters.a1;
a2 = Parameters.a2;
l = [a1;a2];
%
X0_pos = Parameters.xstart - [Parameters.sho_x, Parameters.sho_y]; 
Xf_pos = Parameters.xfinish - [Parameters.sho_x, Parameters.sho_y]; 

X0 = [inv_Position(X0_pos,l);
    0;0;0;0];
Xf = [inv_Position(Xf_pos,l);
    0;0;0;0];
%
R = Parameters.R; 
Q = Parameters.Q;
Phi = Parameters.Phi;
%
numpts = Parameters.numpts; 
tf = Parameters.T;
t0 = 0;
t = linspace(t0,tf,numpts);
%
Parameters.n = n;
Parameters.rr = rr;
Parameters.X0 = X0;
Parameters.Xf = Xf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BVP

% init = [X0; 0; 0; 0; 0];
% SOLINIT = bvpinit(t,init);

% load a better guess
% D = load('SOLINIT_T.mat');
% SOLINIT.x = t;
% SOLINIT.y = D.Z{block_num};

if isfield (Parameters,'SOLINIT')
    D = Parameters.SOLINIT;
    SOLINIT.x = t;
    SOLINIT.y = D;
else
    init = [X0; zeros(n,1)];
    SOLINIT = bvpinit(t,init);
end

Sol = bvp4c(@(t,Z)Sys(t,Z,Parameters), @(z0,zf)BCFUN(z0,zf,Parameters), SOLINIT);

Z = deval(Sol,t);

X = Z(1:n,:)';
Lam = Z(n+1:2*n,:)';

for i = 1:numpts
    
    [~,~, df_du] = TwoLink_Dynamics(t(i),Z(:,i),zeros(6,1),Parameters);
    Parameters.df_du = df_du;
%     u(:,i) = Input_Uineqcon_YZ(t(i),Z(:,i),Parameters);
    u(:,i) = Input(t(i),Z(:,i),Parameters);
    
end

for i = 1:length(X)
    
    temp = forward_Position(X(i,1:2),l);
    %
    x_pos(i) = temp(1); % + Parameters.sho_x;
    y_pos(i) = temp(2); % + Parameters.sho_y;
    
    J = [-l(1)*sin(X(i,1))-l(2)*sin(X(i,1)+X(i,2)) -l(2)*sin(X(i,1)+X(i,2));
    l(1)*cos(X(i,1))+l(2)*cos(X(i,1)+X(i,2)) l(2)*cos(X(i,1)+X(i,2))];
    
    vel_pos(:,i) = J*X(i,3:4)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot

if Parameters.PLOT
% 
figure
plot(u')
title('command: rate of torque')
legend('dTx','dTy')

figure
hold on

plot(Parameters.xstart(1), Parameters.xstart(2),'k*')
plot(Parameters.xfinish(1), Parameters.xfinish(2),'k*')
plot(x_pos,y_pos,'k')
plot(Parameters.sho_x,Parameters.sho_y,'ro')

grid
title('Trajectory in XY cordinate')
axis equal

figure
plot(t,X')
legend('theta1','theta2','theta-dot1','theta-dot2','T1','T2')

end


function [dXdt, df_dX, df_du] = TwoLink_Dynamics(t,Z,u,Param)

m = [Param.m1; Param.m2];
n = Param.n;
l = [Param.a1; Param.a2];
r = [Param.a1_cm; Param.a2_cm];
II = [Param.I_1; Param.I_2];
%
X = Z(1:n);
%
gama = m(1)*r(1)^2 + m(2)*r(2)^2 + m(2)*l(1)^2 + sum(II);
alpha = 2*m(2)*r(2)*l(1);

I = [ gama+ alpha*cos(X(2)), m(2)*r(2)^2+ II(2)+ 1/2*alpha*cos(X(2));
    m(2)*r(2)^2+ II(2)+ 1/2*alpha*cos(X(2)), m(2)*r(2)^2+ II(2)];


sigma = -l(1)^2*m(2)^2*r(2)^2*cos(X(2))^2 +...
l(1)^2*m(2)^2*r(2)^2+...
II(2)*l(1)^2*m(2)+ ...
m(1)*m(2)*r(1)^2*r(2)^2+...
II(1)*m(2)*r(2)^2+...
II(2)*m(1)*r(1)^2+...
II(1)*II(2);

I_inv = 1/sigma*[m(2)*r(2)^2+ II(2), -(m(2)*r(2)^2+ II(2)+ 1/2*alpha*cos(X(2)));
    -(m(2)*r(2)^2+ II(2)+ 1/2*alpha*cos(X(2))), gama+ alpha*cos(X(2))];

% Compute dI_inv_dx2

sigma2 = -l(1)^2*m(2)^2*r(2)^2*cos(X(2))^2 +...
    l(1)^2*m(2)^2*r(2)^2+...
    II(2)*l(1)^2*m(2)+ ...
    m(1)*m(2)*r(1)^2*r(2)^2+...
    II(1)*m(2)*r(2)^2+...
    II(2)*m(1)*r(1)^2+...
    II(1)*II(2);

sigma1 = 1/sigma2*(l(1)*m(2)*r(2)*sin(X(2))) +...
    1/(sigma2^2)*( 2*l(1)^2*m(2)^2*r(2)^2*cos(X(2))*sin(X(2))*( m(2)*r(2)^2 + l(1)*m(2)*r(2)*cos(X(2))+ II(2)) );


dI_inv_dx2 = [-1/sigma2^2*(l(1)^2*m(2)^2*r(2)^2)*cos(X(2))*sin(X(2))*2*(m(2)*r(2)^2 + II(2)),...
    sigma1;
    sigma1,...
    -1/sigma2*(2*l(1)*m(2)*r(2)*sin(X(2)))-1/(sigma2)^2*( l(1)^2*m(2)^2*r(2)^2*cos(X(2))*sin(X(2))*2*(m(2)*l(1)^2+2*m(2)*cos(X(2))*l(1)*r(2)+ m(1)*r(1)^2+ m(2)*r(2)^2+ II(1)+II(2)) )];

%%%%%%%%%%%%

J = [-l(1)*sin(X(1))-l(2)*sin(X(1)+X(2)) -l(2)*sin(X(1)+X(2));
    l(1)*cos(X(1))+l(2)*cos(X(1)+X(2)) l(2)*cos(X(1)+X(2))];

h = -m(2)*l(1)*r(2)*sin(X(2));

C = [h*X(4) h*(X(3)+X(4));
    -h*X(3) 0];

M = X(5:6) -C*[X(3);X(4)];

%%%%%%%%%%% PARTIALS

dh_dx2 = -m(2)*l(1)*r(2)*cos(X(2));

dC_dx2 = [dh_dx2*X(4), dh_dx2*(X(4)+X(3));
    -dh_dx2*X(3), 0];

dC_dx3 = [0, h;
    -h, 0];

dC_dx4 = [h, h;
    0, 0];

%%%%%%%%%%%%%
dM_dx1 =  zeros(2,1); 

dM_dx2 =  - dC_dx2*[X(3);X(4)];

dM_dx3 =  - dC_dx3*[X(3);X(4)] - C*[1; 0];

dM_dx4 = - dC_dx4*[X(3);X(4)] - C*[0; 1]; 

dM_dx5 = [1;0];

dM_dx6 = [0;1];

%%%%%%%%%%%%% 

Si = [inv(I)*dM_dx1, dI_inv_dx2*M+inv(I)*dM_dx2, inv(I)*dM_dx3, inv(I)*dM_dx4, inv(I)*dM_dx5, inv(I)*dM_dx6];

% Compute the partials
df_du = [zeros(n-2,2)
    eye(2)];

df_dX = [0 0 1 0 0 0;
    0 0 0 1 0 0
    Si;
    zeros(2,n)];

dXdt = [X(3);X(4);
    inv(I)*M;
    u(1); u(2)];

function dZdt = Sys(t,Z,Param)

Q = Param.Q;
n = Param.n;
rr = Param.rr;

X = Z(1:n);
Lam = Z(n+1:2*n);

[~,~, df_du] = TwoLink_Dynamics(t,Z,zeros(rr,1),Param);
Param.df_du = df_du;
% u = Input_Uineqcon_YZ(t,Z,Param);
u = Input(t,Z,Param);

[dXdt,df_dx,~] = TwoLink_Dynamics(t,X,u,Param);

dLam = -Q'*X -df_dx'*Lam;

dZdt = [dXdt;dLam];

function u = Input(t,Z,Param)

R = Param.R;
n = Param.n;

X = Z(1:n);
Lam = Z(n+1:end);

[~, ~,df_du] = TwoLink_Dynamics(t,Z,[0;0],Param);

u = -inv(R)*df_du'*Lam;

function X_joint = inv_Position(X_car,a)

x = X_car(1);
y = X_car(2);

a1 = a(1);
a2 = a(2);

D = (x^2 + y^2 - a1^2 - a2^2)/(2*a1*a2);

X_joint(2) = atan2(abs(sqrt(1-D^2)),D);
X_joint(1) = atan2(y,x) - atan2(a2*sin(X_joint(2)),(a1 + a2*cos(X_joint(2))));

X_joint = X_joint';

function X_car = forward_Position( X_joint, a )

a1 = a(1);
a2 = a(2);

X_car(1) = a1 * cos(X_joint(1)) + a2 * cos(X_joint(1) + X_joint(2));

X_car(2) = a1 * sin(X_joint(1)) + a2 * sin(X_joint(1) + X_joint(2));

X_car = X_car';

function Error = BCFUN(Z0,Zf,Param)

n = Param.n;
Phi = Param.Phi;
X0 = Param.X0;
Xf = Param.Xf;

Error(1:n) = X0 - Z0(1:n);
Error(n+1:2*n) = Xf - Zf(1:n);
% Error(n+1:2*n) = -Phi*(Zf(1:n) - Xf) + Zf(n+1:2*n);

