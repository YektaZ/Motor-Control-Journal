
% compute the Minimum joint jerk model for all the subject parameters I have
% compute the Minimum torque model for all the subject parameters I have

%by YZ

clear all
clc

global sho_x sho_y

sho_x = 0;
sho_y = -0.4;

addpath('C:\GoogleDrive\Synced Folder\Matlab\Simple Two-Link Arm')
% addpath('C:\Users\Yeki\Google Drive\Synced Folder\Matlab\Simple Two-Link Arm')

addpath('C:\GoogleDrive\Synced Folder\Matlab\Minimum Jerk Functions');
% addpath('C:\Users\Yeki\Google Drive\Synced Folder\Matlab\Minimum Jerk Functions')

% Group = 'LimbFeedback'; 
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets'; 
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets';
% Subjects = {'NF04', 'NF05', 'NF06', 'NF07', 'NF08', 'NF15', 'NF16', 'NF30', 'NF31', 'NF32'}; 

% Group = 'CursorFeedback';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets_control'; 
% % % % Subjects = {'NF20', 'NF21', 'NF22', 'NF23', 'NF24', 'NF25', 'NF26', 'NF27', 'NF33', 'NF34'};
% Subjects = {'NF20', 'NF21', 'NF22', 'NF23', 'NF24', 'NF25', 'NF26', 'NF33', 'NF34'};

% Group = 'MixFeedback';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Target_mix';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Target_mix';
% Subjects = {'NF42', 'NF43', 'NF44', 'NF45', 'NF46', 'NF47', 'NF48', 'NF49', 'NF50', 'NF51'}; 

Nominal_dF = load('nominal_solution_dF_1sec.mat'); % 1s --> to get the targets

% a2 NF07 = ?
% a1_LFB = [25, 28, 26, 33, 25, 28, 29, 28, 32, 31]*0.01;
% a2_LFB = [28, 30, 28, 34, 26, 31, 32, 26, 28, 33]*0.01;

% a1_MFB = [34, 25, 27, 33, 24, 28, 30, 25, 31, 31]*0.01;
% a2_MFB = [27, 26, 28, 33, 26, 30, 31, 30, 33, 33]*0.01;

%%%%%%%%%%%%%%%%%
a1_allFB = [25, 28, 26, 33, 25, 28, 29, 28, 32, 31, 34, 25, 27, 33, 25, 28, 30, 25, 31, 31]*0.01;
a2_allFB = [28, 30, 28, 34, 26, 31, 32, 26, 28, 33, 27, 26, 28, 33, 26, 30, 31, 30, 33, 33]*0.01;

m1_allFB = [61, 72, 80, 54, 80, 104, 63, 60, 54, 70, 81, 51, 72, 100, 57, 56, 61, 63, 95, 80]*0.033;
m2_allFB = [61, 72, 80, 54, 80, 104, 63, 60, 54, 70, 81, 51, 72, 100, 57, 56, 61, 63, 95, 80]*0.021;


%% MJJ model

N = length(Nominal_dF.ParamTrack{1}.Sol.X);
t = Nominal_dF.ParamTrack{1}.Sol.t;

fig1 = figure;
hold on

for counter = 1:4
    for sub = 1:length(a1_allFB)
        
        a1 = a1_allFB(sub);
        a2 = a2_allFB(sub);
        %
        Xo = Nominal_dF.ParamTrack{counter}.Parameters.xstart(1:2);
        Xf = Nominal_dF.ParamTrack{counter}.Parameters.xfinish(1:2);
        %
        Xo_joint = inv_Position_2(Xo, a1, a2);
        Xf_joint = inv_Position_2(Xf, a1, a2);
        for i = 1:N
            temp_joint(:,i) = min_Jerk_Position( t(i), t(1), t(end), Xo_joint, Xf_joint);
            temp_vel_joint(:,i) = min_Jerk_Velocity( t(i), t(1), t(end), Xo_joint, Xf_joint); %?
            X_joint{counter}(:,i) = forward_Position_2(temp_joint(:,i), a1, a2);
            dX_joint{counter}(:,i) = forward_Velocity( temp_joint, temp_vel_joint(:,i), a1, a2 );
        end
        
        MJ_joint{sub, counter}.Sol.X = [X_joint{counter}; dX_joint{counter}]';
        MJ_joint{sub, counter}.Sol.t = t;
        
        figure(fig1)
        plot(X_joint{counter}(1,:),X_joint{counter}(2,:))
        pause(0.5)
        
    end
end

axis equal
title('MJJ with individual link lengths')
axis([-0.15, 0.15, 0.28, 0.52])

%% plot the old answer

for counter = 1:4
    
    a1 = Nominal_dF.ParamTrack{counter}.Parameters.a1;
    a2 = Nominal_dF.ParamTrack{counter}.Parameters.a2;
    
    Xo = Nominal_dF.ParamTrack{counter}.Parameters.xstart(1:2);
        Xf = Nominal_dF.ParamTrack{counter}.Parameters.xfinish(1:2);
        %
        Xo_joint = inv_Position_2(Xo, a1, a2);
        Xf_joint = inv_Position_2(Xf, a1, a2);
        for i = 1:N
            temp_joint(:,i) = min_Jerk_Position( t(i), t(1), t(end), Xo_joint, Xf_joint);
            temp_vel_joint(:,i) = min_Jerk_Velocity( t(i), t(1), t(end), Xo_joint, Xf_joint); %?
            X_joint{counter}(:,i) = forward_Position_2(temp_joint(:,i), a1, a2);
            dX_joint{counter}(:,i) = forward_Velocity( temp_joint, temp_vel_joint(:,i), a1, a2 );
        end
        
%         MJ_joint{counter}.Sol.X = [X_joint{counter}; dX_joint{counter}]';
%         MJ_joint{counter}.Sol.t = t;
        
        figure(fig1)
        plot(X_joint{counter}(1,:),X_joint{counter}(2,:),'--', 'LineWidth', 4)
        pause(0.5)
        
end
%% save the simulations

save('MJJ_20subs.mat','MJ_joint');

%% average across 20 subjects

X_sum = cell(4,1); 
Y_sum = cell(4,1); 
dX_sum = cell(4,1); 
dY_sum = cell(4,1);
X_avg = cell(4,1);
Y_avg = cell(4,1);
dX_avg = cell(4,1);
dY_avg = cell(4,1);

figure
hold on

for i = 1:4
    
    X_sum{i} = zeros(30,1);
    Y_sum{i} = zeros(30,1);
    dX_sum{i} = zeros(30,1);
    dY_sum{i} = zeros(30,1);
    
    for sub = 1:length(a1_allFB)
        
        X_sum{i} = X_sum{i} + MJ_joint{sub,i}.Sol.X(:,1);
        Y_sum{i} = Y_sum{i} + MJ_joint{sub,i}.Sol.X(:,2);
        dX_sum{i} = dX_sum{i} + MJ_joint{sub,i}.Sol.X(:,3);
        dY_sum{i} = dY_sum{i} + MJ_joint{sub,i}.Sol.X(:,4);
        
    end

    X_avg{i} = X_sum{i}./length(a1_allFB);
    Y_avg{i} = Y_sum{i}./length(a1_allFB);
    dX_avg{i} = dX_sum{i}./length(a1_allFB);
    dY_avg{i} = dY_sum{i}./length(a1_allFB);

    MJJ_AVG{i}.Sol.X = [X_avg{i}, Y_avg{i}, dX_avg{i}, dY_avg{i}];
    MJJ_AVG{i}.Sol.t = MJ_joint{1,1}.Sol.t;
    
    plot(X_avg{i}, Y_avg{i})
    
end

% MJJ_AVG.Sol.X = Sol;
% MJJ_AVG.Y_avg = Y_avg;

%% save the average

save('MJJ_AVG_20subs.mat','MJJ_AVG');


%% (TR) Min rate of torque solution

figure
hold on

for i = 1:4
    
    for sub = 1:length(a1_allFB)
        
        Parameters = Nominal_dF.ParamTrack{i}.Parameters;
        Parameters.a1 = a1_allFB(sub);
        Parameters.a2 = a2_allFB(sub);
        Parameters.m1 = m1_allFB(sub);
        Parameters.m2 = m2_allFB(sub);
        Parameters.a1_cm = 0.47*a1_allFB(sub);
        Parameters.a2_cm = 0.42*a2_allFB(sub);
        Parameters.I_1 = 1/12*m1_allFB(sub)*(a1_allFB(sub))^2;
        Parameters.I_2 = 1/12*m2_allFB(sub)*(a2_allFB(sub))^2;
        Parameters.sho_x = sho_x;
        Parameters.sho_y = sho_y;
        Parameters.R = eye(2); %two inputs to this system
        Parameters.Q = zeros(6);
        Parameters.Phi = zeros(6);
        
        [X,u,x_pos,y_pos, vel_pos] = TwoLinkArm_min_dT(Parameters);
        
        MT{i,sub}.Sol.X = [x_pos;y_pos;vel_pos]';
        MT{i,sub}.Sol.t = Nominal_dF.ParamTrack{i}.Sol.t';
        
        plot(x_pos,y_pos)
        pause(0.05)
        
        sub
    end
    i
end

axis equal
title('min Rate Torque for individual parameters')
axis([-0.2, 0.2, 0.2, 0.55])

%% plot the old answer


for i = 1:4
    
        Parameters = Nominal_dF.ParamTrack{i}.Parameters;
        Parameters.sho_x = sho_x;
        Parameters.sho_y = sho_y;
        Parameters.R = eye(2); %two inputs to this system
        Parameters.Q = zeros(6);
        Parameters.Phi = zeros(6);
        
        [X,u,x_pos,y_pos, vel_pos] = TwoLinkArm_min_dT(Parameters);
        
        plot(x_pos,y_pos,'--','LineWidth',3)
        
    i
end


%% save the simulations

save('TR_20subs.mat','MT');

%% average across them

X_sum = cell(4,1); 
Y_sum = cell(4,1); 
X_avg = cell(4,1);
Y_avg = cell(4,1);
dX_sum = cell(4,1); 
dY_sum = cell(4,1); 
dX_avg = cell(4,1);
dY_avg = cell(4,1);

figure
hold on

for i = 1:4
    
    X_sum{i} = zeros(30,1);
    Y_sum{i} = zeros(30,1);
    dX_sum{i} = zeros(30,1);
    dY_sum{i} = zeros(30,1);
    
    for sub = 1:length(a1_allFB)
        
        X_sum{i} = X_sum{i} + MT{i,sub}.Sol.X(:,1);
        Y_sum{i} = Y_sum{i} + MT{i,sub}.Sol.X(:,2);
        dX_sum{i} = dX_sum{i} + MT{i,sub}.Sol.X(:,3);
        dY_sum{i} = dY_sum{i} + MT{i,sub}.Sol.X(:,4);
        
    end

    X_avg{i} = X_sum{i}./length(a1_allFB);
    Y_avg{i} = Y_sum{i}./length(a1_allFB);
    dX_avg{i} = dX_sum{i}./length(a1_allFB);
    dY_avg{i} = dY_sum{i}./length(a1_allFB);
    
    TR_AVG{i}.Sol.X = [X_avg{i}, Y_avg{i}, dX_avg{i}, dY_avg{i}];
    TR_AVG{i}.Sol.t = MT{1,1}.Sol.t;
    
    plot(X_avg{i}, Y_avg{i})
    
end


%% save the average

save('TR_AVG_20subs.mat','TR_AVG');

%% plot individual paths vs computed path

%% %%%%%%%%%%%%%%%%%%%%%%  muscle model

figure
hold on
counter = 1;

for i = 1:4
    for sub = 1:length(a1_allFB)
        
        Parameters = Nominal_dF.ParamTrack{i}.Parameters;
        Parameters.a1 = a1_allFB(sub);
        Parameters.a2 = a2_allFB(sub);
        Parameters.m1 = m1_allFB(sub);
        Parameters.m2 = m2_allFB(sub);
        Parameters.a1_cm = 0.47*a1_allFB(sub);
        Parameters.a2_cm = 0.42*a2_allFB(sub);
        Parameters.I_1 = 1/12*m1_allFB(sub)*(a1_allFB(sub))^2;
        Parameters.I_2 = 1/12*m2_allFB(sub)*(a2_allFB(sub))^2;
        
        % this solution
        Sol = getOpt2LinkMuscleFFDirectMethod_RateOfForce_ineq_YZ_2(Parameters);
        
        ParamTrack{counter}.Sol = Sol;
        ParamTrack{counter}.Parameters = Parameters;
       
        plot(Sol.X(:,1), Sol.X(:,2))
        
        sub
        i
        counter = counter+1;
        
    end
    
end
axis([-0.2, 0.2, 0.2, 0.55])

%% just plot

load ParamTrack_20subs.mat

figure
hold on

for i = 1:length(ParamTrack)
    
    plot(ParamTrack{i}.Sol.X(:,1), ParamTrack{i}.Sol.X(:,2))
    
end

axis equal
title('FR')
axis([-0.2, 0.2, 0.2, 0.55])

%% saving

save('ParamTrack_20subs.mat','ParamTrack');

%% average across 20 subs

X_sum = cell(4,1); 
Y_sum = cell(4,1); 
X_avg = cell(4,1);
Y_avg = cell(4,1);
dX_sum = cell(4,1); 
dY_sum = cell(4,1); 
dX_avg = cell(4,1);
dY_avg = cell(4,1);

figure
hold on

for i = 1:4
    X_sum{i} = zeros(30,1);
    Y_sum{i} = zeros(30,1);
    dX_sum{i} = zeros(30,1);
    dY_sum{i} = zeros(30,1);
    
    for sub = 1:length(a1_allFB)
        
        X_sum{i} = X_sum{i} + ParamTrack{(i-1)*length(a1_allFB) + sub}.Sol.X(:,1);
        Y_sum{i} = Y_sum{i} + ParamTrack{(i-1)*length(a1_allFB) + sub}.Sol.X(:,2);
        dX_sum{i} = dX_sum{i} + ParamTrack{(i-1)*length(a1_allFB) + sub}.Sol.X(:,3);
        dY_sum{i} = dY_sum{i} + ParamTrack{(i-1)*length(a1_allFB) + sub}.Sol.X(:,4);
        
    end 
    X_avg{i} = X_sum{i}./length(a1_allFB);
    Y_avg{i} = Y_sum{i}./length(a1_allFB);
    dX_avg{i} = dX_sum{i}./length(a1_allFB);
    dY_avg{i} = dY_sum{i}./length(a1_allFB);
    
    FR_AVG{i}.Sol.X = [X_avg{i}, Y_avg{i}, dX_avg{i}, dY_avg{i}];
    FR_AVG{i}.Sol.t = ParamTrack{1}.Sol.t;
    
    plot(X_avg{i}, Y_avg{i})
end

axis equal
axis([-0.15, 0.15, 0.25, 0.55])
title('FR model average 20 subs')


%% save the average

save('FR_AVG_20subs.mat','FR_AVG');

%% use the average across sub parameters to compute one simulation

clc

a1_avg = mean(a1_allFB);
a2_avg = mean(a2_allFB); 

m1_avg = mean(m1_allFB);
m2_avg = mean(m2_allFB); 

N = length(Nominal_dF.ParamTrack{1}.Sol.X);
t = Nominal_dF.ParamTrack{1}.Sol.t;

%%%%%%%%%%%%%%%%% MJJ

fig1 = figure;
hold on


for counter = 1:4
    
    a1 = a1_avg;
    a2 = a2_avg;
    
    Xo = Nominal_dF.ParamTrack{counter}.Parameters.xstart(1:2);
    Xf = Nominal_dF.ParamTrack{counter}.Parameters.xfinish(1:2);
    %
    Xo_joint = inv_Position_2(Xo, a1, a2);
    Xf_joint = inv_Position_2(Xf, a1, a2);
    for i = 1:N
        temp_joint(:,i) = min_Jerk_Position( t(i), t(1), t(end), Xo_joint, Xf_joint);
        temp_vel_joint(:,i) = min_Jerk_Velocity( t(i), t(1), t(end), Xo_joint, Xf_joint); %?
        X_joint{counter}(:,i) = forward_Position_2(temp_joint(:,i), a1, a2);
        dX_joint{counter}(:,i) = forward_Velocity( temp_joint, temp_vel_joint(:,i), a1, a2 );
    end
    
    MJJ_paramAVG{counter}.Sol.X = [X_joint{counter}', dX_joint{counter}']; %Y_avg{i}, dX_avg{i}, dY_avg{i}];
    MJJ_paramAVG{counter}.Sol.t = t;
        
    figure(fig1)
    plot(X_joint{counter}(1,:),X_joint{counter}(2,:),'--', 'LineWidth', 4)
    pause(0.5)
    
end

axis equal
title('MJJ with avg sub link lengths')
axis([-0.2, 0.2, 0.2, 0.55])


%%%%%%%%%%%%%%%%%%% TR
figure
hold on

for i = 1:4
    
    Parameters = Nominal_dF.ParamTrack{i}.Parameters;
    Parameters.a1 = a1_avg;
    Parameters.a2 = a2_avg;
    Parameters.m1 = m1_avg;
    Parameters.m2 = m2_avg;
    Parameters.a1_cm = 0.47*a1_avg;
    Parameters.a2_cm = 0.42*a2_avg;
    Parameters.I_1 = 1/12*m1_avg*(a1_avg)^2;
    Parameters.I_2 = 1/12*m2_avg*(a2_avg)^2;
    Parameters.sho_x = sho_x;
    Parameters.sho_y = sho_y;
    Parameters.R = eye(2); %two inputs to this system
    Parameters.Q = zeros(6);
    Parameters.Phi = zeros(6);
    
    [X,u,x_pos,y_pos, vel_pos] = TwoLinkArm_min_dT(Parameters);
    
    TR_paramAVG{i}.Sol.X = [x_pos;y_pos;vel_pos]';
    TR_paramAVG{i}.Sol.t = Nominal_dF.ParamTrack{i}.Sol.t;
    
    plot(x_pos,y_pos,'--','LineWidth',3)
    
    i
end

axis equal
title('TR with avg subs link lengths')
axis([-0.2, 0.2, 0.2, 0.55])

%%%%%%%%%%%%%%%%%% FR
figure
hold on

for i = 1:4
        
        Parameters = Nominal_dF.ParamTrack{i}.Parameters;
        Parameters.a1 = a1_avg;
        Parameters.a2 = a2_avg;
        Parameters.m1 = m1_avg;
        Parameters.m2 = m2_avg;
        Parameters.a1_cm = 0.47*a1_avg;
        Parameters.a2_cm = 0.42*a2_avg;
        Parameters.I_1 = 1/12*m1_avg*(a1_avg)^2;
        Parameters.I_2 = 1/12*m2_avg*(a2_avg)^2;
        
        % this solution
        Sol = getOpt2LinkMuscleFFDirectMethod_RateOfForce_ineq_YZ_2(Parameters);
        
        FR_paramAVG{i}.Sol.X = Sol.X;
        FR_paramAVG{i}.Sol.t = Sol.t;
    
        plot(Sol.X(:,1), Sol.X(:,2))

        
end

axis equal
axis([-0.2, 0.2, 0.2, 0.55])
title('FR with avg sub link lenghts')

%% saving

%% THE END

