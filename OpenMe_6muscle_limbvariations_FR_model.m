% by YZ

%
% code to manage and run the direct method optimizer code for FF, but with
% a model that uses commands to a series of muscles, with constant moment
% arms
%
% tau = R*F, where F is a vector of muscle forces, and R is a moment arm
% matrix
%

% LOOK AT BOX REACHES W/ 4 MUSCLES, NO NEGATIVE COMMANDS

% 6 muscles
% 1 - monoartiular shoulder flexor
% 2 - biarticular shoulder/elbow flexor
% 3 - monoarticular elbow flexor
% 4 - monoarticular shoulder extensor
% 5 - biarticular shoulder/elbow extensor
% 6 - monoarticular elbow extensor



close all
clear all
clc

%% BOX TARGETS

global a1 a2 m1 m2 a1_cm a2_cm I_1 I_2

BoxTargets = [.1 -.1; .1 .2; -.1 .2; -.1 -.1; .1 -.1];
numTargs = length(BoxTargets) - 1;
    
panels = [6,2,4,8];

%% limb parameters

% a1 = .35;
% a2 = .36;
% m1 = 2.1;
% m2 = 1.65;
% a1_cm = 0.152;
% a2_cm = 0.2465;
I_1 = .0264;
I_2 = .0472;

steps = 2;
m_range = [linspace(1.123,2.580,steps); linspace(1.017,2.370,steps)]; % from the book by Charles Clauser
l_range = [linspace(0.27, 0.32,steps); linspace(0.27, 0.33,steps)];
% r_range = [linspace(42.7/100*l_range(1,1),50.9/100*l_range(1,2),steps);...
%     linspace(44.4/100*l_range(2,1), 67.7/100*l_range(2,2),steps)]; % from the same book
% II_range = [linspace(0.020, 0.038,steps); linspace(0.022, 0.048, steps)]; % from numbers I have seen

%% LOOK AT BOX REACHES W/ 2 MUSCLES, NO BOUNDS



%% LOOK AT BOX REACHES W/ 4/6 MUSCLES, NO NEGATIVE COMMANDS?

% clear params
% 
% % suppose we have 3 muscles, one for shoulder, one for elbow, and one for
% % eblow-shoulder forces
% 
% if 0
%     comm_dim = 4;
%     Marm = [1, -1, 0, 0; 0, 0, 1, -1];
%     %
% else
%     % positive shoulder, negative shoulder, positive cross, neg cross, pos
%     % elbow, negative elbow
%     comm_dim = 6;
%     Marm = [1, -1, 1, -1, 0, 0;
%             0,  0, 1, -1, 1, -1];
% end

comm_dim = 6;
Marm0 = [1, -1, 1, -1, 0, 0;
    0,  0, 1, -1, 1, -1];
Marm = Marm0;

params.T = .75;
params.xstart = zeros(4,1);
params.xfinish = zeros(4,1);
params.comm_dim = comm_dim;
% we can either specify the number of points, or the dt, but not both
params.numpts = 75;
%params.dt = 0.01;
%
params.Q = zeros(4);
params.R = .1*eye(comm_dim);
params.Phi = 100*eye(4);
params.PLOT = false;
params.Curl = zeros(2);
%
% command bounds
params.Umin = 0; %-1;
params.Umax = 5;

xo = 0; %.2;
yo = .2;

N = length(m_range);

counter = 1;

for j1 = 1:N
    m1 = m_range(1,j1);
    for j2 = 1:N
        m2 = m_range(2,j2);
        for j3 = 1:N
            a1 = l_range(1,j3);
            a1_cm = 0.47*a1;
            for j4 = 1:N
                a2 = l_range(1,j4);
                a2_cm = 0.56*a2;
                
                %                 for i = 1:length(Marm0)
                %
                %                     Marm = Marm0;
                %
                %                     temp = Marm(:,i);
                %                     if temp(1)~=0
                %                         Marm(1,i) = sign(Marm0(1,i))*3;
                %                     end
                %                     if temp(2)~=0
                %                         Marm(2,i) = sign(Marm0(2,i))*3;
                %                     end
                for j5 = [1,4]
                    Marm(1,1) = j5;
                    for j6 =[1,4]
                        Marm(1,2) = -j6;
                        for j7 = [1,4]
                            Marm(:,3) = j7;
                            for j8 = [1,4]
                                Marm(:,4) = -j8;
                                for j9 = [1,4]
                                    Marm(2,5) = j9;
                                    for j10 = [1,4];
                                        Marm(2,6) = -j10;
                                        
                                        
                                        params.MomentArms = Marm;
                                        
                                        for targ = 1:length(BoxTargets)-1
                                            
                                            Xo = BoxTargets(targ,:)'+[xo; yo];
                                            Xf = BoxTargets(targ+1,:)'+[xo; yo];
                                            params.xstart = [Xo; 0; 0];
                                            % this target
                                            params.xfinish = [Xf; 0; 0];
                                            % this solution
                                            Sol =   getOpt2LinkMuscleFFDirectMethod_ineq(params);
                                            
                                            ParamTrack{counter}.Marm = Marm;
                                            ParamTrack{counter}.m1 = m1;
                                            ParamTrack{counter}.m2 = m2;
                                            ParamTrack{counter}.a1 = a1;
                                            ParamTrack{counter}.a2 = a2;
                                            ParamTrack{counter}.Xo = Xo;
                                            ParamTrack{counter}.Xf = Xf;
                                            ParamTrack{counter}.Sol = Sol;
                                            
                                            counter = counter + 1
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
            end
            
        end
    end
end

%% saving

save('ParamTrack2.mat','ParamTrack');

%% plot
close all

num_boxes = length(ParamTrack);

counter = 1;
disp('Click on the figure to exit')
FIG = figure;

for fig = 1:(num_boxes/16)
    
    FIG;
    clf;
    
    for box = 1:4
        
        subplot(2,2,box)
        hold on
        title([num2str(counter),' - ',num2str(counter+3)])
        
        plot(ParamTrack{counter}.Sol.X(:,1),ParamTrack{counter}.Sol.X(:,2))
        counter = counter + 1;
        
        plot(ParamTrack{counter}.Sol.X(:,1),ParamTrack{counter}.Sol.X(:,2))
        counter = counter + 1;
        
        plot(ParamTrack{counter}.Sol.X(:,1),ParamTrack{counter}.Sol.X(:,2))
        counter = counter + 1;
        
        plot(ParamTrack{counter}.Sol.X(:,1),ParamTrack{counter}.Sol.X(:,2))
        counter = counter + 1;
        
    end
    
    w = waitforbuttonpress;
    if w == 0
        break
    end
    
end


%% END CODE
