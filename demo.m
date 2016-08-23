clear all;
close all;
warning off;

%% load data:

load('synthetic.mat');
external_data.data=folds;
external_data.fs = folds(1).fs;     
external_data.sample_frequency=256;
external_data.eps = 0.125;
external_data.minlen = 0.5;
clear folds; 


%% settings
options.w_start=0.6;                    %Value of velocity Weight at the begining
options.w_end=0.3;                      %Value of velocity Weight at the end of the pso iterations
options.c1=2;                           %Cognitive Acceleration
options.c2=2;                           %Social Acceleration
options.c3=0;                           %Neighborhood Acceleration
options.w_varyfor=0.9;                  %The fraction of maximum iterations, for which w is linearly varied
options.Dim=6;                          %Dimensions of the problem
options.Chi=1  ;                        %Constriction factor
options.Nhood=2 ;                       %Neighborhood size (nhood=1 ==> 2 neighbors, one on each side)
options.Neighbor=0;                     %Use neighborhood acceleration (in addition to global acceleration)
options.SwarmSize=20;                   %Swarm Size 
options.Iterations=150;                 %Maximum Iterations`
options.ErrGoal=1e-10 ;                 %Error goal
options.f2eval='objectivefunction';     %Function/System to optimize
options.xlimits=[5 0.25 0 -100 -100 2; 56 3 1 100 100 20];
options.lb=options.xlimits(1,:);        %Lower bounds of Initialization
options.ub=options.xlimits(2,:);        %Upper bounds of Initialization
options.Vmax=options.ub-options.lb;
options.DispInt=1;                      %Display interval


results=mm_pso(options, external_data)




