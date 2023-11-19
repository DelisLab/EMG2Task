%%To begin we must pre-process the EMG data, removing noise and making it
%%ready for analysis

%A cell array of EMG recordings from each trial with muscles as columns and
%timepoints as rows in each cell
EMG=%....insert cell array here

%Low-pass filter and rectify the EMG recordings within each trial
order=%insert the filter order
fs=%...insert the sampling frequency
lp_freq=%...insert the frequency which all activity above this will be removed

[b,a]=butter(order,lp_freq/(fs/2),'low');
for i=1:length(EMG)
    EMG{i}=filtfilt(b,a,abs(EMG{i}));

end

%To quantify spatiotemporal modules across trials, we need to align the
%trial durations

len=%...insert a representative time-series length of all trials (e.g. 1000)
t_ = linspace(0, 1, len);
for i=1:length(EMG)
    t = linspace(0,1,length(EMG{i}));
    EMG{i}=interp1(t, EMG{i}, t_, 'spline');
end

%Integrate the signal in uniform steps to make it a more manageble length
step=%...desired number of steps to integrate (e.g. if the timeseries is 1000 timepoints long,
    %a suitable step length could be 20 steps, resulting in a 50 point
    %timeseries

EMG_final=[];
for i=1:length(EMG)
    x=[];
    for ii=1:step:len
        x=[x;trapz(EMG(ii:ii+(step-1),:))];
    end
    EMG_final=cat(3,EMG_final,x);
end

%%The EMG pre-processing is now complete, you now need to specify some
%%discrete variables that describe some feature of the task across trials

Task=%...insert a discrete task variable equal in length to the number of trials (i.e. the 3rd dimension of EMG_final)
    %The variable must consist of integer values starting from the value 0

%Input the EMG and discrete task variables into the following function to
%get task-relevant muscle couplings anc charecterise them as either
%redundant or synergistic

[Task_MI,ST_JOINT,ST_R,ST_S,ST_JOINT_netM,ST_JOINT_netT,ST_R_netM,...
    ST_R_netT,ST_S_netM,ST_S_netT]=SpaceTime_II(EMG_final,Task);

%We are interested in the following outputs: 
%ST_R: The matrix of task-redundant muscle couplings
%ST_S: The matrix of task-synergistic muscle couplings
%ST_R_netM & ST_R_netT: The task-redundant spatial and temporal multiplex
                        %networks
%ST_S_netM & ST_S_netT:The task-synergistic spatial and temporal multiplex
                        %networks

%To identify the number of components to extract, we use a community
%detection protocol:

[Opt_rank_R_SP]=Consensus(ST_R_netM,1);
[Opt_rank_R_TM]=Consensus(ST_R_netT,1);

[Opt_rank_S_SP]=Consensus(ST_S_netM,1);
[Opt_rank_S_TM]=Consensus(ST_S_netT,1);

%Then extract out the Space-Time modules from the task-redundant
%interactions where Wi_R, Wb_R and Acal_R are the temporal and spatial
%modules and activation coefficients respectively.

[Wi_R,Acal_R,Wb_R,TS,VAF,E2]=sNM3F_basic(reshape(ST_R,[size(ST_R,2),size(ST_R,1)*size(ST_R,3)])',Opt_rank_R_SP,Opt_rank_R_TM,size(ST_R,3));


%Then extract out the Space-Time modules from the task-synergistic
%interactions where Wi_S, Wb_S and Acal_S are the temporal and spatial
%modules and activation coefficients respectively.

[Wi_S,Acal_S,Wb_S,TS,VAF,E2]=sNM3F_basic(reshape(ST_S,[size(ST_S,2),size(ST_S,1)*size(ST_S,3)])',Opt_rank_S_SP,Opt_rank_S_TM,size(ST_S,3));

 
%%Now, to determine muscle interactions in continuous task parameter
%%spaces, we must first define some suitable variables

Task=%.....insert a cell array of continuous task variables here (e.g. kinematics, dynamics etc.)
    % that correspond with the EMG trials

%We need to align the time series length of the task variables with the EMG
%data

len=%...insert a representative time-series length of all trials (e.g. 1000)
t_ = linspace(0, 1, size(EMG_final,1));
for i=1:size(Task,3)

    t = linspace(0,1,size(Task,1));
    Task{i}=interp1(t, Task{i}, t_, 'spline');
end

%The task parameters should now be formatted as [Timepoints x Features x
%Trials]
Task=cat(3,Task{:});

%Extract the task redundant and synergistic muscle couplings

[ST_JOINT,ST_COND,ST_S,ST_R,MIs_z_single,...
   netS_joint,netT_joint,netS_cond,netT_cond,netS_S,netT_S,netS_R,netT_R]=SpaceTime_II_Continuous(EMG_final,Task);

%You are interested in the following outputs:
%ST_S: Task synergistic muscle couplings
%ST_R: Task redundant muscle couplings
%netS_S & netT_S: Spatial and temporal task-synergistic multiplex networks
%netS_R & netT_R: Spatial and temporal task-redundant multiplex networks

[Opt_rank_R_SP]=Consensus(netS_R,1);
[Opt_rank_R_TM]=Consensus(netT_R,1);

[Opt_rank_S_SP]=Consensus(netS_S,1);
[Opt_rank_S_TM]=Consensus(netT_S,1);



%Then extract out the Space-Time modules from the task-redundant
%interactions where Wi_R, Wb_R and Acal_R are the temporal and spatial
%modules and activation coefficients respectively.

[Wi_R,Acal_R,Wb_R,TS,VAF,E2]=sNM3F_basic(reshape(ST_R,[size(ST_R,1),size(ST_R,2)*size(ST_R,3)])',Opt_rank_R_SP,Opt_rank_R_TM,size(ST_R,3));


%Then extract out the Space-Time modules from the task-synergistic
%interactions where Wi_S, Wb_S and Acal_S are the temporal and spatial
%modules and activation coefficients respectively.

[Wi_S,Acal_S,Wb_S,TS,VAF,E2]=sNM3F_basic(reshape(ST_S,[size(ST_S,1),size(ST_S,2)*size(ST_S,3)])',Opt_rank_S_SP,Opt_rank_S_TM,size(ST_S,3));

