%%Processing script for manuscript: Functional muscle networks reveal the mechanistic effects 
    %of post-stroke rehabilitation on motor impairment and therapeutic responsiveness.


clear all;clc;
currentDir = cd;
addpath(currentDir);
myFolder_left= append(currentDir,'\EMG_left');
filePatterna_left = fullfile(myFolder_left,'*_sx.mat');
matFiles_left = dir(filePatterna_left);
myFolder_right= append(currentDir,'\EMG_right');
filePatterna_right = fullfile(myFolder_right,'*_dx.mat');
matFiles_right = dir(filePatterna_right);

Affected_side='R';

EMG_left={};
for i=1:size(matFiles_left,1)
    load(fullfile(myFolder_left, matFiles_left(i).name));
    [b,a]=butter(4,20/(SamplingFrequency/2),'low');
    if isequal(Affected_side,'L')
        for ii=2:2:10
            try
                e=Data(1:16,find(Data(17,:)==ii))';
                e=filtfilt(b,a,abs(e));
                %            e=interp1([1:length(e)]',e,[1:1000]');
                %             ef=[];
                %             for iii=1:20:1000
                %                 ef=[ef;trapz(e(iii:iii+19,:))];
                %             end
                EMG_left=cat(1,EMG_left,[e,repelem(i,length(e))']);
            catch message
            end
        end
    else
        for ii=1:2:10
            try
                e=Data(1:16,find(Data(17,:)==ii))';
                e=filtfilt(b,a,abs(e));
                %e=interp1([1:length(e)]',e,[1:1000]');
                %             ef=[];
                %             for iii=1:20:1000
                %                 ef=[ef;trapz(e(iii:iii+19,:))];
                %             end
                EMG_left=cat(1,EMG_left,[e,repelem(i,length(e))']);
            catch message
            end
        end
    end
end

EMG_right={};
for i=1:size(matFiles_right,1)
    load(fullfile(myFolder_right, matFiles_right(i).name));
    [b,a]=butter(4,20/(SamplingFrequency/2),'low');
    if isequal(Affected_side,'R')
        for ii=2:2:10
            try
                e=Data(1:16,find(Data(17,:)==ii))';
                e=filtfilt(b,a,abs(e));
                %e=interp1([1:length(e)]',e,[1:1000]');
                %             ef=[];
                %             for iii=1:20:1000
                %                 ef=[ef;trapz(e(iii:iii+19,:))];
                %             end
                EMG_right=cat(1,EMG_right,[e,repelem(i,length(e))']);
            catch message
            end
        end
    else
        for ii=1:2:10
            try
                e=Data(1:16,find(Data(17,:)==ii))';
                e=filtfilt(b,a,abs(e));
                %e=interp1([1:length(e)]',e,[1:1000]');
                %             ef=[];
                %             for iii=1:20:1000
                %                 ef=[ef;trapz(e(iii:iii+19,:))];
                %             end
                EMG_right=cat(1,EMG_right,[e,repelem(i,length(e))']);
            catch message
            end
        end
    end
end

% EMG_processed=cat(2,reshape(EMG_left,[50,16,10,7]),reshape(EMG_right,[50,16,10,7]));
% Task_parameter=repelem(1:size(EMG_processed,4),size(EMG_processed,3))';

EMG_processed_left=cell2mat(EMG_left);
EMG_processed_right=cell2mat(EMG_right);

combos=nchoosek(1:16,2);
R_Left=[];S_Left=[];
R_Right=[];S_Right=[];
for i=1:length(combos)
    try
        mi1=mi_mixture_gd(copnorm(EMG_processed_left(:,combos(i,1))), EMG_processed_left(:,17)-1, max(EMG_processed_left(:,17)));
    catch message
        mi1=0;
    end
    try
        mi2=mi_mixture_gd(copnorm(EMG_processed_left(:,combos(i,2))), EMG_processed_left(:,17)-1, max(EMG_processed_left(:,17)));
    catch message
        mi2=0;
    end
    
    try
        mi12=mi_mixture_gd(copnorm([EMG_processed_left(:,combos(i,1)),EMG_processed_left(:,combos(i,2))]),EMG_processed_left(:,17)-1, max(EMG_processed_left(:,17)));
    catch message
        mi12=0;
    end
    
    II=-[(mi1+mi2)-mi12];
    
    if II<0
        R_Left=[R_Left;abs(II)];
        S_Left=[S_Left;0];
    else
        S_Left=[S_Left;II];
        R_Left=[R_Left;0];
    end

    try
        mi1=mi_mixture_gd(copnorm(EMG_processed_right(:,combos(i,1))), EMG_processed_right(:,17)-1, max(EMG_processed_right(:,17)));
    catch message
        mi1=0;
    end

    try    
        mi2=mi_mixture_gd(copnorm(EMG_processed_right(:,combos(i,2))), EMG_processed_right(:,17)-1, max(EMG_processed_right(:,17)));
    catch messahe
        mi2=0;
    end
    
    try
        mi12=mi_mixture_gd(copnorm([EMG_processed_right(:,combos(i,1)),EMG_processed_right(:,combos(i,2))]),EMG_processed_right(:,17)-1, max(EMG_processed_left(:,17)));
    catch message
        mi12=0;
    end
    
    II=-[(mi1+mi2)-mi12];
    
    if II<0
        R_Right=[R_Right;abs(II)];
        S_Right=[S_Right;0];
    else
        S_Right=[S_Right;II];
        R_Right=[R_Right;0];
    end


end

net_R_left={};net_R_right={};
net_S_left={};net_S_right={};
for i=1:4
    net_R_left=squareform(R_Left);
    threshold = modified_percolation_analysis(net_R_left);
    net_R_left(net_R_left<threshold)=0;net_R_left(net_R_left<0)=0;

    net_R_right=squareform(R_Right);
    threshold = modified_percolation_analysis(net_R_right);
    net_R_right(net_R_right<threshold)=0;net_R_right(net_R_right<0)=0;

    net_S_left=squareform(S_Left);
    threshold = modified_percolation_analysis(net_S_left);
    net_S_left(net_S_left<threshold)=0;net_S_left(net_S_left<0)=0;

    net_S_right=squareform(S_Right);
    threshold = modified_percolation_analysis(net_S_right);
    net_S_right(net_S_right<threshold)=0;net_S_right(net_S_right<0)=0;
end

mask=tril(true(size(zeros(16,16))),-1);


if isequal(Affected_side,'L')
    Output.Affected.EMG_processed=EMG_processed_left;
    Output.Unaffected.EMG_processed=EMG_processed_right;

    Output.Affected.R=R_Left;Output.Unaffected.R=R_Right;Output.Affected.R=net_R_left(mask);Output.Unaffected.R=net_R_right(mask);
    Output.Affected.S=S_Left;Output.Unaffected.S=S_Right;Output.Affected.S=net_S_left(mask);Output.Unaffected.S=net_S_right(mask);

    Output.Affected.Nets.R=net_R_left;Output.Unaffected.Nets.R=net_R_right;
    Output.Affected.Nets.S=net_S_left;Output.Unaffected.Nets.S=net_S_right;
else
    Output.Affected.EMG_processed=EMG_processed_right;
    Output.Unaffected.EMG_processed=EMG_processed_left;

    Output.Affected.R=R_Right;Output.Unaffected.R=R_Left;Output.Affected.R=net_R_right(mask);Output.Unaffected.R=net_R_left(mask);
    Output.Affected.S=S_Right;Output.Unaffected.S=S_Left;Output.Affected.S=net_S_right(mask);Output.Unaffected.S=net_S_left(mask);

    Output.Affected.Nets.R=net_R_right;Output.Unaffected.Nets.R=net_R_left;
    Output.Affected.Nets.S=net_S_right;Output.Unaffected.Nets.S=net_S_left;
end



save(append(cd,'\Output.mat'),'Output');
