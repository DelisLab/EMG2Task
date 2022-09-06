function [Task_MI,ST_JOINT,ST_R,ST_S,ST_JOINT_netM,ST_JOINT_netT,ST_R_netM,...
    ST_R_netT,ST_S_netM,ST_S_netT]=SpaceTime_II(X,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Co-information: -(II) = I(Mx;T) +I(My;T) - I(Mx,My;T)

%Input 
%X=3D tensor [Timepoints x Muscles x Trials]
%Z=Matrix of discrete task variables equal in length to size(X,3)...
    %Values must be between [0 - Ym-1] (Ym=Max value in Y)

%Output
%Task_MI: Individual Muscle-task encodings
%ST_JOINT: Task-relevant muscle activations
%ST_R: Task-redundant muscle couplings
%ST_S: Task-synergistic muscle couplings
%#_netM: Multiplex network of muscle-wise couplings across temporal scales
%#_netT: Multiplex network of timepoint-wise couplings across spatial
        %scales



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combos=nchoosek(1:size(X,2),2);
combos_time=nchoosek(1:size(X,1),2);

%Individual-level task information
Task_MI=[];
Task_MI_P=[];
for zi=1:size(Z,2)
    z=Z(:,zi);
    mis=[];
    mis_p=[];
    for i=1:size(X,2)
        for ii=1:size(X,1)
            x_var=X(ii,i,:);
            x_var=copnorm(x_var(:));
            
            try
                I=mi_mixture_gd(x_var, z, max(z)+1);
                mis=[mis;I];
                
%Uncomment to significance test individual muscle-task encodings                
%                 perms=[];
%                 for iter=1:iterations
%                     In=mi_mixture_gd(x_var(randperm(length(x_var))), z, max(z)+1);
%                     perms=[perms;In];
%                 end
%                 mu=mean(perms);
%                 sd=std(perms);
%                 zscore=norminv(0.99);
%                 perm=mu+(zscore*sd);
%                 if I<perm
%                     mis_p=[mis_p;0];
%                 else
%                     mis_p=[mis_p;I];
%                 end
            catch message
            end
            
        end
    end
    Task_MI=cat(3,Task_MI,reshape(mis,[size(X,1),size(X,2)]));
    %Task_MI_P=cat(3,Task_MI_P,reshape(mis_p,[size(X,1),size(X,2)]));
    
end




%Pairwise task information
Joint_Task_MI=[];
for zi=1:size(Z,2)
    z=Z(:,zi);
    mis=[];
    for i=1:length(combos)
        for ii=1:length(combos_time)
            x_var1=X(combos_time(ii,1),combos(i,1),:);
            x_var1=x_var1(:);
            
            x_var2=X(combos_time(ii,2),combos(i,2),:);
            x_var2=x_var2(:);
            
            try
                I=mi_mixture_gd(copnorm([x_var1,x_var2]), z, max(z)+1);
                mis=[mis;I];
            catch message
                mis=[mis;ent_g(copnorm([x_var1,x_var2]),true)];
            end
        end
        
        for ii=1:size(X,1)
            x_var1=X(ii,combos(i,1),:);
            x_var1=x_var1(:);
            
            x_var2=X(ii,combos(i,2),:);
            x_var2=x_var2(:);
            
            try
                I=mi_mixture_gd(copnorm([x_var1,x_var2]), z, max(z)+1);
                mis=[mis;I];
            catch message
                mis=[mis;ent_g(copnorm([x_var1,x_var2]),true)];
            end
        end
        
    end
    Joint_Task_MI=cat(3,Joint_Task_MI,reshape(mis,[length(combos_time)+size(X,1),length(combos)]));
end



%Summing individual muscle-task encodings
ST=[];

for i=1:size(Task_MI,3)
    task=reshape(Task_MI(:,:,i),[size(X,1),size(X,2)]);
    st=[];
    for ii=1:length(combos)
        for iii=1:length(combos_time)
            st=[st,sum([task(combos_time(iii,1),combos(ii,1)),task(combos_time(iii,2),combos(ii,2))],2)];
        end
        
        for iii=1:size(X,1)
            st=[st,sum([task(iii,combos(ii,1))',task(iii,combos(ii,1))'],2)];
        end
    end
    ST=cat(3,ST,reshape(st,[length(combos_time)+size(X,1),length(combos)]));
end




%II calculation
ST_II=[];
for i=1:size(ST,3)
    st_mi=reshape(ST(:,:,i),[length(combos_time)+size(X,1),length(combos)]);
    j_mi=reshape(Joint_Task_MI(:,:,i),[length(combos_time)+size(X,1),length(combos)]);
    st_II=[];
    
    for ii=1:size(st_mi,2)
        for iii=1:size(st_mi,1)
            
            st_II=[st_II,st_mi(iii,ii)-j_mi(iii,ii)];
        end
    end
    
    
    ST_II=cat(3,ST_II,reshape(st_II,[size(st_mi,1),size(st_mi,2)]));
    
end






ST_JOINT_M=[];
ST_JOINT_T=[];
ST_JOINT_netM={};
ST_JOINT_netT={};
for i=1:size(Joint_Task_MI,3)
    rep=[];
    for ii=1:size(Joint_Task_MI,1)
        
        A=squareform(reshape(Joint_Task_MI(ii,:,i),[1,length(combos)]));
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_JOINT_netM=cat(2,ST_JOINT_netM,A);
        rep=cat(2,rep,A(mask));
    end
    ST_JOINT_M=cat(3,ST_JOINT_M,rep);
    
    rep=[];
    for ii=1:size(Joint_Task_MI,2)
        A=squareform(reshape(Joint_Task_MI(1:length(combos_time),ii,i),[length(combos_time),1]));
        d=diag(reshape(Joint_Task_MI(length(combos_time)+1:end,ii,i),[size(X,1),1]));
        A=A+d;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_JOINT_netT=cat(2,ST_JOINT_netT,A);
        rep=cat(2,rep,A(mask));
    end
    ST_JOINT_T=cat(3,ST_JOINT_T,rep);

end

%Space-Time heuristic
for t=1:size(ST_JOINT_T,3)
    for col=1:size(ST_JOINT_T,2)
        for row=1:size(ST_JOINT_T,1)
            if ~ST_JOINT_T(row,col,t)>0 || ~ST_JOINT_M(col,row,t)>0
                ST_JOINT_T(row,col,t)=0;
                ST_JOINT_M(col,row,t)=0;
            end
        end
    end
end
ST_JOINT=ST_JOINT_M;



ST_Redundancy=ST_II*-1;
ST_Redundancy(ST_Redundancy>0)=0;
ST_Redundancy=abs(ST_Redundancy);


ST_Synergy=ST_II*-1;
ST_Synergy(ST_Synergy<0)=0;





ST_R_M=[];
ST_S_M=[];

ST_R_netM={};
ST_S_netM={};

ST_R_T=[];
ST_S_T=[];

ST_R_netT={};
ST_S_netT={};

for t=1:size(ST_Redundancy,3)
    rep=[];
    for i=1:size(ST_Redundancy,1)
        A=squareform(reshape(ST_Redundancy(i,:,t),[1,size(ST_Redundancy,2)]));
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_R_netM=cat(2,ST_R_netM,A);
        rep=cat(2,rep,A(mask));
    end
    ST_R_M=cat(3,ST_R_M,rep);
    
    rep=[];
    for i=1:size(ST_Synergy,1)
        A=squareform(reshape(ST_Synergy(i,:,t),[1,size(ST_Synergy,2)]));
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_S_netM=cat(2,ST_S_netM,A);
        rep=cat(2,rep,A(mask));
    end
    ST_S_M=cat(3,ST_S_M,rep);
    
    
    rep=[];
    for i=1:size(ST_Redundancy,2)
        A=squareform(reshape(ST_Redundancy(1:length(combos_time),i,t),[1,length(combos_time)]));
        d=diag(reshape(ST_Redundancy(length(combos_time)+1:end,i,t),[1,size(X,1)]));
        A=A+d;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_R_netT=cat(2,ST_R_netT,A);
        rep=cat(2,rep,A(mask));
    end
    ST_R_T=cat(3,ST_R_T,rep);
    
    rep=[];
    for i=1:size(ST_Synergy,2)
        A=squareform(reshape(ST_Synergy(1:length(combos_time),i,t),[1,length(combos_time)]));
        d=diag(reshape(ST_Synergy(length(combos_time)+1:end,i,t),[1,size(X,1)]));
        A=A+d;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        ST_S_netT=cat(2,ST_S_netT,A);
        rep=cat(2,rep,A(mask));
    end
    ST_S_T=cat(3,ST_S_T,rep);
    
end
    
%Space-Time Heurisitic
for t=1:size(ST_S_T,3)
    for col=1:size(ST_S_T,2)
        for row=1:size(ST_S_T,1)
            if ~ST_S_T(row,col,t)>0 || ~ST_S_M(col,row,t)>0
                ST_S_T(row,col,t)=0;
                ST_S_M(col,row,t)=0;
            end
            
            if ~ST_R_T(row,col,t)>0 || ~ST_R_M(col,row,t)>0
                ST_R_T(row,col,t)=0;
                ST_R_M(col,row,t)=0;
            end
            
            
        end
    end
end
ST_S=ST_S_M;
ST_R=ST_R_M;



end