function [ST_JOINT,ST_COND,ST_S,ST_R,MIs_z_single,...
   netS_joint,netT_joint,netS_cond,netT_cond,netS_S,netT_S,netS_R,netT_R]=SpaceTime_II_Continuous(X,Z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
    %X(EMG) & Z(Contious task variables) are both 3D tensors [Timepoints, Channels,
    %Trials]

%Output in the shape [Spatial interactions x [All Timepoint A and B interactions, All Timepoint A and A interactions]]
    %ST_JOINT: Task-relevant muscle activiations across trials
    %ST_COND: Task-irrelevant muscle couplings across trials
    %ST_S: Task-synergistic muscle couplings across trials
    %ST_R: Task-redundant muscle couplings across trials
    
    %netS_#: Multiplex network of Muscle-wise couplings across temporal
            %scales
    %netT_#: Multiplex network of Timepoint-wise couplings across spatial
            %scales


            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
%Initialise muscle and timepoint pairwise combinations
combos_t=nchoosek(1:size(X,1),2);
combos_m=nchoosek(1:size(X,2),2);



%Muscle-task encodings, task-dependent and -independent space-time models

MIs_z_joint=[];
MIs_z_single=[];
MIs_z_cond=[];
for t=1:size(Z,2)
    
    MIs_joint=[];
    MIs_single=[];
    MIs_cond=[];
    
    
     for ii=1:length(combos_t)
        mis_single=[];
        for iii=1:size(X,2)
            
            vars1=X(combos_t(ii,1),iii,:);
            zi=Z(combos_t(ii,1),t,:);
            try
                I_single=mi_gg(copnorm(vars1(:)),copnorm(zi(:)),true,true);
                mis_single=[mis_single;I_single];
%  Uncomment to significance test individual muscle-task encodings
%
%                 In=mi_gg(copnorm(vars1(:)),copnorm(zi(:)),false,true);
%                 perms=[];
%                 vars1_perm=vars1(:);
%                 for ps=1:50
%                     r=vars1_perm(randperm(length(vars1_perm)));
%                     perms=[perms;mi_gg(copnorm(r),copnorm(zi(:)),false,true)];
%                 end
%                 mu=mean(perms);
%                 sd=std(perms);
%                 zscore=norminv(0.95);
%                 perm=mu+(zscore*sd);
%                 if In<perm
%                     mis_single=[mis_single;0];
%                 else
%                     mis_single=[mis_single;I_single];
%                 end
            catch message
                mis_single=[mis_single;0];
            end
        end
%         
        
        mis_joint=[];
        mis_cond=[];
        for iii=1:length(combos_m)
            vars1=X(combos_t(ii,1),combos_m(iii,1),:);
            vars2=X(combos_t(ii,2),combos_m(iii,2),:);
            zi=Z(combos_t(ii,1),t,:);
            try
                I_joint=mi_gg(copnorm([vars1(:),vars2(:)]),copnorm(zi(:)),true,true);%/min(ent_g(copnorm([vars1(:),zi(:)]),true),ent_g(copnorm(vars2(:)),true));
                mis_joint=[mis_joint;I_joint];
            catch message
                mis_joint=[mis_joint;ent_g([copnorm(vars1(:)),copnorm(zi(:))],true)];
            end
            
            try
                I_cond=cmi_ggg(copnorm(vars1(:)),copnorm(vars2(:)),copnorm(zi(:)),true,true);
                mis_cond=[mis_cond;I_cond];
            catch message
                mis_cond=[mis_cond;0];
                
            end
            
        end
        
        
        
        
        MIs_single=cat(2,MIs_single,reshape(mis_single,[size(X,2),1]));
        MIs_joint=cat(2,MIs_joint,reshape(mis_joint,[length(combos_m),1]));
        MIs_cond=cat(2,MIs_cond,reshape(mis_cond,[length(combos_m),1]));
    end
    
    MIs_linear_cond=[];
    MIs_linear_joint=[];
    MIs_linear_single=[];
    for i=1:size(X,1)
        mis_single=[];
        for ii=1:size(X,2)
            vars1=X(i,ii,:);
            zi=Z(i,t,:);
            try
                I_single=mi_gg(copnorm(vars1(:)),copnorm(zi(:)),true,true);
                mis_single=[mis_single;I_single];
                
%Uncomment to significance test individual muscle-task encodings
%                 In=mi_gg(copnorm(vars1(:)),copnorm(zi(:)),false,true);
%                 perms=[];
%                 vars1_perm=vars1(:);
%                 for ps=1:50
%                     r=vars1_perm(randperm(length(vars1_perm)));
%                     perms=[perms;mi_gg(copnorm(r),copnorm(zi(:)),false,true)];
%                 end
%                 mu=mean(perms);
%                 sd=std(perms);
%                 zscore=norminv(0.95);
%                 perm=mu+(zscore*sd);
%                 if In<perm
%                     mis_single=[mis_single;0];
%                 else
%                     mis_single=[mis_single;I_single];
%                 end
             catch message
                 mis_single=[mis_single;0];
             end
        end
        MIs_linear_single=cat(2,MIs_linear_single,reshape(mis_single,[size(X,2),1]));
        
        mis_joint_linear=[];
        mis_cond_linear=[];
        for ii=1:length(combos_m)
            vars1=X(i,combos_m(ii,1),:);
            vars2=X(i,combos_m(ii,2),:);
            zi=Z(i,t,:);
            try
                I_joint=mi_gg([copnorm(vars1(:)),copnorm(vars2(:))],copnorm(zi(:)),true,true);%/min(ent_g(copnorm([vars1(:),vars2(:)],true)),ent_g(copnorm(zi(:)),true));
                mis_joint_linear=[mis_joint_linear;I_joint];
            catch message
                mis_joint_linear=[mis_joint_linear;ent_g([copnorm(vars1(:)),copnorm(zi(:))],true)];
            end
            
            try
                I_cond=cmi_ggg(copnorm(vars1(:)),copnorm(vars2(:)),copnorm(zi(:)),true,true);
                mis_cond_linear=[mis_cond_linear;I_cond];
            catch message
                mis_cond_linear=[mis_cond_linear;0];
            end
        end
        MIs_linear_joint=cat(2,MIs_linear_joint,reshape(mis_joint_linear,[length(combos_m),1]));
        MIs_linear_cond=cat(2,MIs_linear_cond,reshape(mis_cond_linear,[length(combos_m),1]));
    end
    
    
    MIs_z_single=cat(3,MIs_z_single,cat(2,MIs_single,MIs_linear_single));
    MIs_z_joint=cat(3,MIs_z_joint,cat(2,MIs_joint,MIs_linear_joint));
    MIs_z_cond=cat(3,MIs_z_cond,cat(2,MIs_cond,MIs_linear_cond));
end

%II
%%Redundancy and synergy


II=[];
for t=1:size(MIs_z_single,3)
    ii2=[];
    for i=1:length(combos_t)+size(X,1)
        for ii=1:length(combos_m)
            try
                xs1=MIs_z_single(combos_m(ii,1),combos_t(i,1),t);
                xs2=MIs_z_single(combos_m(ii,2),combos_t(i,2),t);
            catch message
                xs1=MIs_z_single(combos_m(ii,1),i-length(combos_t),t);
                xs2=MIs_z_single(combos_m(ii,2),i-length(combos_t),t);
            end
           
            xjoint=MIs_z_joint(ii,i,t);
            ii2=[ii2;(xs1+xs2)-xjoint];
        end
    end
    
    %Change signs for co-information
    II=cat(3,II,reshape(-(ii2),[length(combos_m),length(combos_t)+size(X,1)]));
end
% 
% 
% 
% 
% 
% 
ST_jointS=[];
ST_condS=[];
ST_jointT=[];
ST_condT=[];
% 

ST_rS=[];
ST_rT=[];
ST_sS=[];
ST_sT=[];

netS_joint={};
netS_cond={};
netT_joint={};
netT_cond={};

netS_R={};
netT_R={};

netS_S={};
netT_S={};

%Sparsify networks using modified percolation analysis
 for t=1:size(MIs_z_joint,3)
     st_jointS=[];
     st_condS=[];
     st_jointT=[];
     st_condT=[];
%     
    st_rS=[];
    st_rT=[];
    st_sS=[];
    st_sT=[];
    for i=1:size(MIs_z_joint,2)
        A=squareform(MIs_z_joint(:,i,t));
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        netS_joint=cat(2,netS_joint,A);
        st_jointS=cat(2,st_jointS,A(mask));
%         
        A=squareform(MIs_z_cond(:,i,t));
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        netS_cond=cat(2,netS_cond,A);
        st_condS=cat(2,st_condS,A(mask));


        A=squareform(II(:,i,t));
        A(A>0)=0;
        A=abs(A);
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        netS_R=cat(2,netS_R,A);
        st_rS=cat(2,st_rS,A(mask));
        


        A=squareform(II(:,i,t));
        A(A<0)=0;
        mask = tril(true(size(A,1)),-1);
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        netS_S=cat(2,netS_S,A);
        st_sS=cat(2,st_sS,A(mask));
        
     end
%     
    for i=1:size(MIs_z_joint,1)
        A=squareform(MIs_z_joint(i,1:length(combos_t),t));
        d=diag(MIs_z_joint(i,length(combos_t)+1:end,t));
        A=A+d;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        netT_joint=cat(2,netT_joint,A);
        st_jointT=cat(2,st_jointT,A(mask));
%         
        A=squareform(MIs_z_cond(i,1:length(combos_t),t));
        d=diag(MIs_z_cond(i,length(combos_t)+1:end,t));
        A=A+d;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        A(A<0)=0;
        netT_cond=cat(2,netT_cond,A);
        st_condT=cat(2,st_condT,A(mask));
        
%         

%         
%      
        A=squareform(II(i,1:length(combos_t),t));
        d=diag(II(i,length(combos_t)+1:end,t));
        A=A+d;
        A(A>0)=0;
        A=abs(A);
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        netT_R=cat(2,netT_R,A);
        st_rT=cat(2,st_rT,A(mask));
        

%         
       
        A=squareform(II(i,1:length(combos_t),t));
        d=diag(II(i,length(combos_t)+1:end,t));
        A=A+d;
        A(A<0)=0;
        mask = tril(true(size(A,1)));
        [threshold] = modified_percolation_analysis(A);
        A(A<threshold)=0;
        netT_S=cat(2,netT_S,A);
        st_s2T=cat(2,st_s2T,A(mask));
        
        
        
     end
     ST_jointS=cat(3,ST_jointS,st_jointS);
     ST_condS=cat(3,ST_condS,st_condS);
     ST_jointT=cat(3,ST_jointT,st_jointT);
%     
     ST_condT=cat(3,ST_condT,st_condT);
%     
     
     ST_rS=cat(3,ST_rS,st_rS);
     
     ST_rT=cat(3,ST_rT,st_rT);
     ST_sS=cat(3,ST_sS,st_sS);
     ST_sT=cat(3,ST_sT,st_s2T);
%     
 end
% 
% 

% 
% 
% 
ST_JOINT=[];
ST_COND=[];
ST_R=[];
ST_S=[];
% 

% Space-Time Heuristic 
for t=1:size(ST_jointT,3)
    for col=1:size(ST_jointT,2)
        for row=1:size(ST_jointT,1)
            if ~ST_jointT(row,col,t)>0 || ~ST_jointS(col,row,t)>0
                ST_jointT(row,col,t)=0;
                ST_jointS(col,row,t)=0;
            end
%             
            if ~ST_condT(row,col,t)>0 || ~ST_condS(col,row,t)>0
                ST_condT(row,col,t)=0;
                ST_condS(col,row,t)=0;
            end
            
%             
            if ~ST_rT(row,col,t)>0 || ~ST_rS(col,row,t)>0
                ST_rT(row,col,t)=0;
                ST_rS(col,row,t)=0;
            end
%             

% %             
            if ~ST_sT(row,col,t)>0 || ~ST_sS(col,row,t)>0
                ST_sT(row,col,t)=0;
                ST_sS(col,row,t)=0;
            end
            
        end
    end
 end
%

%Final Input matrices for dimensionality reduction
ST_JOINT=ST_jointS;
ST_COND=ST_condS;
ST_S=ST_sS;
ST_R=ST_rS;
% 
% 
% 
% 
 end
