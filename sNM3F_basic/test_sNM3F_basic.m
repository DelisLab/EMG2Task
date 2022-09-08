% EXAMPLE OF MULTIPLE RUNS IMPLEMENTATION
% in practice, you may want to run the algorithm multiple times
% because of local minima

clc; close all; clear all;

T=20; % temporal dimensions
M=6; % spatial dimensions (must be divisible by 3 in this example)
S=10; % total number of samples

% Select how many components to extract
P=3; % 3 temporal components
N=3; % 3 spatial components

% Load EMG-like data from known modules and coefficients
load testdata_N=3_P=3.mat

% Run the sNM3F algorithm - in practice, multiple runs should be performed

REP=50; % number of repetitions

st.temporalModules=[]; st.spatialModules=[];
st.combinators=[]; st.VAF=[]; st.ERR=[];
restmp=repmat(st,1,REP); ERR=nan(1,REP);
for r=1:REP       
    fprintf(['Run #' num2str(r) '\n']);
    [Wi,Acal,Wb,VAF,E2]=sNM3F_basic(Mb,P,N,S);  
    restmp(r).temporalModules=Wi;
    restmp(r).spatialModules=Wb;
    restmp(r).combinators=Acal;
    restmp(r).VAF=VAF;
    restmp(r).ERR=sqrt(E2(end)/(size(Mb,1)*size(Mb,2)));
    restmp(r).E2=E2;
    ERR(r)=E2(end);
end

% Return the best run (with respect to reconstruction error)
[ignr,indr]=min(ERR);
RESULT=restmp(indr);