function [Wi,Acal,Wb,VAF,E2]=sNM3F_basic(Mb,P,N,S)
%------------------------------------------------------------------------
% Basic version of the sample-based non-negative matrix tri-factorization
%------------------------------------------------------------------------

%--- Description
% Input params:
%   - Mb = input matrix composed: vertical concatenation of the recorded
%          data across multiple episodes (size T*S x M)
%   - P  = number of temporal (or row) modules to be extracted
%   - N  = number of spatial (or column) modules to be extracted
%   - S  = number of samples (i.e. episodes, trials...)
% Output:
%   - Wi   = P temporal modules
%   - Wb   = N spatial modules
%   - Acal = activation coefficients
%   - VAF  = variance accounted for
%   - E2   = total reconstruction error (across iterations)

%--- GNU GPL Licence
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>.

% If you use this software for your research, please refer to the paper entitled 
% "A unifying model of concurrent spatial and temporal modularity in muscle activity", 
%  by I. Delis, S. Panzeri, T. Pozzo and B. Berret

%--- Contact the authors
%  Bastien Berret (bastien.berret@u-psud.fr - https://sites.google.com/site/bberret) 
%  Ioannis Delis (ioannis.delis@glasgow.ac.uk)


%-------------------------------------------------------------------------
% sNM3F_basic ALGORITHM PARAMETERS - YOU CAN EDIT HERE
%-------------------------------------------------------------------------
MAXITER=1000; % Maximum iteration before stopping is enforced
ERRTOL=1e-8; % Tolerance on the recontruction error changes
ORTHOGONALIZE=1; % Impose or not the constraints Wi'*Wi=Id and Wb*Wb'=Id
CLEANOUTPUT=0; % Normalize the output at each iteration
DISPLAYITER=0; % Display or not some information at each iteration
NSTOPS=5; % Number of steps for the stopping criterion

%-------------------------------------------------------------------------
% DATASET CHARACTERISTIC
%-------------------------------------------------------------------------
if rem(size(Mb,1),S)==0
    T=size(Mb,1)/S; % number of temporal dimensions (e.g. time frames)
else
   error(['Please check your input matrix... It should be [T*S, M] with ' ...
          'T=number of time steps, S=number of samples and M=number of muscles']);  
end
M=size(Mb,2); % number of spatial dimensions (e.g. EMG channels)
% S is the number of episodes (total number of repetitions/trials)

%---- Get the block transpose of Mb
Mi=blockTranspose(Mb,'r',S);

%-------------------------------------------------------------------------
% BASIC sNM3F ALGORITHM
%-------------------------------------------------------------------------
disp(['Start extracting ' num2str(P) ' temporal modules and ' ...
       num2str(N) ' spatial modules']);

%---- Initialization of arrays - YOU CAN EDIT THE INITIAL GUESS
Wb=rand(N,M);
Wi=rand(T,P);         
A=rand(P,N*S);  

%--- Error container
err=NaN(1,MAXITER);

%--- Initialize some variables
count=0; it=0;
telapsed=zeros(1,MAXITER);

%---- Main iterative loop of the algorithm  
%     N.B.: we decompose the data as Mi=Wi*A*Wb (for each sample)
while (count<NSTOPS && it<MAXITER) 
    tic; 
    it=it+1; 
    %---- Clean the Output, by reordering or rescaling etc.
    %     May improve the robustness & convergence
    if CLEANOUTPUT,
     % Norm rows and colums to one
     [A,Wi,Wb] = normalizeOutput(A,Wi,Wb,N,P,S);
     if ORTHOGONALIZE, 
      % Ensure we get matrices of the order of identity for Wi'*Wi and Wb*Wb'
      [A,Wi,Wb] = rescaleOutput(A,Wi,Wb);  %#ok<UNRCH>
     end
     % Then, order the elements
     [A,Wi,Wb] = orderOutput(A,Wi,Wb,P,N,S);
    end
    
    % 1st STEP ----------------------------------------------------------            
    %--- UPDATE Wb 
    %--- We approximate Mb=Mi^{\prime}=Cb^{\prime}*Wb=Cbb*Wb
    Cb=Wi*A; 
    Cbb=blockTranspose(Cb,'c',S);
    numer=Cbb'*Mb;

    if ORTHOGONALIZE,
    	denom=Cbb'*Mb*(Wb'*Wb); %#ok<UNRCH>
    	Wb= Wb.*numer./(denom+eps(numer));
    else
    	denom=(Cbb'*Cbb)*Wb;
    	Wb= Wb.*numer./(denom+eps(numer));
    end
    
    % 2nd STEP ----------------------------------------------------------            
    %--- UPDATE Wi 
    Atil=blockTranspose(A,'c',S);
    Ci=Atil*Wb;
    Cii=blockTranspose(Ci,'r',S);
    numer=Mi*Cii';

    if ORTHOGONALIZE,
    	denom=(Wi*Wi')*Mi*Cii';  %#ok<UNRCH>
    	Wi = Wi.*numer./(denom+eps(numer)); 
    else
    	denom=Wi*(Cii*Cii');
    	Wi= Wi.*numer./(denom+eps(numer));
    end

    % 3rd STEP ----------------------------------------------------------            
    %--- UPDATE A for all samples s=1..S and compute the error err(it)
    err(it)=0;
    for s=1:S
    	denom=(Wi'*Wi)*Atil(P*(s-1)+1:P*s,:)*(Wb*Wb');
        numer=Wi'*Mi(:,M*(s-1)+1:M*s)*Wb';
    	A(:,N*(s-1)+1:N*s)=A(:,N*(s-1)+1:N*s).*(numer./(denom+eps(numer)));
    	err(it)=err(it)+norm(Mi(:,M*(s-1)+1:M*s)-Wi*A(:,N*(s-1)+1:N*s)*Wb,'fro')^2;
    end
    
    telapsed(it)=toc;
	if DISPLAYITER,
    	disp(['iter #' num2str(it) ' | Error=' num2str( err(it)) ' | Time (s)=' num2str(telapsed(it))]);
	end
                 
    %---- Implement convergence criterion
    if it>2
        if (abs(err(it)-err(it-1))<ERRTOL)
            count=count+1;
            if DISPLAYITER,
                disp(['stop counter #' num2str(count)]);
            end
        elseif err(it)>err(it-1)
            %warning('sNM3F:errorIncreased','Error increased during this iteration...');
            break
        else
            count=0;
        end
    end 
                
end


%---- Constant useful to compute the VAF
SST=0;
Mmean=zeros(T,M);
for s=1:S, Mmean=Mmean+Mi(:,M*(s-1)+1:M*s)/S; end
for s=1:S
 SST=SST+norm(Mi(:,M*(s-1)+1:M*s)-Mmean,'fro')^2;
end

E2=err(1:it);
VAF=1-E2(end)/SST;

%---- Clean the Output, by reordering or rescaling etc.
if CLEANOUTPUT,
 % Norm rows and colums to one
 [A,Wi,Wb] = normalizeOutput(A,Wi,Wb,N,P,S);
 if ORTHOGONALIZE, 
  % Ensure we get matrices of the order of identity for Wi'*Wi and Wb*Wb'
  [A,Wi,Wb] = rescaleOutput(A,Wi,Wb);  %#ok<UNRCH>
 end
 % Then, order the elements
 [A,Wi,Wb] = orderOutput(A,Wi,Wb,P,N,S);
end

%---- Build the full Acal cell array of activations
Acal=cell(1,S);
for s=1:S
    Acal{s}=A(:,N*(s-1)+1:s*N);
end

disp(['Finished! VAF is ' num2str(VAF) ' | Err is ' num2str(E2(end))]);
disp(['Total time elapsed ' num2str(sum(telapsed)) ' seconds']);     
fprintf('\n');

end %#EoF sNM3F_basic
