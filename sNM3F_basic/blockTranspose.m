function Mat_out=blockTranspose(Mat_in,type,S)

if strcmpi(type,'r') % in rows, Mat_in is [T*S,M]
    [r,c]=size(Mat_in);
    D=r/S; % equal to T
    Mat_out=zeros(D,c*S); % We are going to get a [T,M*S] matrix
    for s=1:S
        Mat_out(:,c*(s-1)+1:c*s)=Mat_in(D*(s-1)+1:D*s,:);
    end
elseif strcmpi(type,'c') % in columns, Mat_in is [T,M*S]
    [r,c]=size(Mat_in);
    D=c/S; % equal to M
    Mat_out=zeros(r*S,D); % We are going to get a [T*S,M] matrix
    for s=1:S
    	Mat_out(r*(s-1)+1:r*s,:)=Mat_in(:,D*(s-1)+1:D*s); 
    end
else
    error('Unknown option. It is either "r" or "c".')
end

end %#EoF blockTranspose
