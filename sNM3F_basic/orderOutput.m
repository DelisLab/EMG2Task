function [Anew,Wi,Wb] = orderOutput(A,Wi,Wb,P,N,S)

Anew=zeros(size(A));

% Order elements by the occurence of the maximal element
% in columns for Wi
iP=zeros(1,P);
for i=1:P
    [~,iP(i)]=max(Wi(:,i));
end
[~,idxP] = sort(iP,'ascend');

% in rows for Wb
iN=zeros(1,N);
for j=1:N
    [~,iN(j)]=max(Wb(j,:));
end
[~,idxN] = sort(iN,'ascend');

Wi=Wi(:,idxP);
Wb=Wb(idxN,:);

for s=1:S
    for j=1:N
        Anew(:,N*(s-1)+j)=A(:,N*(s-1)+idxN(j));
    end
end
Anew(:,:)=Anew(idxP,:);

end %#EoF orderOutput
