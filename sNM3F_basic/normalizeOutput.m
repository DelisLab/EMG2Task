function [A,Wi,Wb] = normalizeOutput(A,Wi,Wb,N,P,S)

% Row-wise normalization
sP=zeros(1,P);
for i=1:P
   sP(i)=norm(Wi(:,i));
   Wi(:,i)=Wi(:,i)./sP(i);
end
% Column-wise normalization
sN=zeros(1,N);
for j=1:N
   sN(j)=norm(Wb(j,:));
   Wb(j,:)=Wb(j,:)./sN(j);
end
% Element-wise normalization of A to make the error unchanged
for s=1:S
    for i=1:P
       for j=1:N  
        A(i,N*(s-1)+j)=A(i,N*(s-1)+j).*(sP(i)*sN(j));
       end
    end
end

end %#EoF normalizeOutput

