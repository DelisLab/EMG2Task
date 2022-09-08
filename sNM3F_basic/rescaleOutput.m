function [A,Wi,Wb] = rescaleOutput(A,Wi,Wb)

% Make sure that Wi'*Wi and Wb*Wb' scales to identity
sPm=max(max(Wi'*Wi));
Wi=Wi./sqrt(sPm);
sNm=max(max(Wb*Wb'));
Wb=Wb./sqrt(sNm);
A=A.*sqrt(sPm)*sqrt(sNm);

end %#EoF rescaleOutput
