function b = trace(a)
% function B=trace(A)
%
% DESCRIPTION 
%   Sum of diagonal elements
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: sum of the diagonal elements of A. This is
%    also the sum of the eigenvalues of A.
%  
% SYNTAX 
%   B = trace(A);

% 6/7/2002: PJS  Initial Coding  

  
mnsz = min(size(a));  
if isempty(a)
  
  b=polynomial;
  return;
  
elseif mnsz==1
  
  b=a;
  
else
  
  nta = size(a.degmat,1);  
  [nra,nca]=size(a);

  % Grab diagonal elements
  acoef = reshape(a.coefficient,nta,nra*nca);
  idx = find(eye(nra,nca));
  acoef = acoef(:,idx);

  % Sum diagonal elements
  coefficient = sum(acoef,2);
  b = a;
  b.coefficient = coefficient(:);
  b.matdim = [1 1];
  b = combine(b);
  
end
   