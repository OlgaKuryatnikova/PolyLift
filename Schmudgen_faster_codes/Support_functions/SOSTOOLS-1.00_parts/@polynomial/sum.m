function varargout=sum(A,dim)
% function C=sum(A,Dim)
%
% DESCRIPTION 
%   Sum the elements of a polynomial vector or matrix.
%   
% INPUTS 
%   A: polynomial vector or matrix
%   Dim: dimension along which summation is to be performed
%      (when A is a matrix)
%
% OUTPUTS  
%   C: polynomial or polynomial vector, the result of summation 
%  
% SYNTAX 
%   C = sum(A)
%   C = sum(A,1)
%   C = sum(A,2)
%

% Initial coding: 01/31/03 -- SP
% 3/25/03 PJS   If x is a column vector, sum(x,2) should return x.
%               If nargout = 0, display the result.

% Promote A
A = polynomial(A);

if any(size(A) == 1) & nargin==1
   C = polynomial(0);
   for i1 = 1:prod(size(A))
     s.type = '()';
     s.subs = {i1};
     C = C + subsref(A,s); 
   end
elseif nargin == 1 | dim == 1
   C = polynomial(zeros(1,size(A,2)));
   for i1 = 1:size(A,1)
     s.type = '()';
     s.subs = {i1,':'};
     C = C + subsref(A,s);
   end
else
   C = polynomial(zeros(size(A,1),1));
   for i1 = 1:size(A,2);
     s.type = '()';
     s.subs = {':', i1};
     C = C + subsref(A,s);
   end
end

varargout{1} = C;
