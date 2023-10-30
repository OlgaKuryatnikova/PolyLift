function b = diff(a,x)
% function B=diff(A,X)
%
% DESCRIPTION 
%   Differentiate a polynomial.
%   
% INPUTS 
%   A: polynomial 
%   X: differentiate with respect to X.
%
% OUTPUTS  
%   B: polynomial
%  
% SYNTAX 
%   B = diff(A,X);
%     Differentiate the polynomial, A, with to X.  A should be a
%     polynomial and X should be a polynomial variable or a string.  
%     Differentiation is done element-by-element if A is a matrix.

% 11/5/2002: PJS  Initial Coding  

% Error Checking
if nargin~=2 
  error('Error in calling diff');
end

if isa(x,'polynomial')
  x = x.varname{1};
elseif ~ischar(x)
  error('X must be a polynomial variable or a string');
end

% Get polynomial info about a
a = polynomial(a);
a = combine(a);
adeg = a.degmat;
avar = a.varname;
nta = size(adeg,1);
[nra,nca]=size(a);
acoef = reshape(a.coefficient,nta,nra*nca);

% Find variable we are differentiating with respect to.
varnumb = 0;
for i1 = 1:length(avar)
  if strcmp(avar{i1},x)
    varnumb = i1;
  end
end

% Differentiate
if ~varnumb
  b = polynomial(0);
else
  for i1 = 1:nta
    idx = (i1-1)*nca+(1:nca);
    acoef(i1,:) = acoef(i1,:)*adeg(i1,varnumb);
  end
  adeg(:,varnumb) = max( adeg(:,varnumb)-1 , 0);
  b = combine(polynomial(acoef,adeg,avar,[nra nca]));
end



