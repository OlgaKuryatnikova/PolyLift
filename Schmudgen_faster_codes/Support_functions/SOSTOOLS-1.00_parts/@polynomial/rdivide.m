function c = rdivide(a,b)
% function C=rdivide(A,B)
%
% DESCRIPTION 
%   Element-by-element right division
%   
% INPUTS 
%   A: polynomial
%   B: matrix of constants
%
% OUTPUTS  
%   C: polynomial, the result of division.
%  
% SYNTAX 
%   C= A./B
%     A and B must have the same dimensions unless one
%     is a scalar.  Scalars can be divided with anything.     
%   C = rdivide(A,B)
%     Function-call form 

% 10/22/2002: PJS  Initial Coding  
  
% Promote a to polynomial 
a = polynomial(a);
sza = size(a);

% Promote b to polynomial 
b = polynomial(b);
szb = size(b);

if isempty(a) | isempty(b)

  if isempty(a) & isempty(b)
    % empty./empty = empty
    c=polynomial;
    return;
  elseif isempty(a) & max(size(b))==1
    % empty./scalar = empty
    c=polynomial;
    return;
  elseif isempty(b) & max(size(a))==1
    % scalar./empty = empty
    c=polynomial;
    return;
  else
    % matrix./empty  or empty./matrix     
    error('Matrix dimensions must agree.');
    return
  end

elseif all(sza==szb) | all(sza==[1 1]) | ...
      all(szb==[1 1])

  % Matrix./Matrix  
  degidx = find(b.degmat);
  if ~isempty(degidx)  
    error('B must be a constant.');
  end
  b = combine(b);
  bcoef = reshape(b.coefficient,szb(1),szb(2));
  c = times(a,1./bcoef);  
    
else
  
  error('Matrix dimensions must agree');

end









