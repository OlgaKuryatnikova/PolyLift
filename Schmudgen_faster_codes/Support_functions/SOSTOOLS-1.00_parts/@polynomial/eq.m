function c = eq(a,b)
% function C=eq(A,B)
%
% DESCRIPTION 
%  Equal.  
%  Element-by-element comparison between 2 polynomial
%   matrices, A and B, where the result, C, is a matrix of the 
%   same size as A and B.  An entry of C is one if the
%   the corresponding entries of A and B are the same
%   polynomial and zero otherwise.  
%   
% INPUTS 
%   A,B: polynomials
%
% OUTPUTS  
%   C:  matrix of 1's and 0's
%    
% SYNTAX 
%   A==B  
%    A and B must have the same dimensions unless one is a
%    scalar. A scalar can be compared with anything.
%   C=eq(A,B);
%    Functional form 
  
% 6/14/2002: PJS  Initial Coding  

a = polynomial(a);
sza = size(a);

b = polynomial(b);
szb = size(b);

if isempty(a) & isempty(b)
  % Both terms empty:  []==[] returns [] 
  
  c=polynomial;
  return;

elseif isempty(a) & all(szb==[1 1])
  % []==scalar returns empty

  c=polynomial;
  return;

elseif isempty(b) & all(sza==[1 1])
  % scalar==[] returns empty

  c=polynomial;
  return;
  
elseif all(sza==szb)

  d = a-b;
  
  [nrd,ncd] =size(d);
  ntd = size(d.degmat,1);
  dcoef = reshape(d.coefficient,ntd,nrd*ncd);
  dcoef = sum(abs(dcoef),1);
  c = full(dcoef)==0;
  c = reshape(c,nrd,ncd);
  
elseif all(sza==[1 1]) | all(szb==[1 1])
  % Scalar == Matrix:   
  
  % Make first term be the scalar.  
  if ~all(sza==[1 1])
    temp = a;
    a = b;
    
    b = temp;
    szb = sza;
  end

  % Turn into matrix==matrix
  acoef = a.coefficient;
  acoef = acoef(:,ones(1,szb(1)*szb(2)));  
  a.coefficient = acoef(:);
  a.matdim = szb;
  c = eq(a,b);

else  
  error('Matrix dimensions must agree');
end



