function b = power(a,n)
% function B=power(A,N)
%
% DESCRIPTION 
%   Element-by-element powers
%   
% INPUTS 
%   A: polynomial 
%   N: matrix of natural numbers
%
% OUTPUTS  
%   B: polynomial, result of matrix power
%  
% SYNTAX 
%   B = A.^N
%     Element-by-element powers. A and B must have 
%     the same dimensions unless one is a scalar. A scalar 
%     can operate into anything.
%   B=power(A,N)
%     Function-call form 

% 6/7/2002: PJS  Initial Coding  

% Check number of inputs  
error(nargchk(2,2,nargin));
  
sza=size(a);
szn=size(n);

if isempty(a) | isempty(n)

  if isempty(a) & isempty(n)
    % empty.^empty = empty     
    b=polynomial;
    return;    
  elseif isempty(a) & max(size(n))==1
    % empty.^scalar = empty
    b=polynomial;
    return;
  elseif isempty(n) & max(size(a))==1
    % scalar.^empty = empty
    b=polynomial;
    return;
  else
    % matrix.^empty  or empty.^matrix     
    error('Matrix dimensions must agree.');
  end
    
elseif sza==szn
  
  % Matrix .^ Matrix
  % XXX This is pretty slow...
  b = a; 
  L.type = '()';
  for i1 = 1:sza(1);
    for i2 = 1:sza(2);
      L.subs = {[i1] [i2]};
      aij = subsref(a,L);
      bij = mpower(aij,n(i1,i2));
      b = subsasgn(b,L,bij);
    end
  end

  
elseif all(szn==[1 1])
  
  % Matrix.^Scalar: 
  b = polynomial(1);
  for i1 = 1:n
    b = b.*a;
  end
  
elseif all(sza==[1 1])
  
  % Scalar .^ Matrix:     
  % XXX This is pretty slow...
  b = polynomial(zeros(szn));; 
  L.type = '()';
  for i1 = 1:szn(1);
    for i2 = 1:szn(2);
      L.subs = {[i1] [i2]};
      bij = mpower(a,n(i1,i2));
      b = subsasgn(b,L,bij);
    end
  end

else
  
  error('Matrix dimensions must agree');

end


  
  
  
   