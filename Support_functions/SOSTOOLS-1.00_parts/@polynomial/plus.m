function c = plus(a,b)
% function C=plus(A,B)
%
% DESCRIPTION 
%   Add two polynomials and combine common terms.
%   
% INPUTS 
%   A,B: polynomials
%
% OUTPUTS  
%   C: polynomial, the result of addition
%  
% SYNTAX 
%   C= A+B
%     Adds A and B.  A and B must have the same dimensions
%     unless one is a scalar.  A scalar can be added to 
%     anything.
%   C = plus(A,B)
%     Function-call form 
  
% 6/7/2002: PJS  Initial Coding  
% 6/8/2002: PJS  Allow matrices of polynomials

% Promote a to polynomial 
a = polynomial(a);
sza = size(a);

% Promote b to polynomial 
b = polynomial(b);
szb = size(b);

if isempty(a) | isempty(b)

  if isempty(a) & isempty(b)
    % empty+empty = empty
    c=polynomial;
    return;
  elseif isempty(a) & max(size(b))==1
    % scalar+empty = empty
    c=polynomial;
    return;
  elseif isempty(b) & max(size(a))==1
    % empty+scalar = empty
    c=polynomial;
    return;
  else
    % matrix+empty  or empty+matrix     
    error('Matrix dimensions must agree.');
  end
    
elseif all(sza==szb) 
  
  % Matrix + Matrix

  % Get Dimensions
  nta = size(a.degmat,1);  
  nva = length(a.varname); 
  ntb = size(b.degmat,1);  
  nvb = length(b.varname);    
  
  % Stack up coefficients
  acoef = reshape(a.coefficient,nta,sza(1)*sza(2));
  bcoef = reshape(b.coefficient,ntb,szb(1)*szb(2));  
  coefficient = [acoef; bcoef];
  coefficient = coefficient(:);
  
  if nva==0 & nvb==0 
    degmat = zeros(nta+ntb,1);
    varname = {};
  else
    adeg = a.degmat(:,1:nva);
    bdeg = b.degmat(:,1:nvb);
    degmat = spblkdiag(adeg,bdeg); 
    varname = [a.varname(:); b.varname(:)];
  end
  
  % Form polynomial and combine terms
  c = polynomial(coefficient,degmat,varname,sza);
  c = combine(c);
  
elseif all(sza==[1 1]) | all(szb==[1 1])  
  
  % Scalar + Matrix / Matrix + Scalar
  
  % Make first term of sum be the scalar.  
  if ~all(sza==[1 1])
    temp = a;
    a = b;
    
    b = temp;
    szb = sza;
  end

  % Turn into matrix+matrix
  %   (Due to a bug in Matlab 6.5, the following yields an error:
  %       >>b(20) =sparse(1);
  %    Consequently, repmat fails to work for sparse matrices....
  
  acoef = sprepmat(a.coefficient,1,szb(1)*szb(2));
  a.coefficient = acoef(:);
  a.matdim = szb;
  c = plus(a,b);
  
else
  
  error('Matrix dimensions must agree');

end
