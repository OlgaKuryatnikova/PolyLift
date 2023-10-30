function c = times(a,b)
% function C=times(A,B)
%
% DESCRIPTION 
%   Element-by-element multiply of two polynomial matrices.  
%   
% INPUTS 
%   A,B: polynomials
%
% OUTPUTS  
%   C: polynomial, the result of multiplication.
%  
% SYNTAX 
%   C= A.*B
%     Multiplies A and B element-by-element.  A and B must 
%     have the same dimensions unless one is a scalar.  A 
%     scalar can be multiplied to anything.
%   C = times(A,B)
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
    % empty.*empty = empty
    c=polynomial;
    return;
  elseif isempty(a) & max(size(b))==1
    % empty.*scalar = empty
    c=polynomial;
    return;
  elseif isempty(b) & max(size(a))==1
    % scalar.*empty = empty
    c=polynomial;
    return;
  else
    % matrix.*empty  or empty.*matrix     
    error('Matrix dimensions must agree.');
    return
  end
    
elseif sza==szb
  
  % Matrix .* Matrix

  % Get Dimensions
  nta = size(a.degmat,1); 
  nva = length(a.varname); 
  ntb = size(b.degmat,1);  
  nvb = length(b.varname);    
  
  if nva==0 & nvb==0 
    a = combine(a);
    acoef = reshape(a.coefficient,sza(1),sza(2));
    b = combine(b);
    bcoef = reshape(b.coefficient,szb(1),szb(2));
    
    coefficient = acoef.*bcoef;
    degmat = 0;
    varname = {};
  else    
    
    % Get Coefficients 
    acoef = reshape(a.coefficient,nta,sza(1)*sza(2));
    acoef = sprepmat(acoef,ntb,1);

    temp = reshape(b.coefficient,ntb,szb(1)*szb(2));    
    bcoef =spalloc(ntb*nta,szb(1)*szb(2),nta*nnz(temp));
    for i1 = 1:ntb
      idx = (1:nta)+(i1-1)*nta;
      bcoef(idx,:) = sprepmat(temp(i1,:),nta,1);
    end
    coefficient = acoef.*bcoef;
    
    % Form Degmat and Varname
    adeg = a.degmat(:,1:nva);
    bdeg = b.degmat(:,1:nvb);
    bdegcol = [];
    for i1 = 1:ntb
      bdegcol = [bdegcol; sprepmat(bdeg(i1,:),nta,1)];
    end
    degmat = [sprepmat(adeg,ntb,1) bdegcol];
    varname = [a.varname(:); b.varname(:)]; 
  end
  
  % Form polynomial and combine terms
  c = polynomial(coefficient,degmat,varname,sza);
  c = combine(c); 
    
elseif all(sza==[1 1]) | all(szb==[1 1])
  
  % Scalar .* Matrix / Matrix .* Scalar
  
  % Make first term be the scalar.  
  if ~all(sza==[1 1])
    temp = a;
    a = b;
    
    b = temp;
    szb = sza;
  end

  % Turn into matrix.*matrix
  acoef = a.coefficient;
  acoef = acoef(:,ones(1,szb(1)*szb(2)));  
  a.coefficient = acoef(:);
  a.matdim = szb;
  c = times(a,b);
    
else
  
  error('Inner matrix dimensions must agree');

end

  










