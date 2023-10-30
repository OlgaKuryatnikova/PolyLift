function c = mtimes(a,b)
% function C=mtimes(A,B)
%
% DESCRIPTION 
%   Multiply two polynomial matrices and combine terms.
%   
% INPUTS 
%   A,B: polynomials
%
% OUTPUTS  
%   C: polynomial, the result of multiplication
%  
% SYNTAX 
%   C= A*B
%     Multiplies A and B.  A and B must have the same dimensions
%     unless one is a scalar.  A scalar can be multiplied to 
%     anything.
%   C = mtimes(A,B)
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
    % empty*empty = empty
    c=polynomial;
    return;
  elseif isempty(a) & max(size(b))==1
    % empty*scalar = empty
    c=polynomial;
    return;
  elseif isempty(b) & max(size(a))==1
    % scalar*empty = empty
    c=polynomial;
    return;
  else
    % matrix*empty  or empty*matrix     
    disp(['Warning: Empty matrix multiplication with unequal' ...
	     ' inner dimensions.']);
    c=polynomial;
    return
  end
    
elseif sza(2)==szb(1) 
  
  % Matrix * Matrix

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
    
    coefficient = acoef*bcoef;
    degmat = 0;
    varname = {};
  else    
    % Get Coefficients: 
    %  acoefcol: coef. matrices stacked vertically
    %  bcoefcol: coef. matrices stacked horizontally
    acoef = reshape(a.coefficient,nta,sza(1)*sza(2));    
    acoefcol = spalloc(nta*sza(1),sza(2),nnz(acoef));    
    for i1 = 1:nta
      ridx = (i1-1)*sza(1) + (1:sza(1));
      acoefcol(ridx,:) = reshape(acoef(i1,:),sza(1),sza(2));
    end
    bcoef = reshape(b.coefficient,ntb,szb(1)*szb(2));
    bcoefcol = spalloc(szb(1),ntb*szb(2),nnz(bcoef));    
    for i1 = 1:ntb
      cidx = (i1-1)*szb(2) + (1:szb(2));
      bcoefcol(:,cidx) = reshape(bcoef(i1,:),szb(1),szb(2));
    end   
    tempcoef = acoefcol*bcoefcol;

    coefficient = spalloc(sza(1),nta*ntb*szb(2),nnz(tempcoef));
    for i1 = 1:nta 
      ridx = (i1-1)*sza(1) + (1:sza(1));
      cidx = (i1-1)*ntb*szb(2) + (1:ntb*szb(2));
      coefficient(:,cidx) = tempcoef(ridx,:);
    end
    coefficient = reshape(coefficient,sza(1)*szb(2),nta*ntb);
    coefficient = coefficient';
    coefficient = coefficient(:);

    % Form Degmat and Varname
    tempdeg = a.degmat(:,1:nva);
    adeg = spalloc(ntb*nta,nva,ntb*nnz(tempdeg));    
    for i1 = 1:nta
      ridx = (i1-1)*ntb + (1:ntb);
      adeg(ridx,:) = tempdeg(i1*ones(ntb,1),:);
    end
    bdeg = b.degmat(:,1:nvb);
    bdeg = sprepmat(bdeg,nta,1);
    degmat = [adeg bdeg];
    varname = [a.varname(:); b.varname(:)]; 
  end
  
  % Form polynomial and combine terms
  c = polynomial(coefficient,degmat,varname,[sza(1) szb(2)]);
  c = combine(c); 
    
elseif all(sza==[1 1]) | all(szb==[1 1])
  
  % Scalar * Matrix / Matrix * Scalar
  
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










