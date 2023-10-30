function b = PVsubsasgn_1idx(a,L,RHS)
% function B =  PVsubsasgn_1idx(A,L,RHS)
%
% DESCRIPTION 
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects. 
%   
% INPUTS 
%   A: polynomial 
%   L: a 1x1 structure with the fields:
%    type -- string containing '()' 
%    subs -- 1x1 cell array containing indices
%   RHS: Value to be assigned. 
%
% OUTPUTS  
%   B: object after subsassignment
%  
% SYNTAX 
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding

% Check Indices
subsidx = L(1).subs{1};    
ls = length(subsidx);
if min(subsidx)<1
  
  error('Index into matrix is negative or zero.');

end
  
% Get info about polynomials
b = a;
[nrb,ncb]=size(b);    
ntb = size(b.degmat,1);

temp = RHS;
[nrt,nct]=size(temp);        
ntt = size(temp.degmat,1);

% Fan out scalar RHS to correct dimension
if [nrt nct]==[1 1]
  temp.coefficient = sprepmat(temp.coefficient,1,ls);
  temp.coefficient = temp.coefficient(:);
  nrt = 1; 
  nct = ls;
  temp.matdim = [nrt nct];
end
  
% Check dimensions
if ls~=nrt*nct
  error(['In an assignment  A(I) = B, the number of' ...
	 ' elements in B and I must be the same.']);
end
if (min([nrb ncb])>1) &  (max(subsidx) > nrb*ncb)
  error(['In an assignment  A(I) = B with A a matrix, the maximum' ...
	 ' index must be less than the number of elements of A.']);
end

% Embed temp in a polynomial of right dimension
if isempty(b)
  nc = max(subsidx(:));  
  nr = 1;
elseif min([nrb ncb])==1 & ncb >= nrb
  nc = max([ncb; subsidx(:)]);  
  nr = 1;
elseif min([nrb ncb])==1 
  nc = 1;
  nr = max([nrb; subsidx(:)]);  
else
  nc = ncb;
  nr = nrb;
end  
tempcoef = [];
temp.coefficient = reshape(temp.coefficient,ntt,nrt*nct);
for i1 = 1:ntt
  z = sparse(1,nr*nc);
  z(subsidx) = temp.coefficient(i1,:);
  tempcoef = [tempcoef; z];
end
temp.coefficient = tempcoef;
temp.matdim = [nr nc];

if isempty(b)
  b=temp;
else
  % Zero out terms of b and add back in temp
  b.coefficient = reshape(b.coefficient,ntb,nrb*ncb);
  bcoef = [];
  for i1 = 1:ntb    
    tempcoef = sparse(1,nr*nc);
    tempcoef = b.coefficient(i1,:);
    tempcoef(subsidx) = zeros(ls,1);
    bcoef = [bcoef; tempcoef];
  end
  b.coefficient = bcoef;
  b.matdim = [nr nc];
  b = b+temp;
  b = combine(b);
end

  
  





