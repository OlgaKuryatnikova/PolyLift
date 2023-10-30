function b = PVsubsasgn_2idx(a,L,RHS)
% function B =  PVsubsasgn_2idx(A,L,RHS)
%
% DESCRIPTION 
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects. 
%   
% INPUTS 
%   A: polynomial 
%   L: a 1x1 structure with the fields:
%    type -- string containing '()' 
%    subs -- 1x2 cell array containing indices
%   RHS: Value to be assigned. 
%
% OUTPUTS  
%   B: object after subsassignment
%  
% SYNTAX 
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding

% Check Indices
ridx = L(1).subs{1};    
lr = length(ridx);
cidx = L(1).subs{2};
lc = length(cidx);
if min(ridx)<1 | min(cidx)<1  
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
  temp.coefficient = sprepmat(temp.coefficient,1,lr*lc);
  temp.coefficient = temp.coefficient(:);
  nrt = lr; 
  nct = lc;
  temp.matdim = [nrt nct];
end

% Check dimensions
if lr*lc~=nrt*nct
  error(['In an assignment  A(I) = B, the number of' ...
	  ' elements in B and I must be the same.']);
end
  
% Embed temp in a polynomial of right dimension
nr = max([nrb; ridx(:)]);
nc = max([ncb; cidx(:)]);  
tempcoef = [];
temp.coefficient = reshape(temp.coefficient,ntt,nrt*nct);
for i1 = 1:ntt
  z = sparse(nr,nc);
  z(ridx,cidx) = reshape(temp.coefficient(i1,:),nrt,nct);
  tempcoef = [tempcoef; (z(:))'];
end
temp.coefficient = tempcoef;
temp.matdim = [nr nc];

if isempty(b)  
  b = temp;
else
  % Zero out terms of b and add back in temp
  b.coefficient = reshape(b.coefficient,ntb,nrb*ncb);
  bcoef = [];
  for i1 = 1:ntb    
    tempcoef = sparse(nr,nc);
    tempcoef(1:nrb,1:ncb) = reshape(b.coefficient(i1,:),nrb,ncb);
    tempcoef(ridx,cidx) = zeros(lr,lc);
    bcoef = [bcoef; (tempcoef(:))'];
  end
  b.coefficient = bcoef;
  b.matdim = [nr nc];
  b = b+temp;
  b = combine(b);  
end

 
  
  





