function b = PVsubsasgn_colon(a,L,RHS)
% function B =  PVsubsasgn_colon(A,L,RHS)
%
% DESCRIPTION 
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects. 
%   
% INPUTS 
%   A: polynomial 
%   L: a structure array with the fields:
%    type -- string containing '()' 
%    subs -- Cell array or string containing ':'
%   RHS: Value to be assigned. 
%
% OUTPUTS  
%   B: object after subsassignment
%  
% SYNTAX 
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding

  
% Get info about polynomials
b = a;
temp = RHS;
[nrb,ncb]=size(b);    
ntb = size(b.degmat,1);
[nrt,nct]=size(temp);        
ntt = size(temp.degmat,1);

if nrt*nct == nrb*ncb
  b = temp;
elseif [nrt nct]==[1 1] 
  if ~isempty(b)
    b = temp;
    bcoef = reshape(b.coefficient,ntt,1);
    bcoef = sprepmat(bcoef,1,nrb*ncb);
    b.coefficient = bcoef(:);
    b.matdim = [nrb ncb];
  end
else
  error(['In an assignment  A(:) = B, the number of' ...
	 ' elements in A and B must be the same.']);
end

 