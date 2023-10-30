function varargout = subsref(a,L)
% function B=subsref(A,L);
%
% DESCRIPTION 
%   Subsreference for polynomial objects. 
%   
% INPUTS 
%   A: polynomial 
%   L: a structure array with the fields:
%    type -- string containing '()'  or '.' specifying the
%              subscript type.
%    subs -- Cell array or string containing the actual subscripts.
%
% OUTPUTS  
%   B: object after referencing
%  
% SYNTAX 
%   B=subsref(A,L)
  
% 6/7/2002: PJS  Initial Coding  

switch L(1).type
 case '.'

  b = get(a,L(1).subs);
  
 case '()' 

  nta = size(a.degmat,1);  
  [nra,nca]=size(a);
  acoef = a.coefficient;
  b = a;
    
  if length(L(1).subs)==1        
    % Check Indices
    subsidx = L(1).subs{1};
    if ischar(subsidx) & subsidx == ':'
      subsidx = (1:nra*nca)';
    end
    
    if max(subsidx)>nra*nca | min(subsidx)<1
      error('Index exceeds matrix dimensions.');
    end

    % Do subsref
    bcoef = reshape(acoef,nta,nra*nca);
    bcoef = bcoef(:,subsidx);
    b.coefficient = bcoef(:);
    b.matdim = size(subsidx);
    
  elseif length(L(1).subs)==2
    ridx = L(1).subs{1};
    if ischar(ridx) & ridx==':'
     ridx = 1:nra;
    end
    cidx = L(1).subs{2};
    if ischar(cidx) & cidx==':'
      cidx = 1:nca;
    end
    
    % Check Indices
    if max(ridx)>nra | min(ridx)<1 | ...
	  max(cidx)>nca | min(cidx)<1
      error('Index exceeds matrix dimensions.');
    end

    % Do subsref
    bcoef = reshape(acoef,nta,nra*nca);
    idx = [];
    for i1 = 1:length(cidx)
      idx = [idx ridx+(cidx(i1)-1)*nra];
    end
    
    bcoef = bcoef(:,idx);
    b.coefficient = bcoef(:);
    b.matdim = [length(ridx) length(cidx)];
      
  else
      error('Invalid subsref for polynomials');
  end
  
  % Combine down
  b = combine(b);  

 case '{}'
    
  error('{}- like subsreference is not supported for polynomial objects.');
  
end


if length(L)>1  
  b = subsref(b,L(2:end));
end

varargout{1}=b;


