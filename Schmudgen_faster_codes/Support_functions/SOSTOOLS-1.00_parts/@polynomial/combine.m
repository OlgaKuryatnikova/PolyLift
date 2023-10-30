function b = combine(a)
% function B = combine(A)
%
% DESCRIPTION 
%   Combine terms and eliminate unnecessary variables.
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: polynomial, the result of reduction.
%  
% SYNTAX 
%   B=combine(A);
  
% 10/29/2002: PJS  Wrote this wrapper for PVcombine

% Promote a to polynomial 
a = polynomial(a);

% Get polynomial fields/dimensions
adeg = a.degmat;
avar = a.varname;
nta = size(adeg,1);
nva = size(adeg,2);
[nra,nca] = size(a);
acoef = reshape(a.coefficient,nta,nra*nca);

if isempty(a)
  
  b = a;
  return;
  
elseif nva==0
      
  % Add terms of a constant polynomial
  coefficient = sum(acoef,1);
  coefficient = coefficient(:);
  degmat = sparse(0);
  varname = {};
  matdim = [nra nca];
  b = polynomial(coefficient,degmat,varname,matdim);
  return
    
else
  
  
  % Find unique set of variables
  if nva>1    
    [newdeg,uniquevar] = PVuniquevar(adeg,avar);
  else
    newdeg = adeg;
    uniquevar = avar;
  end

  % Find unique set of monomials
  if nta>1
    [newcoef,uniquedeg] = PVuniqueterm(acoef,newdeg,[nra nca]);
  else
    newcoef = acoef;
    uniquedeg = newdeg;
  end
    
  % Eliminate terms with zero coefficient
  [ridx,cidx] = find( newcoef );
  ridx = unique(ridx);
  coefficient = newcoef(ridx,:);
  degmat = uniquedeg(ridx,:);
  if isempty(coefficient)
    coefficient = sparse(1,nra*nca);
    degmat = sparse(0);
  end
  
  % Eliminate variables that no longer exist
  %[ridx,cidx] = find( degmat );
  %cidx = unique(cidx);
  cidx=1:size(uniquevar,1);
  cidx=cidx';
  varname = uniquevar(cidx);
  degmat = degmat(:,cidx);
  nut = size(degmat,1);
  nuv = size(degmat,2);
  
  % Check if polynomial is constant
  if nuv==0
    coefficient = sum(coefficient,1);
    degmat = sparse(0);
    varname = {};
  end

  % Create the polynomial
  coefficient = coefficient(:);
  matdim = [nra nca];
  b = polynomial(coefficient,degmat,varname,matdim);
  return
  
end







