function c = horzcat(varargin)
% function C=horzcat(A,B);
%
% DESCRIPTION 
%   Horizontal concatenation of polynomial objects.
%   
% INPUTS 
%   A,B: polynomials 
%
% OUTPUTS  
%   C:  horizontal concatenation of input matrices.
%  
% SYNTAX 
%   [A B] 
%       Horizontal concatenation of polynomial matrices A 
%       and B.  A and B must have the same number of rows.
%   [A1 A2 A3 ...]
%       Horizontal concatenation of several polynomial matrices.  
%   C = horzcat(A1,A2,...);
%       Function-call form of horizontal concatenation.
  
  
% 6/8/2002: PJS  Initial Coding
  
if nargin==1

  c = varargin{1};

else
  
  % Promote a to polynomial
  a = polynomial(varargin{1});
  [nra,nca] = size(a); 
  
  % Promote b to polynomial
  b = polynomial(varargin{2});
  [nrb,ncb] = size(b);
  
  if isempty(b);
  
    c = a;
    
  elseif isempty(a);
    
    c = b;
    
  elseif nra==nrb

    % Get Dimensions
    nta = size(a.degmat,1);
    nva = length(a.varname);
    ntb = size(b.degmat,1);
    nvb = length(b.varname);

    if nva==0 & nvb==0
      % Combine constant terms
      ar = combine(a); 
      coef1 = reshape(ar.coefficient,nra,nca);
      br = combine(b); 
      coef2 = reshape(br.coefficient,nrb,ncb);
            
      % Stack up Coefficients and Form Polynomial   
      coefficient = [coef1 coef2];          
      c = polynomial(coefficient);
    else      
      % Form Degmat, Varname, and Matdim
      adeg = a.degmat(:,1:nva);
      bdeg = b.degmat(:,1:nvb);
      degmat = spblkdiag(adeg,bdeg); 
      varname = [a.varname(:); b.varname(:)];
      matdim = [nra nca+ncb];
      
      % Stack up Coefficients
      coef1 = [reshape(a.coefficient,nta,nra*nca) zeros(nta,nra*ncb)];
      coef2 = [zeros(ntb,nrb*nca) reshape(b.coefficient,ntb,nrb*ncb)];
      coefficient = [coef1; coef2];    
      coefficient = coefficient(:);
            
      % Form Polynomial and combine terms
      c = polynomial(coefficient,degmat,varname,matdim);
      c = combine(c);
    end    
               
  else
    
    error('All columns must have the same row dimension')
    
  end

  if nargin>2
    c = horzcat(c,varargin{3:end});
  end
  
end

