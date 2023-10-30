function varargout = mpvar(n,m,var);
% function P = mpvar(N,M,Var);
%
% DESCRIPTION 
%   Create a polynomial matrix or vector variable
%   
% INPUTS 
%   N,M: row and column dimensions of polynomial matrix.
%   var: 
%
% OUTPUTS  
%   P: polynomial matrix 
%  
% SYNTAX 
%   P = mpvar(N,M,var)
%       Creates an NxM polynomial matrix with entries
%       var_i_j.  For example, P = mpvar(2,3,'p') creates:
%           P = [ p_1_1 , p_1_3 , p_2_2 ]
%               [ p_1_2 , p_2_1 , p_2_3 ]
%   mpvar(N,M,var)
%      Equivlavent to calling eval([var '=mpvar(N,M,var);']).
  
% 11/10/2002 PP    Initial Coding
% 11/10/2002 PJS   Minor modifications for speed and allow 
%                  calling sequence with no outputs
  
  
% Argument checking
error(nargchk(3,3,nargin))
error(nargchk(0,1,nargout))

% Form Coefficients
nt = n*m;
coefficient = speye(nt);
coefficient = coefficient(:);

% The degrees
degmat = speye(nt);

% Create the names
stridx = int2str( (1:max([n m]))' );
stridx = strjust( stridx ,'left');
stridx = cellstr(stridx);
varname = cell(nt,1);
for i1=1:n ;
  tempvar = [var '_' stridx{i1} '_'];
  for i2 =1:m;
    varname{i1+(i2-1)*n} = [tempvar stridx{i2}];    
  end
end

% Create matrix dimension
matdim = [n m];

P = polynomial(coefficient,degmat,varname,matdim);
if nargout==0
  assignin('caller', var, P);
else
  varargout{1} = P;
end
