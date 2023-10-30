function p = polynomial(coefficient,degmat,varname,matdim);
% function P = polynomial(Coefficient,Degmat,Varname,Matdim);
%
% DESCRIPTION 
%   Creates a polynomial or a matrix of polynomials.  Redundant
%   terms will be combined.
%   
% INPUTS 
%   Coefficient: coefficients of each monomial.
%   Degmat: degrees of each monomial
%   Varname: names of variables
%   Matdim: dimensions of the polynomial matrix  
%
% OUTPUT 
%   P: polynomial object
%
% SYNTAX 
%   P=polynomial  
%     Creates an empty polynomial object.
%   P=polynomial(Coefficient)
%     If Coefficient is a real matrix of dimension NxM, then P is
%     an NxM constant polynomial. 
%   P=polynomial(Coefficient)
%     If Coefficient is a polynomial object, then P=Coefficient.
%   P=polynomial(Coefficient,Degmat,Varname,Matdim,Nterms,Nvars)  
%     If P is an NxM polynomial that is the sum of T terms in V
%     variables, the inputs should be specified as: 
%       Coefficient is a (T*N*M)x1 sparse matrix.
%          The coefficient of the t^th term in the (i,j) entry of
%          the polynomial matrix is stored in the 
%          t + T*(i-1)  + T*N*(j-1)  entry of Coefficient.
%          i.e. [(1,1) coefs; (2,1) coefs; ...; (N,M) coefs]
%       Degmat is a TxV sparse matrix of natural numbers.  Row t 
%          gives the degrees of each variable for the t^th term.
%       Varname is a Vx1 cell array with entry v giving the name of
%          variable v. For a constant polynomial, varname is an
%          empty 1x1 cell.  
%       Matdim is the vector of matrix dimensions, [N M].
  
% 6/7/2002:   PJS  Initial Coding
% 6/8/2002:   PJS  Allow Matrices of Polynomials
% 10/30/2002: PJS  Input is structure not cell arrays
% 11/5/2002:  PJS  Coefficients stored as sparse 3D array.
  
superiorto('double');

%--------------------------------------%
%    Fill in empty arguments
%--------------------------------------%
if nargin==1 & isa(coefficient,'polynomial')  
  % If input is already a polynomial, pass it back  
  p = coefficient;
  return;
elseif nargin==1 & isreal(coefficient) & ndims(coefficient)==2
  p.coefficient = sparse(coefficient(:));
  p.degmat = sparse(0);
  p.varname = {};
  p.matdim = size(coefficient);
  p=class(p,'polynomial');
  return  
elseif nargin == 0
  p.coefficient = sparse([]);
  p.degmat = sparse([]);
  p.varname = {};
  p.matdim = [0 0];
  p=class(p,'polynomial');
  return
elseif nargin~=4
  error(['Not enough inputs: P = polynomial(Coefficient,Degmat,'...
	 'Varname,Matdim)']);
end

%--------------------------------------%
%       Error Checking
%--------------------------------------%
p.coefficient = coefficient(:);
p.degmat = degmat;
p.varname = varname;
p.matdim = matdim(:)';
[flag,errormsg] = PVisvalid(p);
if flag==1
  % p is a  polynomial object
  p.coefficient = sparse(p.coefficient);
  p.degmat = sparse(p.degmat);
  p.varname = p.varname(:);
  p.matdim = matdim(:)';
  p=class(p,'polynomial');
  %p=combine(p);
else
  % p is not a valid polynomial object.
  error(errormsg);
end








