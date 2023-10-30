function pout = newpvar(varargin)
%
% DESCRIPTION 
%   Create variables (i.e. monomials of degree 1).
%   
% INPUTS 
%   X1,X2,...: variable names
%
% OUTPUTS  
%   None
%  
% SYNTAX 
%   newpvar('x1','x2','x3',...)  
%   newpvar x1 x2 x3  
%     Both of these function calls create monomials of degree 1 with
%     the given names.
  
% 15/08/2009 DONE by Juan and Bissan
% now has output!
 
n=length(varargin);
pout = [];
for i = 1:n
	p.coef= sparse(n,1);
	p.coef(i,1) = 1;
	p.degmat= sparse(n,n);
	p.degmat(i,i) = 1;
	p.varname = varargin;
	p.matdim=[1,1]; 
	pp = polynomial(p.coef,p.degmat,p.varname,p.matdim);
	assignin('caller',varargin{i},pp);
	pout = [pout,pp];
%	assignin('caller',varargin(i), p);
end
