function out=fieldnames(a)
% function out=fieldnames(a)
%
% DESCRIPTION 
%   Returns the public data fields of a polynomial.
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   Out: public data fields of a polynomial
%  
% SYNTAX 
%   display(A);
  
% 11/12/2002: PJS  Initial Coding
%    This function allows tab completion of a public data fields.
  
aprop = PVget(a,'PropNames');
out = aprop.GPropNames;
