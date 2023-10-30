function b = isempty(a)
% function B = isempty(A)
%
% DESCRIPTION 
%   Returns true for empty polynomials.
%   
% INPUTS 
%   A: polynomial
%
% OUTPUTS  
%   B: 1 if A is an empty polynomial and 0 otherwise.
%    
% SYNTAX 
%   B = isempty(A);
  
% 6/14/2002: PJS  Initial Coding  

if size(a)==[0 0]
  b = 1;
else
  b = 0;
end

