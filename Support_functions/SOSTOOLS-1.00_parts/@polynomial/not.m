function b = not(a)
% function B=not(A)
%
% DESCRIPTION 
%  Logical Not.
%  Returns a matrix, B, of the same size as A.  An entry of 
%   B is one if the corresponding entry of A is the zero 
%   polynomial and zero otherwise.  The logical not of an 
%   empty polynomial returns the empty polynomial.
%   
% INPUTS 
%   A: polynomial
%
% OUTPUTS  
%   B:  matrix of 1's and 0's.
%    
% SYNTAX 
%   B=~A  
%   B=not(A);
%    Functional form 
  
% 6/14/2002: PJS  Initial Coding  

  
if isempty(a)
  b =[];
else
  sza = size(a);
  z =polynomial(zeros(sza));
  b = (a==z);
end


