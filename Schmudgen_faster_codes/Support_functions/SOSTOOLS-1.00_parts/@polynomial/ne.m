function c = ne(a,b)
% function C=ne(A,B)
%
% DESCRIPTION 
%  Not Equal.
%  Element-by-element comparison between 2 polynomial
%   matrices, A and B, where the result, C, is a matrix of the 
%   same size as A and B.  An entry of C is one if the
%   the corresponding entries of A and B are not the same
%   polynomial and zero otherwise.  
%   
% INPUTS 
%   A,B: polynomials
%
% OUTPUTS  
%   C:  matrix of 1's and 0's.
%    
% SYNTAX 
%   A~=B  
%    A and B must have the same dimensions unless one is a
%    scalar. A scalar can be compared with anything.
%   C=ne(A,B);
%    Functional form 
  
% 6/14/2002: PJS  Initial Coding  

c = (a==b);
c = ~c;


