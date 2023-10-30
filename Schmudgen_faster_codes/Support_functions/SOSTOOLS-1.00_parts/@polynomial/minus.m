function c = minus(a,b)
% function C=minus(A,B)
%
% DESCRIPTION 
%   Subtract two polynomials and combine common terms.
%   
% INPUTS 
%   A,B: polynomials 
%
% OUTPUTS  
%   C: Result of subtraction, A-B
%  
% SYNTAX 
%   C= A-B
%     Subtracts B from A.  A and B must have the same dimensions
%     unless one is a scalar.  A scalar can be subtracted from
%     anything.  
%   C = minus(A,B)
%     Function-call form for substraction.
  
% 6/7/2002: PJS  Initial Coding  
% 6/8/2002: PJS  Allow matrices of polynomials
% 6/9/2002: PJS  Simplified code to use plus/uminus, probably need to
%                rewrite for speed.
  
c = a+(-b);
