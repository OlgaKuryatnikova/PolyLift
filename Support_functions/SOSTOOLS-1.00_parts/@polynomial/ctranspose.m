function b = ctranspose(a)
% function B=ctranspose(A)
%
% DESCRIPTION 
%   Returns the complex conjugate transpose of a polynomial. 
%   This currently assumes that the polynomial variables are real.
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: transpose of A.
%  
% SYNTAX 
%   B = ctranspose(A);
%   B = A';

% 6/14/2002: PJS  Initial Coding  

b=transpose(a);


