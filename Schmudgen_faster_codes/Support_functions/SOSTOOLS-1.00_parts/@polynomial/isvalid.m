function [flag,errormsg] = isvalid(a)
% function [Flag,Errormsg] = isvalid(A)
%
% DESCRIPTION 
%   Check if A is a valid polynomia object.
%   
% INPUTS 
%   A: a polynomial
%
% OUTPUTS  
%   Flag: Flag=1 if polynomial is valid and 0 otherwise.
%   ErrorMsg: if Flag=0, this is a message indicating the 
%      type of error in the polynomial construction.
%  
% SYNTAX 
%   [Flag,Errormsg] = isvalid(A);

% 6/9/2002: PJS  Initial Coding

astruct = struct(a);
[flag,errormsg] = PVisvalid(astruct);




