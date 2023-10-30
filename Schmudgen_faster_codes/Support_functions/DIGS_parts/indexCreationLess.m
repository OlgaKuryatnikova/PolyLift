% APPS Solver
% Version 1.00, June 26, 2011
%
% Written By: Bissan Ghaddar and Juan Vera

function [alpha, numAlpha] = indexCreationLess(n,d,maxDegree)
%Iterates indexCreationEqu for each degree between 0 and d. See what
%indexCreationEqu does for more details.
%
%[alpha, numAlpha] = indexCreationLess(n,d,maxDegree)

if nargin == 2
    maxDegree = d*ones(1,n);
end

numAlpha = 0;
alpha = [];
for i = 0:d
    [alphai, ni] = indexCreationEqu(n,i,maxDegree);
    alpha = [alpha; alphai];
    numAlpha = numAlpha+ni; 
end
    
