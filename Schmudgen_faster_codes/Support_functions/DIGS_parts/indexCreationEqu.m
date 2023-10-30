% APPS Solver
% Version 1.00, June 26, 2011
%
% Written By: Bissan Ghaddar and Juan Vera


function [alpha, numAlpha] = indexCreationEqu(n,d,maxDegree)
%alpha is a matrix where every column is a variable and every row is a
%possible vector such that the entries in that vector gives
%the degree of each variable such that the total degree for the vector is d
%and each variable can't exceed maxDegree (that is a vector of maximum
%degrees). 
%Examples: indexCreationEqu(2,1,[1,1]) = [1,0; 0,1]
%          indexCreationEqu(2,2,[1,1]) = [1,1]
%[alpha, numAlpha] = indexCreationEqu(n,d,maxDegree)

if nargin == 2
%     maxDegree = d;
    maxDegree = d*ones(1,n); %Davide
end

if d ==0
    numAlpha = 1;
    % alpha = zeros(1,n); %Juan
    alpha = sparse(1,n); %Bissan
elseif n == 1
    if maxDegree < d
        numAlpha = 0;
        alpha = [];
    else
        numAlpha =1;
        alpha = d;
    end
else   
    alpha = [];
    numAlpha = 0;
    for i = 0:min(d,maxDegree(1))
        [alphai,ni] = indexCreationEqu(n-1,d-i,maxDegree(2:end));
        if ~isempty(alphai)
            alpha = [ i*ones(size(alphai,1),1) alphai; alpha];
            numAlpha = numAlpha + ni;
        end
    end
end
