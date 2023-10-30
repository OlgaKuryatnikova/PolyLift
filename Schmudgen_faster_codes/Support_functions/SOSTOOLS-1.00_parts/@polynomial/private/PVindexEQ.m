function [alpha, numAlpha] = PVindexEQ(n,d,maxDegree)
%outputs vector of indexes of n vars of degree = d with maxDegree for each
%variable
%[alpha, numAlpha] = indexCreationEqu(n,d,maxDegree)

if nargin == 2
    maxDegree = d*ones(1,n);
end

if d ==0
    numAlpha = 1;
    alpha = zeros(1,n);
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
        [alphai,ni] = PVindexEQ(n-1,d-i,maxDegree(2:end));
        if ~isempty(alphai)
            alpha = [ i*ones(size(alphai,1),1) alphai; alpha];
            numAlpha = numAlpha + ni;
        end
    end
end
