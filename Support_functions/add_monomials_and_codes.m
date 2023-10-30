function [monVec, numMon, codes] = add_monomials_and_codes(numVars,powers,poly,dmax)
% Add the matrix of monimials up to degree dmax and codes to match them with the main vector
% of monomials
%
% Outputs:
% monVec: matrix of monimials
% numMon: number of monimials
% codes: codes for those monomials
% Inputs:
% numVars: number of variables
% powers: codes used as a basis
% poly: if we need the codes for a given polynomial, provide the polynomial
% dmax: if we need the codes for a new polynomial, provide its max degree

if isempty(poly)
    % Monomials
    numer = 0;
    deg = 1;
    [monVec, numMon]= indexCreationLess(numVars,dmax);
else
    numer=0;
    dmax = 1; % redundant here, just need to avoid dmax==0 condition
    % account for possibly numeric entries for poly
    if isnumeric(poly)
        polyTemp=poly;
        clearvars poly
        poly.coef=polyTemp;
        poly.degmat=sparse(1,numVars);
        numer=1;
        deg=0;
    else
        deg=poly.maxdeg; % we need deg since it can happen that a constant polynomial is considered a polynomial, not a number
    end
    monVec=[poly.degmat(end,:);poly.degmat(1:end-1,:)]; % move the all-zeros line to the first row to avoid losing it in further matching
    numMon = size(monVec,1);
end

% Codes for the monomials
if numer == 1 || deg == 0 || dmax == 0
    codes=uint64(0); % for a constant polynomial the resulting codes contain only a zero for each "orbit"
else
    [temp1, temp2, temp3]=find(monVec);
    if size(temp1,2)>1
        temp1=temp1';
        temp2=temp2';
        temp3=temp3';
    end
    codes=powers;
    codes=codes(temp2,:);
    codes=uint64(temp3).*codes;
    codes=accumarray(temp1,codes(:),[], @(x) sum(x,'native'));
end

if ~isempty(poly)
    codes = [codes;codes(1)];
    codes(1) = [];
end
end