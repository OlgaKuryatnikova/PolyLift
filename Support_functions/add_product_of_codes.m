function [rowVec,colVec,valVec] = add_product_of_codes(poly,codesP,numMonP,codesS,numMonS,codesFull,col_begin,col_end,svec_ind,svecTransform)
% Create inputs rowVec,colVec,valVec for a sparse matrix M of the size of 
% the total number of monomials in the problem such that we can write 
% M*svecVar = coefficients(p*S,x), where p(x) is an existing polynomial,
% and S(x) is a new added polynomial with coefficients svecVar.  
% The matrix is sparse since only few of all possible monomials will
% appear in p*S while we need to match them with the full monomial vector.
%
% Inputs: 
% poly: polynomial
% codesP: codes of p obtained from the function add_monomials_and_codes
% numMonP: number of monomials in p (this is just length(codesP))
% codesS: codes of S obtained from the function add_monomials_and_codes
% numMonS: number of monomials in S (this is just length(codesS))
% codesFull: the full vector of codes for all possible monomials
% col_begin,col_end: the first and last columns which M will take
% in the total big matris for all p and S
% svec_ind: do we need an svec transformation?
% svecTransform: the vector for the svec transformation if we need it

if isnumeric(poly)
    polyTemp = poly;
    clearvars poly
    poly.coef = polyTemp;
end

temp1  = kron(codesS,uint64(ones(numMonP,1))); % repeat each monomial in the new polynomial for all monomials in the given polynomial
temp2 = repmat(codesP,numMonS,1); % repeat all monomials in the given polynomial as many times as we have new monomials
temp3 = temp1+temp2; % represent multiplication of all monomials for match
[~,rowVec] = ismember(temp3,codesFull); % find the sequential numbers of rows in the basic "unique" monomial matrix that correspond to the new monomials of the product of polynomials

colVec = kron(col_begin:col_end,ones(1,numMonP)); % repeat sequential number of each monomial of the new polynomial for all monomials in the given polynomial (the number of new monomials)
if svec_ind == 1 % if we deal with an SOS and use the svec transformation
    valVec = repmat(poly.coef',1,numMonS).*kron(svecTransform,ones(1,numMonP));% create the vector of values
else % otherwise
    valVec = repmat(poly.coef',1,numMonS);% create the vector of values
end
