function  [M,numMonS] = add_coef_matr(poly,dmax,degS,numVars,powers,codesFull,sos_ind,var_subsets)
% Create a sparse matrix M of the size
% "the total number of monomials in the problem" x numMonS.
% M*svecVar = coefficients(poly*S,x), where poly(x) is an existing polynomial,
% and S(x) is a newly added polynomial with coefficients svecVar.  
% The matrix is sparse since only few of all possible monomials will
% appear in poly*S while we need to match them with the full monomial vector.
%
% Inputs: 
% poly: polynomial
% dmax: maximal degre in the problem
% degS: degree of the added polynomial, if does not follow immediately from
% the hierarchy
% numVars: the number of variables in the problem
% powers: codes for each variable, used as a basis
% codesFull: the full vector of codes for all possible monomials
% sos_ind: are we adding an sos coefficient or a free one? 1 means sos
% var_subsets: if the newly added polynomial depends only on a subset of
% variables, provide the vector of the indices of this subset
% Olga Kuryatnikova, 2023

if isempty(degS) && ~isnumeric(poly) % if the degree is not filled in, fill the degree of the input polynomial
    degS = dmax-poly.maxdeg;
elseif isempty(degS) && isnumeric(poly)
    degS = dmax;
end

% Add the codes for the existing polynomial
[~, numMonP, codesP] = add_monomials_and_codes(numVars,powers,poly);

% New monomials and their codes, correspond to the polynomial coefficient 
if nargin <= 7
    var_subsets = 1:numVars;
end
numVarsS = length(var_subsets);
powersS = powers(var_subsets);

if sos_ind == 1 % is the added polynomial coefficient an sos?
    degS=floor(degS/2);
    [~, numMonS, codesS] = add_monomials_and_codes(numVarsS,powersS,[],degS);
    ind_triu = triu(true(numMonS),1);
    [codesS1,codesS2] = meshgrid(codesS,codesS);
    codesS =  codesS1+codesS2;
    codesS = [codesS(ind_triu);codesS(1:numMonS+1:numMonS^2)'];
    svec_ind = 1;
    numOffDiag = numMonS*(numMonS-1)/2;
    svecTransform = [sqrt(2)*ones(1,numOffDiag),ones(1,numMonS)];
    numMonS = numMonS*(numMonS+1)/2;
else
    [~, numMonS, codesS] = add_monomials_and_codes(numVarsS,powersS,[],degS);
    svec_ind = 0;
    svecTransform = [];
end

% create the matruix of constraints
[rowVec,colVec,valVec] = add_product_of_codes(poly,codesP,numMonP,codesS,numMonS,codesFull,1,numMonS,svec_ind,svecTransform);
numMonFull = length(codesFull);
M = sparse(rowVec,colVec,valVec,numMonFull,numMonS);
end
