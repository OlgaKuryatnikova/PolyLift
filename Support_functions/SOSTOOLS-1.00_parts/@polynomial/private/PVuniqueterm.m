function [newcoef,uniquedeg] = PVuniqueterm(coef,degmat,matdim)
% function [NewCoef,UniqueDeg] = PVuniqueterm(Coef,Degmat,MatDim);
%
% DESCRIPTION 
%   (INTERNAL FUNCTION)  
%   Find the unique set of terms (i.e. unique rows of Degmat
%   and form coef in terms of this unique set.
%   
% INPUTS 
%   Coef: coefficient matrix for a polynomial
%   Degmat: degree matrix for a polynomial
%   MatDim  dimensions of the polynomial matrix
%
% OUTPUTS  
%   NewCoef: coefficient matrix in terms of new degree matrix
%   UniqueDegmat: degree matrix with unique rows
%  
% SYNTAX 
%   [NewCoef,UniqueDeg] = PVuniqueterm(Coef,Degmat)  

% 11/25/2002: PJS  Initial Coding
%  1/30/2002: PJS  Combine rows using sparse multiply

% Find repeated monomials
%   sortidx: index such that degsort = degmat(sortidx)
%   repeatidx: nonzero indices indicate repeats
%   rvec:  map rows of degsort to numbers 1,..,nut
[degsort,sortidx] = sortrows(degmat);
repeatidx = all(degsort(2:end,:) == degsort(1:end-1,:),2);    
rvec = cumsum([1; ~repeatidx]);

% Create unique set of monomials
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniquedeg = degmat(uniqueidx,:);

% Sum repeated terms
%  (Use sparse multiply, summat*coef, to sum the rows)
nt = size(degmat,1);
nut = size(uniquedeg,1);
summat = spalloc(nut,nt,nt);
idx = sub2ind([nut nt],rvec,sortidx);
summat(idx) = 1;
newcoef = summat*coef;

% Put terms in "alphabetical" order
% (not necessary, but might look better)
uniquedeg = flipud(uniquedeg);
newcoef = flipud(newcoef);





