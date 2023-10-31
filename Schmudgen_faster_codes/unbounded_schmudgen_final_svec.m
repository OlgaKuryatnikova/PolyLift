%% Schmudgen-based certificate with coefficients in constants, dsos, sdsos or sos
%% The case of UNbounded sets.
%% Use some functions from SOSTOOLS 1.00 and DIGS to encode polynomials.
%% Important Comments:
% (1) This code works properly when 
% (total degree + 1)^(number of variables -1) <= 2^64, othrewise we cannot
% match the monomials. This restriction is not present in the Yalmip codes.
% (2) Here we assume that the set is formulated as a subset of the non-negative
% orthant; if this is not the case, we need to translate the possibly
% negative variables into the non-negative orthant using the transformation
% x -> y-z, y>=0, z>=0 before using the case, see the paper for more details.

clear


%% Write here your case following the format of the below example
numVars=5; % the number of variables in the problem, after the transformation to R+ mentioned above, if needed
x = varsVector('x',numVars);
f = 7*(2*x(1)-x(2)+x(3)-2*x(4)-2*x(5)); % the objective function
g = [(7*x(1)-2)^2-49*x(2)^2-(7*x(3)-1)^2-(7*x(5)-1)^2; 49*x(1)*x(3)-49*x(4)*x(5)+49*x(1)^2-1;7*x(3)-49*x(2)^2-49*x(4)^2-1;...
    49*x(1)*x(5)-49*x(2)*x(3)-2;2-sum(x);x']; % vector of inequality constraints, >=0 format
h = []; % vector of equality constraints, =0 format
setType = 3; % choose the type of coefficients: constant term (0), dsos (1), sdsos (2), sos (3)
r = 0; % degree of the hierarchy r in the paper
%

%% Start the contruction
tic

% additional degree to control for numerical issues if needed
deg_aux = 0;

% choose the type of coefficients: constant term (0), dsos (1), sdsos (2), sos (3)
if setType == 0
    dcert = 0;
else
    dcert = 1;
end

% create the vector of degrees and the number of monomials in g
lenG = length(g);
lenGinit = lenG;
Dg = zeros(lenG,1);
numMong = zeros(lenG,1);
for j=1:lenG
    Dg(j) = g(j).maxdeg;
    numMong(j) = length(g(j).coef);
end
numMonMaxg = max([1;numMong]);
Df = f.maxdeg;

% compute the needed additional degrees
dmax = max(max(Dg),ceil(Df/2));
if isempty(h)
    dmax_tot = r+2*dmax+deg_aux;
else
    dmax_tot = max(h.maxdeg,r+2*dmax+deg_aux);
end

% denominator D for the unbounded certificate
deg_D = floor(r/dmax);
deg_temp =  indexCreationLess(lenG,deg_D);
D = 0;
denom = 1+sum(x);
for i=1:size(deg_temp,1)
    beta = deg_temp(i,:);
    multinom_coef = factorial(r)/prod(factorial([r-sum(beta*dmax),beta*dmax]));
    D = D + multinom_coef*denom^(r-beta*Dg)*prod(g.^(beta'));
end
D_full = D*(1+sum(x))^max(0,2*dmax-Df);
f_full = D_full*f;

% add variables x to the set description if they were not there initially, to use for the certificate
for i=1:numVars
    if sum(ismember(x(i),g)) <= 0.5
        g = [g;x(i)];
        Dg = [Dg;1];
    end
end
lenG = length(g);

% generate possible vectors of degrees for Schmudgen terms
monVecG = indexCreationLess(lenG,dmax_tot); % sparse, can be crucial, inside the code can be changed to be full
indTotal = monVecG*Dg;
indBelowDmax = indTotal <= dmax_tot;
numMonG = sum(indBelowDmax);
monVecG = monVecG(indBelowDmax,:);
indTotal = indTotal(indBelowDmax);

if numMonG == 1
    fprintf('There is only one monomial for the certificate, something is wrong, execution stopped. \n')
    return
end

% make sure we do not have a row with all zeros
g = [g;1];
monVecG(1,lenG+1) = 1;

% choose the type of coefficients: constant term (0), dsos (1), sdsos (2), sos (3)
if setType == 0
    dcert = 0;
else
    dcert = 1;
    auxVar = cell(dcert+1,1); % auxiliary variables needed for a certificate
end

% % Start constructing the certificate

% auxiliary input needed to match monomials with each other
powers=uint64(dmax_tot+1).^uint64((1:numVars)'-1);  % base to construct indices for comparison

% create all monomials for dmax and numVars variables
[~, numMon, codesFull] = add_monomials_and_codes(numVars,powers,[],dmax_tot);
varsOrder=(1:numVars)';

% add the right-hand side for f-lambda = certificate, lambda
% adjusted objective function in the unbounded certificate, here and next

%% Finish editing starting from here: equalities the final constraint and
ML = add_coef_matr(D_full,dmax_tot,0,numVars,powers,codesFull,0);
% add the right-hand side for f-lambda = certificate, f
MF = add_coef_matr(f_full,dmax_tot,0,numVars,powers,codesFull,0);

% Initialize svec of all matrices for the coefficients of each degree
svecVar=cell(dcert+1,1);

% Add the coefficient-wise equality, so f-lambda - certificate = 0
% Handle the schmudgen-like terms, so the "certificate" above
% declare constraints
F = [];

% lambda: objective function variable
sdpvar lambda

rowTemp = zeros(numMonMaxg*numMonG*nchoosek(dcert+numVars,dcert),1);
colTemp = zeros(numMonMaxg*numMonG*nchoosek(dcert+numVars,dcert),1);
valTemp = zeros(numMonMaxg*numMonG*nchoosek(dcert+numVars,dcert),1);
sizeFull = 1;
ind_coef = 1;
% ind_schm = 1;
for dd = 0:dcert
    
    if dd < dcert
        ind_used = floor((dmax_tot-indTotal)/2) == dd;
    else
        ind_used = floor((dmax_tot-indTotal)/2) >= dd;
    end
    
    % indices of shcmudgen terms used with sos with the given degree dd
    numMonTemp = sum(ind_used);
    % codes of the newly added polynomial, taking the svec operator into
    % account
    [~, numMonS, codesS] = add_monomials_and_codes(numVars,powers,[],dd);
    ind_triu = triu(true(numMonS),1);
    [codesS1,codesS2] = meshgrid(codesS,codesS);
    codesS =  codesS1+codesS2;
    codesS = [codesS(ind_triu);codesS(1:numMonS+1:numMonS^2)'];
    
    numWithDiag = numMonS*(numMonS+1)/2;
    numOffDiag = numMonS*(numMonS-1)/2;
    svecTransform = [sqrt(2)*ones(1,numOffDiag),ones(1,numMonS)];
    monVecGTemp = monVecG(ind_used,:);
    
    % the main decision variable that corresponds to SOS, DSOS, SDSOS, or constant term
    svecVar{dd+1} = sdpvar(1,numWithDiag*numMonTemp);
    
    % Constraints on the coefficients
    if numMonS == 1 % if the coefficients are constant terms, by choice of the certificate or degree of the schmudgen term
        
        F=F+[svecVar{dd+1}>=0];
        % add the elements of the equality constraints matrix
        for i = 1:numMonTemp
            % option 4, faster than prod(g.^monVecGTemp(i,:))
            ind_nonzero = find(monVecGTemp(i,:));
            prod_schm = g(ind_nonzero(1)).^monVecGTemp(i,ind_nonzero(1));
            for pp=ind_nonzero(2:end)
                prod_schm = prod_schm*(g(pp).^monVecGTemp(i,pp));
            end
            
            [~, numMonP, codesP] = add_monomials_and_codes(numVars,powers,prod_schm);
            [~,rowTemp(ind_coef:ind_coef+numMonP-1)] = ismember(codesP,codesFull); % find the sequential numbers of rows in the basic "unique" monomial matrix that correspond to the new monomials of the product of polynomials
            colTemp(ind_coef:ind_coef+numMonP-1) = i*ones(1,numMonP); % repeat sequential number of each monomial of the new polynomial for all monomials in the given polynomial (the number of new monomials)
            valTemp(ind_coef:ind_coef+numMonP-1) = prod_schm.coef; % create the vector of values
            ind_coef = ind_coef+numMonP;
            
        end
        sizeFull = sizeFull + numMonTemp;
        
    elseif setType == 1 % dsos
        
        % auxiliary optimization variables
        auxVar{dd+1}=sdpvar(1,numOffDiag*numMonTemp);
        
        % indices to avoid a loop in the dsos constraint
        ind_matr = 1:numOffDiag;
        matr = zeros(numMonS);
        matr(ind_triu) = ind_matr;
        matr = matr + matr';
        matr(matr == 0) = [];
        matr = reshape(matr,numMonS-1,numMonS); % use rows for indices in dsos
        
        for i = 1:numMonTemp
            
            % generate the shcmudgen term
            % option 4, faster than prod(g.^monVecGTemp(i,:))
            ind_nonzero = find(monVecGTemp(i,:));
            prod_schm = g(ind_nonzero(1)).^monVecGTemp(i,ind_nonzero(1));
            for pp=ind_nonzero(2:end)
                prod_schm = prod_schm*(g(pp).^monVecGTemp(i,pp));
            end
            % diag last in each svec, the same as in smat
            F=F+[sqrt(2)*svecVar{dd+1}((i-1)*numWithDiag + numOffDiag + 1:i*numWithDiag) >= sum(auxVar{dd+1}((i-1)*numOffDiag + matr),1)];
            F=F+[auxVar{dd+1}((i-1)*numOffDiag+1:i*numOffDiag) >= svecVar{dd+1}((i-1)*numWithDiag+1:(i-1)*numWithDiag + numOffDiag)];
            F=F+[-auxVar{dd+1}((i-1)*numOffDiag+1:i*numOffDiag) <= svecVar{dd+1}((i-1)*numWithDiag+1:(i-1)*numWithDiag + numOffDiag)];
            
            [~, numMonP, codesP] = add_monomials_and_codes(numVars,powers,prod_schm);
            num_new_coef = numMonP*numWithDiag;
            [rowTemp(ind_coef:ind_coef+num_new_coef-1),colTemp(ind_coef:ind_coef+num_new_coef-1),valTemp(ind_coef:ind_coef+num_new_coef-1)] = ...
                add_product_of_codes(prod_schm,codesP,numMonP,codesS,numWithDiag,codesFull,sizeFull+(i-1)*numWithDiag,sizeFull+i*numWithDiag-1,1,svecTransform);
            ind_coef = ind_coef+num_new_coef;
        end
        sizeFull = sizeFull + numWithDiag*numMonTemp;
        
    elseif setType == 2 % sdsos
        % auxiliary optimization variables and indices
        auxVar{dd+1} = sdpvar(3*numOffDiag*numMonTemp,1); % 1 is kk, 2 is kj, 3 is jj
        ind_socp = 1;
        [indRow, indCol] = find(ind_triu);
        
        for i = 1:numMonTemp
            % generate the shcmudgen term
            ind_nonzero = find(monVecGTemp(i,:));
            prod_schm = g(ind_nonzero(1)).^monVecGTemp(i,ind_nonzero(1));
            for pp=ind_nonzero(2:end)
                prod_schm = prod_schm*(g(pp).^monVecGTemp(i,pp));
            end
            % SDSOS constraints
            svecVar{dd+1}((i-1)*numWithDiag+numOffDiag+1:i*numWithDiag) = 0;
            for k = 1:numOffDiag
                svecVar{dd+1}((i-1)*numWithDiag+numOffDiag+indRow(k)) = svecVar{dd+1}((i-1)*numWithDiag+numOffDiag+indRow(k)) + auxVar{dd+1}(ind_socp);
                svecVar{dd+1}((i-1)*numWithDiag+numOffDiag+indCol(k)) = svecVar{dd+1}((i-1)*numWithDiag+numOffDiag+indCol(k)) + auxVar{dd+1}(ind_socp+2);
                svecVar{dd+1}((i-1)*numWithDiag+k) =  sqrt(2)*auxVar{dd+1}(ind_socp+1);
                F = F + [auxVar{dd+1}(ind_socp) + auxVar{dd+1}(ind_socp+2)>=0];
                F = F +[norm([2*auxVar{dd+1}(ind_socp+1);auxVar{dd+1}(ind_socp)-auxVar{dd+1}(ind_socp+2)])<=auxVar{dd+1}(ind_socp)+auxVar{dd+1}(ind_socp+2)];
                ind_socp = ind_socp + 3;
            end
            
            % create inputs for the coefficients matrix
            [~, numMonP, codesP] = add_monomials_and_codes(numVars,powers,prod_schm);
            num_new_coef = numMonP*numWithDiag;
            [rowTemp(ind_coef:ind_coef+num_new_coef-1),colTemp(ind_coef:ind_coef+num_new_coef-1),valTemp(ind_coef:ind_coef+num_new_coef-1)] = ...
                add_product_of_codes(prod_schm,codesP,numMonP,codesS,numWithDiag,codesFull,sizeFull+(i-1)*numWithDiag,sizeFull+i*numWithDiag-1,1,svecTransform);
            ind_coef = ind_coef+num_new_coef;
            
        end
        sizeFull = sizeFull + numWithDiag*numMonTemp;
        
    elseif setType == 3 % sos
        
        % auxiliary inpit needed for the SDP constraint because of the svec operator
        vecStemp = lambda*ones(numMonS);
        
        for i = 1:numMonTemp
            
            % the SDP constraint for the SOS
            mTemp = smat_real(svecVar{dd+1}((i-1)*numWithDiag+1:i*numWithDiag),numMonS,vecStemp);
            F=F+[mTemp >= 0];
            
            % generate the schmudgen term
            ind_nonzero = find(monVecGTemp(i,:));
            prod_schm = g(ind_nonzero(1)).^monVecGTemp(i,ind_nonzero(1));
            for pp=ind_nonzero(2:end)
                prod_schm = prod_schm*(g(pp).^monVecGTemp(i,pp));
            end
            
            % create inputs for the coefficients matrix
            [~, numMonP, codesP] = add_monomials_and_codes(numVars,powers,prod_schm);
            num_new_coef = numMonP*numWithDiag;
            [rowTemp(ind_coef:ind_coef+num_new_coef-1),colTemp(ind_coef:ind_coef+num_new_coef-1),valTemp(ind_coef:ind_coef+num_new_coef-1)] = ...
                add_product_of_codes(prod_schm,codesP,numMonP,codesS,numWithDiag,codesFull,sizeFull+(i-1)*numWithDiag,sizeFull+i*numWithDiag-1,1,svecTransform);
            ind_coef = ind_coef+num_new_coef;
            
        end
        sizeFull = sizeFull + numWithDiag*numMonTemp;
        
    end
end

% merge all main decision variables
svecVar = transpose([svecVar{:}]);
% create the constaints matrix
rowTemp(ind_coef:end) = [];
colTemp(ind_coef:end) = [];
valTemp(ind_coef:end) = [];
Mfull = sparse(rowTemp,colTemp,valTemp,numMon,sizeFull-1);
% account for equality constraints
if isempty(h)
    Mfull = [Mfull,ML];
    [Q,R]=qr(Mfull);
    MF=Q'*MF;
    F=F+[MF == R*[svecVar;lambda]];
    %     F=F+[MF-ML*lambda == Mfull*svecVar];
else
    coefeq = cell(length(h),1);
    Meq = cell(length(h),1);
    sizeq = zeros(length(h),1);
    for rr=1:length(h)
        [Meq{rr},~,sizeq(rr)] = add_coef_matr(h(rr),dmax_tot,[],numVars,powers,codesFull,0);
        coefeq{rr} = sdpvar(1,sizeq(rr));
    end
    Meq=[Meq{:}];
    coefeqVec = transpose([coefeq{:}]);
    Mfull = [Mfull,Meq,ML];
    [Q,R]=qr(Mfull);
    MF=Q'*MF;
    F=F+[MF == R*[svecVar;coefeqVec;lambda]];
    %    F=F+[MF - ML*lambda + Meq*coefeqVec == Mfull*svecVar];
end

% objective
obj=-lambda;
construct_time = toc;
% end the construction

% options for optimization
opt=sdpsettings;
opt.verbose = 0;
opt.dualize = 0;
if setType >= 3
    opt.solver='mosek';
    %     opt.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1E-7;
    %     opt.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1E-7;
    %     opt.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1E-7;
    %     opt.mosek.MSK_DPAR_PRESOLVE_TOL_X=1E-7;
    opt.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=1800;
else
    opt.solver='gurobi';
    opt.gurobi.NumericFocus = 2;
    opt.gurobi.TimeLimit = 1800;
end

% optimization
Yout = optimize(F,obj,opt)

% results
lower_bound = double(lambda)
construct_time
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem