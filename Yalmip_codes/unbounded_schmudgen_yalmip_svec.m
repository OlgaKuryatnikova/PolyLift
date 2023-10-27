% generalized Schmudgen certificate with coefficients in constants, dsos, sdsos or sos
% the case of UNbounded sets

clear
% Write here your case following the format of the below example
numVars=3;
x = sdpvar(numVars,1);
f = sum(x(1:numVars-1).^2)/(numVars-1) + x(end)^2; % objective function
g = [x(end)^2 - 1; sum(x(1:numVars-1).^2)-sum(x(1:numVars-1))*x(end)-(numVars-1)]; % vector of inequality constraints, >=0 format, WITHOUT x>=0
h=[]; % vector of equality constraints
setType = 2; % type of coefficients: constant term (0), dsos (1), sos (2)
r = 0; % level of the hierarchy, r from the paper
%

% % Start the contruction
tic

% additional degree to control for numerica issues if needed
deg_aux = 0;

lenG = length(g);
Dg = zeros(lenG,1);
for j=1:lenG
    Dg(j) = degree(g(j));
end
Df = degree(f);
dpoly = max(max(Dg),ceil(Df/2));
dmax = max(degree(h),r+2*dpoly+deg_aux);

% denominator D for the unbounded certificate
deg_D = floor(r/dpoly);
beta =  indexCreationLess_sym(lenG,deg_D);
D = 0;
denom = 1+sum(x);
for i=1:size(beta,1)
    deg_temp = beta(i,:);
    D = D + multinomial(r,[r-sum(deg_temp*dpoly),deg_temp*dpoly])*denom^(r-deg_temp*Dg)*prod(g.^(deg_temp'));
end

% lambda-variable
sdpvar lambda

% adjusted objective function in the unbounded certificate
f = D*(1+sum(x))^max(0,2*dpoly-Df)*(f - lambda);

% Add the coefficient-wise equality, so Denominator(f-lambda) - certificate = 0, start with Denominator(f-lambda)
sum_poly = f;

% Add variables x to the set description to use for the certificate
g = [g;x];
Dg = [Dg;ones(numVars,1)];
lenG = length(g);

% Generate possible vectors of degrees for Schmudgen terms
monVecG = indexCreationLess_sym(lenG,dmax);
indTotal = monVecG*Dg;
% option 1
indBelowDmax = indTotal <= dmax; % this option provides the same vector of degrees as in the paper (option 2), but it is simpler
% % option 2
% indTotalG = sum(monVecG(:,1:lenGinit),2);
% indTotalx = sum(monVecG(:,lenGinit+1:end),2);
% indBelowDmax = indTotalx + dpoly*indTotalG <= dmax;
monVecG = monVecG(indBelowDmax,:);
indTotal = indTotal(indBelowDmax);
numMonG = sum(indBelowDmax);

if numMonG == 1
    fprintf('There is only one monomial for the certificate, something is wrong, execution stopped. \n')
    return
end

% Choose the type of coefficients: constant term (0), dsos (1), sos (2)
if setType == 0
    dcert = 0;
else
    dcert = 1;
    auxVar = cell(dcert+1,1); % auxiliary variables needed for a certificate
end

% Initialize svec of all matrices for the coefficients of each degree
svecVar=cell(dcert+1,1);

% Handle the schmudgen-like terms, so the "certificate" above
F = [];
for dd = 0:dcert
    if dd < dcert
        ind_used = floor((dmax-indTotal)/2) == dd;
    else
        ind_used = floor((dmax-indTotal)/2) >= dd;
    end
    % indices of shcmudgen terms used with sos with the given degree dd
    numMonTemp = sum(ind_used);
    monS = monolist(x,dd);
    tmonS = transpose(monS);
    numMonS = length(monS);
    numWithDiag = numMonS*(numMonS+1)/2;
    numOffDiag = numMonS*(numMonS-1)/2;
    
    monVecGTemp = monVecG(ind_used,:);
    [ind_row,ind_col, deg] = find(monVecGTemp);
    if sum(sum(monVecGTemp)) == 0 || ind_row(1) == 2 % this does not slow things down because we have to do it only twice for dcert = 1
        ind_row = [1;ind_row];
        ind_col = [1;ind_col];
        deg = [0;deg];
    end
    
    % Constraints on the coefficients
    if numMonS == 1 % if the coefficients are constant terms, by choice of the certificate or degree of the schmudgen term
        svecVar{dd+1} = sdpvar(1,numWithDiag*numMonTemp);
        prod_schm = x(1)*ones(numMonTemp,1);
        for i = 1:numMonTemp
            ind_i = ind_row == i;
            ind_nonzero = ind_col(ind_i);
            deg_nonzero = deg(ind_i);
            prod_schm(i) = prod(g(ind_nonzero).^deg_nonzero);
        end
        sum_poly = sum_poly - svecVar{dd+1}*prod_schm;
        % Non-negativity constraint
        F=F+[svecVar{dd+1}>=0];
        
    elseif setType == 1 % dsos
        ind_triu = triu(true(numMonS),1);
        ind_matr = 1:numOffDiag;
        matr = zeros(numMonS);
        matr(ind_triu) = ind_matr;
        matr = matr + matr';
        matr(matr == 0) = [];
        matr = reshape(matr,numMonS-1,numMonS)'; % use rows for indices in dsos
        
        svecVar{dd+1} = sdpvar(numWithDiag*numMonTemp,1);
        vecStemp = x(1)*ones(numMonS);
        auxVar{dd+1}=sdpvar(numOffDiag*numMonTemp,1);
        
        for i = 1:numMonTemp
            % option 3
            ind_i = ind_row == i;
            ind_nonzero = ind_col(ind_i);
            deg_nonzero = deg(ind_i);
            prod_schm = prod(g(ind_nonzero).^deg_nonzero);
            sum_poly = sum_poly - tmonS*smat_real(svecVar{dd+1}((i-1)*numWithDiag+1:i*numWithDiag),numMonS,vecStemp)*monS*prod_schm;
            % DSOS constraints
            F=F+[sqrt(2)*svecVar{dd+1}((i-1)*numWithDiag + numOffDiag + 1:i*numWithDiag) >= sum(auxVar{dd+1}((i-1)*numOffDiag + matr),2)];
            F=F+[auxVar{dd+1}((i-1)*numOffDiag+1:i*numOffDiag) >= svecVar{dd+1}((i-1)*numWithDiag+1:(i-1)*numWithDiag + numOffDiag)];
            F=F+[-auxVar{dd+1}((i-1)*numOffDiag+1:i*numOffDiag) <= svecVar{dd+1}((i-1)*numWithDiag+1:(i-1)*numWithDiag + numOffDiag)];
        end
        
    elseif setType == 2 % sos
        svecVar{dd+1} = sdpvar(numWithDiag*numMonTemp,1);
        vecStemp = x(1)*ones(numMonS);
        mTemp = x(1)*ones(numMonS);
        for i = 1:numMonTemp
            mTemp = smat_real(svecVar{dd+1}((i-1)*numWithDiag+1:i*numWithDiag),numMonS,vecStemp);
            ind_i = ind_row == i;
            ind_nonzero = ind_col(ind_i);
            deg_nonzero = deg(ind_i);
            prod_schm = prod(g(ind_nonzero).^deg_nonzero);
            sum_poly = sum_poly - tmonS*mTemp*monS*prod_schm;
            % SOS constraint
            F=F+[mTemp >= 0];
        end

    end
end

% equality constraints
if ~isempty(h)
    numeq = length(h);
    Seq = cell(numeq,1);
    sizeq = zeros(numeq,1);
    for rr=1:numeq
        moneq = monolist(x,dmax - degree(h(rr)));
        sizeq(rr) = length(moneq);
        Seq{rr} = sdpvar(1,sizeq(rr));
        sum_poly = sum_poly - Seq{rr}*moneq*h(rr);
    end
end

% total final constraint
F = F + [coefficients(sum_poly,x) == 0];

% objective
obj=-lambda;

construct_time = toc;
% end the construction

% % Options and optimization
opt=sdpsettings;
opt.verbose = 0;
opt.dualize = 0;
if setType >= 2
    opt.solver='mosek';
    opt.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1E-7;
    opt.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1E-7;
    opt.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1E-7;
    opt.mosek.MSK_DPAR_PRESOLVE_TOL_X=1E-7;
    opt.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=1800;
else
    opt.solver='gurobi';
    opt.gurobi.NumericFocus = 2;
    opt.gurobi.TimeLimit = 1800;
end

Yout = optimize(F,obj,opt);

lower_bound = double(lambda)
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem