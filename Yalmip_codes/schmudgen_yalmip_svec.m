%% Schmudgen-based certificate with coefficients in constants, dsos, sdsos or sos
%% The case of compact sets.
%% Use Yalmip to encode polynomials.
clear

%% Write here your case following the format of the below example
numVars=5; % the number of variables in the problem
x = sdpvar(numVars,1);
f = 7*(2*x(1)-x(2)+x(3)-2*x(4)-2*x(5)); % the objective function
g = [(7*x(1)-2)^2-49*x(2)^2-(7*x(3)-1)^2-(7*x(5)-1)^2; 49*x(1)*x(3)-49*x(4)*x(5)+49*x(1)^2-1;7*x(3)-49*x(2)^2-49*x(4)^2-1;...
    49*x(1)*x(5)-49*x(2)*x(3)-2;2-sum(x);x]; % vector of inequality constraints, >=0 format
M = 2;
L = zeros(numVars,1);
U = M*ones(numVars,1);
h = []; % vector of equality constraints, =0 format
dmax = 4; % maximal degreee of the hierarchy
setType = 2; % type of coefficients: constant term (0), dsos (1), sos (2)
%

%% Start the contruction
tic

lenG = length(g);
Dg = zeros(lenG,1);
vec_bounds = zeros(numVars,2); % lb (1 column) ub (2 column)
for j=1:lenG
    Dg(j) = degree(g(j));
    [coefg,mong] = coefficients(g(j));
    numMonGj = length(coefg);
    if numMonGj <=2 && Dg(j)==1
        ind_bound = find(ismember(x,mong(end)));
        if ~isempty(ind_bound)
            if numMonGj == 1 && coefg(end) > 0 % ax >= 0
                L(ind_bound) = 0;
                vec_bounds(ind_bound,1) = 1;
            elseif numMonGj == 1 && coefg(end) < 0 % -ax >= 0
                U(ind_bound) = 0;
                vec_bounds(ind_bound,2) = 1;
            elseif numMonGj == 2 && coefg(end) > 0 % ax - L >= 0
                L(ind_bound) = -coefg(1)/coefg(end);
                vec_bounds(ind_bound,1) = 1;
            else % numMong(j) == 2 && coefg(end) < 0: U - ax >= 0
                U(ind_bound) = -coefg(1)/coefg(end);
                vec_bounds(ind_bound,2) = 1;
            end
        end
    end
end
Df = degree(f);


% add bounds on x to the set description if they were not there initially, to use for the certificate
for i=1:numVars
    if sum(vec_bounds(i,:)) == 0 % if we have no bounds at all
        g = [g;x(i)-L(i);U(i)-x(i)];
        Dg = [Dg;1;1];
    elseif vec_bounds(i,1) == 1 && vec_bounds(i,2) == 0 % we only have a lower bound
        g = [g;U(i)-x(i)];
        Dg = [Dg;1];
    elseif vec_bounds(i,1) == 0 && vec_bounds(i,2) == 1 % we only have an upper bound
        g = [g;x(i)-L(i)];
        Dg = [Dg;1];
    end
end
lenG = length(g);

% lambda-variable
sdpvar lambda

% Add the coefficient-wise equality, so Denominator(f-lambda) - certificate = 0, start with f-lambda
sum_poly = f - lambda;

% Generate possible vectors of degrees for Schmudgen terms
monVecG = indexCreationLess(lenG,dmax);
indTotal = monVecG*Dg;
% option 1
indBelowDmax = indTotal <= dmax;
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

Yout = optimize(F,obj,opt)

lower_bound = double(lambda)
construct_time
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem