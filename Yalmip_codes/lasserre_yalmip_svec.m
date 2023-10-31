%% The Lasserre hierarchy with svec operator instead of the full matrix.
% The svec operator can be useful when the matrices become very large.
%% Use Yalmip to encode polynomials.
clear

% Write here your case following the format of the below example
numVars=5; % the number of variables in the problem
x = sdpvar(numVars,1);
f = sum(x(1:numVars-1).^2)/(numVars-1) + x(end)^2; % objective function
g = [x(end)^2 - 1; sum(x(1:numVars-1).^2)-sum(x(1:numVars-1))*x(end)-(numVars-1);x]; % vector of inequality constraints, >=0 format
h=[]; % vector of equality constraints
dmax = 2; % the maximum degree of the relaxation
%%

%% Start the contruction
tic

% add the constant term to the non-negativity constraints
g = [1;g];
lenG = length(g);

% lambda-variable
sdpvar lambda

% Add the coefficient-wise equality, so f-lambda - certificate = 0
sum_poly = f - lambda;

% Start adding the Putinar terms
F = [];
svecVar = cell(lenG,1);
for i=1:lenG
    dd = floor((dmax - degree(g(i)))/2);
    monS = monolist(x,dd);
    tmonS = transpose(monS);
    numMonS = length(monS);
    numWithDiag = numMonS*(numMonS+1)/2;
    numOffDiag = numMonS*(numMonS-1)/2;
    svecVar{i} = sdpvar(numWithDiag,1);
    vecStemp = x(1)*ones(numMonS);
    if numMonS>1
        mTemp = smat_real(svecVar{i},numMonS,vecStemp);
        sum_poly = sum_poly - tmonS*mTemp*monS*g(i);
        F = F+[mTemp >= 0];
    else
        sum_poly = sum_poly - svecVar{i}*g(i);
        F = F+[svecVar{i} >= 0];
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

F = F + [coefficients(sum_poly,x) == 0];

% objective
obj=-lambda;
construct_time = toc;
% end the construction

% % Options and optimization
opt=sdpsettings;
opt.verbose = 0;
opt.dualize = 0;

opt.solver='mosek';
% opt.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1E-7;
% opt.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1E-7;
% opt.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1E-7;
% opt.mosek.MSK_DPAR_PRESOLVE_TOL_X=1E-7;
opt.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=1800;

Yout = optimize(F,obj,opt)%

lower_bound = double(lambda)
construct_time
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem
