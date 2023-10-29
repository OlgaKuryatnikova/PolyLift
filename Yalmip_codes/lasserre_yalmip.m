% The Lasserre hierarchy implemented in Yalmip
clear

%% Write here your case following the format of the below example
numVars=5;
x = sdpvar(numVars,1);
f = sum(x(1:numVars-1).^2)/(numVars-1) + x(end)^2; % objective function
g = [x(end)^2 - 1; sum(x(1:numVars-1).^2)-sum(x(1:numVars-1))*x(end)-(numVars-1)]; % vector of inequality constraints, >=0 format
h=[]; % vector of equality constraints
dmax = 4; % the maximum degree of the relaxation
%%

% Begin the certificate
g = [1;g];
lenG = length(g);

% % Start the contruction
tic

% lambda-variable
sdpvar lambda

% Add the coefficient-wise equality, so f-lambda - certificate = 0
sum_poly = f - lambda;

% Handle the schmudgen-like terms, so the "certificate" above
F = [];
svecVar = cell(lenG,1);

for i=1:lenG
    dd = floor((dmax - degree(g(i)))/2);
    monS = monolist(x,dd);
    tmonS = transpose(monS);
    numMonS = length(monS);
    svecVar{i} = sdpvar(numMonS);
    sum_poly = sum_poly - tmonS*svecVar{i}*monS*g(i);
    F = F+[svecVar{i} >= 0];
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

% End the certificate

%% Options and optimization
opt=sdpsettings;
opt.verbose = 0;
opt.dualize = 0;

opt.solver='mosek';
opt.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1E-7;
opt.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1E-7;
opt.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1E-7;
opt.mosek.MSK_DPAR_PRESOLVE_TOL_X=1E-7;
opt.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=1800;

% optimize
Yout = optimize(F,obj,opt)%

% outputs
lower_bound = double(lambda)
construct_time
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem
