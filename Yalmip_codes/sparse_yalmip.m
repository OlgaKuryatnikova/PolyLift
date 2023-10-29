% Sparse certificate on a simplex or on a box
% the paper presents the version on a simplex only
clear

%% Write here your case following the format of the below example
numVars=3;
x = sdpvar(numVars,1);
f = sum(x(1:numVars-1).^2)/(numVars-1) + x(end)^2; % objective function
g = [x(end)^2 - 1; sum(x(1:numVars-1).^2)-sum(x(1:numVars-1))*x(end)-(numVars-1)]; % vector of inequality constraints, >=0 format, WITHOUT upper and lower bounds
L = zeros(numVars,1); % lower bounds
U = (numVars-1)*ones(numVars,1); % upper bounds
M = sum(U); % simplex upper bound from the certificate in the paper
h=[]; % vector of equality constraints
dmax = 8; % maximal degreee of the hierarchy
setType = 0; % 0 means simplex, just what is presented in the paper, 1 means box
n_vars_univ = 5; % starting from how many variables we start using univariate sos in the certificate
%%

%% Construct the certificate
tic

degf = degree(f); % degree of the objective
lenG = length(g);
Mhat = M;
Ug = zeros(lenG,1);
Dg = zeros(lenG,1);
varsets = cell(lenG,1);
vec_bounds = zeros(numVars,1); % ub (-1) or lb (1)

for i=1:lenG
    Dg(i)=degree(g(i)); % max degree of g(i)
    varsets{i} = find(degree(g(i),x)); % variables of g(i)
    [coefg,mong] = coefficients(g(i));
    numMonGi = length(coefg);
    
    dFull=[];
    for j=1:length(coefg)
        dFull = [dFull;[Dg(i)- degree(mong(j)),degree(mong(j),x)]]; 
    end
    coefMultinom = abs(coefg)./(factorial(Dg(i))./prod(factorial(dFull),2));
    maxcoef=max(coefMultinom);
    Ug(i) = (1+M+sum(abs(L)-L))^Dg(i)*maxcoef;
    Mhat=Mhat+Ug(i);
    % add redundant constraints if we work with a box
    if setType == 1 % box full
        g(end+1)=Ug(i)-g(i);
        varsets{end+1} = varsets{i};
        Dg(end+1) = Dg(i);
        Ug(end+1) = Ug(i);
    end
end

% add additional terms for x and the large polynomial
if setType == 0 % simplex
    gLast = Mhat-sum(x)-sum(g);
    g = [g;x-L;gLast];
    varsets = [varsets;num2cell(1:numVars)';{1:numVars}];
    Dg = [Dg;ones(numVars,1);degree(gLast)];
elseif setType == 1 % box full
    Ug = [Ug;U-L];
    Dg = [Dg;ones(numVars,1)];
    g = [g;U-x'];
    Mhat = M + sum(Ug); % adjust Mhat to account for the upper bounds on g and x
    g = [g;x-L];
    % in this case Mhat-sum(x)-sum(g) = constant term, so we skip it
    varsets = [varsets;num2cell(1:numVars)';num2cell(1:numVars)'];
    Dg = [Dg;ones(numVars,1)];
else    % The case to use if we do not add upper bounds on g since those add to the total upper bound
    Ug = [Ug;U-L];
    Dg = [Dg;ones(numVars,1)];
    g = [g;U-x];
    Mhat = M + sum(Ug); % adjust Mhat to account for the upper bounds on g and x
    g = [g;x-L;Mhat-sum(U-L)-sum(g)];
    varsets = [varsets;num2cell(1:numVars)';num2cell(1:numVars)';1:numVars];
    Dg = [Dg;ones(numVars,1);max(sum(g(end).degmat,2))];
end
lenG=length(g);

% add the full-variables polynomial
if dmax>=2*max(max(Dg,degf))
    gbig=Mhat^2-sum(g.^2);
else
    gbig=[];
end
gbig=[gbig;1];
lenGbig=length(gbig);
%         save(['nnrsp' num2str(numVars) 'v6d' num2str(numg) 'con' num2str(degf(jj)) 'obj' num2str(ii) 'exp.mat'])

%% Start adding polynomial coefficients

% lambda-variable
sdpvar lambda

% Add the coefficient-wise equality, so Denominator(f-lambda) - certificate = 0, start with f-lambda
sum_poly = f - lambda;
sosBig = cell(lenGbig,1);
F = [];
% add large SOS
for i=1:lenGbig
    dd = floor((dmax - degree(g(i)))/2);
    monS = monolist(x,dd);
    tmonS = transpose(monS);
    numMonS = length(monS);
    sosBig{i} = sdpvar(numMonS);
    vecStemp = x(1)*ones(numMonS);
    sum_poly = sum_poly - tmonS*sosBig{i}*monS*g(i);
    F = F+[sosBig{i} >= 0];
end

% add small SOS
dsmall=floor((dmax-Dg)./(2*Dg));
MG=cell(lenG,1);
sizeSDPsmall=cell(lenG,1);
SOS=cell(lenG,1);

for i=1:lenG
    if length(varsets{i}) >= n_vars_univ % use univariate sos if we have more than n_vars_univ vars
        degs=(0:dsmall(i))';
        monS=g(i).^degs;
        tmonS = transpose(monS);
    else
        dd = floor((dmax - degree(g(i)))/2);
        monS = monolist(x(varsets{i}),dd);
        tmonS = transpose(monS);
    end
    numMonS = length(monS);
    SOS{i} = sdpvar(numMonS);
    F=F+[SOS{i}>=0];
    sum_poly = sum_poly - tmonS*SOS{i}*monS*g(i);
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
construct_time=toc;
% end the construction

%% Options and optimization
opt=sdpsettings;
opt.solver='mosek';
opt.verbose=0;
opt.dualize=1;
opt.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=1800;
Yout = optimize(F,obj,opt)

% Output the results
lower_bound = double(lambda)
construct_time
sol_time = Yout.solvertime
yalmip_time = Yout.yalmiptime
status = Yout.problem
