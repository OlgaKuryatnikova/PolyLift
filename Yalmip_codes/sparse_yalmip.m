%% Sparse certificate on a simplex or on a box.
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
setType = 1; % 0 means simplex, 1 means box
n_vars_univ = 5; % if we have this many variables or more, we start using univariate sos in the certificate
%%

%% Start the contruction
tic

degf = degree(f); % degree of the objective
lenG = length(g);
Mhat = M;
Ug = zeros(lenG,1);
Dg = zeros(lenG,1);
varsets = cell(lenG,1);
numMong = zeros(lenG,1);
% compute the necessary bounds and extend the set of polynomials if we work
% on the full box
for i=1:lenG
    Dg(i)=degree(g(i)); % max degree of g(i)
    varsets{i} = find(degree(g(i),x)); % variables of g(i)
    [coefg,mong] = coefficients(g(i));
    numMong(i) = length(coefg);
    % scale with the corresponding multinomial coefficient
    dFull = zeros(numMong(i),numVars+1);
    for j=1:numMong(i)
        dFull(j,:) = [Dg(i)- degree(mong(j)),degree(mong(j),x)];
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
        coefg = coefficients(g(end));
        numMong(end+1) = length(coefg);
    end
end

% check which bounds on x are present
vec_bounds = zeros(numVars,2); % lb (1 column) ub (2 column)
for i=1:lenG
    if numMong(i) <=2 && Dg(i)==1
        [coefg,mong] = coefficients(g(i));
        ind_bound = find(ismember(x,mong(end)));
        if ~isempty(ind_bound)
            if numMong(i) == 1 && coefg(end) > 0 % ax >= 0
                L(ind_bound) = 0;
                vec_bounds(ind_bound,1) = 1;
            elseif numMong(i) == 1 && coefg(end) < 0 % -ax >= 0
                U(ind_bound) = 0;
                vec_bounds(ind_bound,2) = 1;
            elseif numMong(i) == 2 && coefg(end) > 0 % ax - L >= 0
                L(ind_bound) = -coefg(1)/coefg(end);
                vec_bounds(ind_bound,1) = 1;
            else % numMong(j) == 2 && coefg(end) < 0: U - ax >= 0
                U(ind_bound) = -coefg(1)/coefg(end);
                vec_bounds(ind_bound,2) = 1;
            end
        end
    end
end

% add additional terms for x and the large polynomial
if setType == 0 % simplex
    gLast = Mhat-sum(x)-sum(g);
    for i=1:numVars
        if vec_bounds(i,1) == 0 % we have no lower bound on x
            g = [g;x(i)-L(i)];
            Dg = [Dg;1];
            varsets = [varsets;i];
        end
    end
    g = [g;gLast];
    Dg = [Dg;degree(gLast)];
    varsets = [varsets;{1:numVars}];
elseif setType == 1 % box full
    Ug = [Ug;U-L];
    Mhat = M + sum(Ug); % adjust Mhat to account for the upper bounds on g and x
    for i=1:numVars
        if sum(vec_bounds(i,:)) == 0 % if we have no bounds at all
            g = [g;x(i)-L(i);U(i)-x(i)];
            Dg = [Dg;1;1];
            varsets = [varsets;i;i];
        elseif vec_bounds(i,1) == 1 && vec_bounds(i,2) == 0 % we only have a lower bound
            g = [g;U(i)-x(i)];
            Dg = [Dg;1];
            varsets = [varsets;i];
        elseif vec_bounds(i,1) == 0 && vec_bounds(i,2) == 1 % we only have an upper bound
            g = [g;x(i)-L(i)];
            Dg = [Dg;1];
            varsets = [varsets;i];
        end
    end
    % in this case Mhat-sum(x)-sum(g) = constant term, so we skip it
else    % The case to use if we do not add upper bounds on g since those add to the total upper bound
    Ug = [Ug;U-L];
    Mhat = M + sum(Ug); % adjust Mhat to account for the upper bounds on g and x
    for i=1:numVars
        if sum(vec_bounds(i,:)) == 0 % if we have no bounds at all
            g = [g;x(i)-L(i);U(i)-x(i)];
            Dg = [Dg;1;1];
            varsets = [varsets;i;i];
        elseif vec_bounds(i,1) == 1 && vec_bounds(i,2) == 0 % we only have a lower bound
            g = [g;U(i)-x(i)];
            Dg = [Dg;1];
            varsets = [varsets;i];
        elseif vec_bounds(i,1) == 0 && vec_bounds(i,2) == 1 % we only have an upper bound
            g = [g;x(i)-L(i)];
            Dg = [Dg;1];
            varsets = [varsets;i];
        end
    end
    g = [g;Mhat-sum(U-L)-sum(g)];
    varsets = [varsets;{1:numVars}];
    Dg = [Dg;degree(g(end))];
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
