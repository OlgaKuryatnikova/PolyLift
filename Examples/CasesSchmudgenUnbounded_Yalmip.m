% these UNBOUNDED cases made to run the Schmudgen certificate, in SOSTOOLS

% Modification of Example 3 from the non-archimedean example: multivar,
% unbounded
numVars=3;
c = 0;
x = sdpvar(numVars,1);
f = c - x(1)^2 - x(2)^2  + sum(x(3:end).^3);
g = [x(1)-1/2; x(2)-1/2; 1-x(1)*x(2); x'];
h=[];

% Modification of example 4.5 in Nie et al: multivar, reduced, nonneg
numVars=6;
M = 1;
x = sdpvar(numVars,1);
f = sum(x(1:numVars-1).^2)/(numVars-1) + x(end)^2;
g = [x(end)^2 - 1; sum(x(1:numVars-1).^2)-M*sum(x(1:numVars-1))*x(end)-(numVars-1); x'];
h=[];
