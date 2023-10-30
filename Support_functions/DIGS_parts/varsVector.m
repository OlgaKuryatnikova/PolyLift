% APPS Solver
% Version 1.00, June 26, 2011
%
% Written By: Bissan Ghaddar and Juan Vera

function x = varsVector(name,n)
%returns vector x of length n such that x(i) is variables 'name'0i
%Usage: x = varsVector(name,n) 
%create vars in polynomial using newpvar and append zeros at the beginning
%Input name of the variable and the size n
%Output a vector x

N = (1:n)';
Nasc = int2str(N);
Iblanks = find(Nasc ==' ');                                                 
Nasc(Iblanks) = '0';

St =[kron(['''' name],ones(n,1)), Nasc, kron(['''' ','],ones(n,1))];
St =reshape(St',1,[]);

St(end) = [];

eval(['x = newpvar(',St,');']); 
 
