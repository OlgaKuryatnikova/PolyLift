function Y = smat_real(x,sz,Y)
% Y = smat_real(x,sz,Y)
% Placing the elements of the vector x in the symmetric real-valued matrix Y of the size sz x sz 
% using the adjoint of the standard svec operator. 
% If Y is not given, it is set to be all-zero sparse.

ind=triu(true(sz),1); % faster than "find"

if nargin==2
    Y=sparse(sz,sz);
else
    Y(~ind)=0;
end

Y(ind)=x(1:end-sz)/sqrt(2);
Y=Y+transpose(Y);
Y(1:sz+1:sz^2)=x(end-sz+1:end);

end

