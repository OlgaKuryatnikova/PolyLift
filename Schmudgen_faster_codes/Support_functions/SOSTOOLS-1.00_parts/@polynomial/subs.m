function b = subs(a,old,new)
% function B = subs(A,Old,New);
%
% DESCRIPTION 
%  Symbolic Substitution.
%  This function differs from the command "subs" in the symbolic toolbox in
%  three ways.  First, the output of this function is always a polynomial.
%  Use the command "double" to depromote a constant polynomial to a
%  double. Second, the input "Old" cannot be an arbitrary polynomial
%  expression.  Third, the order of "Old" and "New" cannot be switched.
%   
% INPUTS 
%   A: polynomial
%   Old: polynomial variable or a cell array of polynomial variables
%   New: polynomial or double or a cell array of polynomials and/or 
%          doubles.
%   (Either A and/or New must be 1x1)
%
% OUTPUTS  
%   B:  polynomial
%    
% SYNTAX 
%   B = subs(A);
%       -Replaces all variables in A with values in the workspace.
%   B = subs(A,Old,New);
%       -Replaces Old with New in the symbolic expression A.
%       -If Old and New are cell arrays of the same size, each element 
%       of Old is replaced by the corresponding element of New. 
%       -If A and Old are scalars and New is an array or cell array, 
%       the scalars are expanded to produce an array result.
%   B = subs(A,New);
%       -Equivalent to B=subs(A,A.varname{1},New);
  
% 3/25/2003: PJS  Initial Coding  
  
% Promote a to polynomial 
a = polynomial(a);

% Convert old to a cell array of strings and new to a cell 
% array of doubles and/or polynomials
if nargin == 1
  
  new = cell(0);
  nva = length(a.varname);
  cnt = 1;
  for i1 = 1:nva
    
    if evalin('caller', ['isa(' a.varname{i1} ',''double'')'] );
      old{cnt} = a.varname{i1};
      new{cnt} = evalin('caller',a.varname{i1});
      cnt = cnt +1;
    end
    
    if evalin('caller', ['isa(' a.varname{i1} ',''polynomial'')'] );
      temp = evalin('caller',a.varname{i1});
      if temp~=polynomial(1,1,a.varname(i1),[1 1])
 	old{cnt} = a.varname{i1};
 	new{cnt} = temp;
 	cnt = cnt+1;
      end
    end
    
  end

  if length(new)==0
    b = a;
    return;
  end
  
elseif nargin == 2
  
  b=subs(a,a.varname{1},old);
  return;
        
elseif nargin == 3

  if ~iscell(old)
    old = {old};
  end  
  szo = size(old);
  for i1 = 1:szo(1)
    for i2 = 1:szo(2)
      if isa(old{i1,i2},'polynomial');
	old{i1,i2} = old{i1,i2}.varname{1};
      elseif ~ischar(old{i1,i2})
	error(['Old should be a polynomial variable or a '  ...
	       'cell array of polynomial variables']);
      end
    end
  end    

  if ~iscell(new)
    new = {new};
  end
  szn = size(new);
  for i1 = 1:szn(1)
    for i2 = 1:szn(2)
      if ~isa(new{i1,i2},'polynomial') & ...
	    ~isa(new{i1,i2},'double') 
	error(['New should be a polynomial or double or a ' ... 
	      'cell array of polynomials and/or doubles.']);
      end
    end
  end    
  
end

% Some quick error checking
sza = size(a);
szn = size(new{1});
if any(cellfun('size',new,1)~=szn(1)) | ...
      any(cellfun('size',new,2)~=szn(2)) 
  error('All entries of new must have the same dimensions');
end
if any(szn~=1) & any(sza~=1)
  error('Either A or New must be 1x1');
end

% Make New a single array if A and Old are 1x1
b=a;
szb=size(b);
szo = size(old);
if all(szo==[1 1]) & all(szb==[1 1]) & ...
      all(cellfun('isclass',new,'double'))

  new = {cell2mat(new)};
  szn = size(new{1});
  
elseif all(szo==[1 1]) & all(szb==[1 1])

  tempnew = polynomial;
  for i1 = 1:size(new,1);
    temppoly = polynomial;
    for i2 = 1:size(new,2);
      temppoly = [temppoly new{i1,i2}];    
    end
    tempnew = [tempnew; temppoly];
  end  
  new = {tempnew};
  szn = size(new{1});
  
end


old = old(:);
new = new(:);
if all(cellfun('isclass',new,'double')) % Numeric Substitution
    
  ntb = size(b.degmat,1);
  if all(szb==[1 1])
    b.matdim = szn;
    tempsz = szn(1)*szn(2);
    bcoef = sprepmat(b.coefficient,1,tempsz);    
  else
    tempsz = szb(1)*szb(2);
    bcoef = reshape(b.coefficient,ntb,tempsz);
  end

  for i1 = 1:length(old)  
    idx = strmatch(old{i1},b.varname); 
    if isempty(idx)
      error(['Undefined variable ' old{i1}]);
    end
    tempdeg = b.degmat(:,idx);
    b.degmat(:,idx) = [];
    b.varname(idx)=[];  
    
    tempdeg = sprepmat(tempdeg,1,tempsz);
    if all(szb==[1 1])
      tempnew = sprepmat(new{i1}(:)',ntb,1);
    else      
      tempnew = sprepmat(new{i1}(:)',ntb,tempsz);
    end
    bcoef = bcoef.*(tempnew.^tempdeg); 
    b.coefficient=bcoef;
  end

  
else % Symbolic Substitution  
    
  ntb = size(b.degmat,1);
  nvb = size(b.degmat,2);
  [nrb,ncb]=size(b);
  bcoef = reshape(b.coefficient,ntb,nrb*ncb);
  bdeg = b.degmat; 
  bvar = b.varname;
  bdim = b.matdim;   

  for i1 = 1:length(old);
    temp = strmatch(old{i1},b.varname); 
    if isempty(temp)
      error(['Undefined variable ' old{i1}]);
    else
      idx(i1) = temp;
    end
  end
  idx2=setdiff(1:nvb,idx);
  
  b=polynomial(0);
  for i1 = 1:ntb        
    coef = bcoef(i1,:);
    if isempty(idx2)
      temppoly = polynomial(coef);
    else      
      coef = coef(:);
      temppoly = polynomial(coef,bdeg(i1,idx2),bvar(idx2),bdim);
    end
    
    for i2 = 1:length(old)      
      temppoly = temppoly*power(new{i2},bdeg(i1,idx(i2)));
    end
    b=b+temppoly;
  end
    
end

% Combine down
b = combine(b);










