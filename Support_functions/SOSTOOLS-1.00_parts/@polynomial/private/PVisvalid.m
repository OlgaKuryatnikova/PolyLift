function [flag,errormsg] = PVisvalid(astruct)
% function [Flag,Errormsg] = PVisvalid(Astruct)
%
% DESCRIPTION 
%   (INTERNAL FUNCTION)  
%   Check if Astruct is a valid structure for a polynomial object.
%   
% INPUTS 
%   Astruct: structure from a polynomial object. 
%       The validity of this structure will be checked.
%
% OUTPUTS  
%   Flag: Flag=1 if structure is valid and 0 otherwise.
%   ErrorMsg: if Flag=0, this is a message indicating the 
%      location and type of error.
%  
% SYNTAX 
%   [Flag,Errormsg] = PVisvalid(Astruct);

% 6/9/2002: PJS  Initial Coding

% Initialize Outputs
flag = 1;
errormsg = [];
  
% Check for errors in the structure
if ~isstruct(astruct)
  flag = 0;
  errormsg = 'Input must be a structure';
  return;
end
afields = fieldnames(astruct);
polyfields = {'coefficient'; 'degmat'; 'varname'; 'matdim'};
if ~isempty( setxor(afields,polyfields) )
  flag = 0;
  errormsg = ['Structure fields must be: coefficient, degmat,' ...
	      ' varname, and matdim.'];
  return;
end  
  
% Grab fields
adim = astruct.matdim;
acoef = astruct.coefficient;
adeg = astruct.degmat;
avar = astruct.varname;

% Empty polynomials are valid
if isempty(astruct.coefficient) & isempty(astruct.degmat) & ...
	isempty(astruct.varname) & all(adim==[0 0])
    flag = 1;
    errormsg = [];
    return;
end

% Check matdim: Should be 1x2 vector of integers
if ~isreal(adim) | ndims(adim)~=2 | floor(adim)~=ceil(adim) | ...
       size(adim)~=[1 2] | any(adim <0)
  flag = 0;
  errormsg=['Matdim should be a 1x2 matrix of nautral numbers'];
  return;        
end

% Check coefficient: Should be a real array.
if isempty(acoef)
  flag = 0;
  errormsg=['Coefficient cannot be empty if Degmat or Varname' ... 
	    ' are nonempty.'];
  return;      
else
  if ~isreal(acoef) | ndims(acoef)~=2  | min(size(acoef))~=1
    flag = 0;
    errormsg = ['Coefficient must be a real vector.'];
    return;
  end
end

% Check degmat:  Should be a matrix of natural numbers.
% (Note that acoef is nonempty if we reached this point)
if isempty(adeg)
  flag = 0;
  errormsg=['Degmat cannot be empty if Coefficient is nomempty.'];
  return;
else
  if isreal(adeg) & ndims(adeg)==2    
    adegnz = nonzeros(adeg);
    tempflg = isempty(adegnz) | (adegnz>=0 & floor(adegnz)==ceil(adegnz));
    if ~tempflg
      flag = 0;
      errormsg=['Degmat should be a matrix of nautral numbers'];
      return;
    end	      
  else
    flag = 0;
    errormsg=['Degmat should be a matrix of nautral numbers'];
    return;
  end	    
end

% Check length of Coefficient
if length(acoef)~=adim(1)*adim(2)*size(adeg,1) 
  flag = 0;
  errormsg=['Length of Coefficient must be N*M*T.'];
  return;     
end

% Number of variables
nv = size(adeg,2);
constantpoly = (nv==1) & isempty(adegnz);

% Check varname: Should be column cell arrays of strings with
% length == nv.
% (Note that acoef and adeg are nonempty if we reached this point)
if isempty(avar) 
  
  if ~constantpoly
    flag = 0;
    errormsg=['Varname cannot be empty if polynomial is not constant.'];
    return;
  end
  
else
  
  tempflg = ndims(avar)==2 & iscellstr(avar) & min(size(avar)) == 1;
  if ~tempflg
    flag = 0;
    errormsg=['Varname should be a column cell array of strings'];
    return;
  end	

  if length(avar) ~= nv;
    flag = 0;
    errormsg=['Number of variables implied by Degmat and Varname are ' ... 
	      'inconsistent.'];
    return;     
  end
      	    
end
      




