% INPUT: boundmat expresses low/up bounds: is boundary matrix
%       cond1/2 lo/up.bound on interval size
%
% RETURN: k x 2 matrix of pairs i,j denoting k entries
%       satisfying bounds ie. entries to be pert'd
%	and 3rd entry per row indicates default basis
%
function Entr = permat(boundmat,cond1,cond2)
 
N = size(boundmat,1);
intlen = boundmat-boundmat' ;
Entr = [];

if cond1<eps;
  cond1 = eps;
  fprintf('rectified lo.bound := %.1e on intervals for perturbable entries\n',cond1);
end;

nrm = norm(intlen,1)/N ;
fprintf('1norm / matdim = %.1e\n',nrm);
 
for i=2:N for j=i+1:N
  if (intlen(i,j)>cond1) & (intlen(j,i)<cond2)
        Entr = [ Entr; i, j, min(1,nrm) ];
  end;
end; end;
