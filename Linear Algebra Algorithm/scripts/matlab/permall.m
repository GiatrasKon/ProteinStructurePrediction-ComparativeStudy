% Input: dim of bordered matrix
%
% Return: kx3 matrix: each row has x,y expressing entry to be pert'd
%       These are all non-border entries above diagonal: x<y
%	3rd entry is basis element, default=1: too crude: must improve!
%
function pertEntr = permall (dim)
basis_default = 1;
pertEntr = [];
for i=2:dim for j=i+1:dim
        pertEntr = [pertEntr; i,j, basis_default];
end; end;
