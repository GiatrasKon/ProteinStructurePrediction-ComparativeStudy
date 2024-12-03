% used to be called distmat()
% Input: matrices D,D1
%	kx2 matrix ro contains k entries
%	scalar step indicates amount of perturbation
% Return: symmetric matrix D w/entry ro[1,:] perturbed by step
%       the other entries ro[] copied from D1
%
function D = pert1mat(D,D1,ro,step)
 
  N = size(ro,1);
  x = ro(1,1); y = ro(1,2);
  D(x,y) = D(x,y) + step;
  D(y,x) = D(x,y);
 
  for i=2:N
    x = ro(i,1); y = ro(i,2);
    D(x,y) = D1(x,y);
    D(y,x) = D(x,y);
  end
%
% end pert1mat
