% Input: matrix inpMat
%	 pert bounds the perturbation
%	 Entr entries perturbed by "pert" * rand(0..1) * rand.sgn
%		but should still be non-negative
%
% Return: newMat (symmetric if inpMat symmetric)
% 
function newMat = pertMmat(inpMat,Entr,pert)
 
  if ~isequal(inpMat-inpMat',zeros(size(inpMat)))
        error('input must be symmetric matrix');
        break;
  end;

  newMat = inpMat;
  NumEnt = size(Entr,1);
  for i = 1:NumEnt
    x = Entr(i,1); y = Entr(i,2);
    if rand(1)> .5 sgn=1; else sgn=-1; end;
    newMat(x,y) = newMat(x,y) + pert*rand(1)*sgn*newMat(x,y);
    if newMat(x,y)<0; newMat(x,y)=0; end;
    newMat(y,x) = newMat(x,y);
  end;

  fprintf('Mat.entries perturbed by +-%.2f*(0..1)%%.\n',pert);
