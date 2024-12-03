% I.Emiris: Created 6/01,2002. Last modified: 7/03.
% based on increment/molconf7exp.m by nikitop@csd.uoc.gr
%
% Examples:
%   A = molstruc7exp; svds(A,6); mconf(A,1e-8,0) % approx.symm A
%   B=bounds7; mconf(B,1e-8,1) % up/low bounds B
%   load 'molstruc12.mat' C; svds(C,6); mconf(C,1e-8,1) % symm C
%
% Input: approx.dist.mat A XOR symmetric matrix D XOR bounds B
%	toler = zero threshold for 6th sing.value
%	flag = scalar >= 0
%
% If input is symmetric matrix:
% flag=0: input is symmetric approx.dist.mat.A
% flag>0: perturb away from given dist.matrix C, creating new symm.mat
%
% If input is non-symmetric matrix:
% input matrix contains bounds,
% start from midpoint xor random point in interval iff flag=0 xor ~0 resp.
%
% Action in all cases: perturb matrix until 6th s.val < toler
%	by following gradient that reduces rank
%
function [ newMat, D ] = mconf(D, toler, flag)
 
fprintf('||||||||||||||||||||||||||||||||||||||||||||||||||\n');
fprintf('||||| new mconf() run ||||||||||||||||||||||||||||\n');
fprintf('||||||||||||||||||||||||||||||||||||||||||||||||||\n');

% This part of the function distinguishes between different kinds of input 

if ~isequal(D-D',zeros(size(D)))

  fprintf('Input bounds: NOT symmetric\n');
  boundMat = D;
  PertEntr = perbasis (boundMat,eps,Inf);
  % disp (PertEntr);
  D = bnd2mid(boundMat,flag);
  fprintf ('Starting with symm.mat=\n'); disp (D);

else

  boundMat = maxbound(size(D,1)); 
  PertEntr = permall(size(D,1)); 
  if flag > 0 
  	fprintf('symmetric input to be perturbed\n');
	D = pertMmat(D,PertEntr,flag);
	% printmat(D);
  else
	fprintf('Symmetric input matrix used as start.point\n');
  end;

end;
	
% The rest of the function works the same for all kinds of input 

fprintf('Start S-vals %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e...\n',svds(D,6));
fprintf('Perturb %d entries to make s-val.6 < %1.1e\n',size(PertEntr,1),toler);
solutions = 1;	% #trial directions away from dist.mat: change to generalize function
totaltime = 0;
 
for k=1:solutions

    initsvals = svds(D,7);
 
    t0 = clock;
    % flops(0);
    [newMat, news6, iterCount] = svred(D, PertEntr, toler, boundMat);

    % ff = flops;
    time = etime(clock,t0);

    if (news6 > toler)
	fprintf('cant make sing.val.6 < %1.1e after max#iterations\n',toler);
	if news6 > initsvals(6);
		fprintf('6th sing.value increased: quitting\n');
		break;
	end;
    end; 
    news = svds(newMat,7);

    fprintf('InitSingval(1st..5,6,7th)        FinalSingval(1st..5,6,7th)      #loops sec\n'); 
    fprintf('%1.1e;%1.1e,%1.1e,%1.1e  ',initsvals(1),initsvals(5),initsvals(6),initsvals(7));
    fprintf('%1.1e;%1.1e,%1.1e,%1.1e %6d %2.1f ',news(1),news(5),news6,min(news),iterCount,time);

    totaltime = totaltime + time;
    offBounds = violatebnd (boundMat, newMat, PertEntr);
    fprintf(':violated %d bounds:',size(offBounds,1));
    for i = 1:size(offBounds,1);
      fprintf('%d,%d=[%f',offBounds(i,:),boundMat(offBounds(i,2),offBounds(i,1)));
      fprintf('<%f] ',boundMat(offBounds(i,1),offBounds(i,2)));
    end; 
    fprintf('\n');
 
    % drawing
    % metMat = metric(newMat); xyzMat = embed(metMat); drawmol(xyzMat);
 
end; % for k
 
if 0 & (solutions > 0)
    fprintf('\nNumber of conformations: %2d', solutions);
    fprintf('\nAverage time per conformation = %5.2f sec\n',totaltime/solutions);
end;
