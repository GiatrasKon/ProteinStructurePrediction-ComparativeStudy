% I.Emiris, created 6/02, last modified 7/03.
%
% INPUT: Cand = cand.dist.mat
%	boundMat = non-symmetric matrix of up/low bounds
%	"PertEntr" shows entries which are pert'ble
%	up/low bounds from mconf()
%
% RETURN: new basis for perturbations: divides by 5 (arbitrary): to improve!
%
function [numOff, newEntr] = basis (Cand, boundMat, PertEntr)

  stretch = 1;				% accept beyond interval length

  numEntr = size(PertEntr,1);
  numOff = 0;
  newEntr = [];			% init
  % disp(violatebnd(boundMat, Cand, PertEntr));
  % fprintf(' Candidate matrix (without borders)=\n'); disp(Cand(2:size(Cand,1),2:size(Cand,1)));
  fprintf('Pert.Entries: Basis * 1e+8=\n'); disp(1e+8 * PertEntr(:,3)');

  for i=1:numEntr
      x = PertEntr(i,1); y = PertEntr(i,2);

      if 1 | PertEntr(i,3) > 1e-3;		% option to remove

        off = 2 * abs( Cand(x,y) - (boundMat(x,y)+boundMat(y,x))/2 ) / (boundMat(x,y)-boundMat(y,x)) ; 
        % if ( boundMat(x,y) < Cand(x,y) ) | ( boundMat(y,x) > Cand(x,y) ); 

        if off > stretch ; 
	  newEntr = [ newEntr; x,y, PertEntr(i,3) / 5 ];	% / 5 arbitrary
  	  numOff = numOff + 1;
  	  fprintf(' cand %f violates %.2f<%.2f\n',Cand(x,y),boundMat(y,x),boundMat(x,y)); 
	else
	  newEntr = [ newEntr; x,y, PertEntr(i,3) ];
	end;

      else
	newEntr = [ newEntr; x,y, 0 ];
  	fprintf(' Zero-ing %d-th basis element:',i); % disp(PertEntr');
      end;
	
  end; % for 

  % fprintf(' %d violations out of %d perturbable\n',numOff,size(PertEntr,1)); 

  if 0 & numOff>0 ; 
	fprintf(' basis elements < 1 as follows:\n');
	disp([find(PertEntr(:,3)<1)]');
  end;
%
% end basis()
