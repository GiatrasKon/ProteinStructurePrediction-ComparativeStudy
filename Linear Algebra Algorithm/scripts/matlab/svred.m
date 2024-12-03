% I.Emiris, created 06/02, last modified 07/03
%
% Inputs: input matrix inpMat
%	px3 matrix PertEntr expresses ENTRIES to be PERTURBED and basis
%	toler = threshold of 6th sing.val
%	boundMat = matrix of lower/upper bounds (below/above diagonal)
%
% RETURN: newMat : whose 6-th sing.val = news6
% 	news6 <= toler
%	iterCount = main#iterations needed ie. #actual.perturbations
% 
% perform local optimization for sing.val "which=6"
% bound "iterCount"=total#iterations by "iterBoundTot"
% bound "badpert"=#iterations NOT decreasing 6th s.value by "iterBoundBad"
%
function [newMat, news6, iterCount] = svred(inpMat, PertEntr, toler, boundMat)

%
% INITIALIZATION SECTION
%

   iterBoundTot = 20 ;		% total#perturbations
   if (iterBoundTot < 2) iterBoundTot = 2; end;
   iterBoundBad = 10 ;		% iterations until 6th s.value decreased
   if (iterBoundTot < iterBoundBad) iterBoundTot = iterBoundBad; end;
   maxIterBasis = 0;
   bndIterBasis = 10;		% iterations w/varying basis to respect bounds

   which = 6; 			% Minimize the 6th sing. val.

   dimMat = size(inpMat,1);
   numPer = size(PertEntr,1);

   mmold = inpMat; 		% last best 
   newMat = inpMat;

   pertR = zeros(dimMat,dimMat);
   Delta = zeros(dimMat,dimMat);

   [U S V] = svd(newMat);
   candSigma_k = S(which,which);
   svold = candSigma_k;
 
   iterCount = 0;
   if ( candSigma_k <= toler ) news6 = candSigma_k; end; 
   if (svold < candSigma_k); badpert = 1; else badpert = 0; end;

%
% DOUBLE ITERATION SECTION
%
 
while (iterCount < iterBoundTot) & (badpert < iterBoundBad) & (candSigma_k > toler)

   iterCount = iterCount + 1;

%
% Step : Define Delta such that u_' * Delta * v_ = u_' * newMat * v_
% Iterate over possible perturbations until one respecting bounds
%
  
   numOff = 1;			% exist off bounds perts
   iterBasis = 0;		% count basis-adjustments due to off.bounds

   while (numOff>0 & iterBasis<bndIterBasis)
     % if (numOff<2 & iterBasis>bndIterBasis-2); break; end;

     iterBasis = iterBasis + 1;
     matE = [];			% cols correspond to pert.entrie
     vecF = [];			% col.vector = matrix p x 1
  
     for i = which:dimMat;

       v_i = V(:,i);		% j=i
       u_i = U(:,i);
       vecF = [vecF; v_i' * newMat * u_i];
       vecF = [vecF ; zeros(dimMat-i,1)];
       % for j= i+1:dimMat; vecF = [vecF ; 0]; end;
 
       for j = i:dimMat;
         Erow = [];
         v_j = V(:,j);
         for k=1:numPer;
	 	x = PertEntr(k,1); y = PertEntr(k,2);
	 	Erow = [ Erow, v_j(x)*u_i(y) + v_j(y)*u_i(x) ];
         end;
         matE = [matE ; Erow ];		% ; (Erow .* PertEntr(:,3)') ];
       end; % for j

     end; % for i
 
     % Xi = matE \ vecF;		% mldivide 
     Xi = pinv(matE) * vecF;		% works even w/sing.matE

     % fprintf('Xi = E+(%d %d) F(%d) =\n',size(matE),size(vecF)); disp(Xi');
     % fprintf('svd(E:%dx%d)=',size(matE));disp([svds(matE)]); 
     % fprintf('DIFFERING:\n'); printmat(PertEntr(:,3)'); printmat([Xi - (E0 \ vecF)]);
     % E0 as if basis=1: diff. not at basis elts<>1 so isnt undone by elts<>1
 
     for i=1:numPer;
       x = PertEntr(i,1); y = PertEntr(i,2);
       Delta(x,y) = Xi(i) * PertEntr(i,3);
       Delta(y,x) = Delta(x,y);
     end;
 
     [numOff, newEntr] = basis (inpMat-(pertR+Delta), boundMat, PertEntr);
     PertEntr = newEntr;

     % fprintf(' (%d) kx3 matrix PertEntr.transposed =\n',iterBasis); disp(PertEntr');
     % fprintf('Perturbing by\n'); disp(pertR+Delta);
     fprintf('Perturb by\n'); disp((Xi .* PertEntr(:,3))');
     fprintf('**************************************************\n');

   end; % while numOff>0
     
   if numOff>0;
   	% fprintf('Final perturbation basis =\n'); disp(PertEntr'); 
	error('BAD basis');
   end;
   % fprintf(' BAD'); else fprintf(' OK'); end;
   % fprintf(' basis after %d loops\n',iterBasis);
   if iterBasis>maxIterBasis; maxIterBasis=iterBasis; end;
 
%
% Step : Compute newMat; SVD of newMat - a_k * P to get u_, v_, candSigma_k, f'(a)
%
 
   % newMat = newMat - Delta;		% worse than newMat= inpMat-pertR 

   pertR = pertR + Delta; 
   newMat = inpMat - pertR; 

%
% Step : Compute SVD of newMat - a_k * P to get u_, v_, candSigma_k and f'(a)
%
 
   [U S V] = svd(newMat);
   candSigma_k = S(which,which);
 
   if ( candSigma_k <= toler )
	news6 = candSigma_k;
   	% fprintf('Final perturbation basis =\n'); disp(PertEntr');
   	fprintf('done with <= %d loops per basis correction\n',maxIterBasis);
	break;
   else fprintf('......%d...... candS6=%e > toler=%.2e\n',iterCount,candSigma_k,toler);
   end;
   
   if (svold > candSigma_k);
	fprintf(' svold=%f > cand=%f sets badpert=%d->0\n',svold,candSigma_k,badpert);
	badpert = 0;
   	svold = candSigma_k; 
   	mmold = newMat; 
   else
	badpert = badpert + 1;
   end;
   
   % fprintf('~~~ norm Delta 1:%f, fro=%f\n',norm(Delta,1),norm(Delta,'fro'));
 
end % while iterCount < iterBoundTot
 
%
% If the method doesn't perform well, return input
%
 
if (candSigma_k > toler)
   % fprintf('Final perturbation basis =\n'); disp(PertEntr'); 
   newMat = mmold;
   news6 = svold;
end 
%
% end svred()
