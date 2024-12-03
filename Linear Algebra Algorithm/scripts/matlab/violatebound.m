% INPUT: check BOUNDS in B wrt dist.mat D for pert'd entries R
%
% RETURN: exists VIOLATION iff return NumV>0; 0 iff no violation
%       NumV counts #violated bounds; they are also printed
%
function NumV = violatebound(B,D,R);
 
Np = size(R,1);
NumV = 0;
 
for i=1:Np
   x = R(i,1); y = R(i,2);
   if ( ( B(x,y) < D(x,y) ) | ( B(y,x) > D(y,x) ) )
      NumV = NumV+1;
      % fprintf('violated D(%d,%d)=%f off [%f,%f]\n',x,y,D(x,y),B(y,x),B(x,y));
      % break;
   end;
end;
% if (NumV>0) fprintf('\n'); end;
%
% end violatebound()
