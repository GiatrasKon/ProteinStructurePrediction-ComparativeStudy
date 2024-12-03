% RETURN: k x 2 matrix of k entry pairs to be pert'd
%       = Entr minus entries in "ro"
%
function nEntr = npermat(Entr,ro)
 
N = size(ro,1);
if (N == 0) nEntr = Entr; break; end;
 
Nr = size(Entr,1);
T = Nr-N;
nEntr = zeros(T,2);
step = 0;
 
for i = 1:Nr
    flag = 1;
    for j = 1:N
      if (Entr(i,1) == ro(j,1) & Entr(i,2) == ro(j,2) ) flag = 0; break; end;
    end
    if (flag == 1)
        step = step + 1;
        nEntr(step,1) = Entr(i,1);
        nEntr(step,2) = Entr(i,2);
    end
end
%
% end npermat()
