% T.Nikitopoulos (and I.Emiris) 2002,2003
% INPUT: boundary matrix B w/low,upper bounds
%       P = list of i,j entries to be checked
% RETURN: new boundary matrix satisfying triangle ineq. at P
%
function B = triangle(B, P)
 
Natom = size(B,1);
Nunkn = size(P,1);
 
% fix the Upper and Lower parts of B
Upper = zeros(Natom,Natom);
Lower = zeros(Natom,Natom);
 
for i=1:Natom
   for j=i+1:Natom
      Upper(i,j) = B(i,j);
      Upper(j,i) = B(i,j);
      Lower(i,j) = B(j,i);
      Lower(j,i) = B(j,i);
   end
end
 
% Start a variation of Floyd's algorithm
for l=1:Nunkn
   i = P(l,1);
   j = P(l,2);
   for k=2:Natom
 
% path lengths in left-hand network
      if Upper(i,j) > Upper(i,k) + Upper(k,j)
         Upper(i,j) = Upper(i,k) + Upper(k,j);
         Upper(j,i) = Upper(i,j);
      end
% path lengths from left to right-hand network
      if Lower(i,j) < Lower(i,k) - Upper(k,j)
         Lower(i,j) = Lower(i,k) - Upper(k,j);
         Lower(j,i) = Lower(i,j);
      elseif Lower(i,j) < Lower(j,k) - Upper(k,i)
         Lower(i,j) = Lower(j,k) - Upper(k,i);
         Lower(j,i) = Lower(i,j);
      end
% check for triangle inequality violations
      if Lower(i,j) > Upper(i,j)
         fprintf('  Erroneous bounds %d %d using %d\n'i,j,k);
      end
   end %k
end %l
 
% construct the new boundary matrix B
B = zeros(Natom,Natom);
for i=1:Natom
   for j=i+1:Natom
      B(i,j) = Upper(i,j);
      B(j,i) = Lower(i,j);
    end
end
%
% end of function triangle
