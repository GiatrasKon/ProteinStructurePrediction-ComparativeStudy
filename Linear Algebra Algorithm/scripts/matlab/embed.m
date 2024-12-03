% Outputs: n x 3 matrix of Cartesian coordinates for metric matrix G.
%	It is the best result in a sense of the Frobenius norm.
%	Puts origin in center of mass.
%	n = #points, must be >= 3.
%
function X = embed(G)
 
[U S V] = svd(G);
 
x = sqrt(S(1,1))*U(:,1);
y = sqrt(S(2,2))*U(:,2);
z = sqrt(S(3,3))*U(:,3);
 
X = [x, y, z];
