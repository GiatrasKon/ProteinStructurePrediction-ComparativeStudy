% input = bounds B: upper above, lower below diagonal, equal iff dist.known
%	border of 1's in top row, left column; 0 diagonal
%	move : scalar flag in [0..1]
%
% return = symm matrix, bordered matrix w/1's, 0 diagonal
%	each dist = midpoint iff move=0
%		  = pert'd from midpoint by (+-1)move*rand*interv.size/2
%		    NO perturabtion if up/low bounds equal
%
function R = bnd2mid (B,move)

  R = B;
  dim = size(B,1);
  UplusL = B + B';

  if move == 0 ;

    for i=2:dim for j=i+1:dim
      R(i,j) = UplusL(i,j) / 2;
      R(j,i) = R(i,j);
    end; end; 

  else; % move ~= 0

    if abs(move)>1 | move<0; move=1; end; 
    UminusL = B - B';
    % disp(UminusL);

    for i=2:dim for j=i+1:dim
	if B(i,j) ~= B(j,i);
	  if rand>0.5 ; sgn=1; else sgn=-1; end;
          R(i,j) = (UplusL(i,j) + sgn*rand*move*UminusL(i,j)) / 2;
	end;
	R(j,i) = R(i,j);
    end; end; 

  end; % move
%
% end function
