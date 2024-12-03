% Input: dim of bordered symmetric dist.matrix
% Return: bordered matrix expressing 0,Inf bounds st.
%       all non-border entries below diagonal = 0
%       all non-border entries above diagonal = Inf
%       borders filled in with zeros coz never used
%
function mxBnd = maxbound(dim)
  mxBnd = zeros(dim,dim);
  for i=2:dim for j=i+1:dim
        mxBnd(i,j) = Inf;
  end; end;
% maxbound