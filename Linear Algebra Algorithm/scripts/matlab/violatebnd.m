% INPUT: check whether BOUNDS in Bnds respected by candidate mat Cand
% 	as far as (perturbable) entries in Entr are concerned
%
% RETURN: (number of) violated entries
%
function OffEnt = violatebnd (Bnds, Cand, Entr);
 
Np = size(Entr,1);
OffEnt = [];	
 
for i=1:Np
   x = Entr(i,1); y = Entr(i,2); 
   if ( Bnds(x,y) < Cand(x,y) ) | ( Bnds(y,x) > Cand(y,x) )
	OffEnt = [OffEnt; x,y];
   end; 
end;
%
% end function()
