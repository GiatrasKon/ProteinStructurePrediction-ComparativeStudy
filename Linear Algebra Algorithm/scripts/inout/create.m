% create bndry matrix

% function B = create(B83)
% B = B83(1:sz,1:sz);
% if (B-B') error('submat not symm'); end;
% fd=fopen('./bounds83.m','w');

fd=fopen('./39mol1.m','w');
% fprintf(fd,'function B = bounds19()\n\nB=[%.8f, ',M(1,1));

mdim = length(mat);

fprintf(fd,'sing_vals = [\t% of matrix below\n');
for j=1:mdim-1 fprintf(fd,'%.10f, ',svs(j)); end;
fprintf(fd,'%.10f];\n\n',svs(mdim));

fprintf(fd,'dist39mol1 = [\t% one almost-conformation\n%.10f, ',mat(1,1));
for j=2:mdim-1 fprintf(fd,'%.10f, ',mat(1,j)); end;
fprintf(fd,'%.10f;\n',mat(1,mdim));

for i=2:mdim-1
	for j=1:mdim-1 fprintf(fd,'%.10f, ',mat(i,j)); end;
	fprintf(fd,'%.10f;\n',mat(i,mdim));
end;

for j=1:mdim-1 fprintf(fd,'%.10f, ',mat(mdim,j)); end;
fprintf(fd,'%.10f];\n',mat(mdim,mdim));
fclose(fd);
