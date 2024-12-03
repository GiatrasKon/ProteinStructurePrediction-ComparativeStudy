%   Nikitopoulos, Emiris
%   METRIC(D) finds the metric (Gram) matrix G of a Cayley-Menger matrix D.
%   Computes the squared distances to the center of mass from the
%   squared distances among the points, in the Cayley-Menger matrix.
%
%   See also EMBED, ERRORF.
%
function G  = metric(D)
 
N = size(D,1);
D = sqrt(D);
D = D(2:N,2:N);
N = N-1;
 
for i = 1:N
    tmp1 = 0;
    for j = 1:N
        tmp1 = tmp1 + D(i,j)*D(i,j);
    end
    tmp1 = tmp1/N;
    tmp2 = 0;
    for j = 2:N for k = 1:j-1
            tmp2 = tmp2 + D(j,k)*D(j,k);
    end; end;
    tmp2 = tmp2/(N*N);
    Do(i) = tmp1 - tmp2;
end
 
for i = 1:N
    for j = 1:N
        G(i,j) = ( Do(i) + Do(j) - D(i,j)*D(i,j) ) / 2;
    end;
end
