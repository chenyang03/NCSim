
function [ direction ] = Direct_h( X, Y )
N = length(X);
if (X==Y)
    Z = rand(1,N);
    direction = Z / Dist_h(Z, zeros(1, N));
else
    result = X - Y;
    result(N) = X(N) + Y(N);
    direction = result / Dist_h(X, Y);
end;

