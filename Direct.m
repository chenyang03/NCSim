function [ direction ] = Direct( X, Y )
d = Dist(X,Y);
if (d==0)
    Z = rand(1,length(X));
    direction = Z / (Z*Z').^0.5;
else
    direction = (X - Y) / Dist(X,Y);
end;

