function [ distance ] = Dist_h( X, Y )

    N=length(X);
    distance = ((X(1:N-1) - Y(1:N-1))*(X(1:N-1) - Y(1:N-1))').^0.5+X(N)+Y(N);
%distance = ((X - Y)*(X - Y)').^0.5;

