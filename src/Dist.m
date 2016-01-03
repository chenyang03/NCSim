function [ distance ] = Dist( X, Y )
distance = ((X - Y)*(X - Y)').^0.5;

