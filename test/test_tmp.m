clear; clc; close all;

load ash219-opt.mat

X = data.X;

testoptprecond(data, "suitesparse", "?", "LR");

% [d1, d2] = gettwosidedprecond(X, 5);
% 
% XTX = X' * d2 * X;
% cond(sqrt(d1) \ (sqrt(d1) \ XTX)')