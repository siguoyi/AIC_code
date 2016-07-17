function [Phi]=pnmat(M,m)
% Generates a M x m random matrix whose entries are drawn from {-1,+1} with equal probability.
% Usage: [Phi]=pnmat(M,m)

Phi=sign(rand(M,m)-0.5);