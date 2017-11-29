function [beta] = ols_simple(x,y)
%[beta,bias,stdErr,resid,sigma] = ols(x,y,weights,addBias);
%------------------------------------------------------------
% (weighted) Ordinary Least Squares: Solves
%
%                  y = x*beta + noise
%
%
% by minimizeing
%
%                     ||x*Beta - y||^2
%------------------------------------------------------------
%INPUT:
% <x>:      - (nObs x nDim) matrix of independent input
%             variables
%
% <y>:      - (N x nObs) matrix of N dependent output units

xc = pinv(x'*x); 

% RUN OLS
beta = xc*(x'*y');  % params x cases
