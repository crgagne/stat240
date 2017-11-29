function [beta,bias,stdErr,residuals,sigma] = ols(x,y,weights,addBias);
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
%
% <weights> - (N x 1) vector of observation weights (e.g.
%             standard deviation of input vectors) to use in
%             weighted least squares. If not provided, assumes
%             equal weighting for all inputs.
%
% <addBias> - flag to indicate whether or not to add a acolume
%             of ones to x for bias terms.
%
%OUTPUT:
% <beta>:   - (nDim x N) matrix of mixing weights for input
%
% <bias>:   - (1 x N) vector of dc terms
%
% <stdErr>: - ((nDim+1) x N) matrix standard error of estimated
%             mixing weights (including error for bias term).
%
% <resid>:  - model fit residuals (y - yhat)
%
% <sigma>:  - unbiased estimator of input noise standard
%             deviation
%
%-------------------------------------------------------------
%DES

%if notDefined('addBias')
	addBias = 0;
%end	
% ADD BIAS TERM, IF NECESSARY
if addBias
	x = [ones(size(x,1),1),x];
end

[nMeasure, nDim] = size(x);


%if notDefined('weights')
	W = eye(nMeasure);
%else
%	if numel(weights) == (nMeasure - 1)
%		weights = [1;weights]; % ADD BIAS TERM
%	end
%	W = diag(1./weights.^2);
%end

xc = pinv(x'*W*x); 

% RUN OLS
beta = xc*(x'*(W*y'));  % params x cases

% GET RESIDUALS AND STANDARED ERROR
yhat = x*beta;
residuals = y' - yhat;
sqErr = sum(residuals.^2,1);  % (1 x cases)

dF = nMeasure-nDim; % DEGREES OF FREEDOM
sigma = sqrt(zerodiv(sqErr,dF,NaN)); % NOISE ESTIMATOR
stdErr = sqrt(diag(xc))*sigma;  % (params x cases)

if addBias
	bias = beta(1,:);
	beta = beta(2:end,:);
else
	bias = [];
end
