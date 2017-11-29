function [ beta_threshold, betahat, betahatSig ] = bootstrap_fdr( X, Y, q, nIter)
%BOOTSTRAP_FDR Find beta threshold to achieve FDR of q using bootstrapping
%testing
%   X: design matrix (nObs x p)
%   Y: BOLD responses (nObs x nVox)
%   q: desired FDR
%   nIter: # of iterations to use

Y = Y';
betahat = ols_simple(X,Y);
[n,p] = size(X);
nVox = size(Y,1);
betaStrap = zeros(p,nVox,nIter);
resid = Y'-X*betahat;

for i = 1:nIter
    perm = randsample(n,n,1);
    Ystrap = X*betahat + resid(perm,:);
    betaStrap(:,:,i) = ols_simple(X,Ystrap');
end

testvals = zeros(p,500);
for i = 1:p
    temp = linspace(0,max(abs(betahat(i,:))),501);
    testvals(i,:) = temp(1:500);
end
sigcount = zeros(p,500);
nullcountBoot = zeros(nIter,p,500);
for i = 1:500
    for ptest = 1:p
        testval = testvals(ptest,i);
        sigfind = find(abs(betahat(ptest,:))>testval);
        sigcount(ptest,i) = length(sigfind);
        nullcountBoot(:,ptest,i) = squeeze(sum((betaStrap(ptest,sigfind,:).*repmat(betahat(ptest,sigfind),[1 1 nIter]))<0,2));
    end
end

fdr_mean = squeeze(mean(nullcountBoot))./sigcount;
beta_threshold = zeros(p,1);
for ptest = 1:p
    bthresh_idx = find(fdr_mean(ptest,:)<=q,1);
    if isempty(bthresh_idx)
        beta_threshold(ptest) = NaN;
    else
        beta_threshold(ptest) = testvals(ptest,bthresh_idx);
    end
end

betahatSig = abs(betahat) > repmat(beta_threshold,[1 nVox]);

