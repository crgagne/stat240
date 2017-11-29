function [ beta_threshold, betahat, betahatSig ] = permutation_fdr( X, Y, q, nperm)
%PERMUTATION_FDR Find beta threshold to achieve FDR of q using permutation
%testing
%   X: design matrix (nObs x p)
%   Y: BOLD responses (nObs x nVox)
%   q: desired FDR
%   nperm: # of permutations to use

Y = Y';
betahat = ols_simple(X,Y);
[n,p] = size(X);
nVox = size(Y,1);
betaPerm = zeros(p,nVox,nperm);

for i = 1:nperm
    perm = randperm(n);
    Xperm = X(perm,:);
    betaPerm(:,:,i) = ols_simple(Xperm,Y);
end

testvals = zeros(p,500);
for i = 1:p
    temp = linspace(0,max(abs(betahat(i,:))),501);
    testvals(i,:) = temp(1:500);
end
sigcount = zeros(p,500);
sigcountPerm = zeros(nperm,p,500);
for i = 1:500
    for ptest = 1:p
        testval = testvals(ptest,i);
        sigcount(ptest,i) = sum(abs(betahat(ptest,:))>testval);
        sigcountPerm(:,ptest,i) = squeeze(sum(abs(betaPerm(ptest,:,:))>testval,2));
    end
end

fdr_mean = squeeze(mean(sigcountPerm))./sigcount;
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

