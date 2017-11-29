function [ T_threshold, T_threshold_05, T_threshold_95,betahat, T, betahatSig,betahatSig_05,betahatSig_95 ] = permutation_fdr_t( X, Y, q, nperm)
%PERMUTATION_FDR Find beta threshold to achieve FDR of q using permutation
%testing
%   X: design matrix (nObs x p)
%   Y: BOLD responses (nObs x nVox)
%   q: desired FDR
%   nperm: # of permutations to use

Y = Y';
[betahat,~,stdErr] = ols(X,Y,[],0);
T = betahat./stdErr;
[n,p] = size(X);
nVox = size(Y,1);
betaPerm = zeros(p,nVox,nperm);
TPerm = zeros(p,nVox,nperm);

for i = 1:nperm
    perm = randperm(n);
    Xperm = X(perm,:);
    [betaPerm(:,:,i),~,stdErr] = ols(Xperm,Y);
    TPerm(:,:,i) = betaPerm(:,:,i)./stdErr;
end

testvals = zeros(p,500);
for i = 1:p
    temp = linspace(0,max(abs(T(i,:))),501);
    testvals(i,:) = temp(1:500);
end
sigcount = zeros(p,500);
sigcountPerm = zeros(nperm,p,500);
for i = 1:500
    for ptest = 1:p
        testval = testvals(ptest,i);
        sigcount(ptest,i) = sum(abs(T(ptest,:))>testval);
        sigcountPerm(:,ptest,i) = squeeze(sum(abs(TPerm(ptest,:,:))>testval,2));
    end
end

fdr_mean = squeeze(median(sigcountPerm))./sigcount;

fdr_05 = squeeze(quantile(sigcountPerm,0.05))./sigcount;
fdr_95 = squeeze(quantile(sigcountPerm,0.95))./sigcount;

T_threshold = zeros(p,1);
T_threshold_05 = zeros(p,1);
T_threshold_95 = zeros(p,1);
for ptest = 1:p
    bthresh_idx = find(fdr_mean(ptest,:)<=q,1);
    bthresh_idx_05 = find(fdr_05(ptest,:)<=q,1);
    bthresh_idx_95 = find(fdr_95(ptest,:)<=q,1);
    if isempty(bthresh_idx)
        T_threshold(ptest) = NaN;
    else
        T_threshold(ptest) = testvals(ptest,bthresh_idx);
    end
    if isempty(bthresh_idx_05)
        T_threshold_05(ptest) = NaN;
    else
        T_threshold_05(ptest) = testvals(ptest,bthresh_idx_05);
        end
    if isempty(bthresh_idx_95)
        T_threshold_95(ptest) = NaN;
    else
        T_threshold_95(ptest) = testvals(ptest,bthresh_idx_95);
    end
    
end

betahatSig = abs(T) > repmat(T_threshold,[1 nVox]);
betahatSig_05 = abs(T) > repmat(T_threshold_05,[1 nVox]);
betahatSig_95 = abs(T) > repmat(T_threshold_95,[1 nVox]);
