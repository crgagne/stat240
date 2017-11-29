


load('low.mat')

% test errors
%  e1 = E(:,:,1)
% imshow(e1'*e1)
%

FDR = [];

betahatSigBHall = zeros(2,900,size(B,3));
betahatSigpall = zeros(2,900,size(B,3));

for i = 1:size(B,3)
    
   % BH
   %Bhat=ols_simple(X(:,:,i),Y(:,:,i)');
   %imshow(reshape(Bhat(1,:),30,30))
   %Yhat = X(:,:,i)*Bhat;
   s = size(Y);
   q = s(2);
   n = s(1);
   s = size(B);
   p = s(1);
   
   %sigmahat = sqrt(1/(n-p)*sum((Yhat-Y(:,:,i)).^2)); % 1xq
   
  % covB= inv(X(:,:,i)'*X(:,:,i))*sigmahat
   %diagXX = diag(inv(X(:,:,i)'*X(:,:,i)));
   %P=[];
   %betahatSigBH=[];
   %for bb = 1:p
   %   SE =  diagXX(bb)*sigmahat;
   %   T = Bhat(bb,:)./SE;
   %   P(bb,:) =1-tcdf(abs(T),n-p); % check this
   %end
   
   
   [beta,bias,stdErr,residuals,sigma]=ols(X(:,:,i),Y(:,:,i)');
   T = beta./stdErr;
   P =(1-tcdf(abs(T),n-p))*2;
   Pcat=P(:);
   pthresh=bh(Pcat,0.1);
   betahatSigBH=P<=pthresh;
   %imshow(reshape(P(1,:),30,30)) 
   %imshow(reshape(betahatSigBH(1,:),30,30)) 
   betahatSigBHall(:,:,i)=betahatSigBH;
  
   % Permutation
  [  T_threshold, T_threshold_05, T_threshold_95,betahat, T, betahatSigp,betahatSigp_05,betahatSigp_95]= permutation_fdr_t( X(:,:,i), Y(:,:,i), .1, 100);
  betahatSigpall(:,:,i)=betahatSigp;
  % imshow(reshape(betahatSigp(1,:),30,30))
   %imshow(reshape(T(1,:,idbquit),30,30))
   
   % Bootstrap 
  % [ beta_thresholdb, betahatb, betahatSigb ] = bootstrap_resample_fdr( X(:,:,i), Y(:,:,i), .1, 100);
   % imshow(reshape(betahatSigb(1,:),30,30))
   
   % Assess Performance %
   Bbin = B(:,:,i);
   Bbin(Bbin>0)=1;
   for model =1:2
       switch model
           case 1
               betahatSig=betahatSigBH;
               tp=sum(sum((Bbin==1).*(betahatSig==1)));
               fp=sum(sum((Bbin==0).*(betahatSig==1)));
               fn=sum(sum((Bbin==1).*(betahatSig==0)));
               tn=sum(sum((Bbin==0).*(betahatSig==0)));
               FD = fp;
               FDR(i,model) = fp/(fp+tp);
           case 2
               betahatSig=betahatSigp;
               tp=sum(sum((Bbin==1).*(betahatSig==1)));
               fp=sum(sum((Bbin==0).*(betahatSig==1)));
               fn=sum(sum((Bbin==1).*(betahatSig==0)));
               tn=sum(sum((Bbin==0).*(betahatSig==0)));
               FD = fp;
               FDR(i,model) = fp/(fp+tp);
               
               % 05 and 95 
               betahatSig=betahatSigp_95;
               tp=sum(sum((Bbin==1).*(betahatSig==1)));
               fp=sum(sum((Bbin==0).*(betahatSig==1)));
               fn=sum(sum((Bbin==1).*(betahatSig==0)));
               tn=sum(sum((Bbin==0).*(betahatSig==0)));
               FD = fp;
               FDR_ci(i,1) = fp/(fp+tp);
               betahatSig=betahatSigp_05;
               tp=sum(sum((Bbin==1).*(betahatSig==1)));
               fp=sum(sum((Bbin==0).*(betahatSig==1)));
               fn=sum(sum((Bbin==1).*(betahatSig==0)));
               tn=sum(sum((Bbin==0).*(betahatSig==0)));
               FD = fp;
               FDR_ci(i,2) = fp/(fp+tp);
               if (FDR_ci(i,2)>=FDR(i,model))&&(FDR_ci(i,1)<=FDR(i,model))
                    included(i) = 1;
               end
               
           case 3
               betahatSig=betahatSigb;
       end
       

   end
   
end

figure();
figure('Color','w')
imagesc(reshape(betahatSigBH(1,:),30,30)+reshape(betahatSigBH(2,:),30,30))
caxis([0,1])
colormap('gray')

figure();
figure('Color','w')
imagesc(reshape(betahatSigp_05(1,:),30,30)+reshape(betahatSigp_05(2,:),30,30))
caxis([0,1])
colormap('gray')
c

figure('Color','w')
imagesc(reshape(betahatSigp_05(1,:),30,30)'+reshape(betahatSigp_05(2,:),30,30)')
caxis([0,1])
colormap('gray')


figure('Color','w')
imagesc(reshape(betahatSigp_95(1,:),30,30)'+reshape(betahatSigp_95(2,:),30,30)')
caxis([0,1])
colormap('gray')


figure('Color','w')
imagesc(reshape(betahatSigp(1,:),30,30)'+reshape(betahatSigp(2,:),30,30)')
caxis([0,1])
colormap('gray')



figure('Color','w')
imagesc(reshape(P(1,:,1),30,30))
colormap('jet')
figure('Color','w')
imagesc(reshape(P(2,:,1),30,30))
colormap('jet')


figure();
imagesc(reshape(betahatSigBH(1,:),30,30)+reshape(betahatSigBH(2,:),30,30))
caxis([0,1])
colormap('gray')


figure('Color','w')
imagesc(reshape(B(1,:,1),30,30)+reshape(B(2,:,1),30,30))
caxis([0,1])
colormap('gray')

figure('Color','w')
imagesc(reshape(mean(betahatSigBHall(1,:,:),3),30,30)+reshape(mean(betahatSigBHall(2,:,:),3),30,30))
colormap('jet')
caxis([0,1])

figure('Color','w')
imagesc(reshape(mean(betahatSigpall(1,:,:),3),30,30)+reshape(mean(betahatSigpall(2,:,:),3),30,30))
%caxis([0,max(mean(betahatSigpall(1,:,1:1000)))])
caxis([0,1])
colormap('jet')


figure();
hist(FDR(:,1))
mean(FDR(:,1))

figure();
hist(FDR(:,2))
mean(FDR(:,2))
%figure();
%hist(FDR(:,3))
%mean(FDR(:,3))
mean(included)


%%% BH
figure('Color','w')
hist(FDRlow(:,2),linspace(0,1,20))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.3)
hold on
hist(FDRhigh(:,2),linspace(0,1,20))
h2 = findobj(gca,'Type','patch');
set(h2,'facealpha',0.3)
xlabel('Actual FDR per simulation')
legend('low spatial correlation','high spatial correlation')
xlim([0,.8])


