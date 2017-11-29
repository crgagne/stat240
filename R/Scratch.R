fmri.modelcomparisons<-function(nsims,map,design){
  
  
  
  ##set up storage
  est.mse<-c(0) #probably want per q too..!
  pred.mse<-c(0)
  pred.cor<-c(0)
  modellist=map$modellist
  model=map$modellist[1]
  df<-data.frame(model,est.mse,pred.mse,pred.cor,SNR=map$SNR,p=length(design$conditions),n=design$TRS)
  levels(df$model)=modellist
  
  
  maps.Bhat<-array(0,c(map$dim1,map$dim2,length(design$conditions),length(modellist),nsims))
  maps.pred.cor<-array(0,c(map$dim1,map$dim2,length(modellist),nsims))
  #maps.B<-array(0,c(map$dim1,map$dim2,dim(B)[1],length(modellist)))
  peaks.ts.Y=array(0,c(design$TRS,length(map$which.b),length(modellist)))
  peaks.ts.Yhat=array(0,c(design$TRS,length(map$which.b),length(modellist)))
  
  
  for (subject in seq(nsims)){
    
    
    
    
    ##### Collect Training Data
    design.train<-generate_event_design(design)
    map.train<-generate_regions_signal(design.train,map)
    
    #### Collect Testing Data
    #### maybe change up study parameters.. ### 
    design.test<-generate_event_design(design)
    map.test<-generate_regions_signal(design.test,map)
    
    # get distance matrix for search light
    distrank<-data.matrix(read.table("sim_distrank_masked.txt"))
    #distrank<-distrank[map.train$e.idx,]
    # change the distrank, so the indexes reflect mask
    
    
    # extract data
    X.train<-map.train$X # nxq
    Y.train<-map.train$Y # needs to be nxq
    B<-t(map.train$B)
    X.test<-map.test$X
    Y.test<-map.test$Y # needs to be nxq
    
    p<-length(design$conditions)
    
    
    # loop through models
    mm=1
    for (mname in modellist){
      
      funname = paste(mname,".fit",sep='')
      
      # Fit 
      if (grepl('cw',mname)){
        # Search Light
        fitted<-searchlight.fit(funname,X.train,Y.train,distrank,map$search.size)
      }else{
        # fit on all voxels
        fitted<-do.call(funname,list(X.train=X.train,Y.train=Y.train))
      }
      
      
      # estimation error
      results.e<-estimation.error(B,fitted$B.hat)
      
      # grab Beta maps.. for each model
      for (bb in seq(p)){
        emptymap<-matrix(NA,map$dim1*map$dim2)
        emptymap[map.train$e.idx]<-fitted$B.hat[bb,]
        maps.Bhat[,,bb,mm,subject]=matrix(emptymap,map$dim1,map$dim2)
      }
      #predict test set
      Yhat.test<-predict.me(fitted,X.test)
      
      #prediction error
      results.p<-prediction.error(Y.test,Yhat.test)
      
      emptymap<-matrix(NA,map$dim1*map$dim2)
      emptymap[map.train$e.idx]<-results.p$map.cor
      maps.pred.cor[,,mm,subject]=emptymap
      
      # for (rr in map$which.b){
      #    peaks.ts.Y[,rr,mm]=Y.test[,map$region.centers[rr]]
      #    peaks.ts.Yhat[,rr,mm]=Yhat.test[,map$region.centers[rr]]
      #  } 
      
      ### Store performance ####
      newrow = c(1:dim(df)[2])
      newrow[1]=mname
      newrow[2]=results.e$est.mse2
      newrow[3]=results.p$pred.mse
      newrow[4]=results.p$pred.cor
      newrow[5]=map$SNR
      newrow[6]=p
      newrow[7]=design$TRS
      df = rbind(df,newrow)
      mm=mm+1
      
    }
    
  }
  
  ## make numeric
  df<-df[2:dim(df)[1],]
  for (dd in seq(2,dim(df)[2])){
    df[,dd]<-as.numeric(df[,dd])
  }
  results=list(df=df,maps.Bhat=maps.Bhat,peaks.ts.Y=peaks.ts.Y,peaks.ts.Yhat=peaks.ts.Yhat)
  return(results)
}



##### Writing to matlab 

writeMat() 





