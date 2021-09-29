#Code for implementing BiMM tree: A method for a longitudinal/clustered decision tree for binary outcomes
#and BiMM forest: A method for a longitudinal/clustered random forest for binary outcomes
#Paper by: Jaime Speiser (email: jspeiser@wakehealth.edu), Beth Wolf, Dongjun Chung, Dean Karvellas, David Koch, Valerie Durkalski
#Journal article in Chemometrics and Intelligent Laboraty Systems

#notes: outcome must be binary with values of 0 and 1 (factor variable)

#required packages
library(blme)
library(rpart)
library(randomForest)

#function for BiMM tree model development
bimmTree<-function(formula,random,traindata,method="1iter",h1c=0.5,h2c=1.5,ErrorTolerance=0.006,MaxIterations=1000,verbose=TRUE){
  #rename the dataset
  data=traindata
  #parse formula
  Predictors<-paste(attr(terms(formula),"term.labels"),collapse="+")
  TargetName<-formula[[2]]
  Target<-data[,toString(TargetName)]
  initialRandomEffects=rep(0,length(data[,1]))
  minsize=round(length(data[,1])/10,0)
  tree.control=rpart.control(minbucket=minsize)
  #set up variables for loop
  ContinueCondition<-TRUE
  iterations<-0
  #initial values
  AdjustedTarget<-as.numeric(Target)-initialRandomEffects
  oldlik<- -Inf
  # Make a new data frame to include all the new variables
  newdata <- data
  
  while(ContinueCondition){
    # Current values of variables
    newdata[,"AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations+1
    #build tree
    tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),collapse = "~")),
                  data = newdata, method = "class", control = tree.control)
    if(verbose){
      print(paste(c("Iteration:",iterations)))
      print(tree)
    }
    ## Estimate New Random Effects and Errors using BLMER
    # Get variables that identify the node for each observation
    data[,"nodeInd"] <- 0
    data["nodeInd"] <- tree$where
    # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
    # Check that the fitted tree has at least two nodes.
    if(min(tree$where)==max(tree$where)){
      lmefit <- tryCatch(bglmer(formula(c(paste(paste(c(toString(TargetName),1), collapse="~"), "+(1|random)",sep=""))),
                                data=data,family=binomial,control=glmerControl(optCtrl=list(maxfun=20000))),error=function(cond)"skip")
    } else {
      lmefit <- tryCatch(bglmer(formula(c(paste(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~"), "+(1|random)",sep=""))),
                                data=data, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20000000000))),error=function(cond)"skip")
    }
    if(verbose){
      print(paste(c("Iteration:",iterations)))
      print(lmefit)
    }
    # Get the likelihood to check on convergence
    if(!(class(lmefit)[1]=="character")){
      newlik <- logLik(lmefit)
      ContinueCondition <- (abs(newlik-oldlik)>ErrorTolerance & iterations < MaxIterations)
      oldlik <- newlik
      # Extract random effects to make the new adjusted target
      logit<-predict(tree)[,1]
      logit2<-exp(predict(lmefit))/(1+exp(predict(lmefit)))
      AllEffects <- logit2
      #1 iteration, ignore adjusted target
      if(method=="1iter") {
        ContinueCondition<-FALSE
      }
      if(method=="h1"){
        AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) + AllEffects>h1c,1,0)
      }
      if(method=="h2"){
        AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) + AllEffects<h2c,0,1)
      }
      if(method=="h3"){
        for(k in 1:length(AllEffects)){
          if(as.numeric(Target[k])+AllEffects[k]<.5){AdjustedTarget[k]=0}
          else if(as.numeric(Target[k])+AllEffects[k]>1.5){AdjustedTarget[k]=1}
          else{
            #generate random probability coin flip based on AllEffects (q notation in paper)
            AdjustedTarget[k]<-rbinom(1,1,AllEffects[k])
          }
        }
      }
      #check to see if updated outcomes comprise 2 groups--if not, get out of loop
      if(levels(factor(AdjustedTarget))=="1" &method !="1iter"){
        ContinueCondition<-FALSE
        print("Error: updates are all for one group")
      }
    }
    if((class(lmefit)[1]=="character")){
      ContinueCondition<-FALSE
      print("Error: Bayesian GLMM did not converge")
    }
  }
  
  #return stuff
  return(list(Tree=tree, EffectModel=lmefit,Iterations=iterations,PostLogLike=logLik(lmefit),returndata=data))
  
}

#example
# tree1iter<-bimmTree(formula=comagradelow~Sex+Ethnicity+Age+ALT+AST+Bilirubin+Creat+Phosphate+Lactate+PH+platelets+ammonia+inr+pressors+rrt,random=traindata1$zSubjectCode,
#                     traindata=traindata1,method="1iter")
# #plot tree
# prp(tree1iter$Tree, type=4,extra=2,main="BiMM Tree 1 Iteration")
# #make predictions for new subjects
# preds<-predict(tree1iter$Tree,testdata,type="class")
# table(preds,testdata$comagradelow)



###########################################################################################################################
###########################################################################################################################
#function for BiMM forestmodel development

#note: bimm forest requires complete data (no missing predictors or outcome)

bimmForest<-function(formula,random,traindata,seed=8636,method="1iter",h1c=0.5,h2c=1.5,ErrorTolerance=0.5,MaxIterations=100,verbose=TRUE){
  #rename the dataset
  data=traindata
  #parse formula
  Predictors<-paste(attr(terms(formula),"term.labels"),collapse="+")
  TargetName<-formula[[2]]
  Target<-data[,toString(TargetName)]
  initialRandomEffects=rep(0,length(data[,1]))
  #set up variables for loop
  ContinueCondition<-TRUE
  iterations<-0
  #initial values
  AdjustedTarget<-as.numeric(Target)-initialRandomEffects
  oldlik<- -Inf
  # Make a new data frame to include all the new variables
  newdata <- data
  
  while(ContinueCondition){
    # Current values of variables
    newdata[,"AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations+1
    #build forest
    set.seed(seed)
    forest <- randomForest(formula(paste(c("factor(AdjustedTarget)", Predictors),collapse = "~")),
                           data = data, method = "class")
    forestprob<-predict(forest,type="prob")[,2]
    if(verbose){
      print(paste(c("Iteration:",iterations)))
      print(forest)
    }
    ## Estimate New Random Effects and Errors using BLMER
    options(warn=-1)
    lmefit <- tryCatch(bglmer(formula(c(paste(paste(c(toString(TargetName),"forestprob"), collapse="~"), "+(1|random)",sep=""))),
                              data=data,family=binomial,control=glmerControl(optCtrl=list(maxfun=20000))),error=function(cond)"skip")
    if(verbose){
      print(paste(c("Iteration:",iterations)))
      print(lmefit)
    }
    # Get the likelihood to check on convergence
    if(!(class(lmefit)[1]=="character")){
      newlik <- logLik(lmefit)
      ContinueCondition <- (abs(newlik-oldlik)>ErrorTolerance & iterations < MaxIterations)
      oldlik <- newlik
      # Extract random effects to make the new adjusted target
      logit<-forestprob
      logit2<-exp(predict(lmefit))/(1+exp(predict(lmefit)))
      AllEffects <- logit2
      #1 iteration, ignore adjusted target
      if(method=="1iter") {
        ContinueCondition<-FALSE
      }
      if(method=="h1"){
        AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) + AllEffects-1>h1c,1,0)
      }
      if(method=="h2"){
        AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) + AllEffects-1<h2c,0,1)
      }
      if(method=="h3"){
        for(k in 1:length(AllEffects)){
          if(as.numeric(Target[k])+AllEffects[k]-1<.5){AdjustedTarget[k]=0}
          else if(as.numeric(Target[k])+AllEffects[k]-1>1.5){AdjustedTarget[k]=1}
          else{
            #generate random probability coin flip based on AllEffects (q notation in paper)
            AdjustedTarget[k]<-rbinom(1,1,AllEffects[k])
          }
        }
      }
      #check to see if updated outcomes are the same, if so get out of loop
      if(min(AdjustedTarget)==max(AdjustedTarget)){
        ContinueCondition<-FALSE
        shouldpredict=FALSE
        print("Error: updates are all for one group")
      }
    }
    if((class(lmefit)[1]=="character")){
      ContinueCondition<-FALSE
      print("Error: Bayesian GLMM did not converge")
    }
  }
  
  #return stuff
  return(list(Forest=forest,EffectModel=lmefit,Iterations=iterations,PostLogLike=logLik(lmefit),returndata=data))
  
}

#example
#impute missing values
# set.seed(3829)
# traindata2<-rfImpute(factor(comagradelow)~Sex+Ethnicity+Age+ALT+AST+Bilirubin+Creat+Phosphate+Lactate+PH+platelets+ammonia+inr+pressors+rrt,
#                      data=traindata1)
# attach(traindata1)
# xdata<-cbind(traindata2,comagradelow,zSubjectCode)
# #develop bimm forest model
# forest1iter<-bimmForest(formula=comagradelow~Sex+Ethnicity+Age+ALT+AST+Bilirubin+Creat+Phosphate+Lactate+PH+platelets+ammonia+inr+pressors+rrt,random=xdata$zSubjectCode,traindata=xdata,method="1iter")
# #plot forest
# plot(forest1iter$Forest)
#variable importance plot
# varImpPlot(forest1iter$Forest)