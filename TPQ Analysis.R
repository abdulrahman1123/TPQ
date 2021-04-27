## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#                                                                                               #
#                                       TPQ Analysis                                            #
#                                                                                               #
#     This will be a formal document and a well-written one in terms of understandability       #
#                                                                                               #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #



####################################################################################
#                                 Load Libraries                                   #
####################################################################################
library(ggplot2)
library(extrafont)
library(reshape2)
library(car)
library(ppcor)
library(mlogit)
library(stringr)
library(dplyr)
library(ggthemes)
library(nlme)
library(mvoutlier)
library(devtools)
#library(rstatix)
library(lme4)
library(RGtk2)
library(MANOVA.RM)
library(mvnormtest)
library(ggpubr)
library(MVN)
library(heplots)
library(WRS)
library(rankMANOVA)
#library(ElemStatLearn)


####################################################################################
#                                 Load Functions                                   #
####################################################################################

windowsFonts(Times=windowsFont("Amiri"))
Distribute = function(data,distributed.var,within.var,within.val1,within.val2,common.var="Subject"){
  #distributed is the variable you want to distribute (Genotype or Response for example)
  #The within.var is the variable you want to distribute within (Session for example)
  #within.val1 is the value of the within.var in which there is data to be distributed (e.g. Retest) 
  #within.val2 is the value of the within.var you want to distribute to (e.g. Test)
  #Example, if you want to distribute the values of Response found in the retest session, into the test session, you enter:
  #common: is the variable you want to check its co-occurence between the two groups you distribute into (usually it is the subject or ID variable)
  #Distribute(Data,"Response","Session","Test","Retest","Subject"), and it reads "Distribute the Data of Response in the Test Session into the Retest Session, checking for Subject as a common variable"
  
  #Condition1: is the condition that the data are found in within.var1, Condition2 is abour within.var2
  Condtion1=data[[within.var]]==within.val1
  Condtion2=data[[within.var]]==within.val2
  
  for (i in 1:nrow(data[Condtion1,])){
    if (data[[common.var]][Condtion1][i] %in% data[[common.var]][Condtion2]){
      data[[distributed.var]][data[[common.var]]==data[[common.var]][Condtion1][i]] = data[[distributed.var]][Condtion1][i]
    }#else{
    #data[[distributed.var]][data[[common.var]]==data[[common.var]][Condtion1][i]] = NA
    #}
  }
  return (data)
}
MatchDelete <- function(data,variable,common="Subject"){
  data[[variable]]=as.factor(data[[variable]])
  data$inclusion=""
  Levels=levels(data[[variable]])
  
  #First I define the level of the factor to be matched, and the data for that factor
  #then I define the remaining levels, so that if an item is not found in //any// of the levels it will be deleted
  for (k in 1:length(Levels)){
    DataToMatch=data[[common]][data[[variable]]==Levels[k]]#This is for the subjects in the variable to matched in the specified level
    LevelToMatch = (1:length(Levels))[-k]# This is a variable containing the number of the "other" levels of the variable
    for (item in LevelToMatch){
      DataMatched=data[[common]][data[[variable]]==Levels[item]]#This is one of the levels to be matched to.
      for (i in 1:length(DataToMatch)){
        if (!(DataToMatch[i] %in% DataMatched)){
          data$inclusion[data[[common]]==DataToMatch[i]]=NA
        }
      }
    }
  }
  data=data[complete.cases(data$inclusion),]
  print ("PLEASE CHECK!"); print("#######################")
  for (item in Levels){
    print (paste("#Rows in",item, "is:",nrow(data[data[[variable]]==item,])))
  }
  data$inclusion=NULL
  return(data)
}
UMeans <- function(Data, DV, IV1, IV2=NA_character_,IV3=NA_character_,IV4=NA_character_,GroupBy=NA_character_){
  #remember: you just need to specify the data, and the variables, but note that the variables should be in character format
  #also, there should be no missing data (we might change that)
  #if any of the variables is NA, then make the switch off
  NameIV1=IV1; NameIV2=IV2; NameIV3=IV3; NameIV4=IV4; NameDV=DV
  IV1Switch="on"; IV2Switch="on"; IV3Switch="on"; IV4Switch="on"
  if (is.null(Data[[IV4]])){IV4Switch="off"}; if (is.null(Data[[IV3]])){IV3Switch="off"}; if (is.null(Data[[IV2]])){IV2Switch="off"}
  
  IV1=factor(Data[[IV1]])
  IV2=factor(Data[[IV2]])
  IV3=factor(Data[[IV3]])
  IV4=factor(Data[[IV4]])
  
  NLevel1 <<- nlevels(IV1)
  NLevel2 <<- nlevels(IV2)
  NLevel3 <<- nlevels(IV3)
  NLevel4 <<- nlevels(IV4)
  
  Factor1 <<- levels(IV1)
  Factor2 <<- levels(IV2)
  Factor3 <<- levels(IV3)
  Factor4 <<- levels(IV4)
  
  Mean <<- {}
  SEM <<- {}
  Count <<- {}
  #now NFactor1 means the nams of the factors in the first variable
  NFactor1 <<- {}
  NFactor2 <<- {}
  NFactor3 <<- {}
  NFactor4 <<- {}
  
  if (IV4Switch=="on"){
    #FactorCount is the number of times we should repeat the name of the factor in the final Means-data.frame
    #Some skills: we need to repeat factor 2 (for example) by length(factor3)*length(factor4), for each level, and then we repeat that
    #   by length(factor1)
    Factor1Count <<- length(Factor2)*length(Factor3)*length(Factor4)
    Factor2Count <<- length(Factor3)*length(Factor4)
    Factor3Count <<- length(Factor4)
    Factor4Count <<- 1
    
    NFactor1 <<- append(NFactor1,rep(Factor1,each=Factor1Count))
    NFactor2 <<- append(NFactor2,rep(rep(Factor2,each=Factor2Count),length(Factor1)))
    NFactor3 <<- append(NFactor3,rep(rep(Factor3,each=Factor3Count),length(Factor1)*length(Factor2)))
    NFactor4 <<- append(NFactor4,rep(rep(Factor4,each=Factor4Count),length(Factor1)*length(Factor2)*length(Factor3)))
    for (a in 1:NLevel1){
      for (b in 1:NLevel2){
        for (c in 1:NLevel3){
          for (d in 1:NLevel4){
            Mean <<- append(Mean,mean(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]&IV4==Factor4[d]]))
            SEM <<- append(SEM,sd(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]&IV4==Factor4[d]])/sqrt(length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]&IV4==Factor4[d]])))
            Count <<- append(Count,length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]&IV4==Factor4[d]]))
            
          }
        }
      }
    }
    Final <<- data.frame(Mean,SEM,NFactor1, NFactor2, NFactor3, NFactor4,Count)
    names(Final)=c(NameDV,"SEM",NameIV1,NameIV2,NameIV3,NameIV4,"Count")
    
    
  } else  if (IV3Switch=="on"){
    Factor1Count <<- length(Factor2)*length(Factor3)
    Factor2Count <<- length(Factor3)
    Factor3Count <<- 1
    
    NFactor1 <<- append(NFactor1,rep(Factor1,each=Factor1Count))
    NFactor2 <<- append(NFactor2,rep(rep(Factor2,each=Factor2Count),length(Factor1)))
    NFactor3 <<- append(NFactor3,rep(rep(Factor3,each=Factor3Count),length(Factor1)*length(Factor2)))
    for (a in 1:NLevel1){
      for (b in 1:NLevel2){
        for (c in 1:NLevel3){
          Mean <<- append(Mean,mean(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]]))
          SEM <<- append(SEM,sd(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]])/sqrt(length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]])))
          Count <<- append(Count,length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]&IV3==Factor3[c]]))
        }
      }
    }
    Final <<- data.frame(Mean,SEM,NFactor1, NFactor2, NFactor3,Count)
    names(Final)=c(NameDV,"SEM",NameIV1,NameIV2,NameIV3,"Count")
    
    
  } else  if (IV2Switch=="on"){
    Factor1Count <<- length(Factor2)
    Factor2Count <<- 1
    
    NFactor1 <<- append(NFactor1,rep(Factor1,each=Factor1Count))
    NFactor2 <<- append(NFactor2,rep(rep(Factor2,each=Factor2Count),length(Factor1)))
    for (a in 1:NLevel1){
      for (b in 1:NLevel2){
        Mean <<- append(Mean,mean(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]]))
        SEM <<- append(SEM,sd(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]])/sqrt(length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]])))
        Count <<- append(Count,length(Data[[DV]][IV1==Factor1[a]&IV2==Factor2[b]]))
      }
    }
    Final <<- data.frame(Mean,SEM,NFactor1, NFactor2,Count)
    names(Final)=c(NameDV,"SEM",NameIV1,NameIV2,"Count")
    
    
  } else {
    NFactor1 <<- Factor1
    for (a in 1:NLevel1){
      Mean <<- append(Mean,mean(Data[[DV]][IV1==Factor1[a]]))
      SEM <<- append(SEM,sd(Data[[DV]][IV1==Factor1[a]])/sqrt(length(Data[[DV]][IV1==Factor1[a]])))
      Count <<- append(Count,length(Data[[DV]][IV1==Factor1[a]]))
    }
    Final <<- data.frame(Mean,SEM,NFactor1,Count)
    names(Final)=c(NameDV,"SEM",NameIV1,"Count")
  }
  
  if (!is.na(GroupBy)){
    Final$Groups=paste0(Final[[GroupBy]]," (",Final$Count,")",sep="")
  }
  return (Final)
}
TypicalTheme=theme_bw(base_size = 16,base_family = "Times")+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                                                                  plot.subtitle = element_text(hjust = 0.5,face = "italic"))
FindAnova <- function(EZAnova,Transpose = TRUE){
  #note: the anova must be detailed (detailed=TRUE)
  AnovaModel=EZAnova$ANOVA
  
  NList=NULL #Names list
  DFList=NULL #DF List
  FList=NULL #F.value
  PList=NULL #P.value
  EList=NULL #Eta.value
  AnovaList=list()
  
  for (i in 2:nrow(AnovaModel)){
    NList=append(NList, AnovaModel[i, 1])
  }
  for (i in 2:nrow(AnovaModel)){
    DFValue=paste("(", AnovaModel[i, 2], ", ", AnovaModel[i, 3], ")", sep = "")
    DFList=append(DFList, DFValue)
  }
  for (i in 2:nrow(AnovaModel)){
    FList=append(FList, round(AnovaModel[i, 6], digits = 3))
  }
  for (i in 2:nrow(AnovaModel)){
    PValue=AnovaModel[i, 7]
    if (PValue<0.0001){
      PValue="<0.0001"
    }else{
      PValue=paste("=", round(AnovaModel[i, 7], digits = 3), sep="")
    }
    PList=append(PList, PValue)
  }
  for (i in 2:nrow(AnovaModel)){
    #Remember that the eta squared in the ezAnova is a generalized eta squared, we want the partial eta squared
    #which is the one that only takes into accoutn the SSerror of that variable and not all SSerror
    #So, we need to calculate it by ourselves
    #####EList=append(EList, round(AnovaModel[i, 9], digits = 3))
    PEtaSquared=AnovaModel[i, 4]/(AnovaModel[i, 4]+AnovaModel[i, 5])
    EList=append(EList, round(PEtaSquared, digits = 3))
  }
  AnovaList=c(DFList, FList, PList, EList)
  AnovaList=matrix(AnovaList, nrow = 4, byrow = TRUE)
  
  Eta="EtaSquared"
  if (length(PList)>1){
    Eta="Part.Eta.Sq."
  }
  rownames(AnovaList)=c("df", "F", "p.value", Eta)
  colnames(AnovaList)=NList
  
  if (Transpose){
    return(t(AnovaList))
  }else{
    return(AnovaList)
  }
}
LogisticInfo <- function(Model,PrintList=FALSE, PrintAny=TRUE){
  
  OriginalModel=Model
  Model=summary(Model)
  if (isTRUE(PrintAny)){
    print (Model)
  }
  
  Coefficients=round(Model$coefficients[1,1],digits = 2)
  Names=dimnames(Model$coefficients)[[1]][-1]
  
  for (i in 2:nrow(Model$coefficients)){
    Name = dimnames(Model$coefficients)[[1]][i]
    #if (str_sub(Name,start = 1,end = nchar(Name)/2)==str_sub(Name,start = 1+nchar(Name)/2)){
    #Name = str_sub(Name,start = 1,end = nchar(Name)/2)
    #}
    Zscore=Model$coefficients[i,3]
    Bvalue=Model$coefficients[i,1]
    chi.df=Model$df.null-Model$df.residual
    df=Model$df.null+1
    Deviance=Model$null.deviance
    R=sqrt((Zscore^2 - 2*chi.df)/Deviance)
    
    Coefficients=append(Coefficients,round(Bvalue,digits = 2))
    
    ModelChi <- Model$null.deviance - Model$deviance
    pvalue.chi=1-pchisq(ModelChi,chi.df)
    print (c(ModelChi,chi.df))
    pvalue.model=Model$coefficients[i,4]
    R2.hl = ModelChi/Model$null.deviance
    R2.cs = 1-exp(-1*ModelChi/df)
    R2.n = R2.cs / (1-(exp(-(Model$null.deviance/df))))
    odds.ratio=exp(Model$coefficients)[i,1]
    
    if (isTRUE(PrintAny)){
      print (paste("Results for the variable",Name))
      print (paste("z-score =",round(Zscore,digits = 4)))
      print (paste("B Value =",round(Bvalue,digits = 4)))
      print ("")
      print (paste("Chi square: X^2 (",chi.df,") =  ",round(ModelChi,digits=4),", p-value = ",round(pvalue.chi,digits = 4),sep = ""))
      print (paste("Model p-value =",round(pvalue.model,digits = 4)))
      print ("")
      print (paste("R =",round(R,digits = 4)))
      print (paste("Hosmer and Lemeshow’s R^2 =",round(R2.hl,digits = 4)))
      print (paste("Cox and Snell's R^2 =",round(R2.cs,digits = 4)))
      print (paste("Nagelkerke’s R^2 =",round(R2.n,digits = 4)))
      print ("")
      print (paste("Odds Ratio =",round(odds.ratio,digits = 4)))
      print ("")
    }
  }
  
  Xs=paste(as.character(Coefficients[1]))
  for (i in 2:length(Coefficients)){
    Sign=ifelse((Coefficients[i]>0|Coefficients[i]==0)," + "," - ")
    Xs=paste(Xs,Sign,abs(Coefficients[i])," X",i-1,sep = "")
  }
  
  OriginalData=OriginalModel$data[Names]
  
  Outcome=as.character(OriginalModel$formula[[2]])
  
  OriginalData[Outcome]=OriginalModel$data[Outcome]
  
  PredictedExp=rep(Coefficients[1],nrow(OriginalData))
  
  for (i in length(OriginalData)-1){
    PredictedExp=PredictedExp+OriginalData[Names[i]]*Coefficients[i+1]
  }
  PredictedProb=(1/(1+exp(-1*PredictedExp)))
  PredictedResult=ifelse(PredictedProb<0.5,0,1)
  
  OriginalData["PredictedResult"]=PredictedResult
  OriginalData["Prob"]=PredictedProb
  if (isTRUE(PrintAny)){
    print ("CONFIDENCE INTERVAL OF ODDS RATIO:")
    print (exp(confint(OriginalModel)))
    print (paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,"))",sep = ""))
  }
  
  
  if (isTRUE(PrintList)){
    print (OriginalData)
  }
  PredMatrix=table(OriginalData[,length(OriginalData)-1],OriginalData[,(length(OriginalData)-2)])
  colnames(PredMatrix)=c("0(Predicted)","1(Predicted)"); rownames(PredMatrix)=c("0(Actual)","1(Actual)")
  print (paste("Sensitivity =",round(PredMatrix[4]/(PredMatrix[3]+PredMatrix[4]),digits = 3)*100,"%"))
  print (paste("Specificity =",round(PredMatrix[1]/(PredMatrix[1]+PredMatrix[2]),digits = 3)*100,"%"))
  print (paste("Neg.Pred.Value =",round(PredMatrix[1]/(PredMatrix[1]+PredMatrix[3]),digits = 3)*100,"%"))
  print (paste("Pos.Pred.Value =",round(PredMatrix[4]/(PredMatrix[2]+PredMatrix[4]),digits = 3)*100,"%"))
  return (PredMatrix)
  
}
LogisticFunction = function(Model, Threshold = 0.5, plt_type = "histogram"){
  #The Model is a logistic model (DV~IV1+IV2..., family = binomial())
  Data=Model$data
  Summary=summary(Model)
  n.value=Model$df.null+1
  Coefficients = Summary$coefficients
  
  #The improvement we got (change in deviance) after including predictors (baseline deviance - )
  modelChi = Model$null.deviance-Model$deviance 
  
  #degress of freedom, which are the df for the constant-only model minus the df in the predictors-model
  chidf = Model$df.null - Model$df.residual
  
  #since it has a chisquare distribution, we can calculate significance
  chisq.prob=1-pchisq(modelChi, chidf)
  
  #r.value=sqrt((z.value^2 - 2*chidf)/Model$null.deviance)
  
  R2.hl = modelChi/Model$null.deviance
  R2.cs = 1- exp(-modelChi/n.value)
  R2.n = R2.cs/(1-exp(-(Model$null.deviance/n.value)))
  
  #Odds ratio, upper and lower CI
  odds_ratio=round(exp(Model$coefficients), digits = 3)
  Lower.CI=round(exp(confint(Model))[, 1], digits = 3)
  Upper.CI=round(exp(confint(Model))[, 2], digits = 3)
  
  #Beta and SE
  Beta=round(Coefficients[, 1], digits = 3)
  SE=round(Coefficients[, 2], digits = 3)
  
  #zvalues and p-values
  Zvalues=round(Coefficients[, 3], digits = 3)
  Pvalues=round(Coefficients[, 4], digits = 3)
  
  PredictedVar = Data[[as.character(Model$formula[[2]])]]
  PredictedVarName = as.character(Model$formula[[2]])
  Predictors=as.character(Model$formula[3])
  
  #Sometimes, there are no levels, I will create them
  if (is.null(levels(PredictedVar))){
    PredictedVar = as.factor(PredictedVar)
  }
  BaseLevel=levels(PredictedVar)[1]
  FirstLevel=levels(PredictedVar)[2]
  Data$predicted.probability <-fitted(Model)
  Data$predicted.outcome <-ifelse(fitted(Model)<Threshold, BaseLevel, FirstLevel)
  
  TruePositive=nrow(Data[(PredictedVar==FirstLevel&Data$predicted.outcome==FirstLevel), ])
  FalsePositive=nrow(Data[(PredictedVar==BaseLevel&Data$predicted.outcome==FirstLevel), ])
  TrueNegative=nrow(Data[(PredictedVar==BaseLevel&Data$predicted.outcome==BaseLevel), ])
  FalseNegative=nrow(Data[(PredictedVar==FirstLevel&Data$predicted.outcome==BaseLevel), ])
  
  Sensitivity=round(TruePositive/(TruePositive+FalseNegative), digits = 3)
  Specificity=round(TrueNegative/(TrueNegative+FalsePositive), digits = 3)
  PPV=round(TruePositive/(TruePositive+FalsePositive), digits=3)
  NPV=round(TrueNegative/(TrueNegative+FalseNegative), digits=3)
  
  Xs=round(Coefficients[1,1],digits = 4)
  for (i in 2:length(Coefficients[,1])){
    VAr_name = names(Coefficients[,1])[i]
    factor_value = gsub("-"," - ",round(Coefficients[i,1],digits = 4))
    Xs = paste0(Xs,factor_value,"*",VAr_name)
  }
  #The y.value finder, which is the maximum count on the plot divided by 2, this is for plotting reasons
  TList=NULL
  for (i in 0:20){
    interval=c((i-0.5)*0.05, (i+0.5)*0.05)
    Num=sum((Data$predicted.probability>interval[1]&Data$predicted.probability<interval[2]))
    TList=append(TList, Num)
  }
  y.value=max(TList)/2
  
  #Data $ standardized.residuals <-rstandard(Model)
  #Data $ studentized.residuals <-rstudent(Model)
  #Data $ dfbeta <-dfbeta(Model)
  #Data $ dffit <-dffits(Model)
  #Data $ leverage <-hatvalues(Model)
  #BaseLevel
  #FirstLevel
  
  GText_0=paste("Specificity =", (100*Specificity), "%  || ", "NPV =", (100*NPV), "%")
  GText_1=paste("Sensitivity =", (100*Sensitivity), "%  || ", "PPV =", (100*PPV), "%")
  
  Title = paste("Logistic Regression Function for", as.character(Model$formula[[2]]))
  Subtitle= paste("As Predicted by", Predictors)
  if (plt_type == "histogram"){
    Drawing=ggplot(Data, aes(x=predicted.probability, fill=PredictedVar))+geom_histogram(binwidth=0.02, color="black")+
      scale_fill_manual(name="Group", values=c("#08457E", "#FBEC5D"))+TypicalTheme+geom_vline(xintercept = Threshold)+
      scale_x_continuous("Predicted Probability",limits = c(-0.1, 1.1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+
      annotate(geom = "text", size=5, family="Amiri", x = (Threshold/2), y = y.value*2.1, label = paste("Predicted to be", BaseLevel))+
      annotate(geom = "text", size=5, family="Amiri", x = ((1+Threshold)/2), y = y.value*2.1, label = paste("Predicted to be", FirstLevel))+
      annotate(geom = "text", size=5, family="Amiri", x = 0, y = y.value, label = GText_0, angle = 90)+
      annotate(geom = "text", size=5, family="Amiri", x = 1, y = y.value, label = GText_1, angle = 90)+
      ggtitle(Title, subtitle = Subtitle)+scale_y_continuous("Number of Cases")
  } else if (plt_type == "glm"){
    if (length(strsplit(Predictors,"+",fixed = TRUE)[[1]])==1){
      Data$Values = Data[[Predictors]]
      PredictedLevels = levels(Data[[PredictedVarName]])
      Data$PredictedFactor = as.character(Data[[PredictedVarName]])
      Data$PredictedFactor[Data$PredictedFactor==PredictedLevels[1]]=0
      Data$PredictedFactor[Data$PredictedFactor==PredictedLevels[2]]=1
      Data$PredictedFactor = as.numeric(Data$PredictedFactor)
      Drawing =ggplot(data = Data,mapping= aes(x = Values, y=PredictedFactor))+
        geom_jitter(shape = 21, color = "#7c0a02", fill="white",size = 4,width = 0.1, height = 0)+
        stat_smooth(method = "glm", method.args = list(family = "binomial"),se=F,color="#7c0a02",size=4)+
        scale_y_continuous(name = PredictedVarName,breaks=c(0,1), labels= c(levels(Data[[PredictedVarName]])[1],levels(Data[[PredictedVarName]])[2]))+
        scale_x_continuous(name = Predictors)+
        ggtitle(Title,subtitle = Subtitle)+
        TypicalTheme
    }else{
      Drawing = ("Can't draw glm smoothed figure when there is more than one predictor")
    }
  }
  
  #ggplot(Data, aes(x=predicted.probability, fill=Cured))+geom_bar(color="black")+scale_fill_manual(values=c("#08457E", "#FBEC5D"))+
  #scale_x_continuous(limits = c(-0.1, 1.1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+TypicalTheme+geom_vline(xintercept = 0.5)+
  #annotate(geom = "text", size=5, family="Amiri", x = -0.1, y = 10, label = GText_0, angle = 90)+
  #annotate(geom = "text", size=5, family="Amiri", x = 1.1, y = 10, label = GText_1, angle = 90)
  
  #Create a dataframe that contains only the needed data
  data_frame_essential = data.frame(predicted.probability= Data$predicted.probability,
                                    predicted.outcome = Data$predicted.outcome,
                                    actual.outcome = Data[[PredictedVarName]])
  
  for (i in 2:length(as.character(Model_IC_All$formula[[3]]))){
    local_factor= as.character(Model_IC_All$formula[[3]])[i]
    data_frame_essential[[local_factor]]=Data[[local_factor]]
  }
  
  #Create prediction matrix
  PredMatrix=table(data_frame_essential$predicted.outcome,data_frame_essential$actual.outcome)
  rownames(PredMatrix)=c(paste0(BaseLevel,"(Predicted)"),paste0(FirstLevel,"(Predicted)")); colnames(PredMatrix)=c(paste0(BaseLevel,"(Actual)"),paste0(FirstLevel,"(Actual)"))
  
  LogValues=data.frame(Beta, SE, Lower.CI, odds_ratio, Upper.CI, Zvalues, Pvalues)
  DerivedValues=data.frame(chisq=modelChi, df=chidf, p.value.chi=chisq.prob, r2.hl=R2.hl, r2.cs=R2.cs, r2.n=R2.n, sensitivity=Sensitivity, specificity=Specificity)
  DataList=list(data_frame = data_frame_essential,
                log.values=LogValues, derived.values=DerivedValues,
                prob.formula = paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,")",sep = ""),
                PredMatrix = PredMatrix)
  
  print (Drawing)
  
  return (DataList)
}
SMeans <- function(Data, DV, IVs,GroupBy=NA_character_){
  NamesList=list(IVs)
  FactorList={}
  LevelsCount={}
  for (item in IVs){
    AddedLevels=levels(factor(Data[[item]]))
    NLevels=nlevels(factor(Data[[item]]))
    
    FactorList = append(FactorList,list(AddedLevels))
    LevelsCount = append(LevelsCount, NLevels)
  }
  FactorFrame=data.frame(Dummy=rep(NA,prod(LevelsCount)))
  FactorNumber=length(LevelsCount)
  
  for (i in 1:length(IVs)){
    EachRepetitionAmount = prod(LevelsCount[i:FactorNumber])/LevelsCount[i]
    GeneralRepetitionAmount = prod(LevelsCount[1:i])/LevelsCount[i]
    FactorFrame[[IVs[i]]] = rep(rep(FactorList[[i]],each=EachRepetitionAmount),GeneralRepetitionAmount)
  }
  FactorFrame$Dummy=NULL
  
  Means={}
  SEM={}
  Count={}
  for (i in 1:prod(LevelsCount)){
    Condition=rep(TRUE, nrow(Data))
    for (j in 1:ncol(FactorFrame)){
      Condition=Condition&(as.character(Data[[IVs[j]]])==FactorFrame[[j]][i])
    }
    
    Means = append(Means,mean(Data[[DV]][Condition],na.rm = TRUE))
    SEM = append(SEM,sd(Data[[DV]][Condition],na.rm = TRUE)/sqrt(length(Data[[DV]][Condition])))
    Count = append(Count,length(Data[[DV]][Condition]))
  } 
  
  FinalData=data.frame(Dummy=1:length(Means))
  FinalData[[DV]]=Means
  FinalData$SEM=SEM
  FinalData$Count=Count
  FinalData=cbind(FinalData,FactorFrame)
  FinalData$Dummy=NULL
  
  if (!is.na(GroupBy)){
    FinalData$Groups=paste0(FinalData[[GroupBy]]," (",FinalData$Count,")",sep="")
  }
  return (FinalData)
}
Attach.H <- function(DataFrame1,DataFrame2,SubjectID1="Subject",SubjectID2="Subject"){
  Length=nrow(DataFrame1)
  AvoidedCols=colnames(DataFrame1)#This line is to avoid the repeated columns (e.g. Subject, or diagnosis ...)
  ColumnNames=colnames(DataFrame2)[!(colnames(DataFrame2) %in% AvoidedCols)]
  #The following line creates a data frame with columns with similar names as the 2nd data frame's names, with NA values
  #The number of NA values must be the same as dataframe 1
  CombiningDF=NULL
  for (i in 1:length(ColumnNames)){CombiningDF=cbind(CombiningDF,rep(NA,Length))}; colnames(CombiningDF)=ColumnNames
  CombiningDF=data.frame(CombiningDF,stringsAsFactors = FALSE)
  ResultingDataFrame=DataFrame1
  
  #for (ColN in ColumnNames){
  for (Sub in DataFrame2[[SubjectID2]]){
    if (Sub %in% DataFrame1[[SubjectID1]]){
      for (ColName in (ColumnNames)){
        CombiningDF[[ColName]][DataFrame1[[SubjectID1]]==Sub] = as.character(DataFrame2[[ColName]][DataFrame2[[SubjectID2]]==Sub])
      }
    }
  }
  for (ColName in ColumnNames){
    CombiningDF[[ColName]]=ReturnClass(DataFrame2[[ColName]],CombiningDF[[ColName]])
  }
  ResultingDataFrame=cbind(ResultingDataFrame,CombiningDF)
  return(ResultingDataFrame)
}
ReturnClass <- function(Orig.Var,New.Var){
  #This function will be used inside the other functions
  if (class(Orig.Var)=="factor"){
    New.Var=factor(New.Var)
    #print ("Coercing back to Factor")
  }else if(class(Orig.Var)=="numeric"|class(Orig.Var)=="integer"){
    #print ("Coercing back to Numeric")
    New.Var=as.numeric(New.Var)
  }else if(class(Orig.Var)=="character"){
    #print ("Coercing back to Character")
    New.Var=as.character(New.Var)
  }
  return(New.Var)
}

#Load the Data and keep those with test and retest
####################################################################################
#                                    Load Data                                     #
####################################################################################
if (.Platform$OS.type == "unix"){
  TPQData=read.csv("/home/abdulrahman/Dropbox/TPQ_Paper/TPQ_Analysis_GAD.csv")
}else{
  TPQData=read.csv("file:///C:/Users/jsawa/Dropbox/TPQ_Paper/TPQ_Analysis_GAD.csv")
}
TPQData=MatchDelete(TPQData,"Session")

#Calculate Date of Retest
TPQData$DateRetest=TPQData$TestingDay
TPQData=Distribute(TPQData,"DateRetest","Session","Retest","Test")
TPQData$RetestPeriod=as.numeric(TPQData$DateRetest)-as.numeric(TPQData$TestingDay)


#Distribute GAD and Response States
TPQData=Distribute(TPQData,"Response","Session","Retest","Test")
TPQData=Distribute(TPQData,"GAD","Session","Test","Retest")


#Remove those with incomplete TPQ, diagnosis, response or retest period
TPQData=TPQData[(TPQData$Session=="Retest"|(TPQData$RetestPeriod>27 & TPQData$RetestPeriod<76)),]
TPQData=TPQData[complete.cases(TPQData$GAD=="NoGad"|TPQData$GAD=="GAD"|TPQData$Session=="Retest"|TPQData$Diagnosis=="HC"),]
TPQData=TPQData[(TPQData$Response=="Non-responder"|TPQData$Response=="Responder"|TPQData$Response=="NoResponse"|TPQData$Session=="Retest"|TPQData$Diagnosis=="HC"),]

#Remove those who were on medications
#TPQData=TPQData[((as.character(TPQData$SSRI)=="None" & TPQData$Session=="Test" & TPQData$Diagnosis=="MDD") | TPQData$Session=="Retest" | TPQData$Diagnosis=="HC"),]

#Remove empty TPQs (sometimes they appear as all zeros)
TPQData$NS[(TPQData$NS+TPQData$HA+TPQData$RD)==0]=NA
TPQData=TPQData[(complete.cases(TPQData$NS)),]
TPQData=MatchDelete(TPQData,"Session")

#Order the data according to session, response, GAD and Subject, in order. This is important for the SPSS form of data
TPQData = TPQData[order(TPQData$Session,TPQData$Response,TPQData$GAD,TPQData$Subject),]
TPQData$Subject[TPQData$Session=="Test"] == TPQData$Subject[TPQData$Session=="Retest"]



#Make Age and Education numeric
TPQData$Age=as.numeric(as.character(TPQData$Age))
TPQData$Years.of.Education=as.numeric(as.character(TPQData$Years.of.Education))
TPQData$WAIS.R=as.numeric(as.character(TPQData$WAIS.R))

#Make response a factor
TPQData$Response = factor(TPQData$Response)
TPQData$Subject = factor(TPQData$Subject)



####################################################################################
#                           Data Preparation for  MANOVA                           #
####################################################################################

TPQData_M = data.frame(Subject = factor(TPQData$Subject),ID = TPQData$ID, HA = TPQData$HA, NS = TPQData$NS, RD = TPQData$RD,
                       Session = TPQData$Session, Response = TPQData$Response, GAD = TPQData$GAD,
                       Diagnosis = factor(TPQData$Diagnosis))
TPQData_M$GAD=as.character(TPQData_M$GAD)
TPQData_M$Response=as.character(TPQData_M$Response)

TPQData_M$GAD[TPQData_M$GAD=="HC"]="NoGad"
TPQData_M$GAD=factor(TPQData_M$GAD)

TPQData_M$Response[(TPQData_M$Response=="HC")]="NoResponse"
TPQData_M$Response=factor(TPQData_M$Response)

out_plot = aq.plot(TPQData_M[,3:5])
TPQData_Out = TPQData_M[!out_plot$outliers,]
TPQData_Out=MatchDelete(TPQData_Out,"Session")

TPQData_Out$Subject = factor(TPQData_Out$Subject)

TPQData_Out = TPQData_Out[order(TPQData_Out$Session,TPQData_Out$Response,TPQData_Out$GAD,TPQData_Out$Subject),]
#nrow(TPQData_Out[TPQData_Out$Response=="Non-responder",])
#First, make the SPSS data-equivalent for comparison reasons
SPSSTPQ = data.frame(Subject = TPQData_Out$Subject[TPQData_Out$Session=="Test"],
                     Response = TPQData_Out$Response[TPQData_Out$Session=="Test"],
                     GAD = TPQData_Out$GAD[TPQData_Out$Session=="Test"],
                     HATest = TPQData_Out$HA[TPQData_Out$Session=="Test"],
                     RDTest = TPQData_Out$RD[TPQData_Out$Session=="Test"],
                     NSTest = TPQData_Out$NS[TPQData_Out$Session=="Test"],
                     HARetest = TPQData_Out$HA[TPQData_Out$Session=="Retest"],
                     RDRetest = TPQData_Out$RD[TPQData_Out$Session=="Retest"],
                     NSRetest = TPQData_Out$NS[TPQData_Out$Session=="Retest"])
if (!.Platform$OS.type == "unix"){
  write.table(SPSSTPQ,"file:///C:/Users/jsawa/Desktop/TPQSPSS_GAD.xls", sep="\t",row.names = FALSE)
}

####################################################################################
#                                       MANOVA                                     #
####################################################################################
#Lets try the normal MANOVA using R
# MANOVA test
TPQData_Out$Response = relevel(TPQData_Out$Response, ref = "NoResponse")
TPQData_Out$GAD = relevel(TPQData_Out$GAD, ref = "NoGad")

#Planned contrasts
contrasts(TPQData_Out$Response) = cbind(c(2,-1,-1),c(0,1,-1)) #1: NoResponse compared to all others & 2: responders vs nonresponders
contrasts(TPQData_Out$GAD) = c(1,-1)
TPQScales = cbind(TPQData_Out$HA, TPQData_Out$RD, TPQData_Out$NS)
TPQManova <- manova(TPQScales ~ TPQData_Out$Response*TPQData_Out$GAD*TPQData_Out$Session+Error(TPQData_Out$Subject/TPQData_Out$Session))
summary(TPQManova,type="III")


#Assumption Check, based on "Understanding Statistics Using R", and "https://www.datanovia.com/en/lessons/one-way-manova-in-r/#check-sample-size-assumption",
# I will put anything in the book as critical and if the assumption is only in the website, it will be marked as NotCritical
#NotCritical: Sample size check .. no need, I know that subject numbers are more than factors, but anyway
summary(TPQData_Out)
TPQData_Summary= TPQData_Out %>% group_by(Session,GAD,Response) %>% summarise(N = n()) ; TPQData_Summary

#Critical: Outliers .. 5 of them  ,,, make sure to repeat the MANOVA with and without them
#Outliers were removed in the steps above
outliercount = sum(out_plot$outliers[TPQData_Out$Session=="Test"]); print (outliercount)
TPQData_Out %>% group_by(Session,Response,GAD) %>% summarise(N = n())

#Univariate normality: All normal except  f"NS : Non-responder , GAD , Test  =  0.041"  &  "NS : Responder , GAD , Retest  =  0.011"
ggqqplot(TPQData_Out[TPQData_Out$GAD=="GAD",],x= "NS", facet.by = "Session*Response",
         ylab = "NS", ggtheme = theme_bw(), title = "QQ plot for NS only in patients with GAD")

ResponseL = c("NoResponse","Responder","Non-responder")
GADL = c("GAD","NoGad")
SessionL = c("Test","Retest")
DV = c("HA","NS","RD")
for (i in 1:3){
  for (k in 1:2){
    for (l in 1:2){
      for (n in 1:3){
        if (!(GADL[k]=="GAD"&ResponseL[i]=="HC")){
          Data = TPQData_Out[[DV[n]]][(TPQData_Out$Response == ResponseL[i] & TPQData_Out$GAD == GADL[k]& TPQData_Out$Session == SessionL[l])]
          p_value = shapiro.test(Data)$p.value
          print (paste(DV[n],":",ResponseL[i],",",GADL[k],",",SessionL[l]," = ",round(p_value,3)))
        }
        
      }
    }
  }
}


#Critical: multivariate normality
Test_NoRes_NoGAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD == "NoGad"),2:4]
Test_Res_NoGAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="Responder"&TPQData_Out$GAD == "NoGad"),2:4]
Test_NR_NoGAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD == "NoGad"),2:4]
Test_NoRes_GAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD == "GAD"),2:4]
Test_Res_GAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="Responder"&TPQData_Out$GAD == "GAD"),2:4]
Test_NR_GAD = TPQData_Out[(TPQData_Out$Session=="Test" & TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD == "GAD"),2:4]
Retest_NoRes_NoGAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD == "NoGad"),2:4]
Retest_Res_NoGAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="Responder"&TPQData_Out$GAD == "NoGad"),2:4]
Retest_NR_NoGAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD == "NoGad"),2:4]
Retest_NoRes_GAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD == "GAD"),2:4]
Retest_Res_GAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="Responder"&TPQData_Out$GAD == "GAD"),2:4]
Retest_NR_GAD = TPQData_Out[(TPQData_Out$Session=="Retest" & TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD == "GAD"),2:4]

mshapiro.test(t(Test_NoRes_NoGAD))
mshapiro.test(t(Test_Res_NoGAD))
mshapiro.test(t(Test_NR_NoGAD))
mshapiro.test(t(Test_NoRes_GAD))   #*
mshapiro.test(t(Test_Res_GAD))
mshapiro.test(t(Test_NR_GAD))      #*
mshapiro.test(t(Retest_NoRes_NoGAD))
mshapiro.test(t(Retest_Res_NoGAD))
mshapiro.test(t(Retest_NR_NoGAD))  #*
mshapiro.test(t(Retest_NoRes_GAD)) #*
mshapiro.test(t(Retest_Res_GAD))   #*
mshapiro.test(t(Retest_NR_GAD))    #*

#skewness
mvn(Test_NoRes_NoGAD)
mvn(Test_Res_NoGAD)
mvn(Test_NR_NoGAD)
mvn(Test_NoRes_GAD)
mvn(Test_Res_GAD)
mvn(Test_NR_GAD)      
mvn(Retest_NoRes_NoGAD)
mvn(Retest_Res_NoGAD) 
mvn(Retest_NR_NoGAD)  
mvn(Retest_NoRes_GAD)
mvn(Retest_Res_GAD)   #*
mvn(Retest_NR_GAD)    #*


mvn(Test_NR_GAD, multivariatePlot = "qq")
mvn(Retest_Res_NoGAD, multivariatePlot = "qq")
mvn(Retest_NR_NoGAD, multivariatePlot = "qq")
mvn(Retest_Res_GAD, multivariatePlot = "qq")
mvn(Retest_NR_GAD, multivariatePlot = "qq")

TPQData_Out_2 = TPQData_Out
TPQData_Out_2$FactorSum = factor(paste(TPQData_Out$Response,TPQData_Out$GAD))
TPQData_Out_2$Response=NULL
TPQData_Out_2$GAD=NULL
boxM(TPQData_Out_2[,2:4],TPQData_Out_2[,6],cov=TPQData_Out_2[,5])

#NotCritical: check correlation (you want them to be correlated!)
cor(TPQData_M[,2:4], use = "pairwise.complete.obs") #they are weakly correlated

#Critical: variance-covariance matrix
VC_Mat = by(TPQData_Out[,2:4], list(TPQData_Out$GAD,TPQData_Out$Response,TPQData_Out$Session), cov)
TPQData_Summary[4]
Retest_NR_GAD_B = data.frame(VC_Mat[1][[1]],count = 10)
Retest_NR_NoGAD_B = data.frame(VC_Mat[2][[1]],count = 12)
Retest_NoRes_GAD_B = data.frame(VC_Mat[3][[1]],count = 6)
Retest_NoRes_NoGAD_B = data.frame(VC_Mat[4][[1]],count = 51)
Retest_Res_GAD_B = data.frame(VC_Mat[5][[1]],count = 19)
Retest_Res_NoGAD_B = data.frame(VC_Mat[6][[1]],count = 21)
Test_NR_GAD_B = data.frame(VC_Mat[7][[1]],count = 10)
Test_NR_NoGAD_B = data.frame(VC_Mat[8][[1]],count = 12)
Test_NoRes_GAD_B = data.frame(VC_Mat[9][[1]],count = 6)
Test_NoRes_NoGAD_B = data.frame(VC_Mat[10][[1]],count = 51)
Test_Res_GAD_B = data.frame(VC_Mat[11][[1]],count = 19)
Test_Res_NoGAD_B = data.frame(VC_Mat[12][[1]],count = 21)


#Typically, you should look at variances and covariances and if the larger samples are giving the larger variances and covariances
# then the results of the MANOVA are robust and can be trusted if significant. If the smaller ones are the 
# ones giving larger variances and covariances, then the significant results should be taken with a grain of salt
# I will go further and look into the variances and covariances and correlate the count with the variance 
CorList = list(Retest_NoRes_GAD_B,Retest_NoRes_NoGAD_B, Retest_NR_GAD_B,Retest_NR_NoGAD_B,Retest_Res_GAD_B,Retest_Res_NoGAD_B, Test_NoRes_GAD_B,Test_NoRes_NoGAD_B, Test_NR_GAD_B,Test_NR_NoGAD_B,Test_Res_GAD_B,Test_Res_NoGAD_B)
c(data.frame(S=c(1,2,3),B=c(1,2,33)))
HA_C=NULL
NS_C=NULL
RD_C=NULL
HA_NS_C=NULL
HA_RD_C=NULL
NS_RD_C=NULL
Counts=NULL
for (item in CorList){
  HA_C = append(HA_C,item$HA[1])
  HA_NS_C = append(HA_NS_C,item$HA[2])
  HA_RD_C = append(HA_RD_C,item$HA[3])
  NS_C = append(NS_C,item$NS[2])
  NS_RD_C = append(NS_RD_C,item$NS[3])
  RD_C = append(RD_C,item$RD[3])
  Counts=append(Counts, item$count[1])
}

Var_Cov = data.frame(HA=HA_C,NS = NS_C,RD = RD_C,HA_NS = HA_NS_C,HA_RD = HA_RD_C,NS_RD = NS_RD_C,SampleSize=Counts)
Var_Cov = melt(Var_Cov,id= "SampleSize",variable.name = "Variable",value.name = "Variance")
ggplot(Var_Cov, aes(x= SampleSize, y=Variance,color = Variable,group = Variable))+
  geom_point()+facet_wrap(~Variable)+geom_smooth(method = "lm")+TypicalTheme

Correlate(Var_Cov,"SampleSize","Variance",Factor = "Variable")





####################################################################################
#                                   Robust MANOVA                                  #
####################################################################################

TPQData_Out$Response = relevel(TPQData_Out$Response, ref = "NoResponse")
TPQData_Out$GAD = relevel(TPQData_Out$GAD, ref = "NoGad")

testWildBS_1 = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
testWildBS_2 = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
testWildBS_3 = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
summary(testWildBS_1)
summary(testWildBS_2)
summary(testWildBS_3)

#With contrasts
contrasts(TPQData_Out$Response) = cbind(c(2,-1,-1),c(0,1,-1)) #1: NoResponse compared to all others & 2: responders vs nonresponders
contrasts(TPQData_Out$GAD) = c(1,-1)
testWildBS_1_c = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
testWildBS_2_c = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
testWildBS_3_c = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 100,dec = 4, CPU = 8, resampling = "WildBS")
summary(testWildBS_1_c)
summary(testWildBS_2_c)
summary(testWildBS_3_c)

testWildBS_Final = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 10000,dec = 6, CPU = 8, resampling = "WildBS")
summary(testWildBS_Final)


testBoot = rankMANOVA(cbind(NS,HA,RD) ~ Session*GAD*Response,data = TPQData_Out, iter = 10000,dec = 6, CPU = 8,resampling = "bootstrap")
summary(testBoot)

#For comparison to SPSS, I will do everything without the session
testWildBS_Test = rankMANOVA(cbind(NS,HA,RD) ~ GAD*Response,data = NewTPQData_Out[NewTPQData_Out$Session=="Test",], iter = 10000,dec = 6, CPU = 8, resampling = "WildBS")
testWildBS_Retest = rankMANOVA(cbind(NS,HA,RD) ~ GAD*Response,data = NewTPQData_Out[NewTPQData_Out$Session=="Retest",], iter = 10000,dec = 6, CPU = 8, resampling = "WildBS")
summary(testWildBS_Test)
summary(testWildBS_Retest)

testBoot_Test = rankMANOVA(cbind(NS,HA,RD) ~ GAD*Response,data = NewTPQData_Out[NewTPQData_Out$Session=="Test",], iter = 10000,dec = 6, CPU = 8, resampling = "bootstrap")
testBoot_Retest = rankMANOVA(cbind(NS,HA,RD) ~ GAD*Response,data = NewTPQData_Out[NewTPQData_Out$Session=="Retest",], iter = 10000,dec = 6, CPU = 8, resampling = "bootstrap")
summary(testBoot_Test)
summary(testBoot_Retest)


####################################################################################
#                    Assumption check for Follow up ANOVAs                         #
####################################################################################

#homogeniety of variance for each DV ... GOOD
leveneTest(HA ~ Response*GAD, data = TPQData_Out[TPQData_Out$Session=="Test",])
leveneTest(NS ~ Response*GAD, data = TPQData_Out[TPQData_Out$Session=="Test",])
leveneTest(RD ~ Response*GAD, data = TPQData_Out[TPQData_Out$Session=="Test",])


C1 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD=="GAD"
C2 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD=="NoGad"
C3 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD=="GAD"
C4 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD=="NoGad"
C5 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="Responder"&TPQData_Out$GAD=="GAD"
C6 = TPQData_Out$Session=="Test"&TPQData_Out$Response=="Responder"&TPQData_Out$GAD=="NoGad"
C1R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD=="GAD"
C2R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="NoResponse"&TPQData_Out$GAD=="NoGad"
C3R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD=="GAD"
C4R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="Non-responder"&TPQData_Out$GAD=="NoGad"
C5R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="Responder"&TPQData_Out$GAD=="GAD"
C6R = TPQData_Out$Session=="Retest"&TPQData_Out$Response=="Responder"&TPQData_Out$GAD=="NoGad"

for (item in c("HA","NS","RD")){
  print (paste0("C1",item," = ",shapiro.test(TPQData_Out[[item]][C1])$p.value))
  print (paste0("C2",item," = ",shapiro.test(TPQData_Out[[item]][C2])$p.value))
  print (paste0("C3",item," = ",shapiro.test(TPQData_Out[[item]][C3])$p.value))
  print (paste0("C4",item," = ",shapiro.test(TPQData_Out[[item]][C4])$p.value))
  print (paste0("C5",item," = ",shapiro.test(TPQData_Out[[item]][C5])$p.value))
  print (paste0("C6",item," = ",shapiro.test(TPQData_Out[[item]][C6])$p.value))
  print (paste0("C1R",item," = ",shapiro.test(TPQData_Out[[item]][C1R])$p.value))
  print (paste0("C2R",item," = ",shapiro.test(TPQData_Out[[item]][C2R])$p.value))
  print (paste0("C3R",item," = ",shapiro.test(TPQData_Out[[item]][C3R])$p.value))
  print (paste0("C4R",item," = ",shapiro.test(TPQData_Out[[item]][C4R])$p.value))
  print (paste0("C5R",item," = ",shapiro.test(TPQData_Out[[item]][C5R])$p.value))
  print (paste0("C6R",item," = ",shapiro.test(TPQData_Out[[item]][C6R])$p.value))
}
#All Good, but Not normally distributed are NS: RetestResponderGad  and NS: TestNonresponderGAD


####################################################################################
#                                     ANOVAs                                       #
####################################################################################
#Make the reference level correctly assigned
TPQData_Out$Response = relevel(TPQData_Out$Response, ref = "NoResponse")
TPQData_Out$GAD = relevel(TPQData_Out$GAD, ref = "NoGad")

#Planned contrasts
contrasts(TPQData_Out$Response) = cbind(c(2,-1,-1),c(0,1,-1)) #1: NoResponse compared to all others & 2: responders vs nonresponders
contrasts(TPQData_Out$GAD) = c(1,-1)

Model_NS = aov(NS~Response*GAD*Session + Error(Subject/Session),data = TPQData_Out)
Model_HA = aov(HA~Response*GAD*Session + Error(Subject/Session),data = TPQData_Out)
Model_RD = aov(RD~Response*GAD*Session + Error(Subject/Session),data = TPQData_Out)
summary(Model_NS, type = "III")
summary(Model_HA, type = "III")
summary(Model_RD, type = "III")

AnovaModel_NS = ezANOVA(data = TPQData_Out,dv = NS,within = Session,between = .(Response,GAD),wid = Subject,type = 3,detailed = TRUE)
AnovaModel_HA = ezANOVA(data = TPQData_Out,dv = HA,within = Session,between = .(Response,GAD),wid = Subject,type = 3,detailed = TRUE)
AnovaModel_RD = ezANOVA(data = TPQData_Out,dv = RD,within = Session,between = .(Response,GAD),wid = Subject,type = 3,detailed = TRUE)
Anova_NS = FindAnova(AnovaModel_NS)
Anova_HA = FindAnova(AnovaModel_HA)
Anova_RD = FindAnova(AnovaModel_RD)


AnovaModel_HA_Test = ezANOVA(data = TPQData_Out[TPQData_Out$Session=="Test",],dv = HA,between = .(Response,GAD),wid = Subject,type = 3,detailed = TRUE)

mean(TPQData_Out$HA[TPQData_Out$GAD=="GAD" & TPQData_Out$Session=="Test"])

DataMeans = UMeans(TPQData_Out, DV = "HA", IV1 = "Response", IV2 = "Session")#, IV3 = "GAD")
DataMeans$Response = as.character(DataMeans$Response)
DataMeans$Response[DataMeans$Response=="NoResponse"] = "Non-MDD"
DataMeans$Response = factor (DataMeans$Response,levels = c("Non-MDD","Non-responder","Responder"))
DataMeans$Session = relevel(DataMeans$Session, ref = "Test")
ggplot(DataMeans, aes(x = Session, y= HA, color=Response,group = Response))+geom_line(position = position_dodge(0.15),size=1.2)+
  geom_errorbar(aes(ymax = HA+SEM,ymin = HA-SEM), position = position_dodge(0.15), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white", position = position_dodge(0.15))+
  TypicalTheme+scale_color_manual(values = c("#00A550","#Eb4C42","#0087BD"))+ggtitle("Harm Avoidance Comparison According to Response")

DataMeans = UMeans(TPQData_Out, DV = "HA", IV1 = "Response", IV2 = "Session", IV3 = "GAD")
DataMeans$Session = relevel(DataMeans$Session, ref = "Test")
ggplot(DataMeans, aes(x = Session, y= HA, color=GAD,group = GAD))+geom_line(position = position_dodge(0.15),size=1.2)+
  facet_wrap(~Response)+geom_errorbar(aes(ymax = HA+SEM,ymin = HA-SEM), position = position_dodge(0.15), width = 0.25,color="black")+
  geom_point(shape = 21, color = "black",fill = "white", position = position_dodge(0.15))+
  TypicalTheme+scale_color_manual(values = c("#808CA3","#B9B0AB"))+ggtitle("Harm Avoidance Comparison According to Presence of GAD",subtitle = "Separated According to Response")




####################################################################################
#         Compare HC and MDD using the new ICA factorization by Juergen            #
####################################################################################

#Cloninger factorization
TPQ_Q = Attach.H(TPQData_Out, TPQData,SubjectID1 = "ID",SubjectID2 = "ID")

TPQ_Q$NS1=as.numeric(TPQ_Q$Q2)+as.numeric(TPQ_Q$Q4)+as.numeric(TPQ_Q$Q9)+as.numeric(TPQ_Q$Q11)+as.numeric(TPQ_Q$Q40)+as.numeric(TPQ_Q$Q43)+as.numeric(TPQ_Q$Q85)+as.numeric(TPQ_Q$Q93)+as.numeric(TPQ_Q$Q96)
TPQ_Q$NS2=as.numeric(TPQ_Q$Q30)+as.numeric(TPQ_Q$Q46)+as.numeric(TPQ_Q$Q48)+as.numeric(TPQ_Q$Q50)+as.numeric(TPQ_Q$Q55)+as.numeric(TPQ_Q$Q56)+as.numeric(TPQ_Q$Q81)+as.numeric(TPQ_Q$Q99)
TPQ_Q$NS3=as.numeric(TPQ_Q$Q32)+as.numeric(TPQ_Q$Q66)+as.numeric(TPQ_Q$Q70)+as.numeric(TPQ_Q$Q72)+as.numeric(TPQ_Q$Q76)+as.numeric(TPQ_Q$Q78)+as.numeric(TPQ_Q$Q87)
TPQ_Q$NS4=as.numeric(TPQ_Q$Q13)+as.numeric(TPQ_Q$Q16)+as.numeric(TPQ_Q$Q21)+as.numeric(TPQ_Q$Q22)+as.numeric(TPQ_Q$Q24)+as.numeric(TPQ_Q$Q28)+as.numeric(TPQ_Q$Q35)+as.numeric(TPQ_Q$Q60)+as.numeric(TPQ_Q$Q62)+as.numeric(TPQ_Q$Q65)
TPQ_Q$HA1=as.numeric(TPQ_Q$Q1)+as.numeric(TPQ_Q$Q5)+as.numeric(TPQ_Q$Q8)+as.numeric(TPQ_Q$Q10)+as.numeric(TPQ_Q$Q14)+as.numeric(TPQ_Q$Q82)+as.numeric(TPQ_Q$Q84)+as.numeric(TPQ_Q$Q91)+as.numeric(TPQ_Q$Q95)+as.numeric(TPQ_Q$Q98)
TPQ_Q$HA2=as.numeric(TPQ_Q$Q18)+as.numeric(TPQ_Q$Q19)+as.numeric(TPQ_Q$Q23)+as.numeric(TPQ_Q$Q26)+as.numeric(TPQ_Q$Q29)+as.numeric(TPQ_Q$Q47)+as.numeric(TPQ_Q$Q51)
TPQ_Q$HA3=as.numeric(TPQ_Q$Q33)+as.numeric(TPQ_Q$Q37)+as.numeric(TPQ_Q$Q38)+as.numeric(TPQ_Q$Q42)+as.numeric(TPQ_Q$Q44)+as.numeric(TPQ_Q$Q89)+as.numeric(TPQ_Q$Q100)
TPQ_Q$HA4=as.numeric(TPQ_Q$Q49)+as.numeric(TPQ_Q$Q54)+as.numeric(TPQ_Q$Q57)+as.numeric(TPQ_Q$Q59)+as.numeric(TPQ_Q$Q63)+as.numeric(TPQ_Q$Q68)+as.numeric(TPQ_Q$Q69)+as.numeric(TPQ_Q$Q73)+as.numeric(TPQ_Q$Q75)+as.numeric(TPQ_Q$Q80)
TPQ_Q$RD1=as.numeric(TPQ_Q$Q27)+as.numeric(TPQ_Q$Q31)+as.numeric(TPQ_Q$Q34)+as.numeric(TPQ_Q$Q83)+as.numeric(TPQ_Q$Q94)
TPQ_Q$RD2=as.numeric(TPQ_Q$Q39)+as.numeric(TPQ_Q$Q41)+as.numeric(TPQ_Q$Q45)+as.numeric(TPQ_Q$Q52)+as.numeric(TPQ_Q$Q53)+as.numeric(TPQ_Q$Q77)+as.numeric(TPQ_Q$Q79)+as.numeric(TPQ_Q$Q92)+as.numeric(TPQ_Q$Q97)
TPQ_Q$RD3=as.numeric(TPQ_Q$Q3)+as.numeric(TPQ_Q$Q6)+as.numeric(TPQ_Q$Q7)+as.numeric(TPQ_Q$Q12)+as.numeric(TPQ_Q$Q15)+as.numeric(TPQ_Q$Q64)+as.numeric(TPQ_Q$Q67)+as.numeric(TPQ_Q$Q74)+as.numeric(TPQ_Q$Q86)+as.numeric(TPQ_Q$Q88)+as.numeric(TPQ_Q$Q90)
TPQ_Q$RD4=as.numeric(TPQ_Q$Q17)+as.numeric(TPQ_Q$Q20)+as.numeric(TPQ_Q$Q25)+as.numeric(TPQ_Q$Q36)+as.numeric(TPQ_Q$Q58)

#Old Classification (before 4-May)
TPQ_Q$IC_0 =as.numeric(TPQ_Q$Q6)+as.numeric(TPQ_Q$Q10)+as.numeric(TPQ_Q$Q28)+as.numeric(TPQ_Q$Q37)+as.numeric(TPQ_Q$Q43)+as.numeric(TPQ_Q$Q44)+as.numeric(TPQ_Q$Q49)+as.numeric(TPQ_Q$Q53)+as.numeric(TPQ_Q$Q65)+as.numeric(TPQ_Q$Q68)+as.numeric(TPQ_Q$Q69)+as.numeric(TPQ_Q$Q80)
TPQ_Q$IC_1 =as.numeric(TPQ_Q$Q33)+as.numeric(TPQ_Q$Q43)+as.numeric(TPQ_Q$Q50)+as.numeric(TPQ_Q$Q54)+as.numeric(TPQ_Q$Q57)+as.numeric(TPQ_Q$Q59)+as.numeric(TPQ_Q$Q63)+as.numeric(TPQ_Q$Q68)+as.numeric(TPQ_Q$Q69)+as.numeric(TPQ_Q$Q75)+as.numeric(TPQ_Q$Q80)
TPQ_Q$IC_2 =as.numeric(TPQ_Q$Q32)+as.numeric(TPQ_Q$Q45)+as.numeric(TPQ_Q$Q52)+as.numeric(TPQ_Q$Q53)+as.numeric(TPQ_Q$Q75)+as.numeric(TPQ_Q$Q79)+as.numeric(TPQ_Q$Q80)+as.numeric(TPQ_Q$Q84)+as.numeric(TPQ_Q$Q89)+as.numeric(TPQ_Q$Q90)+as.numeric(TPQ_Q$Q100)
TPQ_Q$IC_3 =as.numeric(TPQ_Q$Q11)+as.numeric(TPQ_Q$Q19)+as.numeric(TPQ_Q$Q20)+as.numeric(TPQ_Q$Q28)+as.numeric(TPQ_Q$Q35)+as.numeric(TPQ_Q$Q46)+as.numeric(TPQ_Q$Q50)+as.numeric(TPQ_Q$Q55)+as.numeric(TPQ_Q$Q100)
TPQ_Q$IC_4 =as.numeric(TPQ_Q$Q16)+as.numeric(TPQ_Q$Q21)+as.numeric(TPQ_Q$Q39)+as.numeric(TPQ_Q$Q41)+as.numeric(TPQ_Q$Q45)+as.numeric(TPQ_Q$Q63)+as.numeric(TPQ_Q$Q92)
TPQ_Q$IC_5 =as.numeric(TPQ_Q$Q11)+as.numeric(TPQ_Q$Q36)+as.numeric(TPQ_Q$Q46)+as.numeric(TPQ_Q$Q48)+as.numeric(TPQ_Q$Q55)+as.numeric(TPQ_Q$Q56)+as.numeric(TPQ_Q$Q77)
TPQ_Q$IC_6 =as.numeric(TPQ_Q$Q3)+as.numeric(TPQ_Q$Q12)+as.numeric(TPQ_Q$Q15)+as.numeric(TPQ_Q$Q74)+as.numeric(TPQ_Q$Q86)+as.numeric(TPQ_Q$Q88)+as.numeric(TPQ_Q$Q90)
TPQ_Q$IC_7 =as.numeric(TPQ_Q$Q1)+as.numeric(TPQ_Q$Q5)+as.numeric(TPQ_Q$Q6)+as.numeric(TPQ_Q$Q8)+as.numeric(TPQ_Q$Q15)+as.numeric(TPQ_Q$Q18)+as.numeric(TPQ_Q$Q19)+as.numeric(TPQ_Q$Q23)+as.numeric(TPQ_Q$Q26)+as.numeric(TPQ_Q$Q28)+as.numeric(TPQ_Q$Q30)+as.numeric(TPQ_Q$Q32)+as.numeric(TPQ_Q$Q34)+as.numeric(TPQ_Q$Q35)+as.numeric(TPQ_Q$Q36)+as.numeric(TPQ_Q$Q63)+as.numeric(TPQ_Q$Q84)+as.numeric(TPQ_Q$Q90)+as.numeric(TPQ_Q$Q91)+as.numeric(TPQ_Q$Q98)
TPQ_Q$IC_8 =as.numeric(TPQ_Q$Q66)+as.numeric(TPQ_Q$Q70)+as.numeric(TPQ_Q$Q72)+as.numeric(TPQ_Q$Q76)+as.numeric(TPQ_Q$Q87)
TPQ_Q$IC_9 =as.numeric(TPQ_Q$Q2)+as.numeric(TPQ_Q$Q4)+as.numeric(TPQ_Q$Q13)+as.numeric(TPQ_Q$Q24)+as.numeric(TPQ_Q$Q26)+as.numeric(TPQ_Q$Q29)+as.numeric(TPQ_Q$Q47)+as.numeric(TPQ_Q$Q51)+as.numeric(TPQ_Q$Q85)
TPQ_Q$IC_10 =as.numeric(TPQ_Q$Q11)+as.numeric(TPQ_Q$Q17)+as.numeric(TPQ_Q$Q19)+as.numeric(TPQ_Q$Q20)+as.numeric(TPQ_Q$Q25)+as.numeric(TPQ_Q$Q65)+as.numeric(TPQ_Q$Q73)+as.numeric(TPQ_Q$Q77)
TPQ_Q$IC_11 =as.numeric(TPQ_Q$Q14)+as.numeric(TPQ_Q$Q20)+as.numeric(TPQ_Q$Q24)+as.numeric(TPQ_Q$Q30)+as.numeric(TPQ_Q$Q33)+as.numeric(TPQ_Q$Q56)+as.numeric(TPQ_Q$Q60)+as.numeric(TPQ_Q$Q62)+as.numeric(TPQ_Q$Q65)
TPQ_Q$IC_12 =as.numeric(TPQ_Q$Q10)+as.numeric(TPQ_Q$Q18)+as.numeric(TPQ_Q$Q19)+as.numeric(TPQ_Q$Q23)+as.numeric(TPQ_Q$Q33)+as.numeric(TPQ_Q$Q35)+as.numeric(TPQ_Q$Q38)+as.numeric(TPQ_Q$Q40)+as.numeric(TPQ_Q$Q50)+as.numeric(TPQ_Q$Q73)+as.numeric(TPQ_Q$Q85)+as.numeric(TPQ_Q$Q89)+as.numeric(TPQ_Q$Q93)


# ICA according to New analysis by Juergin comparing MDD and HC in 4-May
TPQ_Q$IC_0_HC = as.numeric(TPQ_Q$QO43)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO49)+as.numeric(TPQ_Q$QO53)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO77)
TPQ_Q$IC_1_HC = as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO43)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO54)+as.numeric(TPQ_Q$QO57)+as.numeric(TPQ_Q$QO59)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO80)
TPQ_Q$IC_2_HC = as.numeric(TPQ_Q$QO16)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO39)+as.numeric(TPQ_Q$QO41)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO92)
TPQ_Q$IC_3_HC = as.numeric(TPQ_Q$QO19)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO35)+as.numeric(TPQ_Q$QO46)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO100)
TPQ_Q$IC_4_HC = as.numeric(TPQ_Q$QO18)+as.numeric(TPQ_Q$QO19)+as.numeric(TPQ_Q$QO23)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO85)+as.numeric(TPQ_Q$QO98)
TPQ_Q$IC_5_HC = as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO22)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO36)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO79)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_6_HC = as.numeric(TPQ_Q$QO3)+as.numeric(TPQ_Q$QO12)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO74)+as.numeric(TPQ_Q$QO86)+as.numeric(TPQ_Q$QO88)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_7_HC = as.numeric(TPQ_Q$QO1)+as.numeric(TPQ_Q$QO5)+as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO18)+as.numeric(TPQ_Q$QO19)+as.numeric(TPQ_Q$QO23)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO30)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO35)+as.numeric(TPQ_Q$QO44)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO91)+as.numeric(TPQ_Q$QO98)
TPQ_Q$IC_8_HC = as.numeric(TPQ_Q$QO66)+as.numeric(TPQ_Q$QO70)+as.numeric(TPQ_Q$QO72)+as.numeric(TPQ_Q$QO76)+as.numeric(TPQ_Q$QO87)
TPQ_Q$IC_9_HC = as.numeric(TPQ_Q$QO2)+as.numeric(TPQ_Q$QO4)+as.numeric(TPQ_Q$QO13)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO29)+as.numeric(TPQ_Q$QO47)+as.numeric(TPQ_Q$QO51)+as.numeric(TPQ_Q$QO85)
TPQ_Q$IC_10_HC = as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO17)+as.numeric(TPQ_Q$QO21)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO71)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO77)
TPQ_Q$IC_11_HC = as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO30)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO60)+as.numeric(TPQ_Q$QO62)+as.numeric(TPQ_Q$QO65)
TPQ_Q$IC_12_HC = as.numeric(TPQ_Q$QO4)+as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO46)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO71)+as.numeric(TPQ_Q$QO77)

TPQ_Q$IC_0_MDD = as.numeric(TPQ_Q$QO13)+as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO42)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO71)+as.numeric(TPQ_Q$QO77)+as.numeric(TPQ_Q$QO85)
TPQ_Q$IC_1_MDD = as.numeric(TPQ_Q$QO5)+as.numeric(TPQ_Q$QO7)+as.numeric(TPQ_Q$QO10)+as.numeric(TPQ_Q$QO13)+as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO17)+as.numeric(TPQ_Q$QO19)+as.numeric(TPQ_Q$QO23)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO91)+as.numeric(TPQ_Q$QO92)
TPQ_Q$IC_2_MDD = as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO22)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO30)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO35)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO49)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO57)+as.numeric(TPQ_Q$QO66)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO70)+as.numeric(TPQ_Q$QO72)+as.numeric(TPQ_Q$QO74)+as.numeric(TPQ_Q$QO90)+as.numeric(TPQ_Q$QO92)
TPQ_Q$IC_3_MDD = as.numeric(TPQ_Q$QO10)+as.numeric(TPQ_Q$QO17)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO49)+as.numeric(TPQ_Q$QO60)+as.numeric(TPQ_Q$QO77)+as.numeric(TPQ_Q$QO79)+as.numeric(TPQ_Q$QO80)
TPQ_Q$IC_4_MDD = as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO9)+as.numeric(TPQ_Q$QO47)+as.numeric(TPQ_Q$QO51)+as.numeric(TPQ_Q$QO60)+as.numeric(TPQ_Q$QO62)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO66)+as.numeric(TPQ_Q$QO84)
TPQ_Q$IC_5_MDD = as.numeric(TPQ_Q$QO9)+as.numeric(TPQ_Q$QO12)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO21)+as.numeric(TPQ_Q$QO30)+as.numeric(TPQ_Q$QO37)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO59)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO81)+as.numeric(TPQ_Q$QO86)+as.numeric(TPQ_Q$QO88)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_6_MDD = as.numeric(TPQ_Q$QO3)+as.numeric(TPQ_Q$QO5)+as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO49)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO70)+as.numeric(TPQ_Q$QO76)+as.numeric(TPQ_Q$QO78)+as.numeric(TPQ_Q$QO87)+as.numeric(TPQ_Q$QO88)
TPQ_Q$IC_7_MDD = as.numeric(TPQ_Q$QO2)+as.numeric(TPQ_Q$QO3)+as.numeric(TPQ_Q$QO4)+as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO20)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO29)+as.numeric(TPQ_Q$QO34)+as.numeric(TPQ_Q$QO47)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO51)+as.numeric(TPQ_Q$QO58)+as.numeric(TPQ_Q$QO60)+as.numeric(TPQ_Q$QO62)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO86)+as.numeric(TPQ_Q$QO88)
TPQ_Q$IC_8_MDD = as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO12)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO36)+as.numeric(TPQ_Q$QO37)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO57)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO71)+as.numeric(TPQ_Q$QO76)+as.numeric(TPQ_Q$QO78)+as.numeric(TPQ_Q$QO86)+as.numeric(TPQ_Q$QO87)+as.numeric(TPQ_Q$QO96)
TPQ_Q$IC_9_MDD = as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO46)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO59)+as.numeric(TPQ_Q$QO79)+as.numeric(TPQ_Q$QO91)
TPQ_Q$IC_10_MDD = as.numeric(TPQ_Q$QO2)+as.numeric(TPQ_Q$QO3)+as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO29)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO37)+as.numeric(TPQ_Q$QO47)+as.numeric(TPQ_Q$QO51)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO76)+as.numeric(TPQ_Q$QO84)
TPQ_Q$IC_11_MDD = as.numeric(TPQ_Q$QO2)+as.numeric(TPQ_Q$QO10)+as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO62)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO92)+as.numeric(TPQ_Q$QO97)
TPQ_Q$IC_12_MDD = as.numeric(TPQ_Q$QO1)+as.numeric(TPQ_Q$QO4)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO36)+as.numeric(TPQ_Q$QO42)+as.numeric(TPQ_Q$QO44)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO52)+as.numeric(TPQ_Q$QO59)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO77)+as.numeric(TPQ_Q$QO79)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO82)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO95)+as.numeric(TPQ_Q$QO97)+as.numeric(TPQ_Q$QO98)+as.numeric(TPQ_Q$QO100)

TPQ_Q$IC_0_All = as.numeric(TPQ_Q$QO1)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO42)+as.numeric(TPQ_Q$QO44)+as.numeric(TPQ_Q$QO49)+as.numeric(TPQ_Q$QO52)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO78)+as.numeric(TPQ_Q$QO79)+as.numeric(TPQ_Q$QO80)+as.numeric(TPQ_Q$QO82)+as.numeric(TPQ_Q$QO95)+as.numeric(TPQ_Q$QO96)+as.numeric(TPQ_Q$QO97)+as.numeric(TPQ_Q$QO98)
TPQ_Q$IC_1_All = as.numeric(TPQ_Q$QO66)+as.numeric(TPQ_Q$QO70)+as.numeric(TPQ_Q$QO72)+as.numeric(TPQ_Q$QO76)+as.numeric(TPQ_Q$QO87)
TPQ_Q$IC_2_All = as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO32)+as.numeric(TPQ_Q$QO35)+as.numeric(TPQ_Q$QO36)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_3_All = as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO17)+as.numeric(TPQ_Q$QO21)+as.numeric(TPQ_Q$QO46)+as.numeric(TPQ_Q$QO52)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO85)
TPQ_Q$IC_4_All = as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO77)+as.numeric(TPQ_Q$QO80)
TPQ_Q$IC_5_All = as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO60)+as.numeric(TPQ_Q$QO62)+as.numeric(TPQ_Q$QO65)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO100)
TPQ_Q$IC_6_All = as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO39)+as.numeric(TPQ_Q$QO41)+as.numeric(TPQ_Q$QO45)+as.numeric(TPQ_Q$QO92)
TPQ_Q$IC_7_All = as.numeric(TPQ_Q$QO1)+as.numeric(TPQ_Q$QO5)+as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO10)+as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO14)+as.numeric(TPQ_Q$QO18)+as.numeric(TPQ_Q$QO19)+as.numeric(TPQ_Q$QO23)+as.numeric(TPQ_Q$QO30)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO37)+as.numeric(TPQ_Q$QO38)+as.numeric(TPQ_Q$QO40)+as.numeric(TPQ_Q$QO44)+as.numeric(TPQ_Q$QO50)+as.numeric(TPQ_Q$QO73)+as.numeric(TPQ_Q$QO84)+as.numeric(TPQ_Q$QO89)+as.numeric(TPQ_Q$QO91)+as.numeric(TPQ_Q$QO98)
TPQ_Q$IC_8_All = as.numeric(TPQ_Q$QO35)+as.numeric(TPQ_Q$QO46)+as.numeric(TPQ_Q$QO48)+as.numeric(TPQ_Q$QO55)+as.numeric(TPQ_Q$QO56)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO71)+as.numeric(TPQ_Q$QO77)+as.numeric(TPQ_Q$QO91)
TPQ_Q$IC_9_All = as.numeric(TPQ_Q$QO8)+as.numeric(TPQ_Q$QO10)+as.numeric(TPQ_Q$QO11)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO17)+as.numeric(TPQ_Q$QO20)+as.numeric(TPQ_Q$QO25)+as.numeric(TPQ_Q$QO58)+as.numeric(TPQ_Q$QO85)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_10_All = as.numeric(TPQ_Q$QO28)+as.numeric(TPQ_Q$QO33)+as.numeric(TPQ_Q$QO43)+as.numeric(TPQ_Q$QO54)+as.numeric(TPQ_Q$QO57)+as.numeric(TPQ_Q$QO59)+as.numeric(TPQ_Q$QO63)+as.numeric(TPQ_Q$QO68)+as.numeric(TPQ_Q$QO69)+as.numeric(TPQ_Q$QO75)+as.numeric(TPQ_Q$QO80)
TPQ_Q$IC_11_All = as.numeric(TPQ_Q$QO3)+as.numeric(TPQ_Q$QO6)+as.numeric(TPQ_Q$QO12)+as.numeric(TPQ_Q$QO15)+as.numeric(TPQ_Q$QO74)+as.numeric(TPQ_Q$QO86)+as.numeric(TPQ_Q$QO88)+as.numeric(TPQ_Q$QO90)
TPQ_Q$IC_12_All = as.numeric(TPQ_Q$QO2)+as.numeric(TPQ_Q$QO4)+as.numeric(TPQ_Q$QO24)+as.numeric(TPQ_Q$QO26)+as.numeric(TPQ_Q$QO29)+as.numeric(TPQ_Q$QO47)+as.numeric(TPQ_Q$QO51)+as.numeric(TPQ_Q$QO85)

#Finally, the full range of questions
TPQ_Q$Complete_IC_0_HC =  0.048 * as.numeric(TPQ_Q$QO1)+ 0.029 * as.numeric(TPQ_Q$QO2)+ -0.043 * as.numeric(TPQ_Q$QO3)+ -0.132 * as.numeric(TPQ_Q$QO4)+ -0.061 * as.numeric(TPQ_Q$QO5)+ 0.139 * as.numeric(TPQ_Q$QO6)+ -0.072 * as.numeric(TPQ_Q$QO7)+ 0.06 * as.numeric(TPQ_Q$QO8)+ 0.039 * as.numeric(TPQ_Q$QO9)+ -0.194 * as.numeric(TPQ_Q$QO10)+ 0.035 * as.numeric(TPQ_Q$QO11)+ 0.01 * as.numeric(TPQ_Q$QO12)+ -0.085 * as.numeric(TPQ_Q$QO13)+ -0.081 * as.numeric(TPQ_Q$QO14)+ 0.149 * as.numeric(TPQ_Q$QO15)+ -0.119 * as.numeric(TPQ_Q$QO16)+ -0.109 * as.numeric(TPQ_Q$QO17)+ 0.053 * as.numeric(TPQ_Q$QO18)+ 0.211 * as.numeric(TPQ_Q$QO19)+ 0.002 * as.numeric(TPQ_Q$QO20)+ -0.145 * as.numeric(TPQ_Q$QO21)+ -0.065 * as.numeric(TPQ_Q$QO22)+ 0.181 * as.numeric(TPQ_Q$QO23)+ -0.018 * as.numeric(TPQ_Q$QO24)+ 0.07 * as.numeric(TPQ_Q$QO25)+ 0.002 * as.numeric(TPQ_Q$QO26)+ -0.016 * as.numeric(TPQ_Q$QO27)+ -0.159 * as.numeric(TPQ_Q$QO28)+ -0.046 * as.numeric(TPQ_Q$QO29)+ 0.085 * as.numeric(TPQ_Q$QO30)+ 0.027 * as.numeric(TPQ_Q$QO31)+ -0.027 * as.numeric(TPQ_Q$QO32)+ -0.023 * as.numeric(TPQ_Q$QO33)+ -0.066 * as.numeric(TPQ_Q$QO34)+ 0.133 * as.numeric(TPQ_Q$QO35)+ -0.02 * as.numeric(TPQ_Q$QO36)+ -0.126 * as.numeric(TPQ_Q$QO37)+ 0.03 * as.numeric(TPQ_Q$QO38)+ 0.114 * as.numeric(TPQ_Q$QO39)+ 0.026 * as.numeric(TPQ_Q$QO40)+ 0.107 * as.numeric(TPQ_Q$QO41)+ 0.122 * as.numeric(TPQ_Q$QO42)+ -0.326 * as.numeric(TPQ_Q$QO43)+ 0.144 * as.numeric(TPQ_Q$QO44)+ -0.266 * as.numeric(TPQ_Q$QO45)+ -0.016 * as.numeric(TPQ_Q$QO46)+ -0.054 * as.numeric(TPQ_Q$QO47)+ -0.159 * as.numeric(TPQ_Q$QO48)+ -0.253 * as.numeric(TPQ_Q$QO49)+ 0.104 * as.numeric(TPQ_Q$QO50)+ 0.044 * as.numeric(TPQ_Q$QO51)+ -0.084 * as.numeric(TPQ_Q$QO52)+ -0.416 * as.numeric(TPQ_Q$QO53)+ 0.02 * as.numeric(TPQ_Q$QO54)+ 0.024 * as.numeric(TPQ_Q$QO55)+ 0.046 * as.numeric(TPQ_Q$QO56)+ -0.084 * as.numeric(TPQ_Q$QO57)+ -0.01 * as.numeric(TPQ_Q$QO58)+ -0.066 * as.numeric(TPQ_Q$QO59)+ 0.056 * as.numeric(TPQ_Q$QO60)+ -0.056 * as.numeric(TPQ_Q$QO61)+ -0.029 * as.numeric(TPQ_Q$QO62)+ 0.042 * as.numeric(TPQ_Q$QO63)+ -0.069 * as.numeric(TPQ_Q$QO64)+ -0.425 * as.numeric(TPQ_Q$QO65)+ -0.008 * as.numeric(TPQ_Q$QO66)+ -0.001 * as.numeric(TPQ_Q$QO67)+ -0.124 * as.numeric(TPQ_Q$QO68)+ -0.18 * as.numeric(TPQ_Q$QO69)+ -0.069 * as.numeric(TPQ_Q$QO70)+ -0.086 * as.numeric(TPQ_Q$QO71)+ -0.071 * as.numeric(TPQ_Q$QO72)+ -0.005 * as.numeric(TPQ_Q$QO73)+ -0.074 * as.numeric(TPQ_Q$QO74)+ -0.036 * as.numeric(TPQ_Q$QO75)+ 0.02 * as.numeric(TPQ_Q$QO76)+ 0.315 * as.numeric(TPQ_Q$QO77)+ -0.067 * as.numeric(TPQ_Q$QO78)+ -0.13 * as.numeric(TPQ_Q$QO79)+ 0.066 * as.numeric(TPQ_Q$QO80)+ -0.021 * as.numeric(TPQ_Q$QO81)+ -0.018 * as.numeric(TPQ_Q$QO82)+ -0.044 * as.numeric(TPQ_Q$QO83)+ -0.038 * as.numeric(TPQ_Q$QO84)+ -0.069 * as.numeric(TPQ_Q$QO85)+ -0.044 * as.numeric(TPQ_Q$QO86)+ -0.067 * as.numeric(TPQ_Q$QO87)+ -0.002 * as.numeric(TPQ_Q$QO88)+ -0.084 * as.numeric(TPQ_Q$QO89)+ 0.123 * as.numeric(TPQ_Q$QO90)+ -0.117 * as.numeric(TPQ_Q$QO91)+ -0.08 * as.numeric(TPQ_Q$QO92)+ 0.076 * as.numeric(TPQ_Q$QO93)+ 0.019 * as.numeric(TPQ_Q$QO94)+ -0.005 * as.numeric(TPQ_Q$QO95)+ -0.032 * as.numeric(TPQ_Q$QO96)+ 0.121 * as.numeric(TPQ_Q$QO97)+ -0.039 * as.numeric(TPQ_Q$QO98)+ 0.005 * as.numeric(TPQ_Q$QO99)+ 0.031 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_1_HC =  0.117 * as.numeric(TPQ_Q$QO1)+ -0.037 * as.numeric(TPQ_Q$QO2)+ 0.016 * as.numeric(TPQ_Q$QO3)+ 0.051 * as.numeric(TPQ_Q$QO4)+ -0.178 * as.numeric(TPQ_Q$QO5)+ -0.162 * as.numeric(TPQ_Q$QO6)+ -0.05 * as.numeric(TPQ_Q$QO7)+ 0.09 * as.numeric(TPQ_Q$QO8)+ 0.046 * as.numeric(TPQ_Q$QO9)+ -0.182 * as.numeric(TPQ_Q$QO10)+ 0.012 * as.numeric(TPQ_Q$QO11)+ -0.107 * as.numeric(TPQ_Q$QO12)+ -0.128 * as.numeric(TPQ_Q$QO13)+ -0.176 * as.numeric(TPQ_Q$QO14)+ -0.154 * as.numeric(TPQ_Q$QO15)+ 0.097 * as.numeric(TPQ_Q$QO16)+ 0.032 * as.numeric(TPQ_Q$QO17)+ -0.139 * as.numeric(TPQ_Q$QO18)+ -0.168 * as.numeric(TPQ_Q$QO19)+ -0.041 * as.numeric(TPQ_Q$QO20)+ 0.058 * as.numeric(TPQ_Q$QO21)+ -0.103 * as.numeric(TPQ_Q$QO22)+ -0.192 * as.numeric(TPQ_Q$QO23)+ -0.144 * as.numeric(TPQ_Q$QO24)+ 0.044 * as.numeric(TPQ_Q$QO25)+ -0.038 * as.numeric(TPQ_Q$QO26)+ -0.012 * as.numeric(TPQ_Q$QO27)+ -0.182 * as.numeric(TPQ_Q$QO28)+ 0.019 * as.numeric(TPQ_Q$QO29)+ -0.196 * as.numeric(TPQ_Q$QO30)+ 0.019 * as.numeric(TPQ_Q$QO31)+ 0.139 * as.numeric(TPQ_Q$QO32)+ -0.241 * as.numeric(TPQ_Q$QO33)+ -0.133 * as.numeric(TPQ_Q$QO34)+ 0.056 * as.numeric(TPQ_Q$QO35)+ 0.094 * as.numeric(TPQ_Q$QO36)+ -0.121 * as.numeric(TPQ_Q$QO37)+ -0.099 * as.numeric(TPQ_Q$QO38)+ -0.009 * as.numeric(TPQ_Q$QO39)+ -0.16 * as.numeric(TPQ_Q$QO40)+ -0.044 * as.numeric(TPQ_Q$QO41)+ 0.046 * as.numeric(TPQ_Q$QO42)+ -0.301 * as.numeric(TPQ_Q$QO43)+ 0.025 * as.numeric(TPQ_Q$QO44)+ -0.051 * as.numeric(TPQ_Q$QO45)+ -0.021 * as.numeric(TPQ_Q$QO46)+ -0.002 * as.numeric(TPQ_Q$QO47)+ -0.182 * as.numeric(TPQ_Q$QO48)+ -0.161 * as.numeric(TPQ_Q$QO49)+ -0.258 * as.numeric(TPQ_Q$QO50)+ 0.059 * as.numeric(TPQ_Q$QO51)+ 0.068 * as.numeric(TPQ_Q$QO52)+ -0.091 * as.numeric(TPQ_Q$QO53)+ -0.704 * as.numeric(TPQ_Q$QO54)+ 0.026 * as.numeric(TPQ_Q$QO55)+ -0.139 * as.numeric(TPQ_Q$QO56)+ -0.749 * as.numeric(TPQ_Q$QO57)+ 0.071 * as.numeric(TPQ_Q$QO58)+ 0.783 * as.numeric(TPQ_Q$QO59)+ 0.003 * as.numeric(TPQ_Q$QO60)+ -0.167 * as.numeric(TPQ_Q$QO61)+ -0.067 * as.numeric(TPQ_Q$QO62)+ 0.314 * as.numeric(TPQ_Q$QO63)+ -0.046 * as.numeric(TPQ_Q$QO64)+ 0.059 * as.numeric(TPQ_Q$QO65)+ 0.013 * as.numeric(TPQ_Q$QO66)+ 0.025 * as.numeric(TPQ_Q$QO67)+ -0.28 * as.numeric(TPQ_Q$QO68)+ -0.429 * as.numeric(TPQ_Q$QO69)+ -0.113 * as.numeric(TPQ_Q$QO70)+ 0.093 * as.numeric(TPQ_Q$QO71)+ -0.115 * as.numeric(TPQ_Q$QO72)+ -0.122 * as.numeric(TPQ_Q$QO73)+ 0.001 * as.numeric(TPQ_Q$QO74)+ 0.241 * as.numeric(TPQ_Q$QO75)+ -0.026 * as.numeric(TPQ_Q$QO76)+ 0.08 * as.numeric(TPQ_Q$QO77)+ -0.032 * as.numeric(TPQ_Q$QO78)+ 0.021 * as.numeric(TPQ_Q$QO79)+ 0.33 * as.numeric(TPQ_Q$QO80)+ 0.064 * as.numeric(TPQ_Q$QO81)+ 0.086 * as.numeric(TPQ_Q$QO82)+ -0.039 * as.numeric(TPQ_Q$QO83)+ 0.083 * as.numeric(TPQ_Q$QO84)+ -0.117 * as.numeric(TPQ_Q$QO85)+ 0.021 * as.numeric(TPQ_Q$QO86)+ -0.006 * as.numeric(TPQ_Q$QO87)+ -0.001 * as.numeric(TPQ_Q$QO88)+ 0.142 * as.numeric(TPQ_Q$QO89)+ -0.154 * as.numeric(TPQ_Q$QO90)+ 0.078 * as.numeric(TPQ_Q$QO91)+ 0.185 * as.numeric(TPQ_Q$QO92)+ 0.025 * as.numeric(TPQ_Q$QO93)+ 0.017 * as.numeric(TPQ_Q$QO94)+ 0.089 * as.numeric(TPQ_Q$QO95)+ -0.012 * as.numeric(TPQ_Q$QO96)+ 0.154 * as.numeric(TPQ_Q$QO97)+ 0.13 * as.numeric(TPQ_Q$QO98)+ 0.038 * as.numeric(TPQ_Q$QO99)+ 0.105 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_2_HC =  0.035 * as.numeric(TPQ_Q$QO1)+ 0.035 * as.numeric(TPQ_Q$QO2)+ -0.059 * as.numeric(TPQ_Q$QO3)+ -0.102 * as.numeric(TPQ_Q$QO4)+ -0.152 * as.numeric(TPQ_Q$QO5)+ -0.095 * as.numeric(TPQ_Q$QO6)+ 0.07 * as.numeric(TPQ_Q$QO7)+ 0.024 * as.numeric(TPQ_Q$QO8)+ 0.041 * as.numeric(TPQ_Q$QO9)+ 0.042 * as.numeric(TPQ_Q$QO10)+ 0.06 * as.numeric(TPQ_Q$QO11)+ 0.036 * as.numeric(TPQ_Q$QO12)+ 0.058 * as.numeric(TPQ_Q$QO13)+ 0.025 * as.numeric(TPQ_Q$QO14)+ -0.119 * as.numeric(TPQ_Q$QO15)+ -0.244 * as.numeric(TPQ_Q$QO16)+ 0.037 * as.numeric(TPQ_Q$QO17)+ 0.002 * as.numeric(TPQ_Q$QO18)+ 0.042 * as.numeric(TPQ_Q$QO19)+ -0.208 * as.numeric(TPQ_Q$QO20)+ -0.229 * as.numeric(TPQ_Q$QO21)+ 0.024 * as.numeric(TPQ_Q$QO22)+ -0.0 * as.numeric(TPQ_Q$QO23)+ 0.028 * as.numeric(TPQ_Q$QO24)+ 0.059 * as.numeric(TPQ_Q$QO25)+ -0.039 * as.numeric(TPQ_Q$QO26)+ 0.019 * as.numeric(TPQ_Q$QO27)+ -0.264 * as.numeric(TPQ_Q$QO28)+ -0.007 * as.numeric(TPQ_Q$QO29)+ -0.148 * as.numeric(TPQ_Q$QO30)+ -0.008 * as.numeric(TPQ_Q$QO31)+ 0.047 * as.numeric(TPQ_Q$QO32)+ 0.098 * as.numeric(TPQ_Q$QO33)+ -0.053 * as.numeric(TPQ_Q$QO34)+ 0.222 * as.numeric(TPQ_Q$QO35)+ -0.032 * as.numeric(TPQ_Q$QO36)+ 0.074 * as.numeric(TPQ_Q$QO37)+ 0.066 * as.numeric(TPQ_Q$QO38)+ -0.571 * as.numeric(TPQ_Q$QO39)+ 0.162 * as.numeric(TPQ_Q$QO40)+ -0.586 * as.numeric(TPQ_Q$QO41)+ 0.104 * as.numeric(TPQ_Q$QO42)+ 0.106 * as.numeric(TPQ_Q$QO43)+ 0.131 * as.numeric(TPQ_Q$QO44)+ 0.4 * as.numeric(TPQ_Q$QO45)+ -0.026 * as.numeric(TPQ_Q$QO46)+ 0.001 * as.numeric(TPQ_Q$QO47)+ -0.087 * as.numeric(TPQ_Q$QO48)+ 0.078 * as.numeric(TPQ_Q$QO49)+ 0.014 * as.numeric(TPQ_Q$QO50)+ 0.007 * as.numeric(TPQ_Q$QO51)+ -0.004 * as.numeric(TPQ_Q$QO52)+ -0.012 * as.numeric(TPQ_Q$QO53)+ 0.023 * as.numeric(TPQ_Q$QO54)+ -0.05 * as.numeric(TPQ_Q$QO55)+ -0.014 * as.numeric(TPQ_Q$QO56)+ 0.067 * as.numeric(TPQ_Q$QO57)+ 0.04 * as.numeric(TPQ_Q$QO58)+ -0.061 * as.numeric(TPQ_Q$QO59)+ 0.038 * as.numeric(TPQ_Q$QO60)+ -0.036 * as.numeric(TPQ_Q$QO61)+ 0.022 * as.numeric(TPQ_Q$QO62)+ -0.236 * as.numeric(TPQ_Q$QO63)+ -0.029 * as.numeric(TPQ_Q$QO64)+ -0.056 * as.numeric(TPQ_Q$QO65)+ -0.04 * as.numeric(TPQ_Q$QO66)+ 0.079 * as.numeric(TPQ_Q$QO67)+ -0.153 * as.numeric(TPQ_Q$QO68)+ -0.102 * as.numeric(TPQ_Q$QO69)+ -0.07 * as.numeric(TPQ_Q$QO70)+ 0.047 * as.numeric(TPQ_Q$QO71)+ -0.074 * as.numeric(TPQ_Q$QO72)+ 0.085 * as.numeric(TPQ_Q$QO73)+ 0.067 * as.numeric(TPQ_Q$QO74)+ -0.075 * as.numeric(TPQ_Q$QO75)+ -0.041 * as.numeric(TPQ_Q$QO76)+ -0.021 * as.numeric(TPQ_Q$QO77)+ -0.028 * as.numeric(TPQ_Q$QO78)+ 0.097 * as.numeric(TPQ_Q$QO79)+ 0.071 * as.numeric(TPQ_Q$QO80)+ 0.087 * as.numeric(TPQ_Q$QO81)+ -0.017 * as.numeric(TPQ_Q$QO82)+ -0.047 * as.numeric(TPQ_Q$QO83)+ 0.017 * as.numeric(TPQ_Q$QO84)+ 0.095 * as.numeric(TPQ_Q$QO85)+ 0.077 * as.numeric(TPQ_Q$QO86)+ 0.018 * as.numeric(TPQ_Q$QO87)+ -0.047 * as.numeric(TPQ_Q$QO88)+ -0.047 * as.numeric(TPQ_Q$QO89)+ -0.095 * as.numeric(TPQ_Q$QO90)+ 0.097 * as.numeric(TPQ_Q$QO91)+ -0.483 * as.numeric(TPQ_Q$QO92)+ 0.167 * as.numeric(TPQ_Q$QO93)+ -0.033 * as.numeric(TPQ_Q$QO94)+ 0.006 * as.numeric(TPQ_Q$QO95)+ -0.005 * as.numeric(TPQ_Q$QO96)+ -0.127 * as.numeric(TPQ_Q$QO97)+ -0.064 * as.numeric(TPQ_Q$QO98)+ -0.125 * as.numeric(TPQ_Q$QO99)+ 0.065 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_3_HC =  -0.058 * as.numeric(TPQ_Q$QO1)+ 0.052 * as.numeric(TPQ_Q$QO2)+ -0.073 * as.numeric(TPQ_Q$QO3)+ -0.015 * as.numeric(TPQ_Q$QO4)+ 0.173 * as.numeric(TPQ_Q$QO5)+ -0.067 * as.numeric(TPQ_Q$QO6)+ -0.022 * as.numeric(TPQ_Q$QO7)+ -0.107 * as.numeric(TPQ_Q$QO8)+ 0.087 * as.numeric(TPQ_Q$QO9)+ -0.144 * as.numeric(TPQ_Q$QO10)+ 0.222 * as.numeric(TPQ_Q$QO11)+ -0.002 * as.numeric(TPQ_Q$QO12)+ 0.02 * as.numeric(TPQ_Q$QO13)+ -0.188 * as.numeric(TPQ_Q$QO14)+ 0.046 * as.numeric(TPQ_Q$QO15)+ 0.061 * as.numeric(TPQ_Q$QO16)+ 0.128 * as.numeric(TPQ_Q$QO17)+ -0.08 * as.numeric(TPQ_Q$QO18)+ -0.243 * as.numeric(TPQ_Q$QO19)+ 0.197 * as.numeric(TPQ_Q$QO20)+ 0.017 * as.numeric(TPQ_Q$QO21)+ -0.043 * as.numeric(TPQ_Q$QO22)+ -0.194 * as.numeric(TPQ_Q$QO23)+ 0.023 * as.numeric(TPQ_Q$QO24)+ 0.166 * as.numeric(TPQ_Q$QO25)+ -0.064 * as.numeric(TPQ_Q$QO26)+ -0.009 * as.numeric(TPQ_Q$QO27)+ 0.247 * as.numeric(TPQ_Q$QO28)+ 0.111 * as.numeric(TPQ_Q$QO29)+ -0.009 * as.numeric(TPQ_Q$QO30)+ 0.069 * as.numeric(TPQ_Q$QO31)+ 0.01 * as.numeric(TPQ_Q$QO32)+ -0.231 * as.numeric(TPQ_Q$QO33)+ 0.014 * as.numeric(TPQ_Q$QO34)+ 0.263 * as.numeric(TPQ_Q$QO35)+ 0.125 * as.numeric(TPQ_Q$QO36)+ -0.042 * as.numeric(TPQ_Q$QO37)+ -0.088 * as.numeric(TPQ_Q$QO38)+ 0.056 * as.numeric(TPQ_Q$QO39)+ -0.09 * as.numeric(TPQ_Q$QO40)+ 0.026 * as.numeric(TPQ_Q$QO41)+ -0.071 * as.numeric(TPQ_Q$QO42)+ -0.052 * as.numeric(TPQ_Q$QO43)+ 0.096 * as.numeric(TPQ_Q$QO44)+ -0.111 * as.numeric(TPQ_Q$QO45)+ 0.277 * as.numeric(TPQ_Q$QO46)+ -0.033 * as.numeric(TPQ_Q$QO47)+ -0.239 * as.numeric(TPQ_Q$QO48)+ -0.078 * as.numeric(TPQ_Q$QO49)+ -0.272 * as.numeric(TPQ_Q$QO50)+ -0.014 * as.numeric(TPQ_Q$QO51)+ 0.062 * as.numeric(TPQ_Q$QO52)+ 0.078 * as.numeric(TPQ_Q$QO53)+ 0.016 * as.numeric(TPQ_Q$QO54)+ 0.27 * as.numeric(TPQ_Q$QO55)+ 0.238 * as.numeric(TPQ_Q$QO56)+ 0.048 * as.numeric(TPQ_Q$QO57)+ 0.177 * as.numeric(TPQ_Q$QO58)+ -0.082 * as.numeric(TPQ_Q$QO59)+ 0.13 * as.numeric(TPQ_Q$QO60)+ -0.019 * as.numeric(TPQ_Q$QO61)+ 0.011 * as.numeric(TPQ_Q$QO62)+ -0.179 * as.numeric(TPQ_Q$QO63)+ 0.023 * as.numeric(TPQ_Q$QO64)+ -0.04 * as.numeric(TPQ_Q$QO65)+ -0.062 * as.numeric(TPQ_Q$QO66)+ -0.087 * as.numeric(TPQ_Q$QO67)+ -0.044 * as.numeric(TPQ_Q$QO68)+ -0.131 * as.numeric(TPQ_Q$QO69)+ 0.005 * as.numeric(TPQ_Q$QO70)+ -0.221 * as.numeric(TPQ_Q$QO71)+ 0.077 * as.numeric(TPQ_Q$QO72)+ 0.066 * as.numeric(TPQ_Q$QO73)+ -0.086 * as.numeric(TPQ_Q$QO74)+ 0.137 * as.numeric(TPQ_Q$QO75)+ -0.124 * as.numeric(TPQ_Q$QO76)+ -0.062 * as.numeric(TPQ_Q$QO77)+ -0.178 * as.numeric(TPQ_Q$QO78)+ -0.091 * as.numeric(TPQ_Q$QO79)+ 0.072 * as.numeric(TPQ_Q$QO80)+ 0.007 * as.numeric(TPQ_Q$QO81)+ -0.051 * as.numeric(TPQ_Q$QO82)+ -0.086 * as.numeric(TPQ_Q$QO83)+ 0.27 * as.numeric(TPQ_Q$QO84)+ 0.078 * as.numeric(TPQ_Q$QO85)+ -0.006 * as.numeric(TPQ_Q$QO86)+ -0.203 * as.numeric(TPQ_Q$QO87)+ 0.02 * as.numeric(TPQ_Q$QO88)+ 0.243 * as.numeric(TPQ_Q$QO89)+ -0.076 * as.numeric(TPQ_Q$QO90)+ -0.145 * as.numeric(TPQ_Q$QO91)+ -0.007 * as.numeric(TPQ_Q$QO92)+ -0.206 * as.numeric(TPQ_Q$QO93)+ -0.066 * as.numeric(TPQ_Q$QO94)+ 0.019 * as.numeric(TPQ_Q$QO95)+ -0.13 * as.numeric(TPQ_Q$QO96)+ 0.058 * as.numeric(TPQ_Q$QO97)+ -0.106 * as.numeric(TPQ_Q$QO98)+ 0.16 * as.numeric(TPQ_Q$QO99)+ 0.326 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_4_HC =  0.087 * as.numeric(TPQ_Q$QO1)+ -0.068 * as.numeric(TPQ_Q$QO2)+ 0.017 * as.numeric(TPQ_Q$QO3)+ -0.027 * as.numeric(TPQ_Q$QO4)+ 0.187 * as.numeric(TPQ_Q$QO5)+ 0.019 * as.numeric(TPQ_Q$QO6)+ 0.102 * as.numeric(TPQ_Q$QO7)+ 0.107 * as.numeric(TPQ_Q$QO8)+ 0.127 * as.numeric(TPQ_Q$QO9)+ 0.227 * as.numeric(TPQ_Q$QO10)+ 0.036 * as.numeric(TPQ_Q$QO11)+ 0.019 * as.numeric(TPQ_Q$QO12)+ -0.093 * as.numeric(TPQ_Q$QO13)+ 0.105 * as.numeric(TPQ_Q$QO14)+ -0.038 * as.numeric(TPQ_Q$QO15)+ 0.188 * as.numeric(TPQ_Q$QO16)+ -0.221 * as.numeric(TPQ_Q$QO17)+ 0.305 * as.numeric(TPQ_Q$QO18)+ 0.297 * as.numeric(TPQ_Q$QO19)+ -0.129 * as.numeric(TPQ_Q$QO20)+ 0.165 * as.numeric(TPQ_Q$QO21)+ 0.103 * as.numeric(TPQ_Q$QO22)+ 0.333 * as.numeric(TPQ_Q$QO23)+ 0.054 * as.numeric(TPQ_Q$QO24)+ -0.292 * as.numeric(TPQ_Q$QO25)+ 0.101 * as.numeric(TPQ_Q$QO26)+ 0.06 * as.numeric(TPQ_Q$QO27)+ -0.007 * as.numeric(TPQ_Q$QO28)+ 0.054 * as.numeric(TPQ_Q$QO29)+ 0.135 * as.numeric(TPQ_Q$QO30)+ 0.052 * as.numeric(TPQ_Q$QO31)+ 0.178 * as.numeric(TPQ_Q$QO32)+ 0.149 * as.numeric(TPQ_Q$QO33)+ 0.079 * as.numeric(TPQ_Q$QO34)+ 0.18 * as.numeric(TPQ_Q$QO35)+ 0.167 * as.numeric(TPQ_Q$QO36)+ 0.103 * as.numeric(TPQ_Q$QO37)+ 0.013 * as.numeric(TPQ_Q$QO38)+ 0.201 * as.numeric(TPQ_Q$QO39)+ 0.165 * as.numeric(TPQ_Q$QO40)+ 0.154 * as.numeric(TPQ_Q$QO41)+ 0.078 * as.numeric(TPQ_Q$QO42)+ 0.086 * as.numeric(TPQ_Q$QO43)+ 0.072 * as.numeric(TPQ_Q$QO44)+ 0.065 * as.numeric(TPQ_Q$QO45)+ 0.164 * as.numeric(TPQ_Q$QO46)+ 0.01 * as.numeric(TPQ_Q$QO47)+ 0.023 * as.numeric(TPQ_Q$QO48)+ 0.019 * as.numeric(TPQ_Q$QO49)+ 0.233 * as.numeric(TPQ_Q$QO50)+ -0.072 * as.numeric(TPQ_Q$QO51)+ 0.071 * as.numeric(TPQ_Q$QO52)+ 0.213 * as.numeric(TPQ_Q$QO53)+ 0.002 * as.numeric(TPQ_Q$QO54)+ 0.233 * as.numeric(TPQ_Q$QO55)+ 0.255 * as.numeric(TPQ_Q$QO56)+ 0.051 * as.numeric(TPQ_Q$QO57)+ -0.074 * as.numeric(TPQ_Q$QO58)+ 0.086 * as.numeric(TPQ_Q$QO59)+ 0.128 * as.numeric(TPQ_Q$QO60)+ 0.151 * as.numeric(TPQ_Q$QO61)+ 0.113 * as.numeric(TPQ_Q$QO62)+ 0.097 * as.numeric(TPQ_Q$QO63)+ 0.2 * as.numeric(TPQ_Q$QO64)+ 0.061 * as.numeric(TPQ_Q$QO65)+ -0.006 * as.numeric(TPQ_Q$QO66)+ 0.149 * as.numeric(TPQ_Q$QO67)+ -0.177 * as.numeric(TPQ_Q$QO68)+ -0.108 * as.numeric(TPQ_Q$QO69)+ 0.04 * as.numeric(TPQ_Q$QO70)+ -0.022 * as.numeric(TPQ_Q$QO71)+ 0.091 * as.numeric(TPQ_Q$QO72)+ 0.237 * as.numeric(TPQ_Q$QO73)+ 0.154 * as.numeric(TPQ_Q$QO74)+ 0.232 * as.numeric(TPQ_Q$QO75)+ 0.008 * as.numeric(TPQ_Q$QO76)+ -0.027 * as.numeric(TPQ_Q$QO77)+ 0.13 * as.numeric(TPQ_Q$QO78)+ 0.132 * as.numeric(TPQ_Q$QO79)+ 0.3 * as.numeric(TPQ_Q$QO80)+ 0.099 * as.numeric(TPQ_Q$QO81)+ 0.054 * as.numeric(TPQ_Q$QO82)+ 0.086 * as.numeric(TPQ_Q$QO83)+ 0.047 * as.numeric(TPQ_Q$QO84)+ 0.379 * as.numeric(TPQ_Q$QO85)+ 0.024 * as.numeric(TPQ_Q$QO86)+ 0.01 * as.numeric(TPQ_Q$QO87)+ 0.038 * as.numeric(TPQ_Q$QO88)+ -0.026 * as.numeric(TPQ_Q$QO89)+ -0.043 * as.numeric(TPQ_Q$QO90)+ -0.057 * as.numeric(TPQ_Q$QO91)+ 0.167 * as.numeric(TPQ_Q$QO92)+ 0.203 * as.numeric(TPQ_Q$QO93)+ 0.072 * as.numeric(TPQ_Q$QO94)+ -0.002 * as.numeric(TPQ_Q$QO95)+ 0.121 * as.numeric(TPQ_Q$QO96)+ 0.027 * as.numeric(TPQ_Q$QO97)+ 0.253 * as.numeric(TPQ_Q$QO98)+ 0.074 * as.numeric(TPQ_Q$QO99)+ 0.071 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_5_HC =  0.152 * as.numeric(TPQ_Q$QO1)+ 0.163 * as.numeric(TPQ_Q$QO2)+ -0.038 * as.numeric(TPQ_Q$QO3)+ 0.144 * as.numeric(TPQ_Q$QO4)+ -0.141 * as.numeric(TPQ_Q$QO5)+ 0.159 * as.numeric(TPQ_Q$QO6)+ 0.085 * as.numeric(TPQ_Q$QO7)+ -0.017 * as.numeric(TPQ_Q$QO8)+ -0.06 * as.numeric(TPQ_Q$QO9)+ -0.218 * as.numeric(TPQ_Q$QO10)+ -0.011 * as.numeric(TPQ_Q$QO11)+ -0.067 * as.numeric(TPQ_Q$QO12)+ -0.037 * as.numeric(TPQ_Q$QO13)+ -0.241 * as.numeric(TPQ_Q$QO14)+ -0.142 * as.numeric(TPQ_Q$QO15)+ 0.001 * as.numeric(TPQ_Q$QO16)+ 0.013 * as.numeric(TPQ_Q$QO17)+ -0.009 * as.numeric(TPQ_Q$QO18)+ -0.025 * as.numeric(TPQ_Q$QO19)+ -0.16 * as.numeric(TPQ_Q$QO20)+ 0.02 * as.numeric(TPQ_Q$QO21)+ 0.265 * as.numeric(TPQ_Q$QO22)+ -0.128 * as.numeric(TPQ_Q$QO23)+ -0.123 * as.numeric(TPQ_Q$QO24)+ -0.048 * as.numeric(TPQ_Q$QO25)+ -0.24 * as.numeric(TPQ_Q$QO26)+ 0.111 * as.numeric(TPQ_Q$QO27)+ -0.163 * as.numeric(TPQ_Q$QO28)+ -0.125 * as.numeric(TPQ_Q$QO29)+ -0.046 * as.numeric(TPQ_Q$QO30)+ 0.125 * as.numeric(TPQ_Q$QO31)+ 0.077 * as.numeric(TPQ_Q$QO32)+ -0.09 * as.numeric(TPQ_Q$QO33)+ 0.193 * as.numeric(TPQ_Q$QO34)+ -0.22 * as.numeric(TPQ_Q$QO35)+ -0.265 * as.numeric(TPQ_Q$QO36)+ -0.199 * as.numeric(TPQ_Q$QO37)+ -0.342 * as.numeric(TPQ_Q$QO38)+ -0.102 * as.numeric(TPQ_Q$QO39)+ -0.307 * as.numeric(TPQ_Q$QO40)+ -0.152 * as.numeric(TPQ_Q$QO41)+ 0.134 * as.numeric(TPQ_Q$QO42)+ 0.06 * as.numeric(TPQ_Q$QO43)+ 0.134 * as.numeric(TPQ_Q$QO44)+ -0.081 * as.numeric(TPQ_Q$QO45)+ -0.114 * as.numeric(TPQ_Q$QO46)+ -0.039 * as.numeric(TPQ_Q$QO47)+ 0.116 * as.numeric(TPQ_Q$QO48)+ -0.176 * as.numeric(TPQ_Q$QO49)+ 0.043 * as.numeric(TPQ_Q$QO50)+ 0.009 * as.numeric(TPQ_Q$QO51)+ 0.228 * as.numeric(TPQ_Q$QO52)+ -0.226 * as.numeric(TPQ_Q$QO53)+ -0.108 * as.numeric(TPQ_Q$QO54)+ -0.099 * as.numeric(TPQ_Q$QO55)+ -0.052 * as.numeric(TPQ_Q$QO56)+ 0.046 * as.numeric(TPQ_Q$QO57)+ 0.063 * as.numeric(TPQ_Q$QO58)+ -0.012 * as.numeric(TPQ_Q$QO59)+ -0.084 * as.numeric(TPQ_Q$QO60)+ 0.043 * as.numeric(TPQ_Q$QO61)+ -0.024 * as.numeric(TPQ_Q$QO62)+ 0.054 * as.numeric(TPQ_Q$QO63)+ 0.131 * as.numeric(TPQ_Q$QO64)+ -0.037 * as.numeric(TPQ_Q$QO65)+ 0.002 * as.numeric(TPQ_Q$QO66)+ 0.104 * as.numeric(TPQ_Q$QO67)+ -0.168 * as.numeric(TPQ_Q$QO68)+ -0.159 * as.numeric(TPQ_Q$QO69)+ -0.085 * as.numeric(TPQ_Q$QO70)+ 0.071 * as.numeric(TPQ_Q$QO71)+ -0.06 * as.numeric(TPQ_Q$QO72)+ -0.277 * as.numeric(TPQ_Q$QO73)+ 0.1 * as.numeric(TPQ_Q$QO74)+ 0.248 * as.numeric(TPQ_Q$QO75)+ -0.135 * as.numeric(TPQ_Q$QO76)+ -0.064 * as.numeric(TPQ_Q$QO77)+ -0.229 * as.numeric(TPQ_Q$QO78)+ 0.288 * as.numeric(TPQ_Q$QO79)+ 0.268 * as.numeric(TPQ_Q$QO80)+ 0.083 * as.numeric(TPQ_Q$QO81)+ 0.161 * as.numeric(TPQ_Q$QO82)+ 0.116 * as.numeric(TPQ_Q$QO83)+ 0.09 * as.numeric(TPQ_Q$QO84)+ -0.129 * as.numeric(TPQ_Q$QO85)+ -0.0 * as.numeric(TPQ_Q$QO86)+ -0.222 * as.numeric(TPQ_Q$QO87)+ 0.02 * as.numeric(TPQ_Q$QO88)+ 0.231 * as.numeric(TPQ_Q$QO89)+ -0.293 * as.numeric(TPQ_Q$QO90)+ 0.108 * as.numeric(TPQ_Q$QO91)+ 0.057 * as.numeric(TPQ_Q$QO92)+ -0.047 * as.numeric(TPQ_Q$QO93)+ 0.136 * as.numeric(TPQ_Q$QO94)+ 0.14 * as.numeric(TPQ_Q$QO95)+ -0.16 * as.numeric(TPQ_Q$QO96)+ 0.173 * as.numeric(TPQ_Q$QO97)+ 0.187 * as.numeric(TPQ_Q$QO98)+ 0.043 * as.numeric(TPQ_Q$QO99)+ 0.087 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_6_HC =  -0.006 * as.numeric(TPQ_Q$QO1)+ -0.009 * as.numeric(TPQ_Q$QO2)+ -0.684 * as.numeric(TPQ_Q$QO3)+ -0.052 * as.numeric(TPQ_Q$QO4)+ 0.166 * as.numeric(TPQ_Q$QO5)+ -0.214 * as.numeric(TPQ_Q$QO6)+ -0.162 * as.numeric(TPQ_Q$QO7)+ -0.023 * as.numeric(TPQ_Q$QO8)+ 0.012 * as.numeric(TPQ_Q$QO9)+ 0.162 * as.numeric(TPQ_Q$QO10)+ 0.11 * as.numeric(TPQ_Q$QO11)+ 0.685 * as.numeric(TPQ_Q$QO12)+ 0.037 * as.numeric(TPQ_Q$QO13)+ 0.17 * as.numeric(TPQ_Q$QO14)+ 0.387 * as.numeric(TPQ_Q$QO15)+ 0.095 * as.numeric(TPQ_Q$QO16)+ 0.196 * as.numeric(TPQ_Q$QO17)+ 0.064 * as.numeric(TPQ_Q$QO18)+ 0.206 * as.numeric(TPQ_Q$QO19)+ 0.195 * as.numeric(TPQ_Q$QO20)+ 0.031 * as.numeric(TPQ_Q$QO21)+ 0.072 * as.numeric(TPQ_Q$QO22)+ 0.142 * as.numeric(TPQ_Q$QO23)+ 0.048 * as.numeric(TPQ_Q$QO24)+ 0.207 * as.numeric(TPQ_Q$QO25)+ 0.03 * as.numeric(TPQ_Q$QO26)+ -0.069 * as.numeric(TPQ_Q$QO27)+ -0.002 * as.numeric(TPQ_Q$QO28)+ 0.012 * as.numeric(TPQ_Q$QO29)+ -0.057 * as.numeric(TPQ_Q$QO30)+ -0.071 * as.numeric(TPQ_Q$QO31)+ 0.097 * as.numeric(TPQ_Q$QO32)+ 0.159 * as.numeric(TPQ_Q$QO33)+ -0.006 * as.numeric(TPQ_Q$QO34)+ -0.014 * as.numeric(TPQ_Q$QO35)+ 0.131 * as.numeric(TPQ_Q$QO36)+ 0.221 * as.numeric(TPQ_Q$QO37)+ 0.121 * as.numeric(TPQ_Q$QO38)+ 0.093 * as.numeric(TPQ_Q$QO39)+ 0.082 * as.numeric(TPQ_Q$QO40)+ 0.032 * as.numeric(TPQ_Q$QO41)+ -0.016 * as.numeric(TPQ_Q$QO42)+ 0.078 * as.numeric(TPQ_Q$QO43)+ -0.047 * as.numeric(TPQ_Q$QO44)+ 0.007 * as.numeric(TPQ_Q$QO45)+ -0.002 * as.numeric(TPQ_Q$QO46)+ -0.032 * as.numeric(TPQ_Q$QO47)+ -0.015 * as.numeric(TPQ_Q$QO48)+ 0.018 * as.numeric(TPQ_Q$QO49)+ -0.004 * as.numeric(TPQ_Q$QO50)+ -0.019 * as.numeric(TPQ_Q$QO51)+ -0.017 * as.numeric(TPQ_Q$QO52)+ 0.093 * as.numeric(TPQ_Q$QO53)+ 0.049 * as.numeric(TPQ_Q$QO54)+ -0.024 * as.numeric(TPQ_Q$QO55)+ 0.078 * as.numeric(TPQ_Q$QO56)+ -0.001 * as.numeric(TPQ_Q$QO57)+ 0.038 * as.numeric(TPQ_Q$QO58)+ -0.041 * as.numeric(TPQ_Q$QO59)+ 0.017 * as.numeric(TPQ_Q$QO60)+ -0.105 * as.numeric(TPQ_Q$QO61)+ -0.062 * as.numeric(TPQ_Q$QO62)+ 0.066 * as.numeric(TPQ_Q$QO63)+ -0.102 * as.numeric(TPQ_Q$QO64)+ -0.035 * as.numeric(TPQ_Q$QO65)+ 0.037 * as.numeric(TPQ_Q$QO66)+ -0.109 * as.numeric(TPQ_Q$QO67)+ 0.023 * as.numeric(TPQ_Q$QO68)+ 0.022 * as.numeric(TPQ_Q$QO69)+ 0.008 * as.numeric(TPQ_Q$QO70)+ 0.071 * as.numeric(TPQ_Q$QO71)+ 0.067 * as.numeric(TPQ_Q$QO72)+ 0.073 * as.numeric(TPQ_Q$QO73)+ -0.298 * as.numeric(TPQ_Q$QO74)+ -0.007 * as.numeric(TPQ_Q$QO75)+ 0.08 * as.numeric(TPQ_Q$QO76)+ -0.06 * as.numeric(TPQ_Q$QO77)+ 0.047 * as.numeric(TPQ_Q$QO78)+ 0.06 * as.numeric(TPQ_Q$QO79)+ 0.012 * as.numeric(TPQ_Q$QO80)+ -0.025 * as.numeric(TPQ_Q$QO81)+ -0.022 * as.numeric(TPQ_Q$QO82)+ -0.146 * as.numeric(TPQ_Q$QO83)+ -0.086 * as.numeric(TPQ_Q$QO84)+ -0.146 * as.numeric(TPQ_Q$QO85)+ 0.551 * as.numeric(TPQ_Q$QO86)+ 0.003 * as.numeric(TPQ_Q$QO87)+ 0.777 * as.numeric(TPQ_Q$QO88)+ -0.089 * as.numeric(TPQ_Q$QO89)+ 0.4 * as.numeric(TPQ_Q$QO90)+ -0.028 * as.numeric(TPQ_Q$QO91)+ 0.019 * as.numeric(TPQ_Q$QO92)+ 0.162 * as.numeric(TPQ_Q$QO93)+ -0.036 * as.numeric(TPQ_Q$QO94)+ -0.02 * as.numeric(TPQ_Q$QO95)+ 0.087 * as.numeric(TPQ_Q$QO96)+ -0.034 * as.numeric(TPQ_Q$QO97)+ -0.107 * as.numeric(TPQ_Q$QO98)+ 0.016 * as.numeric(TPQ_Q$QO99)+ -0.08 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_7_HC =  -0.294 * as.numeric(TPQ_Q$QO1)+ -0.016 * as.numeric(TPQ_Q$QO2)+ -0.042 * as.numeric(TPQ_Q$QO3)+ 0.021 * as.numeric(TPQ_Q$QO4)+ 0.293 * as.numeric(TPQ_Q$QO5)+ -0.264 * as.numeric(TPQ_Q$QO6)+ 0.159 * as.numeric(TPQ_Q$QO7)+ -0.551 * as.numeric(TPQ_Q$QO8)+ -0.007 * as.numeric(TPQ_Q$QO9)+ 0.208 * as.numeric(TPQ_Q$QO10)+ 0.088 * as.numeric(TPQ_Q$QO11)+ -0.041 * as.numeric(TPQ_Q$QO12)+ 0.003 * as.numeric(TPQ_Q$QO13)+ 0.141 * as.numeric(TPQ_Q$QO14)+ -0.33 * as.numeric(TPQ_Q$QO15)+ -0.01 * as.numeric(TPQ_Q$QO16)+ -0.061 * as.numeric(TPQ_Q$QO17)+ 0.327 * as.numeric(TPQ_Q$QO18)+ 0.388 * as.numeric(TPQ_Q$QO19)+ -0.209 * as.numeric(TPQ_Q$QO20)+ -0.042 * as.numeric(TPQ_Q$QO21)+ 0.035 * as.numeric(TPQ_Q$QO22)+ 0.388 * as.numeric(TPQ_Q$QO23)+ 0.091 * as.numeric(TPQ_Q$QO24)+ -0.202 * as.numeric(TPQ_Q$QO25)+ -0.303 * as.numeric(TPQ_Q$QO26)+ 0.042 * as.numeric(TPQ_Q$QO27)+ 0.277 * as.numeric(TPQ_Q$QO28)+ -0.123 * as.numeric(TPQ_Q$QO29)+ 0.427 * as.numeric(TPQ_Q$QO30)+ -0.028 * as.numeric(TPQ_Q$QO31)+ -0.261 * as.numeric(TPQ_Q$QO32)+ 0.259 * as.numeric(TPQ_Q$QO33)+ 0.215 * as.numeric(TPQ_Q$QO34)+ -0.262 * as.numeric(TPQ_Q$QO35)+ -0.157 * as.numeric(TPQ_Q$QO36)+ 0.184 * as.numeric(TPQ_Q$QO37)+ 0.031 * as.numeric(TPQ_Q$QO38)+ -0.005 * as.numeric(TPQ_Q$QO39)+ 0.018 * as.numeric(TPQ_Q$QO40)+ 0.008 * as.numeric(TPQ_Q$QO41)+ -0.203 * as.numeric(TPQ_Q$QO42)+ 0.185 * as.numeric(TPQ_Q$QO43)+ -0.237 * as.numeric(TPQ_Q$QO44)+ -0.004 * as.numeric(TPQ_Q$QO45)+ 0.027 * as.numeric(TPQ_Q$QO46)+ 0.019 * as.numeric(TPQ_Q$QO47)+ -0.071 * as.numeric(TPQ_Q$QO48)+ -0.032 * as.numeric(TPQ_Q$QO49)+ 0.213 * as.numeric(TPQ_Q$QO50)+ 0.049 * as.numeric(TPQ_Q$QO51)+ 0.017 * as.numeric(TPQ_Q$QO52)+ 0.099 * as.numeric(TPQ_Q$QO53)+ 0.099 * as.numeric(TPQ_Q$QO54)+ -0.047 * as.numeric(TPQ_Q$QO55)+ -0.033 * as.numeric(TPQ_Q$QO56)+ -0.001 * as.numeric(TPQ_Q$QO57)+ -0.083 * as.numeric(TPQ_Q$QO58)+ -0.02 * as.numeric(TPQ_Q$QO59)+ -0.072 * as.numeric(TPQ_Q$QO60)+ 0.216 * as.numeric(TPQ_Q$QO61)+ -0.044 * as.numeric(TPQ_Q$QO62)+ -0.299 * as.numeric(TPQ_Q$QO63)+ 0.175 * as.numeric(TPQ_Q$QO64)+ -0.121 * as.numeric(TPQ_Q$QO65)+ -0.024 * as.numeric(TPQ_Q$QO66)+ 0.018 * as.numeric(TPQ_Q$QO67)+ 0.067 * as.numeric(TPQ_Q$QO68)+ 0.073 * as.numeric(TPQ_Q$QO69)+ -0.031 * as.numeric(TPQ_Q$QO70)+ -0.017 * as.numeric(TPQ_Q$QO71)+ 0.012 * as.numeric(TPQ_Q$QO72)+ 0.292 * as.numeric(TPQ_Q$QO73)+ 0.026 * as.numeric(TPQ_Q$QO74)+ -0.212 * as.numeric(TPQ_Q$QO75)+ 0.091 * as.numeric(TPQ_Q$QO76)+ -0.003 * as.numeric(TPQ_Q$QO77)+ 0.017 * as.numeric(TPQ_Q$QO78)+ -0.169 * as.numeric(TPQ_Q$QO79)+ -0.234 * as.numeric(TPQ_Q$QO80)+ 0.085 * as.numeric(TPQ_Q$QO81)+ -0.097 * as.numeric(TPQ_Q$QO82)+ 0.051 * as.numeric(TPQ_Q$QO83)+ -0.378 * as.numeric(TPQ_Q$QO84)+ 0.082 * as.numeric(TPQ_Q$QO85)+ -0.115 * as.numeric(TPQ_Q$QO86)+ 0.075 * as.numeric(TPQ_Q$QO87)+ -0.039 * as.numeric(TPQ_Q$QO88)+ -0.163 * as.numeric(TPQ_Q$QO89)+ -0.225 * as.numeric(TPQ_Q$QO90)+ -0.418 * as.numeric(TPQ_Q$QO91)+ -0.111 * as.numeric(TPQ_Q$QO92)+ -0.033 * as.numeric(TPQ_Q$QO93)+ 0.035 * as.numeric(TPQ_Q$QO94)+ -0.091 * as.numeric(TPQ_Q$QO95)+ -0.045 * as.numeric(TPQ_Q$QO96)+ -0.108 * as.numeric(TPQ_Q$QO97)+ -0.364 * as.numeric(TPQ_Q$QO98)+ 0.002 * as.numeric(TPQ_Q$QO99)+ -0.185 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_8_HC =  -0.055 * as.numeric(TPQ_Q$QO1)+ 0.017 * as.numeric(TPQ_Q$QO2)+ 0.013 * as.numeric(TPQ_Q$QO3)+ 0.011 * as.numeric(TPQ_Q$QO4)+ 0.016 * as.numeric(TPQ_Q$QO5)+ 0.021 * as.numeric(TPQ_Q$QO6)+ 0.038 * as.numeric(TPQ_Q$QO7)+ -0.038 * as.numeric(TPQ_Q$QO8)+ -0.026 * as.numeric(TPQ_Q$QO9)+ 0.045 * as.numeric(TPQ_Q$QO10)+ -0.073 * as.numeric(TPQ_Q$QO11)+ -0.062 * as.numeric(TPQ_Q$QO12)+ 0.044 * as.numeric(TPQ_Q$QO13)+ 0.043 * as.numeric(TPQ_Q$QO14)+ -0.136 * as.numeric(TPQ_Q$QO15)+ -0.154 * as.numeric(TPQ_Q$QO16)+ -0.023 * as.numeric(TPQ_Q$QO17)+ 0.023 * as.numeric(TPQ_Q$QO18)+ -0.031 * as.numeric(TPQ_Q$QO19)+ -0.008 * as.numeric(TPQ_Q$QO20)+ -0.087 * as.numeric(TPQ_Q$QO21)+ -0.013 * as.numeric(TPQ_Q$QO22)+ -0.05 * as.numeric(TPQ_Q$QO23)+ 0.097 * as.numeric(TPQ_Q$QO24)+ -0.014 * as.numeric(TPQ_Q$QO25)+ -0.038 * as.numeric(TPQ_Q$QO26)+ -0.008 * as.numeric(TPQ_Q$QO27)+ -0.077 * as.numeric(TPQ_Q$QO28)+ -0.061 * as.numeric(TPQ_Q$QO29)+ -0.053 * as.numeric(TPQ_Q$QO30)+ -0.01 * as.numeric(TPQ_Q$QO31)+ -0.101 * as.numeric(TPQ_Q$QO32)+ -0.101 * as.numeric(TPQ_Q$QO33)+ -0.035 * as.numeric(TPQ_Q$QO34)+ -0.023 * as.numeric(TPQ_Q$QO35)+ 0.004 * as.numeric(TPQ_Q$QO36)+ -0.074 * as.numeric(TPQ_Q$QO37)+ 0.03 * as.numeric(TPQ_Q$QO38)+ -0.056 * as.numeric(TPQ_Q$QO39)+ -0.006 * as.numeric(TPQ_Q$QO40)+ 0.003 * as.numeric(TPQ_Q$QO41)+ -0.01 * as.numeric(TPQ_Q$QO42)+ 0.065 * as.numeric(TPQ_Q$QO43)+ 0.028 * as.numeric(TPQ_Q$QO44)+ -0.029 * as.numeric(TPQ_Q$QO45)+ -0.068 * as.numeric(TPQ_Q$QO46)+ -0.054 * as.numeric(TPQ_Q$QO47)+ 0.053 * as.numeric(TPQ_Q$QO48)+ -0.003 * as.numeric(TPQ_Q$QO49)+ 0.106 * as.numeric(TPQ_Q$QO50)+ 0.036 * as.numeric(TPQ_Q$QO51)+ -0.048 * as.numeric(TPQ_Q$QO52)+ -0.076 * as.numeric(TPQ_Q$QO53)+ 0.005 * as.numeric(TPQ_Q$QO54)+ -0.089 * as.numeric(TPQ_Q$QO55)+ -0.12 * as.numeric(TPQ_Q$QO56)+ 0.025 * as.numeric(TPQ_Q$QO57)+ -0.011 * as.numeric(TPQ_Q$QO58)+ 0.027 * as.numeric(TPQ_Q$QO59)+ -0.03 * as.numeric(TPQ_Q$QO60)+ 0.012 * as.numeric(TPQ_Q$QO61)+ -0.037 * as.numeric(TPQ_Q$QO62)+ -0.026 * as.numeric(TPQ_Q$QO63)+ 0.053 * as.numeric(TPQ_Q$QO64)+ 0.012 * as.numeric(TPQ_Q$QO65)+ -0.823 * as.numeric(TPQ_Q$QO66)+ 0.004 * as.numeric(TPQ_Q$QO67)+ 0.021 * as.numeric(TPQ_Q$QO68)+ 0.113 * as.numeric(TPQ_Q$QO69)+ 0.434 * as.numeric(TPQ_Q$QO70)+ 0.059 * as.numeric(TPQ_Q$QO71)+ 0.646 * as.numeric(TPQ_Q$QO72)+ 0.078 * as.numeric(TPQ_Q$QO73)+ 0.052 * as.numeric(TPQ_Q$QO74)+ -0.134 * as.numeric(TPQ_Q$QO75)+ -0.304 * as.numeric(TPQ_Q$QO76)+ -0.03 * as.numeric(TPQ_Q$QO77)+ -0.113 * as.numeric(TPQ_Q$QO78)+ 0.057 * as.numeric(TPQ_Q$QO79)+ -0.159 * as.numeric(TPQ_Q$QO80)+ -0.027 * as.numeric(TPQ_Q$QO81)+ -0.099 * as.numeric(TPQ_Q$QO82)+ 0.017 * as.numeric(TPQ_Q$QO83)+ -0.007 * as.numeric(TPQ_Q$QO84)+ -0.069 * as.numeric(TPQ_Q$QO85)+ -0.118 * as.numeric(TPQ_Q$QO86)+ -0.256 * as.numeric(TPQ_Q$QO87)+ 0.029 * as.numeric(TPQ_Q$QO88)+ -0.016 * as.numeric(TPQ_Q$QO89)+ -0.095 * as.numeric(TPQ_Q$QO90)+ -0.074 * as.numeric(TPQ_Q$QO91)+ -0.217 * as.numeric(TPQ_Q$QO92)+ -0.026 * as.numeric(TPQ_Q$QO93)+ -0.019 * as.numeric(TPQ_Q$QO94)+ -0.034 * as.numeric(TPQ_Q$QO95)+ -0.004 * as.numeric(TPQ_Q$QO96)+ -0.074 * as.numeric(TPQ_Q$QO97)+ -0.025 * as.numeric(TPQ_Q$QO98)+ -0.027 * as.numeric(TPQ_Q$QO99)+ -0.022 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_9_HC =  -0.096 * as.numeric(TPQ_Q$QO1)+ -0.238 * as.numeric(TPQ_Q$QO2)+ -0.047 * as.numeric(TPQ_Q$QO3)+ -0.25 * as.numeric(TPQ_Q$QO4)+ -0.006 * as.numeric(TPQ_Q$QO5)+ 0.015 * as.numeric(TPQ_Q$QO6)+ -0.106 * as.numeric(TPQ_Q$QO7)+ -0.126 * as.numeric(TPQ_Q$QO8)+ 0.047 * as.numeric(TPQ_Q$QO9)+ 0.022 * as.numeric(TPQ_Q$QO10)+ 0.018 * as.numeric(TPQ_Q$QO11)+ -0.001 * as.numeric(TPQ_Q$QO12)+ -0.259 * as.numeric(TPQ_Q$QO13)+ 0.005 * as.numeric(TPQ_Q$QO14)+ -0.045 * as.numeric(TPQ_Q$QO15)+ 0.117 * as.numeric(TPQ_Q$QO16)+ -0.201 * as.numeric(TPQ_Q$QO17)+ 0.161 * as.numeric(TPQ_Q$QO18)+ 0.026 * as.numeric(TPQ_Q$QO19)+ -0.149 * as.numeric(TPQ_Q$QO20)+ 0.092 * as.numeric(TPQ_Q$QO21)+ -0.139 * as.numeric(TPQ_Q$QO22)+ -0.007 * as.numeric(TPQ_Q$QO23)+ -0.239 * as.numeric(TPQ_Q$QO24)+ -0.101 * as.numeric(TPQ_Q$QO25)+ -0.369 * as.numeric(TPQ_Q$QO26)+ 0.029 * as.numeric(TPQ_Q$QO27)+ -0.129 * as.numeric(TPQ_Q$QO28)+ -0.682 * as.numeric(TPQ_Q$QO29)+ -0.072 * as.numeric(TPQ_Q$QO30)+ -0.054 * as.numeric(TPQ_Q$QO31)+ 0.057 * as.numeric(TPQ_Q$QO32)+ 0.096 * as.numeric(TPQ_Q$QO33)+ -0.012 * as.numeric(TPQ_Q$QO34)+ 0.096 * as.numeric(TPQ_Q$QO35)+ -0.086 * as.numeric(TPQ_Q$QO36)+ 0.032 * as.numeric(TPQ_Q$QO37)+ -0.009 * as.numeric(TPQ_Q$QO38)+ -0.056 * as.numeric(TPQ_Q$QO39)+ 0.081 * as.numeric(TPQ_Q$QO40)+ -0.098 * as.numeric(TPQ_Q$QO41)+ 0.078 * as.numeric(TPQ_Q$QO42)+ -0.188 * as.numeric(TPQ_Q$QO43)+ 0.005 * as.numeric(TPQ_Q$QO44)+ -0.001 * as.numeric(TPQ_Q$QO45)+ 0.098 * as.numeric(TPQ_Q$QO46)+ -0.804 * as.numeric(TPQ_Q$QO47)+ -0.215 * as.numeric(TPQ_Q$QO48)+ -0.035 * as.numeric(TPQ_Q$QO49)+ -0.13 * as.numeric(TPQ_Q$QO50)+ -0.783 * as.numeric(TPQ_Q$QO51)+ -0.007 * as.numeric(TPQ_Q$QO52)+ 0.12 * as.numeric(TPQ_Q$QO53)+ 0.009 * as.numeric(TPQ_Q$QO54)+ 0.038 * as.numeric(TPQ_Q$QO55)+ -0.064 * as.numeric(TPQ_Q$QO56)+ 0.049 * as.numeric(TPQ_Q$QO57)+ -0.042 * as.numeric(TPQ_Q$QO58)+ -0.084 * as.numeric(TPQ_Q$QO59)+ -0.141 * as.numeric(TPQ_Q$QO60)+ 0.003 * as.numeric(TPQ_Q$QO61)+ -0.152 * as.numeric(TPQ_Q$QO62)+ -0.108 * as.numeric(TPQ_Q$QO63)+ -0.024 * as.numeric(TPQ_Q$QO64)+ -0.037 * as.numeric(TPQ_Q$QO65)+ 0.027 * as.numeric(TPQ_Q$QO66)+ -0.011 * as.numeric(TPQ_Q$QO67)+ 0.032 * as.numeric(TPQ_Q$QO68)+ 0.082 * as.numeric(TPQ_Q$QO69)+ -0.084 * as.numeric(TPQ_Q$QO70)+ 0.041 * as.numeric(TPQ_Q$QO71)+ -0.086 * as.numeric(TPQ_Q$QO72)+ 0.035 * as.numeric(TPQ_Q$QO73)+ -0.051 * as.numeric(TPQ_Q$QO74)+ -0.195 * as.numeric(TPQ_Q$QO75)+ 0.063 * as.numeric(TPQ_Q$QO76)+ -0.123 * as.numeric(TPQ_Q$QO77)+ 0.034 * as.numeric(TPQ_Q$QO78)+ -0.135 * as.numeric(TPQ_Q$QO79)+ -0.14 * as.numeric(TPQ_Q$QO80)+ -0.017 * as.numeric(TPQ_Q$QO81)+ -0.006 * as.numeric(TPQ_Q$QO82)+ -0.046 * as.numeric(TPQ_Q$QO83)+ -0.082 * as.numeric(TPQ_Q$QO84)+ 0.299 * as.numeric(TPQ_Q$QO85)+ -0.023 * as.numeric(TPQ_Q$QO86)+ 0.016 * as.numeric(TPQ_Q$QO87)+ 0.029 * as.numeric(TPQ_Q$QO88)+ -0.056 * as.numeric(TPQ_Q$QO89)+ 0.054 * as.numeric(TPQ_Q$QO90)+ -0.103 * as.numeric(TPQ_Q$QO91)+ -0.028 * as.numeric(TPQ_Q$QO92)+ 0.094 * as.numeric(TPQ_Q$QO93)+ 0.007 * as.numeric(TPQ_Q$QO94)+ -0.005 * as.numeric(TPQ_Q$QO95)+ 0.101 * as.numeric(TPQ_Q$QO96)+ -0.075 * as.numeric(TPQ_Q$QO97)+ -0.021 * as.numeric(TPQ_Q$QO98)+ 0.046 * as.numeric(TPQ_Q$QO99)+ 0.055 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_10_HC =  -0.11 * as.numeric(TPQ_Q$QO1)+ -0.032 * as.numeric(TPQ_Q$QO2)+ -0.131 * as.numeric(TPQ_Q$QO3)+ 0.056 * as.numeric(TPQ_Q$QO4)+ -0.068 * as.numeric(TPQ_Q$QO5)+ 0.052 * as.numeric(TPQ_Q$QO6)+ 0.003 * as.numeric(TPQ_Q$QO7)+ -0.095 * as.numeric(TPQ_Q$QO8)+ 0.1 * as.numeric(TPQ_Q$QO9)+ -0.094 * as.numeric(TPQ_Q$QO10)+ -0.49 * as.numeric(TPQ_Q$QO11)+ -0.096 * as.numeric(TPQ_Q$QO12)+ -0.216 * as.numeric(TPQ_Q$QO13)+ -0.129 * as.numeric(TPQ_Q$QO14)+ -0.163 * as.numeric(TPQ_Q$QO15)+ -0.14 * as.numeric(TPQ_Q$QO16)+ -0.449 * as.numeric(TPQ_Q$QO17)+ -0.047 * as.numeric(TPQ_Q$QO18)+ -0.214 * as.numeric(TPQ_Q$QO19)+ -0.222 * as.numeric(TPQ_Q$QO20)+ -0.234 * as.numeric(TPQ_Q$QO21)+ -0.119 * as.numeric(TPQ_Q$QO22)+ -0.18 * as.numeric(TPQ_Q$QO23)+ 0.092 * as.numeric(TPQ_Q$QO24)+ -0.458 * as.numeric(TPQ_Q$QO25)+ -0.059 * as.numeric(TPQ_Q$QO26)+ 0.057 * as.numeric(TPQ_Q$QO27)+ -0.073 * as.numeric(TPQ_Q$QO28)+ 0.003 * as.numeric(TPQ_Q$QO29)+ -0.038 * as.numeric(TPQ_Q$QO30)+ -0.005 * as.numeric(TPQ_Q$QO31)+ 0.03 * as.numeric(TPQ_Q$QO32)+ -0.136 * as.numeric(TPQ_Q$QO33)+ -0.107 * as.numeric(TPQ_Q$QO34)+ 0.085 * as.numeric(TPQ_Q$QO35)+ -0.056 * as.numeric(TPQ_Q$QO36)+ -0.217 * as.numeric(TPQ_Q$QO37)+ -0.148 * as.numeric(TPQ_Q$QO38)+ -0.091 * as.numeric(TPQ_Q$QO39)+ -0.197 * as.numeric(TPQ_Q$QO40)+ -0.086 * as.numeric(TPQ_Q$QO41)+ -0.004 * as.numeric(TPQ_Q$QO42)+ -0.039 * as.numeric(TPQ_Q$QO43)+ -0.0 * as.numeric(TPQ_Q$QO44)+ -0.132 * as.numeric(TPQ_Q$QO45)+ -0.089 * as.numeric(TPQ_Q$QO46)+ -0.002 * as.numeric(TPQ_Q$QO47)+ -0.174 * as.numeric(TPQ_Q$QO48)+ -0.138 * as.numeric(TPQ_Q$QO49)+ -0.185 * as.numeric(TPQ_Q$QO50)+ 0.002 * as.numeric(TPQ_Q$QO51)+ -0.119 * as.numeric(TPQ_Q$QO52)+ 0.093 * as.numeric(TPQ_Q$QO53)+ -0.002 * as.numeric(TPQ_Q$QO54)+ -0.056 * as.numeric(TPQ_Q$QO55)+ -0.035 * as.numeric(TPQ_Q$QO56)+ -0.11 * as.numeric(TPQ_Q$QO57)+ -0.205 * as.numeric(TPQ_Q$QO58)+ -0.041 * as.numeric(TPQ_Q$QO59)+ -0.002 * as.numeric(TPQ_Q$QO60)+ 0.043 * as.numeric(TPQ_Q$QO61)+ -0.006 * as.numeric(TPQ_Q$QO62)+ -0.097 * as.numeric(TPQ_Q$QO63)+ -0.053 * as.numeric(TPQ_Q$QO64)+ -0.285 * as.numeric(TPQ_Q$QO65)+ -0.006 * as.numeric(TPQ_Q$QO66)+ -0.02 * as.numeric(TPQ_Q$QO67)+ -0.092 * as.numeric(TPQ_Q$QO68)+ -0.124 * as.numeric(TPQ_Q$QO69)+ -0.04 * as.numeric(TPQ_Q$QO70)+ -0.343 * as.numeric(TPQ_Q$QO71)+ -0.075 * as.numeric(TPQ_Q$QO72)+ -0.367 * as.numeric(TPQ_Q$QO73)+ 0.025 * as.numeric(TPQ_Q$QO74)+ -0.13 * as.numeric(TPQ_Q$QO75)+ 0.021 * as.numeric(TPQ_Q$QO76)+ -0.397 * as.numeric(TPQ_Q$QO77)+ -0.126 * as.numeric(TPQ_Q$QO78)+ -0.029 * as.numeric(TPQ_Q$QO79)+ 0.046 * as.numeric(TPQ_Q$QO80)+ 0.039 * as.numeric(TPQ_Q$QO81)+ 0.011 * as.numeric(TPQ_Q$QO82)+ -0.082 * as.numeric(TPQ_Q$QO83)+ -0.092 * as.numeric(TPQ_Q$QO84)+ -0.094 * as.numeric(TPQ_Q$QO85)+ 0.028 * as.numeric(TPQ_Q$QO86)+ -0.133 * as.numeric(TPQ_Q$QO87)+ 0.027 * as.numeric(TPQ_Q$QO88)+ -0.031 * as.numeric(TPQ_Q$QO89)+ -0.165 * as.numeric(TPQ_Q$QO90)+ -0.141 * as.numeric(TPQ_Q$QO91)+ -0.108 * as.numeric(TPQ_Q$QO92)+ -0.169 * as.numeric(TPQ_Q$QO93)+ -0.02 * as.numeric(TPQ_Q$QO94)+ -0.023 * as.numeric(TPQ_Q$QO95)+ -0.11 * as.numeric(TPQ_Q$QO96)+ -0.172 * as.numeric(TPQ_Q$QO97)+ 0.071 * as.numeric(TPQ_Q$QO98)+ -0.068 * as.numeric(TPQ_Q$QO99)+ -0.116 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_11_HC =  0.035 * as.numeric(TPQ_Q$QO1)+ 0.06 * as.numeric(TPQ_Q$QO2)+ 0.077 * as.numeric(TPQ_Q$QO3)+ 0.157 * as.numeric(TPQ_Q$QO4)+ -0.066 * as.numeric(TPQ_Q$QO5)+ 0.225 * as.numeric(TPQ_Q$QO6)+ 0.042 * as.numeric(TPQ_Q$QO7)+ 0.12 * as.numeric(TPQ_Q$QO8)+ -0.075 * as.numeric(TPQ_Q$QO9)+ 0.116 * as.numeric(TPQ_Q$QO10)+ 0.122 * as.numeric(TPQ_Q$QO11)+ -0.02 * as.numeric(TPQ_Q$QO12)+ 0.128 * as.numeric(TPQ_Q$QO13)+ 0.278 * as.numeric(TPQ_Q$QO14)+ -0.091 * as.numeric(TPQ_Q$QO15)+ -0.072 * as.numeric(TPQ_Q$QO16)+ 0.06 * as.numeric(TPQ_Q$QO17)+ -0.058 * as.numeric(TPQ_Q$QO18)+ -0.077 * as.numeric(TPQ_Q$QO19)+ 0.186 * as.numeric(TPQ_Q$QO20)+ -0.126 * as.numeric(TPQ_Q$QO21)+ 0.041 * as.numeric(TPQ_Q$QO22)+ -0.021 * as.numeric(TPQ_Q$QO23)+ 0.282 * as.numeric(TPQ_Q$QO24)+ 0.086 * as.numeric(TPQ_Q$QO25)+ 0.049 * as.numeric(TPQ_Q$QO26)+ -0.042 * as.numeric(TPQ_Q$QO27)+ -0.035 * as.numeric(TPQ_Q$QO28)+ 0.066 * as.numeric(TPQ_Q$QO29)+ 0.251 * as.numeric(TPQ_Q$QO30)+ 0.02 * as.numeric(TPQ_Q$QO31)+ -0.08 * as.numeric(TPQ_Q$QO32)+ -0.255 * as.numeric(TPQ_Q$QO33)+ 0.032 * as.numeric(TPQ_Q$QO34)+ -0.151 * as.numeric(TPQ_Q$QO35)+ -0.064 * as.numeric(TPQ_Q$QO36)+ -0.069 * as.numeric(TPQ_Q$QO37)+ -0.12 * as.numeric(TPQ_Q$QO38)+ 0.002 * as.numeric(TPQ_Q$QO39)+ 0.046 * as.numeric(TPQ_Q$QO40)+ 0.034 * as.numeric(TPQ_Q$QO41)+ -0.01 * as.numeric(TPQ_Q$QO42)+ 0.014 * as.numeric(TPQ_Q$QO43)+ 0.024 * as.numeric(TPQ_Q$QO44)+ 0.136 * as.numeric(TPQ_Q$QO45)+ -0.223 * as.numeric(TPQ_Q$QO46)+ 0.049 * as.numeric(TPQ_Q$QO47)+ 0.194 * as.numeric(TPQ_Q$QO48)+ 0.13 * as.numeric(TPQ_Q$QO49)+ 0.094 * as.numeric(TPQ_Q$QO50)+ 0.004 * as.numeric(TPQ_Q$QO51)+ 0.043 * as.numeric(TPQ_Q$QO52)+ 0.146 * as.numeric(TPQ_Q$QO53)+ -0.006 * as.numeric(TPQ_Q$QO54)+ -0.206 * as.numeric(TPQ_Q$QO55)+ -0.259 * as.numeric(TPQ_Q$QO56)+ -0.038 * as.numeric(TPQ_Q$QO57)+ -0.032 * as.numeric(TPQ_Q$QO58)+ 0.044 * as.numeric(TPQ_Q$QO59)+ 0.566 * as.numeric(TPQ_Q$QO60)+ 0.05 * as.numeric(TPQ_Q$QO61)+ 0.498 * as.numeric(TPQ_Q$QO62)+ -0.037 * as.numeric(TPQ_Q$QO63)+ -0.019 * as.numeric(TPQ_Q$QO64)+ -0.455 * as.numeric(TPQ_Q$QO65)+ -0.017 * as.numeric(TPQ_Q$QO66)+ -0.001 * as.numeric(TPQ_Q$QO67)+ 0.01 * as.numeric(TPQ_Q$QO68)+ 0.04 * as.numeric(TPQ_Q$QO69)+ 0.125 * as.numeric(TPQ_Q$QO70)+ 0.12 * as.numeric(TPQ_Q$QO71)+ 0.114 * as.numeric(TPQ_Q$QO72)+ -0.023 * as.numeric(TPQ_Q$QO73)+ 0.033 * as.numeric(TPQ_Q$QO74)+ 0.124 * as.numeric(TPQ_Q$QO75)+ 0.12 * as.numeric(TPQ_Q$QO76)+ 0.021 * as.numeric(TPQ_Q$QO77)+ 0.052 * as.numeric(TPQ_Q$QO78)+ 0.086 * as.numeric(TPQ_Q$QO79)+ 0.104 * as.numeric(TPQ_Q$QO80)+ -0.067 * as.numeric(TPQ_Q$QO81)+ -0.024 * as.numeric(TPQ_Q$QO82)+ 0.045 * as.numeric(TPQ_Q$QO83)+ 0.202 * as.numeric(TPQ_Q$QO84)+ -0.012 * as.numeric(TPQ_Q$QO85)+ -0.098 * as.numeric(TPQ_Q$QO86)+ 0.045 * as.numeric(TPQ_Q$QO87)+ -0.005 * as.numeric(TPQ_Q$QO88)+ 0.179 * as.numeric(TPQ_Q$QO89)+ -0.089 * as.numeric(TPQ_Q$QO90)+ 0.093 * as.numeric(TPQ_Q$QO91)+ -0.037 * as.numeric(TPQ_Q$QO92)+ -0.198 * as.numeric(TPQ_Q$QO93)+ 0.011 * as.numeric(TPQ_Q$QO94)+ 0.005 * as.numeric(TPQ_Q$QO95)+ -0.009 * as.numeric(TPQ_Q$QO96)+ -0.037 * as.numeric(TPQ_Q$QO97)+ 0.04 * as.numeric(TPQ_Q$QO98)+ -0.097 * as.numeric(TPQ_Q$QO99)+ 0.133 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_12_HC =  -0.113 * as.numeric(TPQ_Q$QO1)+ -0.154 * as.numeric(TPQ_Q$QO2)+ -0.094 * as.numeric(TPQ_Q$QO3)+ -0.298 * as.numeric(TPQ_Q$QO4)+ -0.057 * as.numeric(TPQ_Q$QO5)+ 0.065 * as.numeric(TPQ_Q$QO6)+ -0.045 * as.numeric(TPQ_Q$QO7)+ 0.167 * as.numeric(TPQ_Q$QO8)+ -0.177 * as.numeric(TPQ_Q$QO9)+ -0.224 * as.numeric(TPQ_Q$QO10)+ 0.261 * as.numeric(TPQ_Q$QO11)+ 0.006 * as.numeric(TPQ_Q$QO12)+ 0.148 * as.numeric(TPQ_Q$QO13)+ -0.211 * as.numeric(TPQ_Q$QO14)+ 0.09 * as.numeric(TPQ_Q$QO15)+ -0.128 * as.numeric(TPQ_Q$QO16)+ 0.09 * as.numeric(TPQ_Q$QO17)+ 0.04 * as.numeric(TPQ_Q$QO18)+ -0.003 * as.numeric(TPQ_Q$QO19)+ -0.004 * as.numeric(TPQ_Q$QO20)+ -0.185 * as.numeric(TPQ_Q$QO21)+ 0.102 * as.numeric(TPQ_Q$QO22)+ -0.059 * as.numeric(TPQ_Q$QO23)+ 0.105 * as.numeric(TPQ_Q$QO24)+ 0.024 * as.numeric(TPQ_Q$QO25)+ 0.076 * as.numeric(TPQ_Q$QO26)+ -0.02 * as.numeric(TPQ_Q$QO27)+ 0.307 * as.numeric(TPQ_Q$QO28)+ 0.07 * as.numeric(TPQ_Q$QO29)+ 0.199 * as.numeric(TPQ_Q$QO30)+ -0.085 * as.numeric(TPQ_Q$QO31)+ -0.236 * as.numeric(TPQ_Q$QO32)+ 0.116 * as.numeric(TPQ_Q$QO33)+ 0.023 * as.numeric(TPQ_Q$QO34)+ 0.078 * as.numeric(TPQ_Q$QO35)+ -0.225 * as.numeric(TPQ_Q$QO36)+ -0.005 * as.numeric(TPQ_Q$QO37)+ 0.053 * as.numeric(TPQ_Q$QO38)+ -0.065 * as.numeric(TPQ_Q$QO39)+ 0.076 * as.numeric(TPQ_Q$QO40)+ -0.012 * as.numeric(TPQ_Q$QO41)+ -0.19 * as.numeric(TPQ_Q$QO42)+ 0.065 * as.numeric(TPQ_Q$QO43)+ -0.167 * as.numeric(TPQ_Q$QO44)+ 0.026 * as.numeric(TPQ_Q$QO45)+ -0.253 * as.numeric(TPQ_Q$QO46)+ 0.018 * as.numeric(TPQ_Q$QO47)+ 0.357 * as.numeric(TPQ_Q$QO48)+ 0.065 * as.numeric(TPQ_Q$QO49)+ 0.115 * as.numeric(TPQ_Q$QO50)+ -0.065 * as.numeric(TPQ_Q$QO51)+ -0.098 * as.numeric(TPQ_Q$QO52)+ 0.037 * as.numeric(TPQ_Q$QO53)+ 0.039 * as.numeric(TPQ_Q$QO54)+ -0.231 * as.numeric(TPQ_Q$QO55)+ -0.249 * as.numeric(TPQ_Q$QO56)+ -0.025 * as.numeric(TPQ_Q$QO57)+ 0.012 * as.numeric(TPQ_Q$QO58)+ -0.02 * as.numeric(TPQ_Q$QO59)+ -0.083 * as.numeric(TPQ_Q$QO60)+ 0.021 * as.numeric(TPQ_Q$QO61)+ -0.11 * as.numeric(TPQ_Q$QO62)+ -0.179 * as.numeric(TPQ_Q$QO63)+ -0.041 * as.numeric(TPQ_Q$QO64)+ -0.02 * as.numeric(TPQ_Q$QO65)+ -0.092 * as.numeric(TPQ_Q$QO66)+ -0.103 * as.numeric(TPQ_Q$QO67)+ -0.124 * as.numeric(TPQ_Q$QO68)+ -0.161 * as.numeric(TPQ_Q$QO69)+ 0.13 * as.numeric(TPQ_Q$QO70)+ -0.277 * as.numeric(TPQ_Q$QO71)+ 0.073 * as.numeric(TPQ_Q$QO72)+ 0.102 * as.numeric(TPQ_Q$QO73)+ -0.057 * as.numeric(TPQ_Q$QO74)+ 0.0 * as.numeric(TPQ_Q$QO75)+ 0.009 * as.numeric(TPQ_Q$QO76)+ -0.231 * as.numeric(TPQ_Q$QO77)+ -0.01 * as.numeric(TPQ_Q$QO78)+ -0.06 * as.numeric(TPQ_Q$QO79)+ 0.111 * as.numeric(TPQ_Q$QO80)+ -0.203 * as.numeric(TPQ_Q$QO81)+ -0.052 * as.numeric(TPQ_Q$QO82)+ -0.101 * as.numeric(TPQ_Q$QO83)+ -0.185 * as.numeric(TPQ_Q$QO84)+ 0.144 * as.numeric(TPQ_Q$QO85)+ -0.07 * as.numeric(TPQ_Q$QO86)+ -0.148 * as.numeric(TPQ_Q$QO87)+ 0.037 * as.numeric(TPQ_Q$QO88)+ -0.189 * as.numeric(TPQ_Q$QO89)+ 0.132 * as.numeric(TPQ_Q$QO90)+ 0.015 * as.numeric(TPQ_Q$QO91)+ -0.023 * as.numeric(TPQ_Q$QO92)+ -0.012 * as.numeric(TPQ_Q$QO93)+ -0.086 * as.numeric(TPQ_Q$QO94)+ -0.017 * as.numeric(TPQ_Q$QO95)+ 0.051 * as.numeric(TPQ_Q$QO96)+ -0.169 * as.numeric(TPQ_Q$QO97)+ -0.01 * as.numeric(TPQ_Q$QO98)+ -0.188 * as.numeric(TPQ_Q$QO99)+ -0.137 * as.numeric(TPQ_Q$QO100)

TPQ_Q$Complete_IC_0_MDD =  0.027 * as.numeric(TPQ_Q$QO1)+ 0.208 * as.numeric(TPQ_Q$QO2)+ -0.139 * as.numeric(TPQ_Q$QO3)+ 0.005 * as.numeric(TPQ_Q$QO4)+ -0.066 * as.numeric(TPQ_Q$QO5)+ 0.038 * as.numeric(TPQ_Q$QO6)+ -0.018 * as.numeric(TPQ_Q$QO7)+ -0.002 * as.numeric(TPQ_Q$QO8)+ 0.111 * as.numeric(TPQ_Q$QO9)+ -0.181 * as.numeric(TPQ_Q$QO10)+ -0.157 * as.numeric(TPQ_Q$QO11)+ 0.185 * as.numeric(TPQ_Q$QO12)+ 0.378 * as.numeric(TPQ_Q$QO13)+ -0.257 * as.numeric(TPQ_Q$QO14)+ 0.21 * as.numeric(TPQ_Q$QO15)+ -0.022 * as.numeric(TPQ_Q$QO16)+ 0.159 * as.numeric(TPQ_Q$QO17)+ -0.011 * as.numeric(TPQ_Q$QO18)+ 0.053 * as.numeric(TPQ_Q$QO19)+ 0.113 * as.numeric(TPQ_Q$QO20)+ 0.003 * as.numeric(TPQ_Q$QO21)+ 0.023 * as.numeric(TPQ_Q$QO22)+ 0.011 * as.numeric(TPQ_Q$QO23)+ 0.22 * as.numeric(TPQ_Q$QO24)+ 0.262 * as.numeric(TPQ_Q$QO25)+ -0.067 * as.numeric(TPQ_Q$QO26)+ -0.072 * as.numeric(TPQ_Q$QO27)+ 0.11 * as.numeric(TPQ_Q$QO28)+ -0.051 * as.numeric(TPQ_Q$QO29)+ -0.006 * as.numeric(TPQ_Q$QO30)+ 0.007 * as.numeric(TPQ_Q$QO31)+ -0.393 * as.numeric(TPQ_Q$QO32)+ -0.267 * as.numeric(TPQ_Q$QO33)+ -0.067 * as.numeric(TPQ_Q$QO34)+ -0.126 * as.numeric(TPQ_Q$QO35)+ 0.036 * as.numeric(TPQ_Q$QO36)+ -0.229 * as.numeric(TPQ_Q$QO37)+ -0.285 * as.numeric(TPQ_Q$QO38)+ -0.029 * as.numeric(TPQ_Q$QO39)+ -0.185 * as.numeric(TPQ_Q$QO40)+ -0.088 * as.numeric(TPQ_Q$QO41)+ -0.32 * as.numeric(TPQ_Q$QO42)+ -0.112 * as.numeric(TPQ_Q$QO43)+ -0.204 * as.numeric(TPQ_Q$QO44)+ 0.057 * as.numeric(TPQ_Q$QO45)+ -0.191 * as.numeric(TPQ_Q$QO46)+ -0.05 * as.numeric(TPQ_Q$QO47)+ 0.011 * as.numeric(TPQ_Q$QO48)+ -0.23 * as.numeric(TPQ_Q$QO49)+ -0.087 * as.numeric(TPQ_Q$QO50)+ 0.002 * as.numeric(TPQ_Q$QO51)+ 0.016 * as.numeric(TPQ_Q$QO52)+ -0.158 * as.numeric(TPQ_Q$QO53)+ -0.1 * as.numeric(TPQ_Q$QO54)+ -0.004 * as.numeric(TPQ_Q$QO55)+ 0.034 * as.numeric(TPQ_Q$QO56)+ -0.104 * as.numeric(TPQ_Q$QO57)+ -0.038 * as.numeric(TPQ_Q$QO58)+ 0.105 * as.numeric(TPQ_Q$QO59)+ 0.0 * as.numeric(TPQ_Q$QO60)+ 0.046 * as.numeric(TPQ_Q$QO61)+ -0.054 * as.numeric(TPQ_Q$QO62)+ 0.093 * as.numeric(TPQ_Q$QO63)+ 0.028 * as.numeric(TPQ_Q$QO64)+ -0.042 * as.numeric(TPQ_Q$QO65)+ -0.08 * as.numeric(TPQ_Q$QO66)+ -0.058 * as.numeric(TPQ_Q$QO67)+ -0.246 * as.numeric(TPQ_Q$QO68)+ -0.006 * as.numeric(TPQ_Q$QO69)+ 0.061 * as.numeric(TPQ_Q$QO70)+ 0.316 * as.numeric(TPQ_Q$QO71)+ 0.047 * as.numeric(TPQ_Q$QO72)+ -0.072 * as.numeric(TPQ_Q$QO73)+ 0.026 * as.numeric(TPQ_Q$QO74)+ 0.006 * as.numeric(TPQ_Q$QO75)+ 0.065 * as.numeric(TPQ_Q$QO76)+ 0.231 * as.numeric(TPQ_Q$QO77)+ -0.112 * as.numeric(TPQ_Q$QO78)+ -0.205 * as.numeric(TPQ_Q$QO79)+ -0.045 * as.numeric(TPQ_Q$QO80)+ -0.184 * as.numeric(TPQ_Q$QO81)+ -0.189 * as.numeric(TPQ_Q$QO82)+ 0.093 * as.numeric(TPQ_Q$QO83)+ 0.038 * as.numeric(TPQ_Q$QO84)+ -0.278 * as.numeric(TPQ_Q$QO85)+ 0.069 * as.numeric(TPQ_Q$QO86)+ -0.135 * as.numeric(TPQ_Q$QO87)+ 0.156 * as.numeric(TPQ_Q$QO88)+ 0.176 * as.numeric(TPQ_Q$QO89)+ 0.174 * as.numeric(TPQ_Q$QO90)+ -0.012 * as.numeric(TPQ_Q$QO91)+ 0.038 * as.numeric(TPQ_Q$QO92)+ -0.222 * as.numeric(TPQ_Q$QO93)+ -0.021 * as.numeric(TPQ_Q$QO94)+ 0.004 * as.numeric(TPQ_Q$QO95)+ -0.216 * as.numeric(TPQ_Q$QO96)+ 0.052 * as.numeric(TPQ_Q$QO97)+ -0.096 * as.numeric(TPQ_Q$QO98)+ 0.043 * as.numeric(TPQ_Q$QO99)+ 0.056 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_1_MDD =  0.088 * as.numeric(TPQ_Q$QO1)+ -0.156 * as.numeric(TPQ_Q$QO2)+ -0.224 * as.numeric(TPQ_Q$QO3)+ -0.131 * as.numeric(TPQ_Q$QO4)+ -0.275 * as.numeric(TPQ_Q$QO5)+ 0.127 * as.numeric(TPQ_Q$QO6)+ -0.245 * as.numeric(TPQ_Q$QO7)+ 0.068 * as.numeric(TPQ_Q$QO8)+ 0.015 * as.numeric(TPQ_Q$QO9)+ -0.316 * as.numeric(TPQ_Q$QO10)+ -0.129 * as.numeric(TPQ_Q$QO11)+ -0.133 * as.numeric(TPQ_Q$QO12)+ -0.278 * as.numeric(TPQ_Q$QO13)+ -0.422 * as.numeric(TPQ_Q$QO14)+ 0.389 * as.numeric(TPQ_Q$QO15)+ 0.002 * as.numeric(TPQ_Q$QO16)+ -0.306 * as.numeric(TPQ_Q$QO17)+ -0.19 * as.numeric(TPQ_Q$QO18)+ -0.28 * as.numeric(TPQ_Q$QO19)+ 0.135 * as.numeric(TPQ_Q$QO20)+ -0.102 * as.numeric(TPQ_Q$QO21)+ -0.043 * as.numeric(TPQ_Q$QO22)+ -0.327 * as.numeric(TPQ_Q$QO23)+ -0.016 * as.numeric(TPQ_Q$QO24)+ 0.156 * as.numeric(TPQ_Q$QO25)+ 0.212 * as.numeric(TPQ_Q$QO26)+ -0.093 * as.numeric(TPQ_Q$QO27)+ -0.031 * as.numeric(TPQ_Q$QO28)+ 0.161 * as.numeric(TPQ_Q$QO29)+ -0.165 * as.numeric(TPQ_Q$QO30)+ 0.018 * as.numeric(TPQ_Q$QO31)+ 0.125 * as.numeric(TPQ_Q$QO32)+ -0.126 * as.numeric(TPQ_Q$QO33)+ -0.1 * as.numeric(TPQ_Q$QO34)+ 0.123 * as.numeric(TPQ_Q$QO35)+ 0.072 * as.numeric(TPQ_Q$QO36)+ -0.132 * as.numeric(TPQ_Q$QO37)+ -0.029 * as.numeric(TPQ_Q$QO38)+ 0.091 * as.numeric(TPQ_Q$QO39)+ 0.071 * as.numeric(TPQ_Q$QO40)+ 0.108 * as.numeric(TPQ_Q$QO41)+ 0.104 * as.numeric(TPQ_Q$QO42)+ -0.132 * as.numeric(TPQ_Q$QO43)+ 0.145 * as.numeric(TPQ_Q$QO44)+ 0.056 * as.numeric(TPQ_Q$QO45)+ -0.066 * as.numeric(TPQ_Q$QO46)+ 0.148 * as.numeric(TPQ_Q$QO47)+ -0.111 * as.numeric(TPQ_Q$QO48)+ -0.167 * as.numeric(TPQ_Q$QO49)+ -0.26 * as.numeric(TPQ_Q$QO50)+ 0.071 * as.numeric(TPQ_Q$QO51)+ -0.216 * as.numeric(TPQ_Q$QO52)+ -0.044 * as.numeric(TPQ_Q$QO53)+ -0.077 * as.numeric(TPQ_Q$QO54)+ -0.008 * as.numeric(TPQ_Q$QO55)+ -0.184 * as.numeric(TPQ_Q$QO56)+ -0.131 * as.numeric(TPQ_Q$QO57)+ 0.158 * as.numeric(TPQ_Q$QO58)+ 0.026 * as.numeric(TPQ_Q$QO59)+ 0.111 * as.numeric(TPQ_Q$QO60)+ -0.157 * as.numeric(TPQ_Q$QO61)+ -0.095 * as.numeric(TPQ_Q$QO62)+ 0.01 * as.numeric(TPQ_Q$QO63)+ -0.199 * as.numeric(TPQ_Q$QO64)+ -0.089 * as.numeric(TPQ_Q$QO65)+ -0.098 * as.numeric(TPQ_Q$QO66)+ -0.143 * as.numeric(TPQ_Q$QO67)+ -0.023 * as.numeric(TPQ_Q$QO68)+ -0.164 * as.numeric(TPQ_Q$QO69)+ 0.076 * as.numeric(TPQ_Q$QO70)+ -0.009 * as.numeric(TPQ_Q$QO71)+ -0.03 * as.numeric(TPQ_Q$QO72)+ -0.049 * as.numeric(TPQ_Q$QO73)+ -0.21 * as.numeric(TPQ_Q$QO74)+ 0.118 * as.numeric(TPQ_Q$QO75)+ -0.044 * as.numeric(TPQ_Q$QO76)+ -0.068 * as.numeric(TPQ_Q$QO77)+ 0.159 * as.numeric(TPQ_Q$QO78)+ 0.065 * as.numeric(TPQ_Q$QO79)+ -0.088 * as.numeric(TPQ_Q$QO80)+ 0.115 * as.numeric(TPQ_Q$QO81)+ -0.131 * as.numeric(TPQ_Q$QO82)+ -0.074 * as.numeric(TPQ_Q$QO83)+ 0.07 * as.numeric(TPQ_Q$QO84)+ -0.208 * as.numeric(TPQ_Q$QO85)+ 0.069 * as.numeric(TPQ_Q$QO86)+ 0.008 * as.numeric(TPQ_Q$QO87)+ 0.036 * as.numeric(TPQ_Q$QO88)+ 0.256 * as.numeric(TPQ_Q$QO89)+ 0.09 * as.numeric(TPQ_Q$QO90)+ 0.25 * as.numeric(TPQ_Q$QO91)+ 0.278 * as.numeric(TPQ_Q$QO92)+ -0.101 * as.numeric(TPQ_Q$QO93)+ -0.031 * as.numeric(TPQ_Q$QO94)+ -0.006 * as.numeric(TPQ_Q$QO95)+ 0.14 * as.numeric(TPQ_Q$QO96)+ 0.111 * as.numeric(TPQ_Q$QO97)+ 0.006 * as.numeric(TPQ_Q$QO98)+ -0.0 * as.numeric(TPQ_Q$QO99)+ 0.058 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_2_MDD =  -0.032 * as.numeric(TPQ_Q$QO1)+ 0.023 * as.numeric(TPQ_Q$QO2)+ -0.076 * as.numeric(TPQ_Q$QO3)+ 0.018 * as.numeric(TPQ_Q$QO4)+ 0.031 * as.numeric(TPQ_Q$QO5)+ -0.051 * as.numeric(TPQ_Q$QO6)+ -0.135 * as.numeric(TPQ_Q$QO7)+ 0.106 * as.numeric(TPQ_Q$QO8)+ -0.08 * as.numeric(TPQ_Q$QO9)+ 0.186 * as.numeric(TPQ_Q$QO10)+ -0.032 * as.numeric(TPQ_Q$QO11)+ 0.226 * as.numeric(TPQ_Q$QO12)+ 0.158 * as.numeric(TPQ_Q$QO13)+ 0.248 * as.numeric(TPQ_Q$QO14)+ 0.125 * as.numeric(TPQ_Q$QO15)+ -0.123 * as.numeric(TPQ_Q$QO16)+ 0.046 * as.numeric(TPQ_Q$QO17)+ 0.128 * as.numeric(TPQ_Q$QO18)+ 0.122 * as.numeric(TPQ_Q$QO19)+ 0.177 * as.numeric(TPQ_Q$QO20)+ 0.026 * as.numeric(TPQ_Q$QO21)+ 0.262 * as.numeric(TPQ_Q$QO22)+ 0.1 * as.numeric(TPQ_Q$QO23)+ 0.29 * as.numeric(TPQ_Q$QO24)+ 0.036 * as.numeric(TPQ_Q$QO25)+ 0.033 * as.numeric(TPQ_Q$QO26)+ -0.077 * as.numeric(TPQ_Q$QO27)+ 0.349 * as.numeric(TPQ_Q$QO28)+ 0.136 * as.numeric(TPQ_Q$QO29)+ 0.285 * as.numeric(TPQ_Q$QO30)+ 0.046 * as.numeric(TPQ_Q$QO31)+ -0.324 * as.numeric(TPQ_Q$QO32)+ 0.229 * as.numeric(TPQ_Q$QO33)+ 0.049 * as.numeric(TPQ_Q$QO34)+ 0.328 * as.numeric(TPQ_Q$QO35)+ 0.145 * as.numeric(TPQ_Q$QO36)+ 0.173 * as.numeric(TPQ_Q$QO37)+ 0.334 * as.numeric(TPQ_Q$QO38)+ 0.082 * as.numeric(TPQ_Q$QO39)+ 0.335 * as.numeric(TPQ_Q$QO40)+ 0.025 * as.numeric(TPQ_Q$QO41)+ -0.225 * as.numeric(TPQ_Q$QO42)+ 0.207 * as.numeric(TPQ_Q$QO43)+ -0.209 * as.numeric(TPQ_Q$QO44)+ 0.084 * as.numeric(TPQ_Q$QO45)+ -0.083 * as.numeric(TPQ_Q$QO46)+ 0.11 * as.numeric(TPQ_Q$QO47)+ 0.433 * as.numeric(TPQ_Q$QO48)+ 0.504 * as.numeric(TPQ_Q$QO49)+ 0.287 * as.numeric(TPQ_Q$QO50)+ 0.008 * as.numeric(TPQ_Q$QO51)+ -0.051 * as.numeric(TPQ_Q$QO52)+ 0.108 * as.numeric(TPQ_Q$QO53)+ 0.142 * as.numeric(TPQ_Q$QO54)+ 0.006 * as.numeric(TPQ_Q$QO55)+ 0.04 * as.numeric(TPQ_Q$QO56)+ 0.342 * as.numeric(TPQ_Q$QO57)+ -0.043 * as.numeric(TPQ_Q$QO58)+ -0.093 * as.numeric(TPQ_Q$QO59)+ 0.091 * as.numeric(TPQ_Q$QO60)+ 0.053 * as.numeric(TPQ_Q$QO61)+ 0.121 * as.numeric(TPQ_Q$QO62)+ -0.09 * as.numeric(TPQ_Q$QO63)+ -0.056 * as.numeric(TPQ_Q$QO64)+ -0.158 * as.numeric(TPQ_Q$QO65)+ -0.407 * as.numeric(TPQ_Q$QO66)+ -0.197 * as.numeric(TPQ_Q$QO67)+ 0.112 * as.numeric(TPQ_Q$QO68)+ 0.252 * as.numeric(TPQ_Q$QO69)+ 0.521 * as.numeric(TPQ_Q$QO70)+ 0.096 * as.numeric(TPQ_Q$QO71)+ 0.593 * as.numeric(TPQ_Q$QO72)+ 0.172 * as.numeric(TPQ_Q$QO73)+ -0.325 * as.numeric(TPQ_Q$QO74)+ -0.146 * as.numeric(TPQ_Q$QO75)+ -0.179 * as.numeric(TPQ_Q$QO76)+ 0.052 * as.numeric(TPQ_Q$QO77)+ 0.195 * as.numeric(TPQ_Q$QO78)+ 0.121 * as.numeric(TPQ_Q$QO79)+ -0.057 * as.numeric(TPQ_Q$QO80)+ -0.141 * as.numeric(TPQ_Q$QO81)+ -0.091 * as.numeric(TPQ_Q$QO82)+ 0.064 * as.numeric(TPQ_Q$QO83)+ -0.029 * as.numeric(TPQ_Q$QO84)+ -0.037 * as.numeric(TPQ_Q$QO85)+ 0.133 * as.numeric(TPQ_Q$QO86)+ -0.212 * as.numeric(TPQ_Q$QO87)+ 0.206 * as.numeric(TPQ_Q$QO88)+ -0.169 * as.numeric(TPQ_Q$QO89)+ 0.268 * as.numeric(TPQ_Q$QO90)+ 0.051 * as.numeric(TPQ_Q$QO91)+ -0.251 * as.numeric(TPQ_Q$QO92)+ 0.145 * as.numeric(TPQ_Q$QO93)+ -0.118 * as.numeric(TPQ_Q$QO94)+ -0.144 * as.numeric(TPQ_Q$QO95)+ 0.103 * as.numeric(TPQ_Q$QO96)+ -0.114 * as.numeric(TPQ_Q$QO97)+ 0.093 * as.numeric(TPQ_Q$QO98)+ -0.057 * as.numeric(TPQ_Q$QO99)+ -0.199 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_3_MDD =  0.003 * as.numeric(TPQ_Q$QO1)+ -0.111 * as.numeric(TPQ_Q$QO2)+ -0.135 * as.numeric(TPQ_Q$QO3)+ 0.0 * as.numeric(TPQ_Q$QO4)+ 0.142 * as.numeric(TPQ_Q$QO5)+ 0.229 * as.numeric(TPQ_Q$QO6)+ -0.121 * as.numeric(TPQ_Q$QO7)+ 0.059 * as.numeric(TPQ_Q$QO8)+ 0.034 * as.numeric(TPQ_Q$QO9)+ 0.247 * as.numeric(TPQ_Q$QO10)+ 0.213 * as.numeric(TPQ_Q$QO11)+ 0.199 * as.numeric(TPQ_Q$QO12)+ 0.19 * as.numeric(TPQ_Q$QO13)+ 0.082 * as.numeric(TPQ_Q$QO14)+ 0.099 * as.numeric(TPQ_Q$QO15)+ -0.014 * as.numeric(TPQ_Q$QO16)+ 0.367 * as.numeric(TPQ_Q$QO17)+ -0.002 * as.numeric(TPQ_Q$QO18)+ -0.146 * as.numeric(TPQ_Q$QO19)+ -0.042 * as.numeric(TPQ_Q$QO20)+ -0.202 * as.numeric(TPQ_Q$QO21)+ 0.008 * as.numeric(TPQ_Q$QO22)+ 0.059 * as.numeric(TPQ_Q$QO23)+ 0.166 * as.numeric(TPQ_Q$QO24)+ 0.23 * as.numeric(TPQ_Q$QO25)+ 0.049 * as.numeric(TPQ_Q$QO26)+ -0.08 * as.numeric(TPQ_Q$QO27)+ 0.075 * as.numeric(TPQ_Q$QO28)+ 0.164 * as.numeric(TPQ_Q$QO29)+ -0.113 * as.numeric(TPQ_Q$QO30)+ -0.069 * as.numeric(TPQ_Q$QO31)+ -0.04 * as.numeric(TPQ_Q$QO32)+ -0.132 * as.numeric(TPQ_Q$QO33)+ -0.148 * as.numeric(TPQ_Q$QO34)+ -0.094 * as.numeric(TPQ_Q$QO35)+ 0.171 * as.numeric(TPQ_Q$QO36)+ 0.085 * as.numeric(TPQ_Q$QO37)+ 0.031 * as.numeric(TPQ_Q$QO38)+ -0.09 * as.numeric(TPQ_Q$QO39)+ 0.09 * as.numeric(TPQ_Q$QO40)+ -0.111 * as.numeric(TPQ_Q$QO41)+ 0.014 * as.numeric(TPQ_Q$QO42)+ 0.037 * as.numeric(TPQ_Q$QO43)+ 0.061 * as.numeric(TPQ_Q$QO44)+ 0.231 * as.numeric(TPQ_Q$QO45)+ -0.125 * as.numeric(TPQ_Q$QO46)+ 0.107 * as.numeric(TPQ_Q$QO47)+ 0.034 * as.numeric(TPQ_Q$QO48)+ 0.303 * as.numeric(TPQ_Q$QO49)+ -0.141 * as.numeric(TPQ_Q$QO50)+ 0.119 * as.numeric(TPQ_Q$QO51)+ 0.141 * as.numeric(TPQ_Q$QO52)+ 0.054 * as.numeric(TPQ_Q$QO53)+ 0.154 * as.numeric(TPQ_Q$QO54)+ -0.188 * as.numeric(TPQ_Q$QO55)+ -0.193 * as.numeric(TPQ_Q$QO56)+ 0.134 * as.numeric(TPQ_Q$QO57)+ 0.121 * as.numeric(TPQ_Q$QO58)+ -0.17 * as.numeric(TPQ_Q$QO59)+ 0.333 * as.numeric(TPQ_Q$QO60)+ -0.167 * as.numeric(TPQ_Q$QO61)+ 0.17 * as.numeric(TPQ_Q$QO62)+ -0.091 * as.numeric(TPQ_Q$QO63)+ -0.096 * as.numeric(TPQ_Q$QO64)+ 0.067 * as.numeric(TPQ_Q$QO65)+ 0.202 * as.numeric(TPQ_Q$QO66)+ -0.122 * as.numeric(TPQ_Q$QO67)+ 0.195 * as.numeric(TPQ_Q$QO68)+ 0.156 * as.numeric(TPQ_Q$QO69)+ -0.202 * as.numeric(TPQ_Q$QO70)+ -0.175 * as.numeric(TPQ_Q$QO71)+ -0.214 * as.numeric(TPQ_Q$QO72)+ 0.092 * as.numeric(TPQ_Q$QO73)+ -0.047 * as.numeric(TPQ_Q$QO74)+ -0.111 * as.numeric(TPQ_Q$QO75)+ 0.064 * as.numeric(TPQ_Q$QO76)+ -0.269 * as.numeric(TPQ_Q$QO77)+ -0.158 * as.numeric(TPQ_Q$QO78)+ -0.232 * as.numeric(TPQ_Q$QO79)+ -0.24 * as.numeric(TPQ_Q$QO80)+ -0.084 * as.numeric(TPQ_Q$QO81)+ 0.047 * as.numeric(TPQ_Q$QO82)+ -0.225 * as.numeric(TPQ_Q$QO83)+ 0.052 * as.numeric(TPQ_Q$QO84)+ -0.009 * as.numeric(TPQ_Q$QO85)+ -0.071 * as.numeric(TPQ_Q$QO86)+ -0.042 * as.numeric(TPQ_Q$QO87)+ 0.082 * as.numeric(TPQ_Q$QO88)+ 0.033 * as.numeric(TPQ_Q$QO89)+ 0.178 * as.numeric(TPQ_Q$QO90)+ 0.196 * as.numeric(TPQ_Q$QO91)+ -0.053 * as.numeric(TPQ_Q$QO92)+ 0.123 * as.numeric(TPQ_Q$QO93)+ -0.104 * as.numeric(TPQ_Q$QO94)+ 0.032 * as.numeric(TPQ_Q$QO95)+ -0.066 * as.numeric(TPQ_Q$QO96)+ -0.152 * as.numeric(TPQ_Q$QO97)+ 0.129 * as.numeric(TPQ_Q$QO98)+ -0.113 * as.numeric(TPQ_Q$QO99)+ -0.091 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_4_MDD =  -0.108 * as.numeric(TPQ_Q$QO1)+ 0.157 * as.numeric(TPQ_Q$QO2)+ 0.158 * as.numeric(TPQ_Q$QO3)+ 0.218 * as.numeric(TPQ_Q$QO4)+ 0.213 * as.numeric(TPQ_Q$QO5)+ 0.289 * as.numeric(TPQ_Q$QO6)+ -0.036 * as.numeric(TPQ_Q$QO7)+ -0.093 * as.numeric(TPQ_Q$QO8)+ 0.244 * as.numeric(TPQ_Q$QO9)+ -0.008 * as.numeric(TPQ_Q$QO10)+ -0.028 * as.numeric(TPQ_Q$QO11)+ -0.191 * as.numeric(TPQ_Q$QO12)+ 0.063 * as.numeric(TPQ_Q$QO13)+ -0.073 * as.numeric(TPQ_Q$QO14)+ 0.135 * as.numeric(TPQ_Q$QO15)+ -0.029 * as.numeric(TPQ_Q$QO16)+ -0.171 * as.numeric(TPQ_Q$QO17)+ 0.035 * as.numeric(TPQ_Q$QO18)+ 0.042 * as.numeric(TPQ_Q$QO19)+ -0.074 * as.numeric(TPQ_Q$QO20)+ -0.152 * as.numeric(TPQ_Q$QO21)+ 0.188 * as.numeric(TPQ_Q$QO22)+ 0.032 * as.numeric(TPQ_Q$QO23)+ -0.022 * as.numeric(TPQ_Q$QO24)+ -0.196 * as.numeric(TPQ_Q$QO25)+ -0.158 * as.numeric(TPQ_Q$QO26)+ -0.074 * as.numeric(TPQ_Q$QO27)+ 0.113 * as.numeric(TPQ_Q$QO28)+ -0.187 * as.numeric(TPQ_Q$QO29)+ 0.18 * as.numeric(TPQ_Q$QO30)+ 0.161 * as.numeric(TPQ_Q$QO31)+ -0.143 * as.numeric(TPQ_Q$QO32)+ 0.08 * as.numeric(TPQ_Q$QO33)+ 0.025 * as.numeric(TPQ_Q$QO34)+ 0.146 * as.numeric(TPQ_Q$QO35)+ -0.04 * as.numeric(TPQ_Q$QO36)+ -0.207 * as.numeric(TPQ_Q$QO37)+ -0.037 * as.numeric(TPQ_Q$QO38)+ 0.161 * as.numeric(TPQ_Q$QO39)+ 0.164 * as.numeric(TPQ_Q$QO40)+ 0.04 * as.numeric(TPQ_Q$QO41)+ 0.136 * as.numeric(TPQ_Q$QO42)+ 0.078 * as.numeric(TPQ_Q$QO43)+ 0.016 * as.numeric(TPQ_Q$QO44)+ 0.077 * as.numeric(TPQ_Q$QO45)+ -0.102 * as.numeric(TPQ_Q$QO46)+ -0.311 * as.numeric(TPQ_Q$QO47)+ 0.18 * as.numeric(TPQ_Q$QO48)+ -0.069 * as.numeric(TPQ_Q$QO49)+ 0.132 * as.numeric(TPQ_Q$QO50)+ -0.259 * as.numeric(TPQ_Q$QO51)+ 0.108 * as.numeric(TPQ_Q$QO52)+ 0.214 * as.numeric(TPQ_Q$QO53)+ 0.013 * as.numeric(TPQ_Q$QO54)+ -0.11 * as.numeric(TPQ_Q$QO55)+ -0.103 * as.numeric(TPQ_Q$QO56)+ -0.017 * as.numeric(TPQ_Q$QO57)+ -0.197 * as.numeric(TPQ_Q$QO58)+ 0.011 * as.numeric(TPQ_Q$QO59)+ 0.292 * as.numeric(TPQ_Q$QO60)+ 0.052 * as.numeric(TPQ_Q$QO61)+ 0.241 * as.numeric(TPQ_Q$QO62)+ 0.06 * as.numeric(TPQ_Q$QO63)+ 0.062 * as.numeric(TPQ_Q$QO64)+ -0.287 * as.numeric(TPQ_Q$QO65)+ 0.343 * as.numeric(TPQ_Q$QO66)+ -0.002 * as.numeric(TPQ_Q$QO67)+ -0.058 * as.numeric(TPQ_Q$QO68)+ 0.001 * as.numeric(TPQ_Q$QO69)+ -0.037 * as.numeric(TPQ_Q$QO70)+ 0.022 * as.numeric(TPQ_Q$QO71)+ -0.08 * as.numeric(TPQ_Q$QO72)+ -0.01 * as.numeric(TPQ_Q$QO73)+ -0.206 * as.numeric(TPQ_Q$QO74)+ -0.038 * as.numeric(TPQ_Q$QO75)+ -0.017 * as.numeric(TPQ_Q$QO76)+ 0.053 * as.numeric(TPQ_Q$QO77)+ -0.049 * as.numeric(TPQ_Q$QO78)+ -0.031 * as.numeric(TPQ_Q$QO79)+ 0.018 * as.numeric(TPQ_Q$QO80)+ 0.097 * as.numeric(TPQ_Q$QO81)+ 0.056 * as.numeric(TPQ_Q$QO82)+ 0.171 * as.numeric(TPQ_Q$QO83)+ -0.313 * as.numeric(TPQ_Q$QO84)+ 0.048 * as.numeric(TPQ_Q$QO85)+ -0.097 * as.numeric(TPQ_Q$QO86)+ 0.157 * as.numeric(TPQ_Q$QO87)+ -0.011 * as.numeric(TPQ_Q$QO88)+ -0.096 * as.numeric(TPQ_Q$QO89)+ 0.084 * as.numeric(TPQ_Q$QO90)+ -0.203 * as.numeric(TPQ_Q$QO91)+ 0.15 * as.numeric(TPQ_Q$QO92)+ -0.039 * as.numeric(TPQ_Q$QO93)+ -0.003 * as.numeric(TPQ_Q$QO94)+ 0.01 * as.numeric(TPQ_Q$QO95)+ 0.007 * as.numeric(TPQ_Q$QO96)+ -0.002 * as.numeric(TPQ_Q$QO97)+ 0.107 * as.numeric(TPQ_Q$QO98)+ 0.019 * as.numeric(TPQ_Q$QO99)+ 0.114 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_5_MDD =  -0.126 * as.numeric(TPQ_Q$QO1)+ 0.079 * as.numeric(TPQ_Q$QO2)+ -0.228 * as.numeric(TPQ_Q$QO3)+ 0.161 * as.numeric(TPQ_Q$QO4)+ 0.008 * as.numeric(TPQ_Q$QO5)+ -0.082 * as.numeric(TPQ_Q$QO6)+ -0.203 * as.numeric(TPQ_Q$QO7)+ 0.125 * as.numeric(TPQ_Q$QO8)+ 0.237 * as.numeric(TPQ_Q$QO9)+ -0.102 * as.numeric(TPQ_Q$QO10)+ 0.023 * as.numeric(TPQ_Q$QO11)+ 0.327 * as.numeric(TPQ_Q$QO12)+ 0.075 * as.numeric(TPQ_Q$QO13)+ 0.076 * as.numeric(TPQ_Q$QO14)+ 0.32 * as.numeric(TPQ_Q$QO15)+ 0.208 * as.numeric(TPQ_Q$QO16)+ 0.099 * as.numeric(TPQ_Q$QO17)+ -0.007 * as.numeric(TPQ_Q$QO18)+ -0.025 * as.numeric(TPQ_Q$QO19)+ 0.051 * as.numeric(TPQ_Q$QO20)+ 0.261 * as.numeric(TPQ_Q$QO21)+ -0.079 * as.numeric(TPQ_Q$QO22)+ 0.106 * as.numeric(TPQ_Q$QO23)+ -0.022 * as.numeric(TPQ_Q$QO24)+ 0.195 * as.numeric(TPQ_Q$QO25)+ 0.186 * as.numeric(TPQ_Q$QO26)+ 0.031 * as.numeric(TPQ_Q$QO27)+ -0.178 * as.numeric(TPQ_Q$QO28)+ 0.067 * as.numeric(TPQ_Q$QO29)+ -0.233 * as.numeric(TPQ_Q$QO30)+ -0.065 * as.numeric(TPQ_Q$QO31)+ 0.214 * as.numeric(TPQ_Q$QO32)+ 0.173 * as.numeric(TPQ_Q$QO33)+ -0.05 * as.numeric(TPQ_Q$QO34)+ -0.001 * as.numeric(TPQ_Q$QO35)+ 0.045 * as.numeric(TPQ_Q$QO36)+ 0.281 * as.numeric(TPQ_Q$QO37)+ 0.362 * as.numeric(TPQ_Q$QO38)+ -0.001 * as.numeric(TPQ_Q$QO39)+ 0.189 * as.numeric(TPQ_Q$QO40)+ -0.15 * as.numeric(TPQ_Q$QO41)+ -0.027 * as.numeric(TPQ_Q$QO42)+ 0.077 * as.numeric(TPQ_Q$QO43)+ 0.018 * as.numeric(TPQ_Q$QO44)+ 0.061 * as.numeric(TPQ_Q$QO45)+ 0.105 * as.numeric(TPQ_Q$QO46)+ 0.061 * as.numeric(TPQ_Q$QO47)+ 0.085 * as.numeric(TPQ_Q$QO48)+ 0.025 * as.numeric(TPQ_Q$QO49)+ 0.037 * as.numeric(TPQ_Q$QO50)+ 0.168 * as.numeric(TPQ_Q$QO51)+ -0.155 * as.numeric(TPQ_Q$QO52)+ 0.069 * as.numeric(TPQ_Q$QO53)+ -0.08 * as.numeric(TPQ_Q$QO54)+ 0.091 * as.numeric(TPQ_Q$QO55)+ 0.27 * as.numeric(TPQ_Q$QO56)+ -0.229 * as.numeric(TPQ_Q$QO57)+ 0.095 * as.numeric(TPQ_Q$QO58)+ 0.238 * as.numeric(TPQ_Q$QO59)+ 0.117 * as.numeric(TPQ_Q$QO60)+ -0.028 * as.numeric(TPQ_Q$QO61)+ 0.078 * as.numeric(TPQ_Q$QO62)+ 0.1 * as.numeric(TPQ_Q$QO63)+ -0.076 * as.numeric(TPQ_Q$QO64)+ 0.066 * as.numeric(TPQ_Q$QO65)+ 0.139 * as.numeric(TPQ_Q$QO66)+ -0.11 * as.numeric(TPQ_Q$QO67)+ -0.365 * as.numeric(TPQ_Q$QO68)+ -0.22 * as.numeric(TPQ_Q$QO69)+ -0.153 * as.numeric(TPQ_Q$QO70)+ -0.148 * as.numeric(TPQ_Q$QO71)+ -0.201 * as.numeric(TPQ_Q$QO72)+ -0.108 * as.numeric(TPQ_Q$QO73)+ -0.146 * as.numeric(TPQ_Q$QO74)+ 0.236 * as.numeric(TPQ_Q$QO75)+ 0.153 * as.numeric(TPQ_Q$QO76)+ 0.011 * as.numeric(TPQ_Q$QO77)+ -0.052 * as.numeric(TPQ_Q$QO78)+ 0.102 * as.numeric(TPQ_Q$QO79)+ 0.334 * as.numeric(TPQ_Q$QO80)+ 0.24 * as.numeric(TPQ_Q$QO81)+ -0.171 * as.numeric(TPQ_Q$QO82)+ -0.038 * as.numeric(TPQ_Q$QO83)+ -0.012 * as.numeric(TPQ_Q$QO84)+ -0.107 * as.numeric(TPQ_Q$QO85)+ 0.239 * as.numeric(TPQ_Q$QO86)+ 0.165 * as.numeric(TPQ_Q$QO87)+ 0.241 * as.numeric(TPQ_Q$QO88)+ -0.247 * as.numeric(TPQ_Q$QO89)+ 0.304 * as.numeric(TPQ_Q$QO90)+ 0.094 * as.numeric(TPQ_Q$QO91)+ -0.215 * as.numeric(TPQ_Q$QO92)+ 0.066 * as.numeric(TPQ_Q$QO93)+ -0.067 * as.numeric(TPQ_Q$QO94)+ -0.031 * as.numeric(TPQ_Q$QO95)+ -0.084 * as.numeric(TPQ_Q$QO96)+ 0.129 * as.numeric(TPQ_Q$QO97)+ -0.068 * as.numeric(TPQ_Q$QO98)+ -0.06 * as.numeric(TPQ_Q$QO99)+ 0.172 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_6_MDD =  0.174 * as.numeric(TPQ_Q$QO1)+ -0.053 * as.numeric(TPQ_Q$QO2)+ 0.365 * as.numeric(TPQ_Q$QO3)+ 0.095 * as.numeric(TPQ_Q$QO4)+ -0.244 * as.numeric(TPQ_Q$QO5)+ 0.322 * as.numeric(TPQ_Q$QO6)+ 0.012 * as.numeric(TPQ_Q$QO7)+ 0.305 * as.numeric(TPQ_Q$QO8)+ -0.187 * as.numeric(TPQ_Q$QO9)+ 0.089 * as.numeric(TPQ_Q$QO10)+ 0.135 * as.numeric(TPQ_Q$QO11)+ -0.198 * as.numeric(TPQ_Q$QO12)+ 0.095 * as.numeric(TPQ_Q$QO13)+ 0.104 * as.numeric(TPQ_Q$QO14)+ -0.125 * as.numeric(TPQ_Q$QO15)+ -0.168 * as.numeric(TPQ_Q$QO16)+ -0.023 * as.numeric(TPQ_Q$QO17)+ -0.036 * as.numeric(TPQ_Q$QO18)+ 0.08 * as.numeric(TPQ_Q$QO19)+ -0.084 * as.numeric(TPQ_Q$QO20)+ -0.101 * as.numeric(TPQ_Q$QO21)+ -0.068 * as.numeric(TPQ_Q$QO22)+ 0.045 * as.numeric(TPQ_Q$QO23)+ 0.191 * as.numeric(TPQ_Q$QO24)+ 0.23 * as.numeric(TPQ_Q$QO25)+ 0.082 * as.numeric(TPQ_Q$QO26)+ -0.08 * as.numeric(TPQ_Q$QO27)+ 0.071 * as.numeric(TPQ_Q$QO28)+ -0.155 * as.numeric(TPQ_Q$QO29)+ 0.026 * as.numeric(TPQ_Q$QO30)+ -0.03 * as.numeric(TPQ_Q$QO31)+ -0.099 * as.numeric(TPQ_Q$QO32)+ 0.125 * as.numeric(TPQ_Q$QO33)+ -0.161 * as.numeric(TPQ_Q$QO34)+ 0.121 * as.numeric(TPQ_Q$QO35)+ -0.143 * as.numeric(TPQ_Q$QO36)+ 0.093 * as.numeric(TPQ_Q$QO37)+ 0.31 * as.numeric(TPQ_Q$QO38)+ -0.203 * as.numeric(TPQ_Q$QO39)+ -0.027 * as.numeric(TPQ_Q$QO40)+ -0.147 * as.numeric(TPQ_Q$QO41)+ 0.035 * as.numeric(TPQ_Q$QO42)+ -0.062 * as.numeric(TPQ_Q$QO43)+ 0.071 * as.numeric(TPQ_Q$QO44)+ 0.076 * as.numeric(TPQ_Q$QO45)+ -0.061 * as.numeric(TPQ_Q$QO46)+ -0.059 * as.numeric(TPQ_Q$QO47)+ -0.106 * as.numeric(TPQ_Q$QO48)+ 0.24 * as.numeric(TPQ_Q$QO49)+ -0.138 * as.numeric(TPQ_Q$QO50)+ -0.003 * as.numeric(TPQ_Q$QO51)+ 0.032 * as.numeric(TPQ_Q$QO52)+ -0.147 * as.numeric(TPQ_Q$QO53)+ -0.103 * as.numeric(TPQ_Q$QO54)+ -0.071 * as.numeric(TPQ_Q$QO55)+ -0.112 * as.numeric(TPQ_Q$QO56)+ -0.079 * as.numeric(TPQ_Q$QO57)+ 0.113 * as.numeric(TPQ_Q$QO58)+ 0.055 * as.numeric(TPQ_Q$QO59)+ 0.039 * as.numeric(TPQ_Q$QO60)+ -0.094 * as.numeric(TPQ_Q$QO61)+ -0.003 * as.numeric(TPQ_Q$QO62)+ 0.275 * as.numeric(TPQ_Q$QO63)+ -0.018 * as.numeric(TPQ_Q$QO64)+ 0.072 * as.numeric(TPQ_Q$QO65)+ 0.174 * as.numeric(TPQ_Q$QO66)+ 0.172 * as.numeric(TPQ_Q$QO67)+ -0.189 * as.numeric(TPQ_Q$QO68)+ -0.042 * as.numeric(TPQ_Q$QO69)+ 0.264 * as.numeric(TPQ_Q$QO70)+ 0.014 * as.numeric(TPQ_Q$QO71)+ 0.193 * as.numeric(TPQ_Q$QO72)+ -0.078 * as.numeric(TPQ_Q$QO73)+ 0.086 * as.numeric(TPQ_Q$QO74)+ 0.138 * as.numeric(TPQ_Q$QO75)+ 0.232 * as.numeric(TPQ_Q$QO76)+ 0.031 * as.numeric(TPQ_Q$QO77)+ 0.234 * as.numeric(TPQ_Q$QO78)+ -0.207 * as.numeric(TPQ_Q$QO79)+ 0.201 * as.numeric(TPQ_Q$QO80)+ -0.102 * as.numeric(TPQ_Q$QO81)+ -0.015 * as.numeric(TPQ_Q$QO82)+ -0.003 * as.numeric(TPQ_Q$QO83)+ -0.15 * as.numeric(TPQ_Q$QO84)+ 0.086 * as.numeric(TPQ_Q$QO85)+ -0.083 * as.numeric(TPQ_Q$QO86)+ 0.284 * as.numeric(TPQ_Q$QO87)+ -0.371 * as.numeric(TPQ_Q$QO88)+ 0.064 * as.numeric(TPQ_Q$QO89)+ 0.014 * as.numeric(TPQ_Q$QO90)+ 0.179 * as.numeric(TPQ_Q$QO91)+ -0.157 * as.numeric(TPQ_Q$QO92)+ -0.219 * as.numeric(TPQ_Q$QO93)+ 0.003 * as.numeric(TPQ_Q$QO94)+ 0.075 * as.numeric(TPQ_Q$QO95)+ 0.103 * as.numeric(TPQ_Q$QO96)+ 0.044 * as.numeric(TPQ_Q$QO97)+ 0.055 * as.numeric(TPQ_Q$QO98)+ -0.011 * as.numeric(TPQ_Q$QO99)+ 0.133 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_7_MDD =  -0.07 * as.numeric(TPQ_Q$QO1)+ -0.338 * as.numeric(TPQ_Q$QO2)+ -0.267 * as.numeric(TPQ_Q$QO3)+ -0.336 * as.numeric(TPQ_Q$QO4)+ -0.028 * as.numeric(TPQ_Q$QO5)+ -0.318 * as.numeric(TPQ_Q$QO6)+ -0.172 * as.numeric(TPQ_Q$QO7)+ -0.059 * as.numeric(TPQ_Q$QO8)+ -0.208 * as.numeric(TPQ_Q$QO9)+ -0.083 * as.numeric(TPQ_Q$QO10)+ -0.002 * as.numeric(TPQ_Q$QO11)+ 0.056 * as.numeric(TPQ_Q$QO12)+ -0.131 * as.numeric(TPQ_Q$QO13)+ -0.229 * as.numeric(TPQ_Q$QO14)+ -0.015 * as.numeric(TPQ_Q$QO15)+ -0.205 * as.numeric(TPQ_Q$QO16)+ -0.213 * as.numeric(TPQ_Q$QO17)+ -0.009 * as.numeric(TPQ_Q$QO18)+ -0.066 * as.numeric(TPQ_Q$QO19)+ -0.276 * as.numeric(TPQ_Q$QO20)+ -0.012 * as.numeric(TPQ_Q$QO21)+ -0.049 * as.numeric(TPQ_Q$QO22)+ -0.044 * as.numeric(TPQ_Q$QO23)+ -0.362 * as.numeric(TPQ_Q$QO24)+ -0.231 * as.numeric(TPQ_Q$QO25)+ -0.273 * as.numeric(TPQ_Q$QO26)+ -0.141 * as.numeric(TPQ_Q$QO27)+ -0.105 * as.numeric(TPQ_Q$QO28)+ -0.443 * as.numeric(TPQ_Q$QO29)+ -0.138 * as.numeric(TPQ_Q$QO30)+ -0.199 * as.numeric(TPQ_Q$QO31)+ -0.028 * as.numeric(TPQ_Q$QO32)+ 0.083 * as.numeric(TPQ_Q$QO33)+ -0.248 * as.numeric(TPQ_Q$QO34)+ 0.014 * as.numeric(TPQ_Q$QO35)+ -0.063 * as.numeric(TPQ_Q$QO36)+ 0.056 * as.numeric(TPQ_Q$QO37)+ -0.018 * as.numeric(TPQ_Q$QO38)+ -0.138 * as.numeric(TPQ_Q$QO39)+ 0.007 * as.numeric(TPQ_Q$QO40)+ -0.184 * as.numeric(TPQ_Q$QO41)+ -0.016 * as.numeric(TPQ_Q$QO42)+ -0.132 * as.numeric(TPQ_Q$QO43)+ -0.168 * as.numeric(TPQ_Q$QO44)+ -0.106 * as.numeric(TPQ_Q$QO45)+ 0.042 * as.numeric(TPQ_Q$QO46)+ -0.432 * as.numeric(TPQ_Q$QO47)+ -0.235 * as.numeric(TPQ_Q$QO48)+ -0.042 * as.numeric(TPQ_Q$QO49)+ -0.163 * as.numeric(TPQ_Q$QO50)+ -0.296 * as.numeric(TPQ_Q$QO51)+ 0.072 * as.numeric(TPQ_Q$QO52)+ -0.039 * as.numeric(TPQ_Q$QO53)+ -0.063 * as.numeric(TPQ_Q$QO54)+ 0.075 * as.numeric(TPQ_Q$QO55)+ -0.018 * as.numeric(TPQ_Q$QO56)+ 0.035 * as.numeric(TPQ_Q$QO57)+ -0.301 * as.numeric(TPQ_Q$QO58)+ 0.061 * as.numeric(TPQ_Q$QO59)+ -0.352 * as.numeric(TPQ_Q$QO60)+ -0.106 * as.numeric(TPQ_Q$QO61)+ -0.386 * as.numeric(TPQ_Q$QO62)+ 0.101 * as.numeric(TPQ_Q$QO63)+ -0.164 * as.numeric(TPQ_Q$QO64)+ 0.027 * as.numeric(TPQ_Q$QO65)+ 0.096 * as.numeric(TPQ_Q$QO66)+ -0.219 * as.numeric(TPQ_Q$QO67)+ -0.022 * as.numeric(TPQ_Q$QO68)+ 0.019 * as.numeric(TPQ_Q$QO69)+ -0.182 * as.numeric(TPQ_Q$QO70)+ -0.194 * as.numeric(TPQ_Q$QO71)+ -0.167 * as.numeric(TPQ_Q$QO72)+ -0.018 * as.numeric(TPQ_Q$QO73)+ -0.205 * as.numeric(TPQ_Q$QO74)+ -0.253 * as.numeric(TPQ_Q$QO75)+ 0.042 * as.numeric(TPQ_Q$QO76)+ -0.104 * as.numeric(TPQ_Q$QO77)+ -0.026 * as.numeric(TPQ_Q$QO78)+ -0.158 * as.numeric(TPQ_Q$QO79)+ 0.066 * as.numeric(TPQ_Q$QO80)+ -0.112 * as.numeric(TPQ_Q$QO81)+ 0.141 * as.numeric(TPQ_Q$QO82)+ -0.12 * as.numeric(TPQ_Q$QO83)+ -0.032 * as.numeric(TPQ_Q$QO84)+ -0.028 * as.numeric(TPQ_Q$QO85)+ 0.357 * as.numeric(TPQ_Q$QO86)+ 0.121 * as.numeric(TPQ_Q$QO87)+ 0.312 * as.numeric(TPQ_Q$QO88)+ -0.145 * as.numeric(TPQ_Q$QO89)+ 0.107 * as.numeric(TPQ_Q$QO90)+ 0.067 * as.numeric(TPQ_Q$QO91)+ -0.07 * as.numeric(TPQ_Q$QO92)+ 0.174 * as.numeric(TPQ_Q$QO93)+ -0.098 * as.numeric(TPQ_Q$QO94)+ 0.103 * as.numeric(TPQ_Q$QO95)+ 0.134 * as.numeric(TPQ_Q$QO96)+ -0.18 * as.numeric(TPQ_Q$QO97)+ 0.111 * as.numeric(TPQ_Q$QO98)+ -0.045 * as.numeric(TPQ_Q$QO99)+ -0.199 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_8_MDD =  0.202 * as.numeric(TPQ_Q$QO1)+ 0.134 * as.numeric(TPQ_Q$QO2)+ 0.115 * as.numeric(TPQ_Q$QO3)+ 0.187 * as.numeric(TPQ_Q$QO4)+ -0.064 * as.numeric(TPQ_Q$QO5)+ 0.051 * as.numeric(TPQ_Q$QO6)+ -0.14 * as.numeric(TPQ_Q$QO7)+ 0.095 * as.numeric(TPQ_Q$QO8)+ -0.042 * as.numeric(TPQ_Q$QO9)+ -0.138 * as.numeric(TPQ_Q$QO10)+ -0.236 * as.numeric(TPQ_Q$QO11)+ -0.283 * as.numeric(TPQ_Q$QO12)+ -0.051 * as.numeric(TPQ_Q$QO13)+ 0.036 * as.numeric(TPQ_Q$QO14)+ -0.28 * as.numeric(TPQ_Q$QO15)+ -0.077 * as.numeric(TPQ_Q$QO16)+ -0.126 * as.numeric(TPQ_Q$QO17)+ -0.123 * as.numeric(TPQ_Q$QO18)+ -0.1 * as.numeric(TPQ_Q$QO19)+ -0.203 * as.numeric(TPQ_Q$QO20)+ -0.157 * as.numeric(TPQ_Q$QO21)+ 0.102 * as.numeric(TPQ_Q$QO22)+ -0.094 * as.numeric(TPQ_Q$QO23)+ 0.154 * as.numeric(TPQ_Q$QO24)+ -0.287 * as.numeric(TPQ_Q$QO25)+ 0.016 * as.numeric(TPQ_Q$QO26)+ -0.016 * as.numeric(TPQ_Q$QO27)+ -0.235 * as.numeric(TPQ_Q$QO28)+ 0.071 * as.numeric(TPQ_Q$QO29)+ -0.113 * as.numeric(TPQ_Q$QO30)+ 0.092 * as.numeric(TPQ_Q$QO31)+ 0.135 * as.numeric(TPQ_Q$QO32)+ -0.157 * as.numeric(TPQ_Q$QO33)+ -0.035 * as.numeric(TPQ_Q$QO34)+ 0.097 * as.numeric(TPQ_Q$QO35)+ -0.298 * as.numeric(TPQ_Q$QO36)+ -0.329 * as.numeric(TPQ_Q$QO37)+ -0.171 * as.numeric(TPQ_Q$QO38)+ -0.119 * as.numeric(TPQ_Q$QO39)+ -0.301 * as.numeric(TPQ_Q$QO40)+ -0.168 * as.numeric(TPQ_Q$QO41)+ -0.031 * as.numeric(TPQ_Q$QO42)+ -0.028 * as.numeric(TPQ_Q$QO43)+ -0.03 * as.numeric(TPQ_Q$QO44)+ -0.265 * as.numeric(TPQ_Q$QO45)+ -0.124 * as.numeric(TPQ_Q$QO46)+ 0.175 * as.numeric(TPQ_Q$QO47)+ 0.052 * as.numeric(TPQ_Q$QO48)+ -0.072 * as.numeric(TPQ_Q$QO49)+ -0.081 * as.numeric(TPQ_Q$QO50)+ 0.171 * as.numeric(TPQ_Q$QO51)+ -0.075 * as.numeric(TPQ_Q$QO52)+ -0.037 * as.numeric(TPQ_Q$QO53)+ -0.178 * as.numeric(TPQ_Q$QO54)+ -0.029 * as.numeric(TPQ_Q$QO55)+ -0.233 * as.numeric(TPQ_Q$QO56)+ -0.283 * as.numeric(TPQ_Q$QO57)+ -0.115 * as.numeric(TPQ_Q$QO58)+ 0.07 * as.numeric(TPQ_Q$QO59)+ 0.067 * as.numeric(TPQ_Q$QO60)+ -0.033 * as.numeric(TPQ_Q$QO61)+ 0.142 * as.numeric(TPQ_Q$QO62)+ 0.01 * as.numeric(TPQ_Q$QO63)+ -0.041 * as.numeric(TPQ_Q$QO64)+ -0.244 * as.numeric(TPQ_Q$QO65)+ -0.213 * as.numeric(TPQ_Q$QO66)+ -0.056 * as.numeric(TPQ_Q$QO67)+ -0.252 * as.numeric(TPQ_Q$QO68)+ -0.275 * as.numeric(TPQ_Q$QO69)+ 0.169 * as.numeric(TPQ_Q$QO70)+ -0.354 * as.numeric(TPQ_Q$QO71)+ 0.131 * as.numeric(TPQ_Q$QO72)+ -0.213 * as.numeric(TPQ_Q$QO73)+ -0.008 * as.numeric(TPQ_Q$QO74)+ 0.202 * as.numeric(TPQ_Q$QO75)+ -0.236 * as.numeric(TPQ_Q$QO76)+ -0.09 * as.numeric(TPQ_Q$QO77)+ -0.499 * as.numeric(TPQ_Q$QO78)+ -0.033 * as.numeric(TPQ_Q$QO79)+ 0.169 * as.numeric(TPQ_Q$QO80)+ 0.009 * as.numeric(TPQ_Q$QO81)+ 0.146 * as.numeric(TPQ_Q$QO82)+ -0.041 * as.numeric(TPQ_Q$QO83)+ 0.042 * as.numeric(TPQ_Q$QO84)+ -0.164 * as.numeric(TPQ_Q$QO85)+ -0.283 * as.numeric(TPQ_Q$QO86)+ -0.355 * as.numeric(TPQ_Q$QO87)+ -0.215 * as.numeric(TPQ_Q$QO88)+ -0.036 * as.numeric(TPQ_Q$QO89)+ -0.121 * as.numeric(TPQ_Q$QO90)+ -0.134 * as.numeric(TPQ_Q$QO91)+ 0.028 * as.numeric(TPQ_Q$QO92)+ -0.05 * as.numeric(TPQ_Q$QO93)+ -0.035 * as.numeric(TPQ_Q$QO94)+ 0.116 * as.numeric(TPQ_Q$QO95)+ -0.569 * as.numeric(TPQ_Q$QO96)+ 0.08 * as.numeric(TPQ_Q$QO97)+ 0.102 * as.numeric(TPQ_Q$QO98)+ -0.078 * as.numeric(TPQ_Q$QO99)+ -0.022 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_9_MDD =  -0.111 * as.numeric(TPQ_Q$QO1)+ 0.045 * as.numeric(TPQ_Q$QO2)+ -0.041 * as.numeric(TPQ_Q$QO3)+ -0.055 * as.numeric(TPQ_Q$QO4)+ 0.036 * as.numeric(TPQ_Q$QO5)+ -0.119 * as.numeric(TPQ_Q$QO6)+ 0.123 * as.numeric(TPQ_Q$QO7)+ 0.077 * as.numeric(TPQ_Q$QO8)+ -0.223 * as.numeric(TPQ_Q$QO9)+ 0.153 * as.numeric(TPQ_Q$QO10)+ 0.049 * as.numeric(TPQ_Q$QO11)+ 0.098 * as.numeric(TPQ_Q$QO12)+ -0.156 * as.numeric(TPQ_Q$QO13)+ 0.106 * as.numeric(TPQ_Q$QO14)+ -0.014 * as.numeric(TPQ_Q$QO15)+ -0.106 * as.numeric(TPQ_Q$QO16)+ -0.155 * as.numeric(TPQ_Q$QO17)+ -0.009 * as.numeric(TPQ_Q$QO18)+ -0.052 * as.numeric(TPQ_Q$QO19)+ -0.215 * as.numeric(TPQ_Q$QO20)+ -0.184 * as.numeric(TPQ_Q$QO21)+ 0.125 * as.numeric(TPQ_Q$QO22)+ 0.025 * as.numeric(TPQ_Q$QO23)+ -0.027 * as.numeric(TPQ_Q$QO24)+ 0.079 * as.numeric(TPQ_Q$QO25)+ -0.055 * as.numeric(TPQ_Q$QO26)+ 0.002 * as.numeric(TPQ_Q$QO27)+ 0.079 * as.numeric(TPQ_Q$QO28)+ 0.013 * as.numeric(TPQ_Q$QO29)+ -0.012 * as.numeric(TPQ_Q$QO30)+ 0.061 * as.numeric(TPQ_Q$QO31)+ -0.248 * as.numeric(TPQ_Q$QO32)+ -0.017 * as.numeric(TPQ_Q$QO33)+ -0.052 * as.numeric(TPQ_Q$QO34)+ 0.019 * as.numeric(TPQ_Q$QO35)+ -0.195 * as.numeric(TPQ_Q$QO36)+ -0.13 * as.numeric(TPQ_Q$QO37)+ 0.012 * as.numeric(TPQ_Q$QO38)+ 0.047 * as.numeric(TPQ_Q$QO39)+ 0.137 * as.numeric(TPQ_Q$QO40)+ 0.019 * as.numeric(TPQ_Q$QO41)+ 0.003 * as.numeric(TPQ_Q$QO42)+ -0.032 * as.numeric(TPQ_Q$QO43)+ 0.044 * as.numeric(TPQ_Q$QO44)+ 0.02 * as.numeric(TPQ_Q$QO45)+ -0.383 * as.numeric(TPQ_Q$QO46)+ 0.098 * as.numeric(TPQ_Q$QO47)+ 0.246 * as.numeric(TPQ_Q$QO48)+ 0.007 * as.numeric(TPQ_Q$QO49)+ 0.267 * as.numeric(TPQ_Q$QO50)+ 0.025 * as.numeric(TPQ_Q$QO51)+ -0.101 * as.numeric(TPQ_Q$QO52)+ 0.029 * as.numeric(TPQ_Q$QO53)+ -0.156 * as.numeric(TPQ_Q$QO54)+ -0.255 * as.numeric(TPQ_Q$QO55)+ -0.245 * as.numeric(TPQ_Q$QO56)+ -0.229 * as.numeric(TPQ_Q$QO57)+ -0.153 * as.numeric(TPQ_Q$QO58)+ 0.346 * as.numeric(TPQ_Q$QO59)+ -0.033 * as.numeric(TPQ_Q$QO60)+ -0.034 * as.numeric(TPQ_Q$QO61)+ 0.085 * as.numeric(TPQ_Q$QO62)+ -0.079 * as.numeric(TPQ_Q$QO63)+ 0.1 * as.numeric(TPQ_Q$QO64)+ -0.08 * as.numeric(TPQ_Q$QO65)+ 0.115 * as.numeric(TPQ_Q$QO66)+ 0.121 * as.numeric(TPQ_Q$QO67)+ 0.001 * as.numeric(TPQ_Q$QO68)+ 0.035 * as.numeric(TPQ_Q$QO69)+ -0.095 * as.numeric(TPQ_Q$QO70)+ 0.194 * as.numeric(TPQ_Q$QO71)+ -0.032 * as.numeric(TPQ_Q$QO72)+ -0.029 * as.numeric(TPQ_Q$QO73)+ 0.174 * as.numeric(TPQ_Q$QO74)+ -0.006 * as.numeric(TPQ_Q$QO75)+ -0.068 * as.numeric(TPQ_Q$QO76)+ 0.01 * as.numeric(TPQ_Q$QO77)+ -0.048 * as.numeric(TPQ_Q$QO78)+ 0.296 * as.numeric(TPQ_Q$QO79)+ 0.101 * as.numeric(TPQ_Q$QO80)+ -0.105 * as.numeric(TPQ_Q$QO81)+ -0.145 * as.numeric(TPQ_Q$QO82)+ 0.014 * as.numeric(TPQ_Q$QO83)+ 0.158 * as.numeric(TPQ_Q$QO84)+ -0.133 * as.numeric(TPQ_Q$QO85)+ -0.089 * as.numeric(TPQ_Q$QO86)+ 0.076 * as.numeric(TPQ_Q$QO87)+ 0.007 * as.numeric(TPQ_Q$QO88)+ 0.018 * as.numeric(TPQ_Q$QO89)+ -0.192 * as.numeric(TPQ_Q$QO90)+ 0.31 * as.numeric(TPQ_Q$QO91)+ -0.144 * as.numeric(TPQ_Q$QO92)+ 0.066 * as.numeric(TPQ_Q$QO93)+ 0.067 * as.numeric(TPQ_Q$QO94)+ 0.019 * as.numeric(TPQ_Q$QO95)+ 0.134 * as.numeric(TPQ_Q$QO96)+ -0.028 * as.numeric(TPQ_Q$QO97)+ 0.009 * as.numeric(TPQ_Q$QO98)+ -0.183 * as.numeric(TPQ_Q$QO99)+ 0.031 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_10_MDD =  0.025 * as.numeric(TPQ_Q$QO1)+ 0.238 * as.numeric(TPQ_Q$QO2)+ -0.279 * as.numeric(TPQ_Q$QO3)+ -0.037 * as.numeric(TPQ_Q$QO4)+ 0.014 * as.numeric(TPQ_Q$QO5)+ -0.277 * as.numeric(TPQ_Q$QO6)+ -0.097 * as.numeric(TPQ_Q$QO7)+ 0.048 * as.numeric(TPQ_Q$QO8)+ -0.041 * as.numeric(TPQ_Q$QO9)+ 0.028 * as.numeric(TPQ_Q$QO10)+ 0.183 * as.numeric(TPQ_Q$QO11)+ -0.001 * as.numeric(TPQ_Q$QO12)+ -0.158 * as.numeric(TPQ_Q$QO13)+ -0.067 * as.numeric(TPQ_Q$QO14)+ -0.056 * as.numeric(TPQ_Q$QO15)+ 0.052 * as.numeric(TPQ_Q$QO16)+ -0.119 * as.numeric(TPQ_Q$QO17)+ 0.011 * as.numeric(TPQ_Q$QO18)+ 0.161 * as.numeric(TPQ_Q$QO19)+ -0.033 * as.numeric(TPQ_Q$QO20)+ -0.169 * as.numeric(TPQ_Q$QO21)+ 0.005 * as.numeric(TPQ_Q$QO22)+ 0.103 * as.numeric(TPQ_Q$QO23)+ 0.211 * as.numeric(TPQ_Q$QO24)+ -0.166 * as.numeric(TPQ_Q$QO25)+ 0.124 * as.numeric(TPQ_Q$QO26)+ -0.028 * as.numeric(TPQ_Q$QO27)+ 0.072 * as.numeric(TPQ_Q$QO28)+ 0.41 * as.numeric(TPQ_Q$QO29)+ 0.136 * as.numeric(TPQ_Q$QO30)+ -0.124 * as.numeric(TPQ_Q$QO31)+ -0.033 * as.numeric(TPQ_Q$QO32)+ 0.29 * as.numeric(TPQ_Q$QO33)+ -0.06 * as.numeric(TPQ_Q$QO34)+ -0.018 * as.numeric(TPQ_Q$QO35)+ -0.013 * as.numeric(TPQ_Q$QO36)+ 0.234 * as.numeric(TPQ_Q$QO37)+ 0.12 * as.numeric(TPQ_Q$QO38)+ -0.016 * as.numeric(TPQ_Q$QO39)+ -0.129 * as.numeric(TPQ_Q$QO40)+ -0.078 * as.numeric(TPQ_Q$QO41)+ -0.062 * as.numeric(TPQ_Q$QO42)+ 0.036 * as.numeric(TPQ_Q$QO43)+ -0.115 * as.numeric(TPQ_Q$QO44)+ 0.205 * as.numeric(TPQ_Q$QO45)+ -0.014 * as.numeric(TPQ_Q$QO46)+ 0.406 * as.numeric(TPQ_Q$QO47)+ 0.056 * as.numeric(TPQ_Q$QO48)+ 0.027 * as.numeric(TPQ_Q$QO49)+ 0.128 * as.numeric(TPQ_Q$QO50)+ 0.339 * as.numeric(TPQ_Q$QO51)+ 0.022 * as.numeric(TPQ_Q$QO52)+ 0.021 * as.numeric(TPQ_Q$QO53)+ -0.066 * as.numeric(TPQ_Q$QO54)+ -0.02 * as.numeric(TPQ_Q$QO55)+ -0.088 * as.numeric(TPQ_Q$QO56)+ -0.075 * as.numeric(TPQ_Q$QO57)+ -0.085 * as.numeric(TPQ_Q$QO58)+ -0.014 * as.numeric(TPQ_Q$QO59)+ 0.075 * as.numeric(TPQ_Q$QO60)+ 0.043 * as.numeric(TPQ_Q$QO61)+ -0.146 * as.numeric(TPQ_Q$QO62)+ 0.014 * as.numeric(TPQ_Q$QO63)+ -0.151 * as.numeric(TPQ_Q$QO64)+ -0.1 * as.numeric(TPQ_Q$QO65)+ 0.225 * as.numeric(TPQ_Q$QO66)+ -0.088 * as.numeric(TPQ_Q$QO67)+ 0.117 * as.numeric(TPQ_Q$QO68)+ -0.113 * as.numeric(TPQ_Q$QO69)+ -0.054 * as.numeric(TPQ_Q$QO70)+ 0.172 * as.numeric(TPQ_Q$QO71)+ -0.095 * as.numeric(TPQ_Q$QO72)+ 0.009 * as.numeric(TPQ_Q$QO73)+ -0.034 * as.numeric(TPQ_Q$QO74)+ -0.25 * as.numeric(TPQ_Q$QO75)+ 0.349 * as.numeric(TPQ_Q$QO76)+ 0.046 * as.numeric(TPQ_Q$QO77)+ -0.018 * as.numeric(TPQ_Q$QO78)+ -0.102 * as.numeric(TPQ_Q$QO79)+ -0.114 * as.numeric(TPQ_Q$QO80)+ -0.128 * as.numeric(TPQ_Q$QO81)+ -0.041 * as.numeric(TPQ_Q$QO82)+ 0.053 * as.numeric(TPQ_Q$QO83)+ -0.292 * as.numeric(TPQ_Q$QO84)+ 0.084 * as.numeric(TPQ_Q$QO85)+ 0.181 * as.numeric(TPQ_Q$QO86)+ 0.218 * as.numeric(TPQ_Q$QO87)+ 0.18 * as.numeric(TPQ_Q$QO88)+ -0.084 * as.numeric(TPQ_Q$QO89)+ 0.098 * as.numeric(TPQ_Q$QO90)+ -0.092 * as.numeric(TPQ_Q$QO91)+ 0.053 * as.numeric(TPQ_Q$QO92)+ 0.102 * as.numeric(TPQ_Q$QO93)+ -0.058 * as.numeric(TPQ_Q$QO94)+ 0.091 * as.numeric(TPQ_Q$QO95)+ -0.009 * as.numeric(TPQ_Q$QO96)+ 0.181 * as.numeric(TPQ_Q$QO97)+ -0.106 * as.numeric(TPQ_Q$QO98)+ 0.06 * as.numeric(TPQ_Q$QO99)+ -0.154 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_11_MDD =  -0.143 * as.numeric(TPQ_Q$QO1)+ 0.296 * as.numeric(TPQ_Q$QO2)+ -0.066 * as.numeric(TPQ_Q$QO3)+ -0.039 * as.numeric(TPQ_Q$QO4)+ -0.022 * as.numeric(TPQ_Q$QO5)+ 0.088 * as.numeric(TPQ_Q$QO6)+ -0.078 * as.numeric(TPQ_Q$QO7)+ 0.089 * as.numeric(TPQ_Q$QO8)+ -0.057 * as.numeric(TPQ_Q$QO9)+ 0.248 * as.numeric(TPQ_Q$QO10)+ 0.3 * as.numeric(TPQ_Q$QO11)+ 0.138 * as.numeric(TPQ_Q$QO12)+ 0.132 * as.numeric(TPQ_Q$QO13)+ 0.288 * as.numeric(TPQ_Q$QO14)+ 0.207 * as.numeric(TPQ_Q$QO15)+ -0.101 * as.numeric(TPQ_Q$QO16)+ 0.106 * as.numeric(TPQ_Q$QO17)+ -0.015 * as.numeric(TPQ_Q$QO18)+ 0.073 * as.numeric(TPQ_Q$QO19)+ 0.067 * as.numeric(TPQ_Q$QO20)+ -0.002 * as.numeric(TPQ_Q$QO21)+ -0.11 * as.numeric(TPQ_Q$QO22)+ 0.135 * as.numeric(TPQ_Q$QO23)+ -0.128 * as.numeric(TPQ_Q$QO24)+ 0.089 * as.numeric(TPQ_Q$QO25)+ 0.02 * as.numeric(TPQ_Q$QO26)+ -0.058 * as.numeric(TPQ_Q$QO27)+ -0.041 * as.numeric(TPQ_Q$QO28)+ -0.161 * as.numeric(TPQ_Q$QO29)+ 0.021 * as.numeric(TPQ_Q$QO30)+ -0.145 * as.numeric(TPQ_Q$QO31)+ 0.128 * as.numeric(TPQ_Q$QO32)+ 0.108 * as.numeric(TPQ_Q$QO33)+ 0.003 * as.numeric(TPQ_Q$QO34)+ 0.144 * as.numeric(TPQ_Q$QO35)+ 0.095 * as.numeric(TPQ_Q$QO36)+ -0.082 * as.numeric(TPQ_Q$QO37)+ 0.103 * as.numeric(TPQ_Q$QO38)+ 0.04 * as.numeric(TPQ_Q$QO39)+ -0.12 * as.numeric(TPQ_Q$QO40)+ 0.174 * as.numeric(TPQ_Q$QO41)+ 0.037 * as.numeric(TPQ_Q$QO42)+ 0.165 * as.numeric(TPQ_Q$QO43)+ -0.176 * as.numeric(TPQ_Q$QO44)+ -0.159 * as.numeric(TPQ_Q$QO45)+ -0.12 * as.numeric(TPQ_Q$QO46)+ -0.068 * as.numeric(TPQ_Q$QO47)+ -0.057 * as.numeric(TPQ_Q$QO48)+ -0.168 * as.numeric(TPQ_Q$QO49)+ -0.122 * as.numeric(TPQ_Q$QO50)+ -0.017 * as.numeric(TPQ_Q$QO51)+ 0.002 * as.numeric(TPQ_Q$QO52)+ -0.017 * as.numeric(TPQ_Q$QO53)+ -0.137 * as.numeric(TPQ_Q$QO54)+ -0.071 * as.numeric(TPQ_Q$QO55)+ -0.124 * as.numeric(TPQ_Q$QO56)+ -0.096 * as.numeric(TPQ_Q$QO57)+ -0.037 * as.numeric(TPQ_Q$QO58)+ 0.099 * as.numeric(TPQ_Q$QO59)+ -0.079 * as.numeric(TPQ_Q$QO60)+ -0.035 * as.numeric(TPQ_Q$QO61)+ -0.254 * as.numeric(TPQ_Q$QO62)+ 0.207 * as.numeric(TPQ_Q$QO63)+ -0.034 * as.numeric(TPQ_Q$QO64)+ 0.329 * as.numeric(TPQ_Q$QO65)+ 0.017 * as.numeric(TPQ_Q$QO66)+ 0.069 * as.numeric(TPQ_Q$QO67)+ 0.226 * as.numeric(TPQ_Q$QO68)+ 0.024 * as.numeric(TPQ_Q$QO69)+ 0.155 * as.numeric(TPQ_Q$QO70)+ -0.062 * as.numeric(TPQ_Q$QO71)+ 0.121 * as.numeric(TPQ_Q$QO72)+ 0.01 * as.numeric(TPQ_Q$QO73)+ -0.214 * as.numeric(TPQ_Q$QO74)+ -0.009 * as.numeric(TPQ_Q$QO75)+ -0.135 * as.numeric(TPQ_Q$QO76)+ -0.031 * as.numeric(TPQ_Q$QO77)+ -0.011 * as.numeric(TPQ_Q$QO78)+ 0.225 * as.numeric(TPQ_Q$QO79)+ -0.11 * as.numeric(TPQ_Q$QO80)+ 0.075 * as.numeric(TPQ_Q$QO81)+ 0.071 * as.numeric(TPQ_Q$QO82)+ 0.073 * as.numeric(TPQ_Q$QO83)+ -0.049 * as.numeric(TPQ_Q$QO84)+ -0.021 * as.numeric(TPQ_Q$QO85)+ 0.038 * as.numeric(TPQ_Q$QO86)+ -0.093 * as.numeric(TPQ_Q$QO87)+ -0.007 * as.numeric(TPQ_Q$QO88)+ 0.143 * as.numeric(TPQ_Q$QO89)+ 0.168 * as.numeric(TPQ_Q$QO90)+ 0.149 * as.numeric(TPQ_Q$QO91)+ 0.383 * as.numeric(TPQ_Q$QO92)+ -0.086 * as.numeric(TPQ_Q$QO93)+ -0.006 * as.numeric(TPQ_Q$QO94)+ -0.071 * as.numeric(TPQ_Q$QO95)+ -0.018 * as.numeric(TPQ_Q$QO96)+ 0.423 * as.numeric(TPQ_Q$QO97)+ -0.048 * as.numeric(TPQ_Q$QO98)+ -0.022 * as.numeric(TPQ_Q$QO99)+ 0.024 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_12_MDD =  0.57 * as.numeric(TPQ_Q$QO1)+ 0.216 * as.numeric(TPQ_Q$QO2)+ 0.119 * as.numeric(TPQ_Q$QO3)+ 0.37 * as.numeric(TPQ_Q$QO4)+ -0.111 * as.numeric(TPQ_Q$QO5)+ -0.009 * as.numeric(TPQ_Q$QO6)+ 0.155 * as.numeric(TPQ_Q$QO7)+ 0.188 * as.numeric(TPQ_Q$QO8)+ 0.116 * as.numeric(TPQ_Q$QO9)+ -0.04 * as.numeric(TPQ_Q$QO10)+ -0.059 * as.numeric(TPQ_Q$QO11)+ 0.075 * as.numeric(TPQ_Q$QO12)+ -0.077 * as.numeric(TPQ_Q$QO13)+ -0.151 * as.numeric(TPQ_Q$QO14)+ -0.057 * as.numeric(TPQ_Q$QO15)+ 0.148 * as.numeric(TPQ_Q$QO16)+ 0.117 * as.numeric(TPQ_Q$QO17)+ -0.061 * as.numeric(TPQ_Q$QO18)+ -0.14 * as.numeric(TPQ_Q$QO19)+ 0.139 * as.numeric(TPQ_Q$QO20)+ 0.106 * as.numeric(TPQ_Q$QO21)+ 0.126 * as.numeric(TPQ_Q$QO22)+ -0.153 * as.numeric(TPQ_Q$QO23)+ -0.051 * as.numeric(TPQ_Q$QO24)+ 0.186 * as.numeric(TPQ_Q$QO25)+ -0.045 * as.numeric(TPQ_Q$QO26)+ 0.021 * as.numeric(TPQ_Q$QO27)+ -0.171 * as.numeric(TPQ_Q$QO28)+ 0.148 * as.numeric(TPQ_Q$QO29)+ -0.102 * as.numeric(TPQ_Q$QO30)+ 0.127 * as.numeric(TPQ_Q$QO31)+ 0.326 * as.numeric(TPQ_Q$QO32)+ -0.145 * as.numeric(TPQ_Q$QO33)+ -0.031 * as.numeric(TPQ_Q$QO34)+ -0.064 * as.numeric(TPQ_Q$QO35)+ 0.318 * as.numeric(TPQ_Q$QO36)+ -0.195 * as.numeric(TPQ_Q$QO37)+ -0.228 * as.numeric(TPQ_Q$QO38)+ -0.057 * as.numeric(TPQ_Q$QO39)+ -0.175 * as.numeric(TPQ_Q$QO40)+ -0.0 * as.numeric(TPQ_Q$QO41)+ 0.301 * as.numeric(TPQ_Q$QO42)+ -0.108 * as.numeric(TPQ_Q$QO43)+ 0.436 * as.numeric(TPQ_Q$QO44)+ 0.303 * as.numeric(TPQ_Q$QO45)+ 0.202 * as.numeric(TPQ_Q$QO46)+ 0.103 * as.numeric(TPQ_Q$QO47)+ -0.059 * as.numeric(TPQ_Q$QO48)+ 0.07 * as.numeric(TPQ_Q$QO49)+ -0.076 * as.numeric(TPQ_Q$QO50)+ 0.179 * as.numeric(TPQ_Q$QO51)+ 0.561 * as.numeric(TPQ_Q$QO52)+ 0.086 * as.numeric(TPQ_Q$QO53)+ -0.092 * as.numeric(TPQ_Q$QO54)+ 0.215 * as.numeric(TPQ_Q$QO55)+ 0.196 * as.numeric(TPQ_Q$QO56)+ -0.051 * as.numeric(TPQ_Q$QO57)+ 0.093 * as.numeric(TPQ_Q$QO58)+ 0.297 * as.numeric(TPQ_Q$QO59)+ -0.01 * as.numeric(TPQ_Q$QO60)+ -0.062 * as.numeric(TPQ_Q$QO61)+ 0.155 * as.numeric(TPQ_Q$QO62)+ 0.402 * as.numeric(TPQ_Q$QO63)+ 0.026 * as.numeric(TPQ_Q$QO64)+ 0.105 * as.numeric(TPQ_Q$QO65)+ 0.228 * as.numeric(TPQ_Q$QO66)+ 0.125 * as.numeric(TPQ_Q$QO67)+ -0.204 * as.numeric(TPQ_Q$QO68)+ -0.093 * as.numeric(TPQ_Q$QO69)+ -0.118 * as.numeric(TPQ_Q$QO70)+ 0.145 * as.numeric(TPQ_Q$QO71)+ -0.111 * as.numeric(TPQ_Q$QO72)+ -0.109 * as.numeric(TPQ_Q$QO73)+ -0.008 * as.numeric(TPQ_Q$QO74)+ 0.435 * as.numeric(TPQ_Q$QO75)+ 0.104 * as.numeric(TPQ_Q$QO76)+ 0.28 * as.numeric(TPQ_Q$QO77)+ -0.056 * as.numeric(TPQ_Q$QO78)+ 0.398 * as.numeric(TPQ_Q$QO79)+ 0.366 * as.numeric(TPQ_Q$QO80)+ -0.03 * as.numeric(TPQ_Q$QO81)+ 0.612 * as.numeric(TPQ_Q$QO82)+ 0.036 * as.numeric(TPQ_Q$QO83)+ 0.505 * as.numeric(TPQ_Q$QO84)+ -0.08 * as.numeric(TPQ_Q$QO85)+ 0.157 * as.numeric(TPQ_Q$QO86)+ 0.104 * as.numeric(TPQ_Q$QO87)+ 0.117 * as.numeric(TPQ_Q$QO88)+ 0.553 * as.numeric(TPQ_Q$QO89)+ -0.169 * as.numeric(TPQ_Q$QO90)+ 0.192 * as.numeric(TPQ_Q$QO91)+ 0.214 * as.numeric(TPQ_Q$QO92)+ 0.148 * as.numeric(TPQ_Q$QO93)+ 0.106 * as.numeric(TPQ_Q$QO94)+ 0.699 * as.numeric(TPQ_Q$QO95)+ -0.091 * as.numeric(TPQ_Q$QO96)+ 0.345 * as.numeric(TPQ_Q$QO97)+ 0.547 * as.numeric(TPQ_Q$QO98)+ 0.12 * as.numeric(TPQ_Q$QO99)+ 0.359 * as.numeric(TPQ_Q$QO100)


TPQ_Q$Complete_IC_0_All =  0.255 * as.numeric(TPQ_Q$QO1)+ 0.227 * as.numeric(TPQ_Q$QO2)+ 0.048 * as.numeric(TPQ_Q$QO3)+ 0.212 * as.numeric(TPQ_Q$QO4)+ -0.11 * as.numeric(TPQ_Q$QO5)+ 0.215 * as.numeric(TPQ_Q$QO6)+ 0.1 * as.numeric(TPQ_Q$QO7)+ 0.118 * as.numeric(TPQ_Q$QO8)+ 0.15 * as.numeric(TPQ_Q$QO9)+ -0.171 * as.numeric(TPQ_Q$QO10)+ -0.195 * as.numeric(TPQ_Q$QO11)+ -0.055 * as.numeric(TPQ_Q$QO12)+ -0.105 * as.numeric(TPQ_Q$QO13)+ -0.186 * as.numeric(TPQ_Q$QO14)+ -0.049 * as.numeric(TPQ_Q$QO15)+ 0.075 * as.numeric(TPQ_Q$QO16)+ -0.035 * as.numeric(TPQ_Q$QO17)+ 0.07 * as.numeric(TPQ_Q$QO18)+ 0.04 * as.numeric(TPQ_Q$QO19)+ -0.077 * as.numeric(TPQ_Q$QO20)+ 0.1 * as.numeric(TPQ_Q$QO21)+ 0.168 * as.numeric(TPQ_Q$QO22)+ 0.011 * as.numeric(TPQ_Q$QO23)+ -0.113 * as.numeric(TPQ_Q$QO24)+ 0.003 * as.numeric(TPQ_Q$QO25)+ -0.095 * as.numeric(TPQ_Q$QO26)+ 0.12 * as.numeric(TPQ_Q$QO27)+ -0.352 * as.numeric(TPQ_Q$QO28)+ -0.128 * as.numeric(TPQ_Q$QO29)+ -0.087 * as.numeric(TPQ_Q$QO30)+ 0.157 * as.numeric(TPQ_Q$QO31)+ 0.296 * as.numeric(TPQ_Q$QO32)+ 0.05 * as.numeric(TPQ_Q$QO33)+ 0.129 * as.numeric(TPQ_Q$QO34)+ -0.015 * as.numeric(TPQ_Q$QO35)+ -0.106 * as.numeric(TPQ_Q$QO36)+ -0.173 * as.numeric(TPQ_Q$QO37)+ -0.261 * as.numeric(TPQ_Q$QO38)+ 0.013 * as.numeric(TPQ_Q$QO39)+ -0.267 * as.numeric(TPQ_Q$QO40)+ -0.145 * as.numeric(TPQ_Q$QO41)+ 0.288 * as.numeric(TPQ_Q$QO42)+ 0.024 * as.numeric(TPQ_Q$QO43)+ 0.281 * as.numeric(TPQ_Q$QO44)+ -0.099 * as.numeric(TPQ_Q$QO45)+ 0.084 * as.numeric(TPQ_Q$QO46)+ -0.076 * as.numeric(TPQ_Q$QO47)+ -0.055 * as.numeric(TPQ_Q$QO48)+ -0.266 * as.numeric(TPQ_Q$QO49)+ 0.019 * as.numeric(TPQ_Q$QO50)+ 0.038 * as.numeric(TPQ_Q$QO51)+ 0.23 * as.numeric(TPQ_Q$QO52)+ -0.126 * as.numeric(TPQ_Q$QO53)+ -0.111 * as.numeric(TPQ_Q$QO54)+ 0.089 * as.numeric(TPQ_Q$QO55)+ 0.1 * as.numeric(TPQ_Q$QO56)+ -0.006 * as.numeric(TPQ_Q$QO57)+ 0.089 * as.numeric(TPQ_Q$QO58)+ 0.016 * as.numeric(TPQ_Q$QO59)+ 0.002 * as.numeric(TPQ_Q$QO60)+ 0.063 * as.numeric(TPQ_Q$QO61)+ 0.068 * as.numeric(TPQ_Q$QO62)+ 0.159 * as.numeric(TPQ_Q$QO63)+ 0.096 * as.numeric(TPQ_Q$QO64)+ -0.092 * as.numeric(TPQ_Q$QO65)+ 0.06 * as.numeric(TPQ_Q$QO66)+ 0.167 * as.numeric(TPQ_Q$QO67)+ -0.305 * as.numeric(TPQ_Q$QO68)+ -0.275 * as.numeric(TPQ_Q$QO69)+ -0.161 * as.numeric(TPQ_Q$QO70)+ 0.008 * as.numeric(TPQ_Q$QO71)+ -0.115 * as.numeric(TPQ_Q$QO72)+ -0.306 * as.numeric(TPQ_Q$QO73)+ 0.147 * as.numeric(TPQ_Q$QO74)+ 0.318 * as.numeric(TPQ_Q$QO75)+ -0.104 * as.numeric(TPQ_Q$QO76)+ 0.025 * as.numeric(TPQ_Q$QO77)+ -0.249 * as.numeric(TPQ_Q$QO78)+ 0.312 * as.numeric(TPQ_Q$QO79)+ 0.331 * as.numeric(TPQ_Q$QO80)+ 0.152 * as.numeric(TPQ_Q$QO81)+ 0.286 * as.numeric(TPQ_Q$QO82)+ 0.123 * as.numeric(TPQ_Q$QO83)+ 0.113 * as.numeric(TPQ_Q$QO84)+ -0.127 * as.numeric(TPQ_Q$QO85)+ 0.029 * as.numeric(TPQ_Q$QO86)+ -0.194 * as.numeric(TPQ_Q$QO87)+ -0.028 * as.numeric(TPQ_Q$QO88)+ 0.188 * as.numeric(TPQ_Q$QO89)+ -0.225 * as.numeric(TPQ_Q$QO90)+ 0.038 * as.numeric(TPQ_Q$QO91)+ 0.062 * as.numeric(TPQ_Q$QO92)+ 0.105 * as.numeric(TPQ_Q$QO93)+ 0.142 * as.numeric(TPQ_Q$QO94)+ 0.249 * as.numeric(TPQ_Q$QO95)+ -0.252 * as.numeric(TPQ_Q$QO96)+ 0.286 * as.numeric(TPQ_Q$QO97)+ 0.285 * as.numeric(TPQ_Q$QO98)+ 0.097 * as.numeric(TPQ_Q$QO99)+ 0.149 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_1_All =  0.077 * as.numeric(TPQ_Q$QO1)+ -0.046 * as.numeric(TPQ_Q$QO2)+ -0.038 * as.numeric(TPQ_Q$QO3)+ -0.013 * as.numeric(TPQ_Q$QO4)+ -0.035 * as.numeric(TPQ_Q$QO5)+ 0.024 * as.numeric(TPQ_Q$QO6)+ -0.041 * as.numeric(TPQ_Q$QO7)+ 0.086 * as.numeric(TPQ_Q$QO8)+ 0.045 * as.numeric(TPQ_Q$QO9)+ -0.045 * as.numeric(TPQ_Q$QO10)+ 0.095 * as.numeric(TPQ_Q$QO11)+ 0.054 * as.numeric(TPQ_Q$QO12)+ -0.067 * as.numeric(TPQ_Q$QO13)+ -0.081 * as.numeric(TPQ_Q$QO14)+ 0.144 * as.numeric(TPQ_Q$QO15)+ 0.144 * as.numeric(TPQ_Q$QO16)+ 0.028 * as.numeric(TPQ_Q$QO17)+ -0.036 * as.numeric(TPQ_Q$QO18)+ -0.01 * as.numeric(TPQ_Q$QO19)+ -0.025 * as.numeric(TPQ_Q$QO20)+ 0.072 * as.numeric(TPQ_Q$QO21)+ -0.014 * as.numeric(TPQ_Q$QO22)+ 0.001 * as.numeric(TPQ_Q$QO23)+ -0.127 * as.numeric(TPQ_Q$QO24)+ 0.031 * as.numeric(TPQ_Q$QO25)+ 0.063 * as.numeric(TPQ_Q$QO26)+ -0.003 * as.numeric(TPQ_Q$QO27)+ -0.006 * as.numeric(TPQ_Q$QO28)+ 0.056 * as.numeric(TPQ_Q$QO29)+ -0.013 * as.numeric(TPQ_Q$QO30)+ -0.02 * as.numeric(TPQ_Q$QO31)+ 0.151 * as.numeric(TPQ_Q$QO32)+ 0.074 * as.numeric(TPQ_Q$QO33)+ -0.02 * as.numeric(TPQ_Q$QO34)+ -0.076 * as.numeric(TPQ_Q$QO35)+ 0.022 * as.numeric(TPQ_Q$QO36)+ 0.071 * as.numeric(TPQ_Q$QO37)+ -0.007 * as.numeric(TPQ_Q$QO38)+ 0.029 * as.numeric(TPQ_Q$QO39)+ 0.003 * as.numeric(TPQ_Q$QO40)+ -0.019 * as.numeric(TPQ_Q$QO41)+ 0.077 * as.numeric(TPQ_Q$QO42)+ -0.097 * as.numeric(TPQ_Q$QO43)+ 0.033 * as.numeric(TPQ_Q$QO44)+ 0.094 * as.numeric(TPQ_Q$QO45)+ 0.048 * as.numeric(TPQ_Q$QO46)+ 0.028 * as.numeric(TPQ_Q$QO47)+ -0.041 * as.numeric(TPQ_Q$QO48)+ -0.009 * as.numeric(TPQ_Q$QO49)+ -0.082 * as.numeric(TPQ_Q$QO50)+ -0.004 * as.numeric(TPQ_Q$QO51)+ 0.076 * as.numeric(TPQ_Q$QO52)+ 0.042 * as.numeric(TPQ_Q$QO53)+ -0.016 * as.numeric(TPQ_Q$QO54)+ 0.038 * as.numeric(TPQ_Q$QO55)+ 0.063 * as.numeric(TPQ_Q$QO56)+ -0.038 * as.numeric(TPQ_Q$QO57)+ 0.011 * as.numeric(TPQ_Q$QO58)+ -0.006 * as.numeric(TPQ_Q$QO59)+ 0.058 * as.numeric(TPQ_Q$QO60)+ -0.04 * as.numeric(TPQ_Q$QO61)+ 0.027 * as.numeric(TPQ_Q$QO62)+ 0.096 * as.numeric(TPQ_Q$QO63)+ -0.077 * as.numeric(TPQ_Q$QO64)+ 0.014 * as.numeric(TPQ_Q$QO65)+ 0.802 * as.numeric(TPQ_Q$QO66)+ -0.007 * as.numeric(TPQ_Q$QO67)+ -0.024 * as.numeric(TPQ_Q$QO68)+ -0.102 * as.numeric(TPQ_Q$QO69)+ -0.463 * as.numeric(TPQ_Q$QO70)+ -0.058 * as.numeric(TPQ_Q$QO71)+ -0.653 * as.numeric(TPQ_Q$QO72)+ -0.09 * as.numeric(TPQ_Q$QO73)+ -0.003 * as.numeric(TPQ_Q$QO74)+ 0.102 * as.numeric(TPQ_Q$QO75)+ 0.351 * as.numeric(TPQ_Q$QO76)+ 0.009 * as.numeric(TPQ_Q$QO77)+ 0.101 * as.numeric(TPQ_Q$QO78)+ -0.095 * as.numeric(TPQ_Q$QO79)+ 0.111 * as.numeric(TPQ_Q$QO80)+ 0.006 * as.numeric(TPQ_Q$QO81)+ 0.137 * as.numeric(TPQ_Q$QO82)+ -0.03 * as.numeric(TPQ_Q$QO83)+ -0.02 * as.numeric(TPQ_Q$QO84)+ 0.092 * as.numeric(TPQ_Q$QO85)+ 0.13 * as.numeric(TPQ_Q$QO86)+ 0.336 * as.numeric(TPQ_Q$QO87)+ -0.002 * as.numeric(TPQ_Q$QO88)+ 0.007 * as.numeric(TPQ_Q$QO89)+ 0.111 * as.numeric(TPQ_Q$QO90)+ 0.097 * as.numeric(TPQ_Q$QO91)+ 0.192 * as.numeric(TPQ_Q$QO92)+ 0.073 * as.numeric(TPQ_Q$QO93)+ 0.013 * as.numeric(TPQ_Q$QO94)+ 0.106 * as.numeric(TPQ_Q$QO95)+ 0.037 * as.numeric(TPQ_Q$QO96)+ 0.07 * as.numeric(TPQ_Q$QO97)+ 0.066 * as.numeric(TPQ_Q$QO98)+ 0.012 * as.numeric(TPQ_Q$QO99)+ 0.013 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_2_All =  0.104 * as.numeric(TPQ_Q$QO1)+ -0.066 * as.numeric(TPQ_Q$QO2)+ 0.023 * as.numeric(TPQ_Q$QO3)+ -0.078 * as.numeric(TPQ_Q$QO4)+ -0.036 * as.numeric(TPQ_Q$QO5)+ 0.017 * as.numeric(TPQ_Q$QO6)+ -0.093 * as.numeric(TPQ_Q$QO7)+ 0.285 * as.numeric(TPQ_Q$QO8)+ 0.159 * as.numeric(TPQ_Q$QO9)+ 0.158 * as.numeric(TPQ_Q$QO10)+ -0.077 * as.numeric(TPQ_Q$QO11)+ 0.118 * as.numeric(TPQ_Q$QO12)+ -0.108 * as.numeric(TPQ_Q$QO13)+ 0.18 * as.numeric(TPQ_Q$QO14)+ 0.21 * as.numeric(TPQ_Q$QO15)+ 0.092 * as.numeric(TPQ_Q$QO16)+ -0.132 * as.numeric(TPQ_Q$QO17)+ 0.056 * as.numeric(TPQ_Q$QO18)+ 0.102 * as.numeric(TPQ_Q$QO19)+ 0.056 * as.numeric(TPQ_Q$QO20)+ 0.12 * as.numeric(TPQ_Q$QO21)+ -0.135 * as.numeric(TPQ_Q$QO22)+ 0.151 * as.numeric(TPQ_Q$QO23)+ -0.012 * as.numeric(TPQ_Q$QO24)+ 0.017 * as.numeric(TPQ_Q$QO25)+ 0.256 * as.numeric(TPQ_Q$QO26)+ -0.032 * as.numeric(TPQ_Q$QO27)+ -0.236 * as.numeric(TPQ_Q$QO28)+ 0.024 * as.numeric(TPQ_Q$QO29)+ -0.213 * as.numeric(TPQ_Q$QO30)+ -0.039 * as.numeric(TPQ_Q$QO31)+ 0.273 * as.numeric(TPQ_Q$QO32)+ 0.21 * as.numeric(TPQ_Q$QO33)+ -0.15 * as.numeric(TPQ_Q$QO34)+ 0.423 * as.numeric(TPQ_Q$QO35)+ 0.256 * as.numeric(TPQ_Q$QO36)+ 0.18 * as.numeric(TPQ_Q$QO37)+ 0.298 * as.numeric(TPQ_Q$QO38)+ 0.042 * as.numeric(TPQ_Q$QO39)+ 0.286 * as.numeric(TPQ_Q$QO40)+ -0.006 * as.numeric(TPQ_Q$QO41)+ 0.197 * as.numeric(TPQ_Q$QO42)+ -0.056 * as.numeric(TPQ_Q$QO43)+ 0.159 * as.numeric(TPQ_Q$QO44)+ 0.096 * as.numeric(TPQ_Q$QO45)+ 0.193 * as.numeric(TPQ_Q$QO46)+ -0.036 * as.numeric(TPQ_Q$QO47)+ -0.138 * as.numeric(TPQ_Q$QO48)+ 0.113 * as.numeric(TPQ_Q$QO49)+ -0.025 * as.numeric(TPQ_Q$QO50)+ -0.051 * as.numeric(TPQ_Q$QO51)+ -0.12 * as.numeric(TPQ_Q$QO52)+ 0.103 * as.numeric(TPQ_Q$QO53)+ 0.0 * as.numeric(TPQ_Q$QO54)+ 0.216 * as.numeric(TPQ_Q$QO55)+ 0.192 * as.numeric(TPQ_Q$QO56)+ -0.012 * as.numeric(TPQ_Q$QO57)+ 0.016 * as.numeric(TPQ_Q$QO58)+ 0.05 * as.numeric(TPQ_Q$QO59)+ 0.095 * as.numeric(TPQ_Q$QO60)+ -0.085 * as.numeric(TPQ_Q$QO61)+ 0.028 * as.numeric(TPQ_Q$QO62)+ 0.151 * as.numeric(TPQ_Q$QO63)+ -0.078 * as.numeric(TPQ_Q$QO64)+ 0.12 * as.numeric(TPQ_Q$QO65)+ -0.014 * as.numeric(TPQ_Q$QO66)+ 0.064 * as.numeric(TPQ_Q$QO67)+ -0.084 * as.numeric(TPQ_Q$QO68)+ -0.036 * as.numeric(TPQ_Q$QO69)+ 0.036 * as.numeric(TPQ_Q$QO70)+ -0.0 * as.numeric(TPQ_Q$QO71)+ 0.049 * as.numeric(TPQ_Q$QO72)+ 0.085 * as.numeric(TPQ_Q$QO73)+ -0.008 * as.numeric(TPQ_Q$QO74)+ 0.038 * as.numeric(TPQ_Q$QO75)+ 0.029 * as.numeric(TPQ_Q$QO76)+ 0.055 * as.numeric(TPQ_Q$QO77)+ 0.183 * as.numeric(TPQ_Q$QO78)+ 0.077 * as.numeric(TPQ_Q$QO79)+ 0.126 * as.numeric(TPQ_Q$QO80)+ 0.091 * as.numeric(TPQ_Q$QO81)+ -0.018 * as.numeric(TPQ_Q$QO82)+ -0.036 * as.numeric(TPQ_Q$QO83)+ 0.057 * as.numeric(TPQ_Q$QO84)+ 0.156 * as.numeric(TPQ_Q$QO85)+ 0.145 * as.numeric(TPQ_Q$QO86)+ 0.142 * as.numeric(TPQ_Q$QO87)+ 0.042 * as.numeric(TPQ_Q$QO88)+ -0.166 * as.numeric(TPQ_Q$QO89)+ 0.269 * as.numeric(TPQ_Q$QO90)+ 0.123 * as.numeric(TPQ_Q$QO91)+ -0.063 * as.numeric(TPQ_Q$QO92)+ 0.227 * as.numeric(TPQ_Q$QO93)+ -0.033 * as.numeric(TPQ_Q$QO94)+ -0.048 * as.numeric(TPQ_Q$QO95)+ 0.167 * as.numeric(TPQ_Q$QO96)+ 0.009 * as.numeric(TPQ_Q$QO97)+ 0.159 * as.numeric(TPQ_Q$QO98)+ 0.026 * as.numeric(TPQ_Q$QO99)+ 0.064 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_3_All =  0.131 * as.numeric(TPQ_Q$QO1)+ -0.014 * as.numeric(TPQ_Q$QO2)+ 0.126 * as.numeric(TPQ_Q$QO3)+ -0.01 * as.numeric(TPQ_Q$QO4)+ 0.123 * as.numeric(TPQ_Q$QO5)+ -0.142 * as.numeric(TPQ_Q$QO6)+ 0.098 * as.numeric(TPQ_Q$QO7)+ 0.019 * as.numeric(TPQ_Q$QO8)+ -0.026 * as.numeric(TPQ_Q$QO9)+ 0.118 * as.numeric(TPQ_Q$QO10)+ 0.321 * as.numeric(TPQ_Q$QO11)+ 0.063 * as.numeric(TPQ_Q$QO12)+ 0.118 * as.numeric(TPQ_Q$QO13)+ 0.0 * as.numeric(TPQ_Q$QO14)+ 0.028 * as.numeric(TPQ_Q$QO15)+ 0.198 * as.numeric(TPQ_Q$QO16)+ 0.319 * as.numeric(TPQ_Q$QO17)+ 0.042 * as.numeric(TPQ_Q$QO18)+ 0.005 * as.numeric(TPQ_Q$QO19)+ 0.123 * as.numeric(TPQ_Q$QO20)+ 0.271 * as.numeric(TPQ_Q$QO21)+ 0.116 * as.numeric(TPQ_Q$QO22)+ -0.003 * as.numeric(TPQ_Q$QO23)+ -0.138 * as.numeric(TPQ_Q$QO24)+ 0.216 * as.numeric(TPQ_Q$QO25)+ -0.0 * as.numeric(TPQ_Q$QO26)+ 0.013 * as.numeric(TPQ_Q$QO27)+ 0.181 * as.numeric(TPQ_Q$QO28)+ 0.039 * as.numeric(TPQ_Q$QO29)+ -0.017 * as.numeric(TPQ_Q$QO30)+ 0.043 * as.numeric(TPQ_Q$QO31)+ 0.074 * as.numeric(TPQ_Q$QO32)+ -0.002 * as.numeric(TPQ_Q$QO33)+ 0.104 * as.numeric(TPQ_Q$QO34)+ -0.032 * as.numeric(TPQ_Q$QO35)+ 0.191 * as.numeric(TPQ_Q$QO36)+ 0.174 * as.numeric(TPQ_Q$QO37)+ 0.019 * as.numeric(TPQ_Q$QO38)+ 0.012 * as.numeric(TPQ_Q$QO39)+ 0.101 * as.numeric(TPQ_Q$QO40)+ 0.053 * as.numeric(TPQ_Q$QO41)+ -0.031 * as.numeric(TPQ_Q$QO42)+ 0.096 * as.numeric(TPQ_Q$QO43)+ 0.045 * as.numeric(TPQ_Q$QO44)+ 0.215 * as.numeric(TPQ_Q$QO45)+ 0.266 * as.numeric(TPQ_Q$QO46)+ -0.032 * as.numeric(TPQ_Q$QO47)+ 0.06 * as.numeric(TPQ_Q$QO48)+ 0.176 * as.numeric(TPQ_Q$QO49)+ 0.038 * as.numeric(TPQ_Q$QO50)+ -0.084 * as.numeric(TPQ_Q$QO51)+ 0.248 * as.numeric(TPQ_Q$QO52)+ 0.106 * as.numeric(TPQ_Q$QO53)+ 0.026 * as.numeric(TPQ_Q$QO54)+ 0.236 * as.numeric(TPQ_Q$QO55)+ 0.264 * as.numeric(TPQ_Q$QO56)+ 0.198 * as.numeric(TPQ_Q$QO57)+ 0.175 * as.numeric(TPQ_Q$QO58)+ 0.036 * as.numeric(TPQ_Q$QO59)+ -0.138 * as.numeric(TPQ_Q$QO60)+ 0.009 * as.numeric(TPQ_Q$QO61)+ -0.06 * as.numeric(TPQ_Q$QO62)+ 0.07 * as.numeric(TPQ_Q$QO63)+ 0.143 * as.numeric(TPQ_Q$QO64)+ 0.471 * as.numeric(TPQ_Q$QO65)+ 0.025 * as.numeric(TPQ_Q$QO66)+ 0.056 * as.numeric(TPQ_Q$QO67)+ 0.07 * as.numeric(TPQ_Q$QO68)+ 0.131 * as.numeric(TPQ_Q$QO69)+ 0.021 * as.numeric(TPQ_Q$QO70)+ 0.219 * as.numeric(TPQ_Q$QO71)+ 0.075 * as.numeric(TPQ_Q$QO72)+ 0.287 * as.numeric(TPQ_Q$QO73)+ 0.031 * as.numeric(TPQ_Q$QO74)+ 0.203 * as.numeric(TPQ_Q$QO75)+ -0.08 * as.numeric(TPQ_Q$QO76)+ 0.135 * as.numeric(TPQ_Q$QO77)+ 0.102 * as.numeric(TPQ_Q$QO78)+ 0.11 * as.numeric(TPQ_Q$QO79)+ 0.058 * as.numeric(TPQ_Q$QO80)+ 0.007 * as.numeric(TPQ_Q$QO81)+ 0.09 * as.numeric(TPQ_Q$QO82)+ 0.06 * as.numeric(TPQ_Q$QO83)+ 0.254 * as.numeric(TPQ_Q$QO84)+ 0.249 * as.numeric(TPQ_Q$QO85)+ 0.018 * as.numeric(TPQ_Q$QO86)+ 0.061 * as.numeric(TPQ_Q$QO87)+ -0.038 * as.numeric(TPQ_Q$QO88)+ 0.211 * as.numeric(TPQ_Q$QO89)+ -0.023 * as.numeric(TPQ_Q$QO90)+ 0.121 * as.numeric(TPQ_Q$QO91)+ 0.151 * as.numeric(TPQ_Q$QO92)+ 0.12 * as.numeric(TPQ_Q$QO93)+ 0.031 * as.numeric(TPQ_Q$QO94)+ 0.093 * as.numeric(TPQ_Q$QO95)+ 0.115 * as.numeric(TPQ_Q$QO96)+ 0.111 * as.numeric(TPQ_Q$QO97)+ 0.107 * as.numeric(TPQ_Q$QO98)+ 0.14 * as.numeric(TPQ_Q$QO99)+ 0.198 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_4_All =  0.014 * as.numeric(TPQ_Q$QO1)+ -0.124 * as.numeric(TPQ_Q$QO2)+ -0.011 * as.numeric(TPQ_Q$QO3)+ -0.209 * as.numeric(TPQ_Q$QO4)+ -0.08 * as.numeric(TPQ_Q$QO5)+ 0.135 * as.numeric(TPQ_Q$QO6)+ -0.012 * as.numeric(TPQ_Q$QO7)+ 0.267 * as.numeric(TPQ_Q$QO8)+ -0.137 * as.numeric(TPQ_Q$QO9)+ -0.099 * as.numeric(TPQ_Q$QO10)+ 0.103 * as.numeric(TPQ_Q$QO11)+ -0.025 * as.numeric(TPQ_Q$QO12)+ 0.105 * as.numeric(TPQ_Q$QO13)+ -0.133 * as.numeric(TPQ_Q$QO14)+ 0.107 * as.numeric(TPQ_Q$QO15)+ -0.125 * as.numeric(TPQ_Q$QO16)+ -0.013 * as.numeric(TPQ_Q$QO17)+ 0.057 * as.numeric(TPQ_Q$QO18)+ -0.024 * as.numeric(TPQ_Q$QO19)+ -0.092 * as.numeric(TPQ_Q$QO20)+ -0.154 * as.numeric(TPQ_Q$QO21)+ 0.219 * as.numeric(TPQ_Q$QO22)+ -0.065 * as.numeric(TPQ_Q$QO23)+ 0.112 * as.numeric(TPQ_Q$QO24)+ -0.056 * as.numeric(TPQ_Q$QO25)+ 0.137 * as.numeric(TPQ_Q$QO26)+ 0.004 * as.numeric(TPQ_Q$QO27)+ 0.222 * as.numeric(TPQ_Q$QO28)+ 0.092 * as.numeric(TPQ_Q$QO29)+ 0.208 * as.numeric(TPQ_Q$QO30)+ -0.012 * as.numeric(TPQ_Q$QO31)+ -0.152 * as.numeric(TPQ_Q$QO32)+ 0.128 * as.numeric(TPQ_Q$QO33)+ 0.037 * as.numeric(TPQ_Q$QO34)+ 0.134 * as.numeric(TPQ_Q$QO35)+ -0.19 * as.numeric(TPQ_Q$QO36)+ -0.018 * as.numeric(TPQ_Q$QO37)+ 0.041 * as.numeric(TPQ_Q$QO38)+ -0.06 * as.numeric(TPQ_Q$QO39)+ 0.103 * as.numeric(TPQ_Q$QO40)+ -0.057 * as.numeric(TPQ_Q$QO41)+ -0.115 * as.numeric(TPQ_Q$QO42)+ 0.12 * as.numeric(TPQ_Q$QO43)+ -0.084 * as.numeric(TPQ_Q$QO44)+ 0.08 * as.numeric(TPQ_Q$QO45)+ -0.221 * as.numeric(TPQ_Q$QO46)+ 0.05 * as.numeric(TPQ_Q$QO47)+ 0.449 * as.numeric(TPQ_Q$QO48)+ 0.126 * as.numeric(TPQ_Q$QO49)+ 0.242 * as.numeric(TPQ_Q$QO50)+ -0.063 * as.numeric(TPQ_Q$QO51)+ -0.037 * as.numeric(TPQ_Q$QO52)+ 0.08 * as.numeric(TPQ_Q$QO53)+ 0.022 * as.numeric(TPQ_Q$QO54)+ -0.163 * as.numeric(TPQ_Q$QO55)+ -0.167 * as.numeric(TPQ_Q$QO56)+ 0.042 * as.numeric(TPQ_Q$QO57)+ -0.028 * as.numeric(TPQ_Q$QO58)+ -0.012 * as.numeric(TPQ_Q$QO59)+ -0.006 * as.numeric(TPQ_Q$QO60)+ 0.038 * as.numeric(TPQ_Q$QO61)+ 0.025 * as.numeric(TPQ_Q$QO62)+ -0.092 * as.numeric(TPQ_Q$QO63)+ -0.002 * as.numeric(TPQ_Q$QO64)+ -0.067 * as.numeric(TPQ_Q$QO65)+ -0.072 * as.numeric(TPQ_Q$QO66)+ -0.048 * as.numeric(TPQ_Q$QO67)+ -0.188 * as.numeric(TPQ_Q$QO68)+ -0.133 * as.numeric(TPQ_Q$QO69)+ 0.133 * as.numeric(TPQ_Q$QO70)+ -0.157 * as.numeric(TPQ_Q$QO71)+ 0.094 * as.numeric(TPQ_Q$QO72)+ 0.016 * as.numeric(TPQ_Q$QO73)+ 0.012 * as.numeric(TPQ_Q$QO74)+ 0.119 * as.numeric(TPQ_Q$QO75)+ -0.075 * as.numeric(TPQ_Q$QO76)+ -0.254 * as.numeric(TPQ_Q$QO77)+ 0.011 * as.numeric(TPQ_Q$QO78)+ 0.088 * as.numeric(TPQ_Q$QO79)+ 0.24 * as.numeric(TPQ_Q$QO80)+ -0.145 * as.numeric(TPQ_Q$QO81)+ 0.019 * as.numeric(TPQ_Q$QO82)+ -0.01 * as.numeric(TPQ_Q$QO83)+ -0.092 * as.numeric(TPQ_Q$QO84)+ 0.189 * as.numeric(TPQ_Q$QO85)+ -0.067 * as.numeric(TPQ_Q$QO86)+ -0.14 * as.numeric(TPQ_Q$QO87)+ -0.017 * as.numeric(TPQ_Q$QO88)+ -0.142 * as.numeric(TPQ_Q$QO89)+ 0.116 * as.numeric(TPQ_Q$QO90)+ 0.107 * as.numeric(TPQ_Q$QO91)+ -0.056 * as.numeric(TPQ_Q$QO92)+ 0.105 * as.numeric(TPQ_Q$QO93)+ -0.044 * as.numeric(TPQ_Q$QO94)+ -0.001 * as.numeric(TPQ_Q$QO95)+ 0.057 * as.numeric(TPQ_Q$QO96)+ -0.139 * as.numeric(TPQ_Q$QO97)+ 0.185 * as.numeric(TPQ_Q$QO98)+ -0.172 * as.numeric(TPQ_Q$QO99)+ -0.135 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_5_All =  0.063 * as.numeric(TPQ_Q$QO1)+ 0.102 * as.numeric(TPQ_Q$QO2)+ 0.08 * as.numeric(TPQ_Q$QO3)+ 0.209 * as.numeric(TPQ_Q$QO4)+ 0.142 * as.numeric(TPQ_Q$QO5)+ 0.179 * as.numeric(TPQ_Q$QO6)+ 0.077 * as.numeric(TPQ_Q$QO7)+ 0.047 * as.numeric(TPQ_Q$QO8)+ 0.061 * as.numeric(TPQ_Q$QO9)+ 0.185 * as.numeric(TPQ_Q$QO10)+ 0.119 * as.numeric(TPQ_Q$QO11)+ 0.007 * as.numeric(TPQ_Q$QO12)+ 0.107 * as.numeric(TPQ_Q$QO13)+ 0.277 * as.numeric(TPQ_Q$QO14)+ -0.065 * as.numeric(TPQ_Q$QO15)+ 0.051 * as.numeric(TPQ_Q$QO16)+ 0.049 * as.numeric(TPQ_Q$QO17)+ 0.023 * as.numeric(TPQ_Q$QO18)+ -0.019 * as.numeric(TPQ_Q$QO19)+ 0.17 * as.numeric(TPQ_Q$QO20)+ -0.052 * as.numeric(TPQ_Q$QO21)+ 0.071 * as.numeric(TPQ_Q$QO22)+ 0.054 * as.numeric(TPQ_Q$QO23)+ 0.299 * as.numeric(TPQ_Q$QO24)+ 0.038 * as.numeric(TPQ_Q$QO25)+ 0.04 * as.numeric(TPQ_Q$QO26)+ -0.023 * as.numeric(TPQ_Q$QO27)+ 0.053 * as.numeric(TPQ_Q$QO28)+ 0.143 * as.numeric(TPQ_Q$QO29)+ 0.216 * as.numeric(TPQ_Q$QO30)+ 0.1 * as.numeric(TPQ_Q$QO31)+ -0.019 * as.numeric(TPQ_Q$QO32)+ -0.215 * as.numeric(TPQ_Q$QO33)+ 0.053 * as.numeric(TPQ_Q$QO34)+ -0.027 * as.numeric(TPQ_Q$QO35)+ 0.099 * as.numeric(TPQ_Q$QO36)+ -0.05 * as.numeric(TPQ_Q$QO37)+ -0.092 * as.numeric(TPQ_Q$QO38)+ 0.088 * as.numeric(TPQ_Q$QO39)+ 0.083 * as.numeric(TPQ_Q$QO40)+ 0.08 * as.numeric(TPQ_Q$QO41)+ 0.01 * as.numeric(TPQ_Q$QO42)+ 0.041 * as.numeric(TPQ_Q$QO43)+ 0.095 * as.numeric(TPQ_Q$QO44)+ 0.147 * as.numeric(TPQ_Q$QO45)+ -0.059 * as.numeric(TPQ_Q$QO46)+ 0.066 * as.numeric(TPQ_Q$QO47)+ 0.114 * as.numeric(TPQ_Q$QO48)+ 0.154 * as.numeric(TPQ_Q$QO49)+ 0.118 * as.numeric(TPQ_Q$QO50)+ 0.004 * as.numeric(TPQ_Q$QO51)+ 0.109 * as.numeric(TPQ_Q$QO52)+ 0.216 * as.numeric(TPQ_Q$QO53)+ 0.031 * as.numeric(TPQ_Q$QO54)+ -0.029 * as.numeric(TPQ_Q$QO55)+ -0.026 * as.numeric(TPQ_Q$QO56)+ 0.007 * as.numeric(TPQ_Q$QO57)+ -0.014 * as.numeric(TPQ_Q$QO58)+ 0.035 * as.numeric(TPQ_Q$QO59)+ 0.613 * as.numeric(TPQ_Q$QO60)+ 0.079 * as.numeric(TPQ_Q$QO61)+ 0.527 * as.numeric(TPQ_Q$QO62)+ 0.002 * as.numeric(TPQ_Q$QO63)+ 0.077 * as.numeric(TPQ_Q$QO64)+ -0.389 * as.numeric(TPQ_Q$QO65)+ 0.029 * as.numeric(TPQ_Q$QO66)+ 0.048 * as.numeric(TPQ_Q$QO67)+ -0.021 * as.numeric(TPQ_Q$QO68)+ 0.04 * as.numeric(TPQ_Q$QO69)+ 0.125 * as.numeric(TPQ_Q$QO70)+ 0.108 * as.numeric(TPQ_Q$QO71)+ 0.137 * as.numeric(TPQ_Q$QO72)+ 0.093 * as.numeric(TPQ_Q$QO73)+ 0.029 * as.numeric(TPQ_Q$QO74)+ 0.198 * as.numeric(TPQ_Q$QO75)+ 0.09 * as.numeric(TPQ_Q$QO76)+ 0.057 * as.numeric(TPQ_Q$QO77)+ 0.066 * as.numeric(TPQ_Q$QO78)+ 0.073 * as.numeric(TPQ_Q$QO79)+ 0.161 * as.numeric(TPQ_Q$QO80)+ -0.025 * as.numeric(TPQ_Q$QO81)+ 0.006 * as.numeric(TPQ_Q$QO82)+ 0.053 * as.numeric(TPQ_Q$QO83)+ 0.231 * as.numeric(TPQ_Q$QO84)+ 0.096 * as.numeric(TPQ_Q$QO85)+ -0.091 * as.numeric(TPQ_Q$QO86)+ 0.048 * as.numeric(TPQ_Q$QO87)+ 0.024 * as.numeric(TPQ_Q$QO88)+ 0.208 * as.numeric(TPQ_Q$QO89)+ -0.09 * as.numeric(TPQ_Q$QO90)+ -0.05 * as.numeric(TPQ_Q$QO91)+ 0.016 * as.numeric(TPQ_Q$QO92)+ -0.132 * as.numeric(TPQ_Q$QO93)+ 0.024 * as.numeric(TPQ_Q$QO94)+ 0.03 * as.numeric(TPQ_Q$QO95)+ -0.005 * as.numeric(TPQ_Q$QO96)+ 0.007 * as.numeric(TPQ_Q$QO97)+ 0.111 * as.numeric(TPQ_Q$QO98)+ 0.008 * as.numeric(TPQ_Q$QO99)+ 0.237 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_6_All =  -0.001 * as.numeric(TPQ_Q$QO1)+ -0.059 * as.numeric(TPQ_Q$QO2)+ 0.052 * as.numeric(TPQ_Q$QO3)+ 0.013 * as.numeric(TPQ_Q$QO4)+ 0.147 * as.numeric(TPQ_Q$QO5)+ 0.127 * as.numeric(TPQ_Q$QO6)+ -0.101 * as.numeric(TPQ_Q$QO7)+ 0.023 * as.numeric(TPQ_Q$QO8)+ -0.024 * as.numeric(TPQ_Q$QO9)+ 0.002 * as.numeric(TPQ_Q$QO10)+ -0.031 * as.numeric(TPQ_Q$QO11)+ 0.008 * as.numeric(TPQ_Q$QO12)+ -0.028 * as.numeric(TPQ_Q$QO13)+ 0.007 * as.numeric(TPQ_Q$QO14)+ 0.253 * as.numeric(TPQ_Q$QO15)+ 0.201 * as.numeric(TPQ_Q$QO16)+ -0.083 * as.numeric(TPQ_Q$QO17)+ 0.021 * as.numeric(TPQ_Q$QO18)+ 0.081 * as.numeric(TPQ_Q$QO19)+ 0.227 * as.numeric(TPQ_Q$QO20)+ 0.205 * as.numeric(TPQ_Q$QO21)+ -0.007 * as.numeric(TPQ_Q$QO22)+ 0.102 * as.numeric(TPQ_Q$QO23)+ -0.055 * as.numeric(TPQ_Q$QO24)+ -0.038 * as.numeric(TPQ_Q$QO25)+ 0.094 * as.numeric(TPQ_Q$QO26)+ -0.035 * as.numeric(TPQ_Q$QO27)+ 0.239 * as.numeric(TPQ_Q$QO28)+ 0.022 * as.numeric(TPQ_Q$QO29)+ 0.153 * as.numeric(TPQ_Q$QO30)+ 0.017 * as.numeric(TPQ_Q$QO31)+ -0.04 * as.numeric(TPQ_Q$QO32)+ -0.054 * as.numeric(TPQ_Q$QO33)+ 0.04 * as.numeric(TPQ_Q$QO34)+ 0.035 * as.numeric(TPQ_Q$QO35)+ 0.1 * as.numeric(TPQ_Q$QO36)+ -0.104 * as.numeric(TPQ_Q$QO37)+ -0.01 * as.numeric(TPQ_Q$QO38)+ 0.589 * as.numeric(TPQ_Q$QO39)+ -0.014 * as.numeric(TPQ_Q$QO40)+ 0.659 * as.numeric(TPQ_Q$QO41)+ -0.064 * as.numeric(TPQ_Q$QO42)+ -0.109 * as.numeric(TPQ_Q$QO43)+ -0.086 * as.numeric(TPQ_Q$QO44)+ -0.393 * as.numeric(TPQ_Q$QO45)+ 0.045 * as.numeric(TPQ_Q$QO46)+ -0.015 * as.numeric(TPQ_Q$QO47)+ 0.083 * as.numeric(TPQ_Q$QO48)+ -0.096 * as.numeric(TPQ_Q$QO49)+ 0.048 * as.numeric(TPQ_Q$QO50)+ -0.058 * as.numeric(TPQ_Q$QO51)+ -0.08 * as.numeric(TPQ_Q$QO52)+ -0.002 * as.numeric(TPQ_Q$QO53)+ 0.032 * as.numeric(TPQ_Q$QO54)+ 0.121 * as.numeric(TPQ_Q$QO55)+ 0.096 * as.numeric(TPQ_Q$QO56)+ -0.049 * as.numeric(TPQ_Q$QO57)+ -0.042 * as.numeric(TPQ_Q$QO58)+ 0.029 * as.numeric(TPQ_Q$QO59)+ -0.07 * as.numeric(TPQ_Q$QO60)+ 0.014 * as.numeric(TPQ_Q$QO61)+ -0.056 * as.numeric(TPQ_Q$QO62)+ 0.224 * as.numeric(TPQ_Q$QO63)+ 0.016 * as.numeric(TPQ_Q$QO64)+ 0.063 * as.numeric(TPQ_Q$QO65)+ 0.033 * as.numeric(TPQ_Q$QO66)+ -0.076 * as.numeric(TPQ_Q$QO67)+ 0.149 * as.numeric(TPQ_Q$QO68)+ 0.105 * as.numeric(TPQ_Q$QO69)+ 0.085 * as.numeric(TPQ_Q$QO70)+ -0.046 * as.numeric(TPQ_Q$QO71)+ 0.08 * as.numeric(TPQ_Q$QO72)+ 0.021 * as.numeric(TPQ_Q$QO73)+ -0.094 * as.numeric(TPQ_Q$QO74)+ 0.099 * as.numeric(TPQ_Q$QO75)+ -0.015 * as.numeric(TPQ_Q$QO76)+ 0.156 * as.numeric(TPQ_Q$QO77)+ 0.113 * as.numeric(TPQ_Q$QO78)+ -0.069 * as.numeric(TPQ_Q$QO79)+ -0.035 * as.numeric(TPQ_Q$QO80)+ -0.056 * as.numeric(TPQ_Q$QO81)+ -0.028 * as.numeric(TPQ_Q$QO82)+ 0.059 * as.numeric(TPQ_Q$QO83)+ -0.001 * as.numeric(TPQ_Q$QO84)+ -0.025 * as.numeric(TPQ_Q$QO85)+ -0.04 * as.numeric(TPQ_Q$QO86)+ -0.023 * as.numeric(TPQ_Q$QO87)+ 0.064 * as.numeric(TPQ_Q$QO88)+ 0.014 * as.numeric(TPQ_Q$QO89)+ 0.215 * as.numeric(TPQ_Q$QO90)+ -0.079 * as.numeric(TPQ_Q$QO91)+ 0.492 * as.numeric(TPQ_Q$QO92)+ -0.089 * as.numeric(TPQ_Q$QO93)+ 0.028 * as.numeric(TPQ_Q$QO94)+ -0.08 * as.numeric(TPQ_Q$QO95)+ 0.105 * as.numeric(TPQ_Q$QO96)+ 0.132 * as.numeric(TPQ_Q$QO97)+ 0.07 * as.numeric(TPQ_Q$QO98)+ 0.125 * as.numeric(TPQ_Q$QO99)+ -0.03 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_7_All =  -0.247 * as.numeric(TPQ_Q$QO1)+ -0.01 * as.numeric(TPQ_Q$QO2)+ -0.002 * as.numeric(TPQ_Q$QO3)+ -0.041 * as.numeric(TPQ_Q$QO4)+ 0.31 * as.numeric(TPQ_Q$QO5)+ -0.147 * as.numeric(TPQ_Q$QO6)+ 0.089 * as.numeric(TPQ_Q$QO7)+ -0.265 * as.numeric(TPQ_Q$QO8)+ -0.036 * as.numeric(TPQ_Q$QO9)+ 0.288 * as.numeric(TPQ_Q$QO10)+ 0.231 * as.numeric(TPQ_Q$QO11)+ 0.135 * as.numeric(TPQ_Q$QO12)+ 0.118 * as.numeric(TPQ_Q$QO13)+ 0.305 * as.numeric(TPQ_Q$QO14)+ -0.097 * as.numeric(TPQ_Q$QO15)+ 0.063 * as.numeric(TPQ_Q$QO16)+ 0.089 * as.numeric(TPQ_Q$QO17)+ 0.384 * as.numeric(TPQ_Q$QO18)+ 0.64 * as.numeric(TPQ_Q$QO19)+ -0.04 * as.numeric(TPQ_Q$QO20)+ 0.074 * as.numeric(TPQ_Q$QO21)+ 0.07 * as.numeric(TPQ_Q$QO22)+ 0.661 * as.numeric(TPQ_Q$QO23)+ 0.103 * as.numeric(TPQ_Q$QO24)+ 0.001 * as.numeric(TPQ_Q$QO25)+ -0.121 * as.numeric(TPQ_Q$QO26)+ 0.014 * as.numeric(TPQ_Q$QO27)+ 0.197 * as.numeric(TPQ_Q$QO28)+ -0.121 * as.numeric(TPQ_Q$QO29)+ 0.399 * as.numeric(TPQ_Q$QO30)+ -0.066 * as.numeric(TPQ_Q$QO31)+ -0.229 * as.numeric(TPQ_Q$QO32)+ 0.509 * as.numeric(TPQ_Q$QO33)+ 0.165 * as.numeric(TPQ_Q$QO34)+ -0.072 * as.numeric(TPQ_Q$QO35)+ -0.125 * as.numeric(TPQ_Q$QO36)+ 0.346 * as.numeric(TPQ_Q$QO37)+ 0.263 * as.numeric(TPQ_Q$QO38)+ 0.142 * as.numeric(TPQ_Q$QO39)+ 0.241 * as.numeric(TPQ_Q$QO40)+ 0.14 * as.numeric(TPQ_Q$QO41)+ -0.133 * as.numeric(TPQ_Q$QO42)+ 0.223 * as.numeric(TPQ_Q$QO43)+ -0.23 * as.numeric(TPQ_Q$QO44)+ 0.043 * as.numeric(TPQ_Q$QO45)+ 0.005 * as.numeric(TPQ_Q$QO46)+ 0.003 * as.numeric(TPQ_Q$QO47)+ 0.101 * as.numeric(TPQ_Q$QO48)+ 0.077 * as.numeric(TPQ_Q$QO49)+ 0.372 * as.numeric(TPQ_Q$QO50)+ 0.005 * as.numeric(TPQ_Q$QO51)+ -0.055 * as.numeric(TPQ_Q$QO52)+ 0.097 * as.numeric(TPQ_Q$QO53)+ 0.15 * as.numeric(TPQ_Q$QO54)+ -0.026 * as.numeric(TPQ_Q$QO55)+ -0.009 * as.numeric(TPQ_Q$QO56)+ 0.024 * as.numeric(TPQ_Q$QO57)+ -0.044 * as.numeric(TPQ_Q$QO58)+ -0.039 * as.numeric(TPQ_Q$QO59)+ -0.056 * as.numeric(TPQ_Q$QO60)+ 0.179 * as.numeric(TPQ_Q$QO61)+ -0.043 * as.numeric(TPQ_Q$QO62)+ -0.14 * as.numeric(TPQ_Q$QO63)+ 0.151 * as.numeric(TPQ_Q$QO64)+ 0.011 * as.numeric(TPQ_Q$QO65)+ -0.006 * as.numeric(TPQ_Q$QO66)+ 0.039 * as.numeric(TPQ_Q$QO67)+ 0.08 * as.numeric(TPQ_Q$QO68)+ 0.129 * as.numeric(TPQ_Q$QO69)+ 0.036 * as.numeric(TPQ_Q$QO70)+ 0.116 * as.numeric(TPQ_Q$QO71)+ 0.068 * as.numeric(TPQ_Q$QO72)+ 0.38 * as.numeric(TPQ_Q$QO73)+ 0.023 * as.numeric(TPQ_Q$QO74)+ -0.185 * as.numeric(TPQ_Q$QO75)+ 0.149 * as.numeric(TPQ_Q$QO76)+ 0.132 * as.numeric(TPQ_Q$QO77)+ 0.144 * as.numeric(TPQ_Q$QO78)+ -0.087 * as.numeric(TPQ_Q$QO79)+ -0.163 * as.numeric(TPQ_Q$QO80)+ 0.059 * as.numeric(TPQ_Q$QO81)+ -0.159 * as.numeric(TPQ_Q$QO82)+ 0.1 * as.numeric(TPQ_Q$QO83)+ -0.465 * as.numeric(TPQ_Q$QO84)+ 0.124 * as.numeric(TPQ_Q$QO85)+ -0.036 * as.numeric(TPQ_Q$QO86)+ 0.164 * as.numeric(TPQ_Q$QO87)+ 0.023 * as.numeric(TPQ_Q$QO88)+ -0.352 * as.numeric(TPQ_Q$QO89)+ 0.086 * as.numeric(TPQ_Q$QO90)+ -0.275 * as.numeric(TPQ_Q$QO91)+ -0.054 * as.numeric(TPQ_Q$QO92)+ 0.127 * as.numeric(TPQ_Q$QO93)+ 0.034 * as.numeric(TPQ_Q$QO94)+ -0.174 * as.numeric(TPQ_Q$QO95)+ 0.117 * as.numeric(TPQ_Q$QO96)+ -0.05 * as.numeric(TPQ_Q$QO97)+ -0.275 * as.numeric(TPQ_Q$QO98)+ -0.007 * as.numeric(TPQ_Q$QO99)+ -0.186 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_8_All =  0.174 * as.numeric(TPQ_Q$QO1)+ 0.046 * as.numeric(TPQ_Q$QO2)+ 0.07 * as.numeric(TPQ_Q$QO3)+ 0.085 * as.numeric(TPQ_Q$QO4)+ -0.215 * as.numeric(TPQ_Q$QO5)+ 0.157 * as.numeric(TPQ_Q$QO6)+ 0.011 * as.numeric(TPQ_Q$QO7)+ 0.196 * as.numeric(TPQ_Q$QO8)+ -0.19 * as.numeric(TPQ_Q$QO9)+ 0.18 * as.numeric(TPQ_Q$QO10)+ 0.011 * as.numeric(TPQ_Q$QO11)+ 0.097 * as.numeric(TPQ_Q$QO12)+ 0.108 * as.numeric(TPQ_Q$QO13)+ 0.223 * as.numeric(TPQ_Q$QO14)+ 0.063 * as.numeric(TPQ_Q$QO15)+ -0.102 * as.numeric(TPQ_Q$QO16)+ 0.048 * as.numeric(TPQ_Q$QO17)+ -0.075 * as.numeric(TPQ_Q$QO18)+ 0.073 * as.numeric(TPQ_Q$QO19)+ 0.013 * as.numeric(TPQ_Q$QO20)+ -0.023 * as.numeric(TPQ_Q$QO21)+ 0.099 * as.numeric(TPQ_Q$QO22)+ 0.031 * as.numeric(TPQ_Q$QO23)+ -0.024 * as.numeric(TPQ_Q$QO24)+ 0.124 * as.numeric(TPQ_Q$QO25)+ 0.036 * as.numeric(TPQ_Q$QO26)+ -0.03 * as.numeric(TPQ_Q$QO27)+ -0.205 * as.numeric(TPQ_Q$QO28)+ -0.112 * as.numeric(TPQ_Q$QO29)+ -0.059 * as.numeric(TPQ_Q$QO30)+ -0.038 * as.numeric(TPQ_Q$QO31)+ -0.038 * as.numeric(TPQ_Q$QO32)+ 0.017 * as.numeric(TPQ_Q$QO33)+ -0.005 * as.numeric(TPQ_Q$QO34)+ -0.249 * as.numeric(TPQ_Q$QO35)+ -0.029 * as.numeric(TPQ_Q$QO36)+ -0.032 * as.numeric(TPQ_Q$QO37)+ 0.029 * as.numeric(TPQ_Q$QO38)+ -0.046 * as.numeric(TPQ_Q$QO39)+ 0.086 * as.numeric(TPQ_Q$QO40)+ 0.021 * as.numeric(TPQ_Q$QO41)+ 0.119 * as.numeric(TPQ_Q$QO42)+ 0.042 * as.numeric(TPQ_Q$QO43)+ 0.007 * as.numeric(TPQ_Q$QO44)+ 0.15 * as.numeric(TPQ_Q$QO45)+ -0.294 * as.numeric(TPQ_Q$QO46)+ -0.012 * as.numeric(TPQ_Q$QO47)+ 0.28 * as.numeric(TPQ_Q$QO48)+ 0.16 * as.numeric(TPQ_Q$QO49)+ 0.185 * as.numeric(TPQ_Q$QO50)+ -0.027 * as.numeric(TPQ_Q$QO51)+ 0.025 * as.numeric(TPQ_Q$QO52)+ -0.11 * as.numeric(TPQ_Q$QO53)+ 0.031 * as.numeric(TPQ_Q$QO54)+ -0.286 * as.numeric(TPQ_Q$QO55)+ -0.23 * as.numeric(TPQ_Q$QO56)+ 0.099 * as.numeric(TPQ_Q$QO57)+ -0.069 * as.numeric(TPQ_Q$QO58)+ 0.024 * as.numeric(TPQ_Q$QO59)+ -0.047 * as.numeric(TPQ_Q$QO60)+ -0.072 * as.numeric(TPQ_Q$QO61)+ 0.06 * as.numeric(TPQ_Q$QO62)+ 0.217 * as.numeric(TPQ_Q$QO63)+ -0.036 * as.numeric(TPQ_Q$QO64)+ 0.071 * as.numeric(TPQ_Q$QO65)+ 0.012 * as.numeric(TPQ_Q$QO66)+ 0.077 * as.numeric(TPQ_Q$QO67)+ 0.168 * as.numeric(TPQ_Q$QO68)+ 0.254 * as.numeric(TPQ_Q$QO69)+ 0.09 * as.numeric(TPQ_Q$QO70)+ 0.354 * as.numeric(TPQ_Q$QO71)+ 0.076 * as.numeric(TPQ_Q$QO72)+ -0.048 * as.numeric(TPQ_Q$QO73)+ 0.029 * as.numeric(TPQ_Q$QO74)+ -0.036 * as.numeric(TPQ_Q$QO75)+ 0.052 * as.numeric(TPQ_Q$QO76)+ 0.233 * as.numeric(TPQ_Q$QO77)+ 0.177 * as.numeric(TPQ_Q$QO78)+ 0.185 * as.numeric(TPQ_Q$QO79)+ -0.088 * as.numeric(TPQ_Q$QO80)+ -0.083 * as.numeric(TPQ_Q$QO81)+ 0.043 * as.numeric(TPQ_Q$QO82)+ 0.093 * as.numeric(TPQ_Q$QO83)+ 0.053 * as.numeric(TPQ_Q$QO84)+ -0.134 * as.numeric(TPQ_Q$QO85)+ 0.066 * as.numeric(TPQ_Q$QO86)+ 0.154 * as.numeric(TPQ_Q$QO87)+ 0.058 * as.numeric(TPQ_Q$QO88)+ 0.032 * as.numeric(TPQ_Q$QO89)+ 0.072 * as.numeric(TPQ_Q$QO90)+ 0.352 * as.numeric(TPQ_Q$QO91)+ 0.029 * as.numeric(TPQ_Q$QO92)+ 0.093 * as.numeric(TPQ_Q$QO93)+ 0.066 * as.numeric(TPQ_Q$QO94)+ 0.033 * as.numeric(TPQ_Q$QO95)+ 0.165 * as.numeric(TPQ_Q$QO96)+ 0.036 * as.numeric(TPQ_Q$QO97)+ 0.101 * as.numeric(TPQ_Q$QO98)+ -0.113 * as.numeric(TPQ_Q$QO99)+ -0.094 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_9_All =  -0.014 * as.numeric(TPQ_Q$QO1)+ -0.103 * as.numeric(TPQ_Q$QO2)+ 0.045 * as.numeric(TPQ_Q$QO3)+ 0.041 * as.numeric(TPQ_Q$QO4)+ 0.09 * as.numeric(TPQ_Q$QO5)+ -0.188 * as.numeric(TPQ_Q$QO6)+ 0.154 * as.numeric(TPQ_Q$QO7)+ -0.233 * as.numeric(TPQ_Q$QO8)+ 0.038 * as.numeric(TPQ_Q$QO9)+ 0.232 * as.numeric(TPQ_Q$QO10)+ -0.277 * as.numeric(TPQ_Q$QO11)+ -0.166 * as.numeric(TPQ_Q$QO12)+ -0.229 * as.numeric(TPQ_Q$QO13)+ 0.081 * as.numeric(TPQ_Q$QO14)+ -0.322 * as.numeric(TPQ_Q$QO15)+ 0.035 * as.numeric(TPQ_Q$QO16)+ -0.466 * as.numeric(TPQ_Q$QO17)+ 0.161 * as.numeric(TPQ_Q$QO18)+ 0.047 * as.numeric(TPQ_Q$QO19)+ -0.417 * as.numeric(TPQ_Q$QO20)+ 0.018 * as.numeric(TPQ_Q$QO21)+ 0.034 * as.numeric(TPQ_Q$QO22)+ 0.019 * as.numeric(TPQ_Q$QO23)+ -0.058 * as.numeric(TPQ_Q$QO24)+ -0.626 * as.numeric(TPQ_Q$QO25)+ -0.116 * as.numeric(TPQ_Q$QO26)+ 0.087 * as.numeric(TPQ_Q$QO27)+ 0.01 * as.numeric(TPQ_Q$QO28)+ -0.034 * as.numeric(TPQ_Q$QO29)+ 0.1 * as.numeric(TPQ_Q$QO30)+ 0.022 * as.numeric(TPQ_Q$QO31)+ 0.054 * as.numeric(TPQ_Q$QO32)+ 0.031 * as.numeric(TPQ_Q$QO33)+ 0.103 * as.numeric(TPQ_Q$QO34)+ -0.185 * as.numeric(TPQ_Q$QO35)+ 0.052 * as.numeric(TPQ_Q$QO36)+ -0.029 * as.numeric(TPQ_Q$QO37)+ -0.125 * as.numeric(TPQ_Q$QO38)+ -0.012 * as.numeric(TPQ_Q$QO39)+ -0.035 * as.numeric(TPQ_Q$QO40)+ 0.033 * as.numeric(TPQ_Q$QO41)+ 0.003 * as.numeric(TPQ_Q$QO42)+ 0.07 * as.numeric(TPQ_Q$QO43)+ -0.064 * as.numeric(TPQ_Q$QO44)+ 0.033 * as.numeric(TPQ_Q$QO45)+ 0.048 * as.numeric(TPQ_Q$QO46)+ 0.003 * as.numeric(TPQ_Q$QO47)+ -0.02 * as.numeric(TPQ_Q$QO48)+ 0.016 * as.numeric(TPQ_Q$QO49)+ 0.177 * as.numeric(TPQ_Q$QO50)+ -0.074 * as.numeric(TPQ_Q$QO51)+ 0.042 * as.numeric(TPQ_Q$QO52)+ 0.173 * as.numeric(TPQ_Q$QO53)+ 0.04 * as.numeric(TPQ_Q$QO54)+ 0.049 * as.numeric(TPQ_Q$QO55)+ 0.138 * as.numeric(TPQ_Q$QO56)+ 0.086 * as.numeric(TPQ_Q$QO57)+ -0.277 * as.numeric(TPQ_Q$QO58)+ 0.035 * as.numeric(TPQ_Q$QO59)+ -0.16 * as.numeric(TPQ_Q$QO60)+ 0.162 * as.numeric(TPQ_Q$QO61)+ 0.0 * as.numeric(TPQ_Q$QO62)+ -0.032 * as.numeric(TPQ_Q$QO63)+ 0.165 * as.numeric(TPQ_Q$QO64)+ 0.031 * as.numeric(TPQ_Q$QO65)+ 0.045 * as.numeric(TPQ_Q$QO66)+ 0.119 * as.numeric(TPQ_Q$QO67)+ 0.039 * as.numeric(TPQ_Q$QO68)+ 0.105 * as.numeric(TPQ_Q$QO69)+ -0.04 * as.numeric(TPQ_Q$QO70)+ 0.0 * as.numeric(TPQ_Q$QO71)+ -0.028 * as.numeric(TPQ_Q$QO72)+ 0.042 * as.numeric(TPQ_Q$QO73)+ 0.15 * as.numeric(TPQ_Q$QO74)+ -0.066 * as.numeric(TPQ_Q$QO75)+ 0.003 * as.numeric(TPQ_Q$QO76)+ -0.104 * as.numeric(TPQ_Q$QO77)+ 0.089 * as.numeric(TPQ_Q$QO78)+ 0.037 * as.numeric(TPQ_Q$QO79)+ 0.009 * as.numeric(TPQ_Q$QO80)+ 0.068 * as.numeric(TPQ_Q$QO81)+ 0.08 * as.numeric(TPQ_Q$QO82)+ 0.084 * as.numeric(TPQ_Q$QO83)+ -0.028 * as.numeric(TPQ_Q$QO84)+ 0.258 * as.numeric(TPQ_Q$QO85)+ -0.044 * as.numeric(TPQ_Q$QO86)+ 0.077 * as.numeric(TPQ_Q$QO87)+ -0.036 * as.numeric(TPQ_Q$QO88)+ -0.04 * as.numeric(TPQ_Q$QO89)+ -0.272 * as.numeric(TPQ_Q$QO90)+ -0.154 * as.numeric(TPQ_Q$QO91)+ 0.048 * as.numeric(TPQ_Q$QO92)+ 0.103 * as.numeric(TPQ_Q$QO93)+ 0.082 * as.numeric(TPQ_Q$QO94)+ -0.001 * as.numeric(TPQ_Q$QO95)+ 0.074 * as.numeric(TPQ_Q$QO96)+ -0.137 * as.numeric(TPQ_Q$QO97)+ 0.099 * as.numeric(TPQ_Q$QO98)+ 0.004 * as.numeric(TPQ_Q$QO99)+ -0.173 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_10_All =  0.198 * as.numeric(TPQ_Q$QO1)+ 0.014 * as.numeric(TPQ_Q$QO2)+ 0.029 * as.numeric(TPQ_Q$QO3)+ 0.07 * as.numeric(TPQ_Q$QO4)+ -0.225 * as.numeric(TPQ_Q$QO5)+ -0.135 * as.numeric(TPQ_Q$QO6)+ -0.012 * as.numeric(TPQ_Q$QO7)+ 0.165 * as.numeric(TPQ_Q$QO8)+ 0.036 * as.numeric(TPQ_Q$QO9)+ -0.206 * as.numeric(TPQ_Q$QO10)+ 0.03 * as.numeric(TPQ_Q$QO11)+ -0.139 * as.numeric(TPQ_Q$QO12)+ -0.16 * as.numeric(TPQ_Q$QO13)+ -0.197 * as.numeric(TPQ_Q$QO14)+ -0.151 * as.numeric(TPQ_Q$QO15)+ 0.089 * as.numeric(TPQ_Q$QO16)+ -0.03 * as.numeric(TPQ_Q$QO17)+ -0.12 * as.numeric(TPQ_Q$QO18)+ -0.171 * as.numeric(TPQ_Q$QO19)+ -0.106 * as.numeric(TPQ_Q$QO20)+ 0.049 * as.numeric(TPQ_Q$QO21)+ -0.046 * as.numeric(TPQ_Q$QO22)+ -0.207 * as.numeric(TPQ_Q$QO23)+ -0.143 * as.numeric(TPQ_Q$QO24)+ 0.029 * as.numeric(TPQ_Q$QO25)+ -0.003 * as.numeric(TPQ_Q$QO26)+ -0.005 * as.numeric(TPQ_Q$QO27)+ -0.237 * as.numeric(TPQ_Q$QO28)+ 0.013 * as.numeric(TPQ_Q$QO29)+ -0.183 * as.numeric(TPQ_Q$QO30)+ 0.032 * as.numeric(TPQ_Q$QO31)+ 0.164 * as.numeric(TPQ_Q$QO32)+ -0.242 * as.numeric(TPQ_Q$QO33)+ -0.135 * as.numeric(TPQ_Q$QO34)+ -0.022 * as.numeric(TPQ_Q$QO35)+ 0.058 * as.numeric(TPQ_Q$QO36)+ -0.151 * as.numeric(TPQ_Q$QO37)+ -0.146 * as.numeric(TPQ_Q$QO38)+ -0.057 * as.numeric(TPQ_Q$QO39)+ -0.218 * as.numeric(TPQ_Q$QO40)+ -0.132 * as.numeric(TPQ_Q$QO41)+ 0.063 * as.numeric(TPQ_Q$QO42)+ -0.325 * as.numeric(TPQ_Q$QO43)+ 0.07 * as.numeric(TPQ_Q$QO44)+ -0.037 * as.numeric(TPQ_Q$QO45)+ -0.07 * as.numeric(TPQ_Q$QO46)+ 0.025 * as.numeric(TPQ_Q$QO47)+ -0.144 * as.numeric(TPQ_Q$QO48)+ -0.22 * as.numeric(TPQ_Q$QO49)+ -0.189 * as.numeric(TPQ_Q$QO50)+ 0.083 * as.numeric(TPQ_Q$QO51)+ 0.124 * as.numeric(TPQ_Q$QO52)+ -0.162 * as.numeric(TPQ_Q$QO53)+ -0.719 * as.numeric(TPQ_Q$QO54)+ -0.012 * as.numeric(TPQ_Q$QO55)+ -0.141 * as.numeric(TPQ_Q$QO56)+ -0.725 * as.numeric(TPQ_Q$QO57)+ 0.044 * as.numeric(TPQ_Q$QO58)+ 0.788 * as.numeric(TPQ_Q$QO59)+ 0.008 * as.numeric(TPQ_Q$QO60)+ -0.144 * as.numeric(TPQ_Q$QO61)+ -0.044 * as.numeric(TPQ_Q$QO62)+ 0.355 * as.numeric(TPQ_Q$QO63)+ -0.023 * as.numeric(TPQ_Q$QO64)+ 0.012 * as.numeric(TPQ_Q$QO65)+ 0.044 * as.numeric(TPQ_Q$QO66)+ 0.073 * as.numeric(TPQ_Q$QO67)+ -0.367 * as.numeric(TPQ_Q$QO68)+ -0.502 * as.numeric(TPQ_Q$QO69)+ -0.11 * as.numeric(TPQ_Q$QO70)+ 0.131 * as.numeric(TPQ_Q$QO71)+ -0.122 * as.numeric(TPQ_Q$QO72)+ -0.168 * as.numeric(TPQ_Q$QO73)+ 0.045 * as.numeric(TPQ_Q$QO74)+ 0.333 * as.numeric(TPQ_Q$QO75)+ -0.011 * as.numeric(TPQ_Q$QO76)+ 0.096 * as.numeric(TPQ_Q$QO77)+ -0.053 * as.numeric(TPQ_Q$QO78)+ 0.1 * as.numeric(TPQ_Q$QO79)+ 0.429 * as.numeric(TPQ_Q$QO80)+ 0.039 * as.numeric(TPQ_Q$QO81)+ 0.149 * as.numeric(TPQ_Q$QO82)+ -0.024 * as.numeric(TPQ_Q$QO83)+ 0.134 * as.numeric(TPQ_Q$QO84)+ -0.109 * as.numeric(TPQ_Q$QO85)+ 0.002 * as.numeric(TPQ_Q$QO86)+ 0.006 * as.numeric(TPQ_Q$QO87)+ -0.035 * as.numeric(TPQ_Q$QO88)+ 0.183 * as.numeric(TPQ_Q$QO89)+ -0.213 * as.numeric(TPQ_Q$QO90)+ 0.13 * as.numeric(TPQ_Q$QO91)+ 0.196 * as.numeric(TPQ_Q$QO92)+ 0.024 * as.numeric(TPQ_Q$QO93)+ 0.033 * as.numeric(TPQ_Q$QO94)+ 0.178 * as.numeric(TPQ_Q$QO95)+ -0.045 * as.numeric(TPQ_Q$QO96)+ 0.213 * as.numeric(TPQ_Q$QO97)+ 0.181 * as.numeric(TPQ_Q$QO98)+ 0.013 * as.numeric(TPQ_Q$QO99)+ 0.144 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_11_All =  -0.022 * as.numeric(TPQ_Q$QO1)+ -0.031 * as.numeric(TPQ_Q$QO2)+ -0.689 * as.numeric(TPQ_Q$QO3)+ -0.088 * as.numeric(TPQ_Q$QO4)+ 0.179 * as.numeric(TPQ_Q$QO5)+ -0.302 * as.numeric(TPQ_Q$QO6)+ -0.154 * as.numeric(TPQ_Q$QO7)+ -0.082 * as.numeric(TPQ_Q$QO8)+ 0.017 * as.numeric(TPQ_Q$QO9)+ 0.117 * as.numeric(TPQ_Q$QO10)+ 0.077 * as.numeric(TPQ_Q$QO11)+ 0.631 * as.numeric(TPQ_Q$QO12)+ 0.016 * as.numeric(TPQ_Q$QO13)+ 0.077 * as.numeric(TPQ_Q$QO14)+ 0.348 * as.numeric(TPQ_Q$QO15)+ 0.087 * as.numeric(TPQ_Q$QO16)+ 0.126 * as.numeric(TPQ_Q$QO17)+ 0.06 * as.numeric(TPQ_Q$QO18)+ 0.116 * as.numeric(TPQ_Q$QO19)+ 0.108 * as.numeric(TPQ_Q$QO20)+ 0.034 * as.numeric(TPQ_Q$QO21)+ 0.072 * as.numeric(TPQ_Q$QO22)+ 0.06 * as.numeric(TPQ_Q$QO23)+ 0.008 * as.numeric(TPQ_Q$QO24)+ 0.094 * as.numeric(TPQ_Q$QO25)+ -0.009 * as.numeric(TPQ_Q$QO26)+ -0.053 * as.numeric(TPQ_Q$QO27)+ 0.051 * as.numeric(TPQ_Q$QO28)+ 0.043 * as.numeric(TPQ_Q$QO29)+ -0.053 * as.numeric(TPQ_Q$QO30)+ -0.067 * as.numeric(TPQ_Q$QO31)+ 0.061 * as.numeric(TPQ_Q$QO32)+ 0.109 * as.numeric(TPQ_Q$QO33)+ 0.002 * as.numeric(TPQ_Q$QO34)+ -0.065 * as.numeric(TPQ_Q$QO35)+ 0.185 * as.numeric(TPQ_Q$QO36)+ 0.187 * as.numeric(TPQ_Q$QO37)+ 0.064 * as.numeric(TPQ_Q$QO38)+ 0.08 * as.numeric(TPQ_Q$QO39)+ 0.064 * as.numeric(TPQ_Q$QO40)+ 0.028 * as.numeric(TPQ_Q$QO41)+ -0.065 * as.numeric(TPQ_Q$QO42)+ 0.072 * as.numeric(TPQ_Q$QO43)+ -0.082 * as.numeric(TPQ_Q$QO44)+ 0.024 * as.numeric(TPQ_Q$QO45)+ 0.013 * as.numeric(TPQ_Q$QO46)+ -0.032 * as.numeric(TPQ_Q$QO47)+ 0.012 * as.numeric(TPQ_Q$QO48)+ 0.017 * as.numeric(TPQ_Q$QO49)+ 0.036 * as.numeric(TPQ_Q$QO50)+ -0.041 * as.numeric(TPQ_Q$QO51)+ 0.008 * as.numeric(TPQ_Q$QO52)+ 0.107 * as.numeric(TPQ_Q$QO53)+ 0.056 * as.numeric(TPQ_Q$QO54)+ 0.012 * as.numeric(TPQ_Q$QO55)+ 0.153 * as.numeric(TPQ_Q$QO56)+ 0.039 * as.numeric(TPQ_Q$QO57)+ -0.022 * as.numeric(TPQ_Q$QO58)+ -0.026 * as.numeric(TPQ_Q$QO59)+ -0.037 * as.numeric(TPQ_Q$QO60)+ -0.072 * as.numeric(TPQ_Q$QO61)+ -0.086 * as.numeric(TPQ_Q$QO62)+ 0.035 * as.numeric(TPQ_Q$QO63)+ -0.085 * as.numeric(TPQ_Q$QO64)+ -0.027 * as.numeric(TPQ_Q$QO65)+ 0.01 * as.numeric(TPQ_Q$QO66)+ -0.139 * as.numeric(TPQ_Q$QO67)+ 0.023 * as.numeric(TPQ_Q$QO68)+ 0.042 * as.numeric(TPQ_Q$QO69)+ -0.008 * as.numeric(TPQ_Q$QO70)+ 0.068 * as.numeric(TPQ_Q$QO71)+ 0.044 * as.numeric(TPQ_Q$QO72)+ 0.096 * as.numeric(TPQ_Q$QO73)+ -0.31 * as.numeric(TPQ_Q$QO74)+ -0.031 * as.numeric(TPQ_Q$QO75)+ 0.05 * as.numeric(TPQ_Q$QO76)+ -0.046 * as.numeric(TPQ_Q$QO77)+ 0.036 * as.numeric(TPQ_Q$QO78)+ 0.053 * as.numeric(TPQ_Q$QO79)+ 0.028 * as.numeric(TPQ_Q$QO80)+ -0.046 * as.numeric(TPQ_Q$QO81)+ -0.016 * as.numeric(TPQ_Q$QO82)+ -0.127 * as.numeric(TPQ_Q$QO83)+ -0.016 * as.numeric(TPQ_Q$QO84)+ -0.092 * as.numeric(TPQ_Q$QO85)+ 0.536 * as.numeric(TPQ_Q$QO86)+ -0.039 * as.numeric(TPQ_Q$QO87)+ 0.78 * as.numeric(TPQ_Q$QO88)+ -0.048 * as.numeric(TPQ_Q$QO89)+ 0.348 * as.numeric(TPQ_Q$QO90)+ -0.053 * as.numeric(TPQ_Q$QO91)+ 0.017 * as.numeric(TPQ_Q$QO92)+ 0.197 * as.numeric(TPQ_Q$QO93)+ -0.041 * as.numeric(TPQ_Q$QO94)+ -0.006 * as.numeric(TPQ_Q$QO95)+ 0.1 * as.numeric(TPQ_Q$QO96)+ -0.042 * as.numeric(TPQ_Q$QO97)+ -0.07 * as.numeric(TPQ_Q$QO98)+ 0.013 * as.numeric(TPQ_Q$QO99)+ -0.094 * as.numeric(TPQ_Q$QO100)
TPQ_Q$Complete_IC_12_All =  -0.143 * as.numeric(TPQ_Q$QO1)+ -0.261 * as.numeric(TPQ_Q$QO2)+ -0.041 * as.numeric(TPQ_Q$QO3)+ -0.286 * as.numeric(TPQ_Q$QO4)+ 0.068 * as.numeric(TPQ_Q$QO5)+ 0.043 * as.numeric(TPQ_Q$QO6)+ -0.084 * as.numeric(TPQ_Q$QO7)+ -0.176 * as.numeric(TPQ_Q$QO8)+ -0.009 * as.numeric(TPQ_Q$QO9)+ 0.009 * as.numeric(TPQ_Q$QO10)+ 0.06 * as.numeric(TPQ_Q$QO11)+ -0.001 * as.numeric(TPQ_Q$QO12)+ -0.166 * as.numeric(TPQ_Q$QO13)+ -0.011 * as.numeric(TPQ_Q$QO14)+ -0.021 * as.numeric(TPQ_Q$QO15)+ 0.025 * as.numeric(TPQ_Q$QO16)+ -0.172 * as.numeric(TPQ_Q$QO17)+ 0.148 * as.numeric(TPQ_Q$QO18)+ 0.075 * as.numeric(TPQ_Q$QO19)+ -0.132 * as.numeric(TPQ_Q$QO20)+ 0.017 * as.numeric(TPQ_Q$QO21)+ -0.111 * as.numeric(TPQ_Q$QO22)+ 0.055 * as.numeric(TPQ_Q$QO23)+ -0.238 * as.numeric(TPQ_Q$QO24)+ -0.087 * as.numeric(TPQ_Q$QO25)+ -0.396 * as.numeric(TPQ_Q$QO26)+ 0.005 * as.numeric(TPQ_Q$QO27)+ -0.056 * as.numeric(TPQ_Q$QO28)+ -0.681 * as.numeric(TPQ_Q$QO29)+ -0.013 * as.numeric(TPQ_Q$QO30)+ -0.05 * as.numeric(TPQ_Q$QO31)+ -0.04 * as.numeric(TPQ_Q$QO32)+ 0.065 * as.numeric(TPQ_Q$QO33)+ -0.031 * as.numeric(TPQ_Q$QO34)+ 0.131 * as.numeric(TPQ_Q$QO35)+ -0.115 * as.numeric(TPQ_Q$QO36)+ -0.004 * as.numeric(TPQ_Q$QO37)+ -0.018 * as.numeric(TPQ_Q$QO38)+ -0.037 * as.numeric(TPQ_Q$QO39)+ 0.128 * as.numeric(TPQ_Q$QO40)+ -0.044 * as.numeric(TPQ_Q$QO41)+ 0.051 * as.numeric(TPQ_Q$QO42)+ -0.163 * as.numeric(TPQ_Q$QO43)+ -0.024 * as.numeric(TPQ_Q$QO44)+ -0.048 * as.numeric(TPQ_Q$QO45)+ 0.056 * as.numeric(TPQ_Q$QO46)+ -0.8 * as.numeric(TPQ_Q$QO47)+ -0.203 * as.numeric(TPQ_Q$QO48)+ -0.044 * as.numeric(TPQ_Q$QO49)+ -0.093 * as.numeric(TPQ_Q$QO50)+ -0.769 * as.numeric(TPQ_Q$QO51)+ -0.023 * as.numeric(TPQ_Q$QO52)+ 0.086 * as.numeric(TPQ_Q$QO53)+ 0.055 * as.numeric(TPQ_Q$QO54)+ 0.018 * as.numeric(TPQ_Q$QO55)+ -0.068 * as.numeric(TPQ_Q$QO56)+ 0.058 * as.numeric(TPQ_Q$QO57)+ -0.079 * as.numeric(TPQ_Q$QO58)+ -0.117 * as.numeric(TPQ_Q$QO59)+ -0.126 * as.numeric(TPQ_Q$QO60)+ -0.008 * as.numeric(TPQ_Q$QO61)+ -0.165 * as.numeric(TPQ_Q$QO62)+ -0.117 * as.numeric(TPQ_Q$QO63)+ -0.0 * as.numeric(TPQ_Q$QO64)+ -0.096 * as.numeric(TPQ_Q$QO65)+ 0.002 * as.numeric(TPQ_Q$QO66)+ -0.027 * as.numeric(TPQ_Q$QO67)+ 0.052 * as.numeric(TPQ_Q$QO68)+ 0.109 * as.numeric(TPQ_Q$QO69)+ -0.059 * as.numeric(TPQ_Q$QO70)+ 0.001 * as.numeric(TPQ_Q$QO71)+ -0.054 * as.numeric(TPQ_Q$QO72)+ 0.085 * as.numeric(TPQ_Q$QO73)+ -0.083 * as.numeric(TPQ_Q$QO74)+ -0.218 * as.numeric(TPQ_Q$QO75)+ 0.026 * as.numeric(TPQ_Q$QO76)+ -0.094 * as.numeric(TPQ_Q$QO77)+ 0.049 * as.numeric(TPQ_Q$QO78)+ -0.162 * as.numeric(TPQ_Q$QO79)+ -0.121 * as.numeric(TPQ_Q$QO80)+ -0.032 * as.numeric(TPQ_Q$QO81)+ -0.041 * as.numeric(TPQ_Q$QO82)+ -0.051 * as.numeric(TPQ_Q$QO83)+ -0.115 * as.numeric(TPQ_Q$QO84)+ 0.258 * as.numeric(TPQ_Q$QO85)+ -0.019 * as.numeric(TPQ_Q$QO86)+ 0.01 * as.numeric(TPQ_Q$QO87)+ 0.038 * as.numeric(TPQ_Q$QO88)+ -0.066 * as.numeric(TPQ_Q$QO89)+ 0.059 * as.numeric(TPQ_Q$QO90)+ -0.138 * as.numeric(TPQ_Q$QO91)+ -0.036 * as.numeric(TPQ_Q$QO92)+ 0.051 * as.numeric(TPQ_Q$QO93)+ 0.002 * as.numeric(TPQ_Q$QO94)+ -0.063 * as.numeric(TPQ_Q$QO95)+ 0.153 * as.numeric(TPQ_Q$QO96)+ -0.123 * as.numeric(TPQ_Q$QO97)+ -0.032 * as.numeric(TPQ_Q$QO98)+ 0.032 * as.numeric(TPQ_Q$QO99)+ 0.055 * as.numeric(TPQ_Q$QO100)


#Log regression
#LogisticInfo(Model1)
#LogisticFunction(Model1)
TPQ_Q_log = TPQ_Q[(TPQ_Q$Diagnosis!="GAD" & TPQ_Q$Session=="Test"),]
TPQ_Q_log$Diagnosis=factor(TPQ_Q_log$Diagnosis)
TPQ_Q_log_MDD = TPQ_Q_log[TPQ_Q_log$Diagnosis=="MDD",]
TPQ_Q_log_MDD$Response = factor (TPQ_Q_log_MDD$Response)


#1- for diagnosis
for (Variable in c("NS1","NS2","NS3","NS4","HA1","HA2","HA3","HA4","RD1","RD2","RD3","RD4")){
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  p_value = summary(Model1)$coefficients[2,4]
  if (p_value<0.05){
    print (Formula)
    print (p_value)
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i)
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  print (Variable)
  print (summary(Model1))
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_HC")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  p_value = summary(Model1)$coefficients[2,4]
  if (p_value<0.05){
    print (Formula)
    print (p_value)
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_MDD")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  p_value = summary(Model1)$coefficients[2,4]
  if (p_value<0.05){
    print (Formula)
    print (p_value)
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_All")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  p_value = summary(Model1)$coefficients[2,4]
  if (p_value<0.05){
    print (Formula)
    print (p_value)
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_HC")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_MDD")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_All")
  Data = TPQ_Q_log[!is.na(TPQ_Q_log[[Variable]]),]
  Formula = as.formula(paste0("Diagnosis ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

#1- for response
for (Variable in c("NS1","NS2","NS3","NS4","HA1","HA2","HA3","HA4","RD1","RD2","RD3","RD4")){
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i)
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_HC")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_MDD")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("IC_",i,"_All")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_HC")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_MDD")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.05){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

for (i in 0:12){
  Variable=paste0("Complete_IC_",i,"_All")
  Data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD[[Variable]]),]
  Formula = as.formula(paste0("Response ~ ",Variable))
  Model1=glm(Formula, data = Data, family = binomial())
  if (summary(Model1)$coefficients[2,4]<0.1){
    print (Variable)
    print (summary(Model1)$coefficients[2,4])
  }
}

# RD4, IC_10, Complete_IC_8_MDD and Complete_IC_10_MDD were significant
#Also, in the 'All' data set, IC_6_All and IC_9_All were significant
Model_RD4 = glm(Response ~ RD4, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$RD4),],family = binomial())
summary(Model_RD4)
LogisticFunction(Model_RD4,0.7,"glm")
LogisticFunction(Model_RD4,0.7)

Model_IC10 = glm(Response ~ IC_10, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$IC_10),],family = binomial())
summary(Model_IC10)
LogisticFunction(Model_IC10,0.5, "glm")

Model_IC_8 = glm(Response ~ Complete_IC_8_MDD, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$Complete_IC_8_MDD),],family = binomial())
summary(Model_IC_8)
LogisticFunction(Model_IC_8,FALSE,0.5,"glm")
LogisticFunction(Model_IC_8,FALSE,0.5)

Model_IC_10 = glm(Response ~ Complete_IC_10_MDD, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$Complete_IC_10_MDD),],family = binomial())
summary(Model_IC_10)
LogisticFunction(Model_IC_10,FALSE,0.7,"glm")
LogisticFunction(Model_IC_10,FALSE,0.7)

Model_IC_6 = glm(Response ~ IC_6_All, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$IC_6_All),],family = binomial())
summary(Model_IC_6)
LogisticFunction(Model_IC_6,FALSE,0.7,"glm")
LogisticFunction(Model_IC_6,FALSE,0.7)

Model_IC_9 = glm(Response ~ IC_9_All, data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$IC_9_All),],family = binomial())
summary(Model_IC_9)
LogisticFunction(Model_IC_9,FALSE,0.7,"glm")
LogisticFunction(Model_IC_9,FALSE,0.7)

TPQ_Q_log_MDD$NumRes = as.character(TPQ_Q_log_MDD$Response)
TPQ_Q_log_MDD$NumRes[TPQ_Q_log_MDD$NumRes=="Responder"]=1
TPQ_Q_log_MDD$NumRes[TPQ_Q_log_MDD$NumRes=="Non-responder"]=0
TPQ_Q_log_MDD$NumRes = as.numeric(TPQ_Q_log_MDD$NumRes)

Model_IC_All = step(glm(Response ~ RD4+Complete_IC_8_MDD+Complete_IC_10_MDD+IC_6_All+IC_9_All, 
                        data = TPQ_Q_log_MDD[!(is.na(TPQ_Q_log_MDD$Complete_IC_8_MDD)|is.na(TPQ_Q_log_MDD$Complete_IC_10_MDD)|
                                                 is.na(TPQ_Q_log_MDD$Complete_IC_10_MDD)|is.na(TPQ_Q_log_MDD$IC_9_All)|is.na(TPQ_Q_log_MDD$IC_6_All)|is.na(TPQ_Q_log_MDD$RD4)),],
                        family=binomial()), direction = "both")
summary(Model_IC_All)
LogisticFunction(Model_IC_All,0.7)#,"glm")
LogisticInfo(Model_IC_All,PrintList = TRUE)
#for plotting to send by email

LogisticFunction(Model_RD4,FALSE,0.7,"glm")
LogisticFunction(Model_RD4,FALSE,0.7)
LogisticFunction(Model_IC_8,FALSE,0.5,"glm")
LogisticFunction(Model_IC_8,FALSE,0.5)
LogisticFunction(Model_IC_10,FALSE,0.5,"glm")
LogisticFunction(Model_IC_10,FALSE,0.5)

TPQ_Q_log_MDD$ResponseNum[TPQ_Q_log_MDD$Response=="Responder"]=1
TPQ_Q_log_MDD$ResponseNum[TPQ_Q_log_MDD$Response=="Non-responder"]=0
ggplot(data = TPQ_Q_log_MDD[!is.na(TPQ_Q_log_MDD$RD4),],mapping = aes(x=RD4,y=ResponseNum))+
  geom_jitter(shape = 21, color = "#7c0a02", fill="white",size = 4,width = 0.1, height = 0)+
  stat_smooth(method = "glm", method.args = list(family = "binomial"),se=F,color="#7c0a02",size=4)+
  TypicalTheme


TPQ_Q_MDD_HC = TPQ_Q_log#[TPQ_Q_log$GAD!="GAD",]
TPQ_Q_MDD_HC$ResponseNum[TPQ_Q_MDD_HC$Response=="Responder"]=2
TPQ_Q_MDD_HC$ResponseNum[TPQ_Q_MDD_HC$Response=="Non-responder"]=1
TPQ_Q_MDD_HC$ResponseNum[TPQ_Q_MDD_HC$Response=="NoResponse"]=0

Model_RD4_All = glm(Response ~ RD4, data = TPQ_Q_MDD_HC[!is.na(TPQ_Q_MDD_HC$RD4),],family = binomial())
summary(Model_RD4_All)
ggplot(data = TPQ_Q_MDD_HC[!is.na(TPQ_Q_MDD_HC$RD4),],mapping = aes(x=RD4,y=ResponseNum))+
  geom_jitter(shape = 21, color = "#7c0a02", fill="white",size = 4,width = 0.1, height = 0)+
  stat_smooth(method = "glm", method.args = list(family = "binomial"),se=F,color="#7c0a02",size=4)+
  TypicalTheme


# Now I will do glm for all variables included ... for Response
TPQ_Q_log_MDD$Response=as.character(TPQ_Q_log_MDD$Response)
TPQ_Q_log_MDD$Response[TPQ_Q_log_MDD$Response=="Non-responder"]="NonResponder"
TPQ_Q_log_MDD$Response = as.factor(TPQ_Q_log_MDD$Response)

fitControl1 = trainControl(method = "cv", number = 5, savePredictions = T,classProbs=TRUE)

IC_Complete_List = NULL
IC_Complete_MDD_List = NULL
IC_List = NULL
CloningerList = c("NS","NS1","NS2","NS2","NS3","NS4","HA","HA1","HA2","HA2","HA3","HA4","RD","RD1","RD2","RD2","RD3","RD4")
for (i in 0:12){
  IC_Complete_List = append(IC_Complete_List, paste0("Complete_IC_",i,"_All"))
  IC_Complete_MDD_List = append(IC_Complete_MDD_List, paste0("Complete_IC_",i,"_MDD"))
  IC_List = append(IC_List, paste0("IC_",i,"_All"))
}


IC_Complete_Formula = paste0("Response ~ ",IC_Complete_List[1])
IC_Complete_MDD_Formula = paste0("Response ~ ",IC_Complete_MDD_List[1])
IC_Formula = paste0("Response ~ ",IC_List[1])
Cloninger_Formula = paste0("Response ~ ",CloningerList[1])
for (i in 2:length(IC_Complete_List)){
  IC_Complete_Formula=paste0(IC_Complete_Formula," + ",IC_Complete_List[i])
  IC_Complete_MDD_Formula=paste0(IC_Complete_MDD_Formula," + ",IC_Complete_MDD_List[i])
  IC_Formula=paste0(IC_Formula," + ",IC_List[i])
  Cloninger_Formula=paste0(Cloninger_Formula," + ",CloningerList[i])
}

No.NA.Data_Complete = na.omit(TPQ_Q_log_MDD[append("Response",IC_Complete_List)])
No.NA.Data_Complete_MDD = na.omit(TPQ_Q_log_MDD[append("Response",IC_Complete_MDD_List)])
No.NA.Data = na.omit(TPQ_Q_log_MDD[append("Response",IC_List)])
No.NA.Data_Cloninger = na.omit(TPQ_Q_log_MDD[append("Response",CloningerList)])

Model_Complete_All = step(glm(as.formula(IC_Complete_Formula), data = No.NA.Data_Complete, family = binomial()), direction = "both")
Model_Complete_MDD = step(glm(as.formula(IC_Complete_MDD_Formula), data = No.NA.Data_Complete_MDD, family = binomial()), direction = "both")
Model_All = step(glm(as.formula(IC_Formula), data = No.NA.Data, family = binomial()), direction = "both")
Model_Cloninger = step(glm(as.formula(Cloninger_Formula), data = No.NA.Data_Cloninger, family = binomial()), direction = "both")
cv_mod_Complete_All = train(Model_Complete_All$formula, data = No.NA.Data_Complete, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_Complete_MDD = train(Model_Complete_MDD$formula, data = No.NA.Data_Complete_MDD, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_All = train(Model_All$formula, data = No.NA.Data, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_Cloninger = train(Model_Cloninger$formula, data = No.NA.Data_Cloninger, method = "glm", family = "binomial",trControl = fitControl1)


summary(Model_Complete_All)
summary(Model_Complete_MDD)
summary(Model_All)
summary(Model_Cloninger)
summary(cv_mod_Complete_All)
summary(cv_mod_Complete_MDD)
summary(cv_mod_All)
summary(cv_mod_Cloninger)

confusionMatrix(table((cv_mod_Complete_All$pred)$pred,(cv_mod_Complete_All$pred)$obs))
confusionMatrix(table((cv_mod_Complete_MDD$pred)$pred,(cv_mod_Complete_MDD$pred)$obs))
confusionMatrix(table((cv_mod_All$pred)$pred,(cv_mod_All$pred)$obs))
confusionMatrix(table((cv_mod_Cloninger$pred)$pred,(cv_mod_Cloninger$pred)$obs))


LogisticFunction(Model_Complete_All,Threshold = 0.75)
LogisticFunction(Model_Complete_MDD,Threshold = 0.6)
LogisticFunction(Model_All,Threshold = 0.75)
LogisticFunction(Model_Cloninger,Threshold = 0.7)


confusionMatrix(table((cv_mod_Complete_All$pred)$pred,(cv_mod_Complete_All$pred)$obs))
confusionMatrix(table((cv_mod_Complete_MDD$pred)$pred,(cv_mod_Complete_MDD$pred)$obs))
confusionMatrix(table((cv_mod_All$pred)$pred,(cv_mod_All$pred)$obs))
confusionMatrix(table((cv_mod_Cloninger$pred)$pred,(cv_mod_Cloninger$pred)$obs))

#ALL
DF = (cv_mod_Complete_All$pred)
DF$pred = ifelse(DF$Responder>0.75,"Responder","NonResponder")
DF$pred=factor(DF$pred,levels = c("Responder","NonResponder"))
DF$obs = factor(DF$obs,levels = c("Responder","NonResponder"))
confusionMatrix(table(DF$pred,DF$obs))

              Responder NonResponder
Responder           15            4
NonResponder        19           14
                   N.R2 : 0.2295154
               Accuracy : 0.5577
            Sensitivity : 0.4412          
            Specificity : 0.7778          
         Pos Pred Value : 0.7895          
         Neg Pred Value : 0.4242  

#Not Complete
DF = (cv_mod_All$pred)
DF$pred = ifelse(DF$Responder>0.75,"Responder","NonResponder")
DF$pred=factor(DF$pred,levels = c("Responder","NonResponder"))
DF$obs = factor(DF$obs,levels = c("Responder","NonResponder"))
confusionMatrix(table(DF$pred,DF$obs))
  Responder NonResponder
Responder 21  4
NonResponder  14  14
N.R2: 0.38867
Accuracy 0.6604
Sensitivity 0.6000          
Specificity 0.7778          
Pos Pred Value 0.8400          
Neg Pred Value 0.5000
#Cloninger
DF = (cv_mod_Cloninger$pred)
DF$pred = ifelse(DF$Responder>0.7,"Responder","NonResponder")
DF$pred=factor(DF$pred,levels = c("Responder","NonResponder"))
DF$obs = factor(DF$obs,levels = c("Responder","NonResponder"))
confusionMatrix(table(DF$pred,DF$obs))


  Responder NonResponder
Responder 19  5
NonResponder  15  13
N.R2 0.2042228
Accuracy 0.6154         
Sensitivity 0.5588         
Specificity 0.7222         
PPV 0.7917         
NPV 0.4643
###################################################################################################
##trail, delete
###################################################################################################
No.NA.Data_Complete$Response2 = as.character(No.NA.Data_Complete$Response)
No.NA.Data_Complete$Response2[No.NA.Data_Complete$Response2=="Non-responder"] = "NonResponder"
cv_mod_Complete_ROC = train(Response2 ~ Complete_IC_7_All + Complete_IC_9_All + Complete_IC_12_All,
                            metric="ROC",data = No.NA.Data_Complete,
                            method = "glm", family = "binomial",
                            trControl = trainControl(method = "cv",
                                                     classProbs = TRUE,savePredictions = T))

roc0 <- roc(No.NA.Data_Complete$Response2, 
            predict(cv_mod_Complete_ROC, No.NA.Data_Complete, type = "prob")[,1], 
            levels = rev(levels(No.NA.Data_Complete$Response2)))


No.NA.Data_Complete2 = No.NA.Data_Complete
No.NA.Data_Complete2$Response=as.character(No.NA.Data_Complete2$Response)
No.NA.Data_Complete2$Response[No.NA.Data_Complete2$Response=="Non-responder"]=0
No.NA.Data_Complete2$Response[No.NA.Data_Complete2$Response=="Responder"]=1
No.NA.Data_Complete2$Response = as.numeric(No.NA.Data_Complete2$Response)
swissmodel = train(as.formula(IC_Complete_Formula),
                   data = No.NA.Data_Complete2,
                   method = "leapSeq",
                   tuneGrid = data.frame(nvmax=1:13),
                   family = "binomial",
                   trControl = trainControl(method = "cv",number=10))

summary(swissmodel$finalModel)
coef(swissmodel$finalModel,7)
summary(swissmodel)
################################################################################################### 

#Same will be done for GAD
TPQ_Q_log_MDD$GAD = relevel(TPQ_Q_log_MDD$GAD,"NoGad")
IC_Complete_List = NULL
IC_Complete_MDD_List = NULL
IC_List = NULL
CloningerList = c("NS","NS1","NS2","NS2","NS3","NS4","HA","HA1","HA2","HA2","HA3","HA4","RD","RD1","RD2","RD2","RD3","RD4")
for (i in 0:12){
  IC_Complete_List = append(IC_Complete_List, paste0("Complete_IC_",i,"_All"))
  IC_Complete_MDD_List = append(IC_Complete_MDD_List, paste0("Complete_IC_",i,"_MDD"))
  IC_List = append(IC_List, paste0("IC_",i,"_All"))
}


IC_Complete_Formula = paste0("GAD ~ ",IC_Complete_List[1])
IC_Complete_MDD_Formula = paste0("GAD ~ ",IC_Complete_MDD_List[1])
IC_Formula = paste0("GAD ~ ",IC_List[1])
Cloninger_Formula = paste0("GAD ~ ",CloningerList[1])
for (i in 2:length(IC_Complete_List)){
  IC_Complete_Formula=paste0(IC_Complete_Formula," + ",IC_Complete_List[i])
  IC_Complete_MDD_Formula=paste0(IC_Complete_MDD_Formula," + ",IC_Complete_MDD_List[i])
  IC_Formula=paste0(IC_Formula," + ",IC_List[i])
  Cloninger_Formula=paste0(Cloninger_Formula," + ",CloningerList[i])
}

No.NA.Data_Complete = na.omit(TPQ_Q_log_MDD[append("GAD",IC_Complete_List)])
No.NA.Data_Complete_MDD = na.omit(TPQ_Q_log_MDD[append("GAD",IC_Complete_MDD_List)])
No.NA.Data = na.omit(TPQ_Q_log_MDD[append("GAD",IC_List)])
No.NA.Data_Cloninger = na.omit(TPQ_Q_log_MDD[append("GAD",CloningerList)])
Model_Complete_All = step(glm(as.formula(IC_Complete_Formula), data = No.NA.Data_Complete, family = binomial()), direction = "back")
Model_Complete_MDD = step(glm(as.formula(IC_Complete_MDD_Formula), data = No.NA.Data_Complete_MDD, family = binomial()), direction = "back")
Model_All = step(glm(as.formula(IC_Formula), data = No.NA.Data, family = binomial()), direction = "back")
Model_Cloninger = step(glm(as.formula(Cloninger_Formula), data = No.NA.Data_Cloninger, family = binomial()), direction = "back")
cv_mod_Complete_All = train(Model_Complete_All$formula, data = No.NA.Data_Complete, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_Complete_MDD = train(Model_Complete_MDD$formula, data = No.NA.Data_Complete_MDD, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_All = train(Model_All$formula, data = No.NA.Data, method = "glm", family = "binomial",trControl = fitControl1)
cv_mod_Cloninger = train(Model_Cloninger$formula, data = No.NA.Data_Cloninger, method = "glm", family = "binomial",trControl = fitControl1)

summary(Model_Complete_All)
summary(Model_Complete_MDD)
summary(Model_All)
summary(Model_Cloninger)

summary(cv_mod_Complete_All)
summary(cv_mod_Complete_MDD)
summary(cv_mod_All)
summary(cv_mod_Cloninger)

confusionMatrix(table((cv_mod_Complete_All$pred)$pred,(cv_mod_Complete_All$pred)$obs))
confusionMatrix(table((cv_mod_Complete_MDD$pred)$pred,(cv_mod_Complete_MDD$pred)$obs))
confusionMatrix(table((cv_mod_All$pred)$pred,(cv_mod_All$pred)$obs))
confusionMatrix(table((cv_mod_Cloninger$pred)$pred,(cv_mod_Cloninger$pred)$obs))


LogisticFunction(Model_All,Threshold = 0.3)
LogisticFunction(Model_Complete_All,Threshold = 0.3)
LogisticFunction(Model_Complete_MDD,Threshold = 0.5)
LogisticFunction(Model_Cloninger,Threshold = 0.275)

#ALL
DF = (cv_mod_All$pred)
DF$pred = ifelse(DF$GAD>0.3,"GAD","NoGad")
DF$pred=factor(DF$pred,levels = c("GAD","NoGad"))
DF$obs = factor(DF$obs,levels = c("GAD","NoGad"))
confusionMatrix(table(DF$pred,DF$obs))



      GAD NoGad
GAD    22    10
NoGad   4    16
N.R2 = 0.4984462
Accuracy : 0.75
Sensitivity : 0.8462          
Specificity : 0.6154          
Pos Pred Value : 0.6875          
Neg Pred Value : 0.8000 

#Not complete
DF = (cv_mod_Complete_All$pred)
DF$pred = ifelse(DF$GAD>0.5,"GAD","NoGad")
DF$pred=factor(DF$pred,levels = c("GAD","NoGad"))
DF$obs = factor(DF$obs,levels = c("GAD","NoGad"))
confusionMatrix(table(DF$pred,DF$obs))
  GAD NoGad
GAD    19     7
NoGad   7    19

N.R2=0.6017277
Accuracy : 0.7308          
Sensitivity : 0.7308          
Specificity : 0.7308          
Pos Pred Value : 0.7308          
Neg Pred Value : 0.7308      
#Cloninger
DF = (cv_mod_Cloninger$pred)
DF$pred = ifelse(DF$GAD>0.275,"GAD","NoGad")
DF$pred=factor(DF$pred,levels = c("GAD","NoGad"))
DF$obs = factor(DF$obs,levels = c("GAD","NoGad"))
confusionMatrix(table(DF$pred,DF$obs))

  GAD NoGad
GAD 25 17
NoGad 1 9
N.R2:0.2082975
Accuracy:0.6154
Sensitivity:0.9615          
Specificity:0.3462          
Pos Pred Value:0.5952          
Neg Pred Value:0.9000

No.NA.Data_RD = na.omit(TPQ_Q_log_MDD[c("GAD","RD4")])
Model_RD4_GAD = glm(GAD ~ RD4, data = No.NA.Data_RD,family = binomial())
summary(Model_RD4_GAD)
LogisticFunction(Model_RD4_GAD)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#
#Now we move on to Mohammad's request to highlight questions with the biggest differences
#In the email send by Juergen, the components IC2, IC4, IC7, IC10 and IC12 had the largest differences between MDD and HCs
#
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

variable_list = mapply(paste0,"Complete_IC_",0:12,"_HC",SIMPLIFY = T) #Just a list to define the required columns
variable_names = mapply(paste0,"IC_",0:12,SIMPLIFY = T)
variable_list = append(c("Subject","Diagnosis","Response","Session"),variable_list)

NewData = TPQ_Q_log[variable_list]

#Take the last 12 columns and rename them
colnames(NewData)[(ncol(NewData)-12):ncol(NewData)] = 0:12
NewData_Melt = melt(NewData,id.vars = c("Subject","Diagnosis","Response","Session"),value.name = "Score", variable.name = "IC_Number")  

DateMeans = SMeans(NewData_Melt,"Score",c("IC_Number","Diagnosis"),GroupBy = "Diagnosis")
DateMeans$IC_Number = factor(DateMeans$IC_Number, levels = 0:12)
ggplot(DateMeans,aes(x=IC_Number,y= Score, fill=Groups))+geom_bar(stat="identity",position = position_dodge(0.75),color="black")+
  TypicalTheme+scale_fill_manual(values = c("#C33E3B","#4EA3DF","#6cBE58"))+
  geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM),width = 0.5, position = position_dodge(0.75))+
  ggtitle("TPQ Scores for the 13 ICA Components",subtitle = "With each Question Multiplied by its Loading")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                                                                                            #
# Now I will make the analysis without the weights ... exact same code, just different columns of inclusion  #
#                                                                                                            #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

variable_list = mapply(paste0,"IC_",0:12,"_MDD",SIMPLIFY = T) #Just a list to define the required columns
variable_names = mapply(paste0,"IC_",0:12,SIMPLIFY = T)
variable_list = append(c("Subject","Diagnosis","Response","Session"),variable_list)

NewData = TPQ_Q_log[variable_list]
colnames(NewData)[(ncol(NewData)-12):ncol(NewData)] = 0:12
NewData_Melt = melt(NewData,id.vars = c("Subject","Diagnosis","Response","Session"),value.name = "Score", variable.name = "IC_Number")  


DateMeans = SMeans(NewData_Melt,"Score",c("IC_Number","Diagnosis"),GroupBy = "Diagnosis")
DateMeans$IC_Number = factor(DateMeans$IC_Number, levels = 0:12)
ggplot(DateMeans,aes(x=IC_Number,y= Score, fill=Groups))+geom_bar(stat="identity",position = position_dodge(0.75),color="black")+
  TypicalTheme+scale_fill_manual(values = c("#C33E3B","#4EA3DF","#6cBE58"))+
  geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM),width = 0.5, position = position_dodge(0.75))+
  ggtitle("TPQ Scores for the 13 ICA Components",subtitle = "Simple Summation of Scores")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                                                                                            #
#                       Now I will make the same analysis, but for the MDD data only                         #
#                                                                                                            #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

variable_list = mapply(paste0,"Complete_IC_",0:12,"_HC",SIMPLIFY = T) #Just a list to define the required columns
variable_names = mapply(paste0,"IC_",0:12,SIMPLIFY = T)
variable_list = append(c("Subject","Diagnosis","Response","Session"),variable_list)

NewData = TPQ_Q_log_MDD[variable_list]
colnames(NewData)[(ncol(NewData)-12):ncol(NewData)] = 0:12
NewData_Melt = melt(NewData,id.vars = c("Subject","Diagnosis","Response","Session"),value.name = "Score", variable.name = "IC_Number")  

DateMeans = SMeans(NewData_Melt,"Score",c("IC_Number","Response"),GroupBy = "Response")
DateMeans$IC_Number = factor(DateMeans$IC_Number, levels = 0:12)
ggplot(DateMeans,aes(x=IC_Number,y= Score, fill=Groups))+geom_bar(stat="identity",position = position_dodge(0.75),color="black")+
  TypicalTheme+scale_fill_manual(values = c("#Eb4C42","#00A550"))+
  geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM),width = 0.5, position = position_dodge(0.75))+
  ggtitle("TPQ Scores for the 13 ICA Components in MDD data Only",subtitle = "With each Question Multiplied by its Weight")


####################################################################################
#                                   Experiment 1                                   #
####################################################################################

Experiment1 = read.csv("/home/abdul-rahman/Dropbox/TPQ_Paper/TPQ_Genetics.csv",header = TRUE)
Experiment1 = Experiment1[Experiment1$Age!=""&Experiment1$Age!="#N/A",]
Experiment1 = Experiment1[Experiment1$WAIS.R!=""&Experiment1$MMSE!="#N/A",]
Experiment1 = Experiment1[(Experiment1$WAIS.R!=""&Experiment1$MMSE!="#N/A"&!is.na(Experiment1$BAI)),]
Experiment1 = Experiment1[Experiment1$BDI.II<15,]
Experiment1 = Experiment1[Experiment1$BAI<22,]

Experiment1$Age = as.numeric(as.character(Experiment1$Age))
Experiment1$MMSE = as.numeric(as.character(Experiment1$MMSE))
Experiment1$WAIS.R = as.numeric(as.character(Experiment1$WAIS.R))
Experiment1$BDI.II = as.numeric(as.character(Experiment1$BDI.II))
Experiment1$BAI = as.numeric(as.character(Experiment1$BAI))
Experiment1$Years.of.Education = as.numeric(as.character(Experiment1$Years.of.Education))

MaCo=Experiment1$Gender=="Male"
FeCo=Experiment1$Gender=="Female"
DataNS1=data.frame(Experiment1$X2, Experiment1$X4, Experiment1$X9, Experiment1$X11, Experiment1$X40, Experiment1$X43, Experiment1$X85, Experiment1$X93, Experiment1$X96)
DataNS1M=DataNS1[MaCo,]
DataNS1F=DataNS1[FeCo,]
DataNS2=data.frame(Experiment1$X30, Experiment1$X46, Experiment1$X48, Experiment1$X50, Experiment1$X55, Experiment1$X56, Experiment1$X81, Experiment1$X99)
DataNS2M=DataNS2[MaCo,]
DataNS2F=DataNS2[FeCo,]
DataNS3=data.frame(Experiment1$X32, Experiment1$X66, Experiment1$X70, Experiment1$X72, Experiment1$X76, Experiment1$X78, Experiment1$X87)
DataNS3M=DataNS3[MaCo,]
DataNS3F=DataNS3[FeCo,]
DataNS4=data.frame(Experiment1$X13, Experiment1$X16, Experiment1$X21, Experiment1$X22, Experiment1$X24, Experiment1$X28, Experiment1$X35, Experiment1$X60, Experiment1$X62, Experiment1$X65)
DataNS4M=DataNS4[MaCo,]
DataNS4F=DataNS4[FeCo,]
DataHA1=data.frame(Experiment1$X1, Experiment1$X5, Experiment1$X8, Experiment1$X10, Experiment1$X14, Experiment1$X82, Experiment1$X84, Experiment1$X91, Experiment1$X95, Experiment1$X98)
DataHA1M=DataHA1[MaCo,]
DataHA1F=DataHA1[FeCo,]
DataHA2=data.frame(Experiment1$X18, Experiment1$X19, Experiment1$X23, Experiment1$X26, Experiment1$X29, Experiment1$X47, Experiment1$X51)
DataHA2M=DataHA2[MaCo,]
DataHA2F=DataHA2[FeCo,]
DataHA3=data.frame(Experiment1$X33, Experiment1$X37, Experiment1$X38, Experiment1$X42, Experiment1$X44, Experiment1$X89, Experiment1$X100)
DataHA3M=DataHA3[MaCo,]
DataHA3F=DataHA3[FeCo,]
DataHA4=data.frame(Experiment1$X49, Experiment1$X54, Experiment1$X57, Experiment1$X59, Experiment1$X63, Experiment1$X68, Experiment1$X69, Experiment1$X73, Experiment1$X75, Experiment1$X80)
DataHA4M=DataHA4[MaCo,]
DataHA4F=DataHA4[FeCo,]
DataRD1=data.frame(Experiment1$X27, Experiment1$X31, Experiment1$X34, Experiment1$X83, Experiment1$X94)
DataRD1M=DataRD1[MaCo,]
DataRD1F=DataRD1[FeCo,]
DataRD2=data.frame(Experiment1$X39, Experiment1$X41, Experiment1$X45, Experiment1$X52, Experiment1$X53, Experiment1$X77, Experiment1$X79, Experiment1$X92, Experiment1$X97)
DataRD2M=DataRD2[MaCo,]
DataRD2F=DataRD2[FeCo,]
DataRD3=data.frame(Experiment1$X3, Experiment1$X6, Experiment1$X7, Experiment1$X12, Experiment1$X15, Experiment1$X64, Experiment1$X67, Experiment1$X74, Experiment1$X86, Experiment1$X88, Experiment1$X90)
DataRD3M=DataRD3[MaCo,]
DataRD3F=DataRD3[FeCo,]
DataRD4=data.frame(Experiment1$X17, Experiment1$X20, Experiment1$X25, Experiment1$X36, Experiment1$X58)
DataRD4M=DataRD4[MaCo,]
DataRD4F=DataRD4[FeCo,]

DataNS=data.frame(Experiment1$X2, Experiment1$X4, Experiment1$X9, Experiment1$X11, Experiment1$X40, Experiment1$X43, Experiment1$X85, Experiment1$X93, Experiment1$X96, Experiment1$X30, Experiment1$X46, Experiment1$X48, Experiment1$X50, Experiment1$X55, Experiment1$X56, Experiment1$X81, Experiment1$X99, Experiment1$X32, Experiment1$X66, Experiment1$X70, Experiment1$X72, Experiment1$X76, Experiment1$X78, Experiment1$X87, Experiment1$X13, Experiment1$X16, Experiment1$X21, Experiment1$X22, Experiment1$X24, Experiment1$X28, Experiment1$X35, Experiment1$X60, Experiment1$X62, Experiment1$X65)
DataNSMale=DataNS[MaCo, ]
DataNSFemale=DataNS[FeCo, ]

DataHA=data.frame(Experiment1$X1, Experiment1$X5, Experiment1$X8, Experiment1$X10, Experiment1$X14, Experiment1$X82, Experiment1$X84, Experiment1$X91, Experiment1$X95, Experiment1$X98, Experiment1$X18, Experiment1$X19, Experiment1$X23, Experiment1$X26, Experiment1$X29, Experiment1$X47, Experiment1$X51, Experiment1$X33, Experiment1$X37, Experiment1$X38, Experiment1$X42, Experiment1$X44, Experiment1$X89, Experiment1$X100,Experiment1$X49, Experiment1$X54, Experiment1$X57, Experiment1$X59, Experiment1$X63, Experiment1$X68, Experiment1$X69, Experiment1$X73, Experiment1$X75, Experiment1$X80)

DataHAMale=DataHA[MaCo, ]
DataHAFemale=DataHA[FeCo, ]

DataRD=data.frame(Experiment1$X27, Experiment1$X31, Experiment1$X34, Experiment1$X83, Experiment1$X94, Experiment1$X39, Experiment1$X41, Experiment1$X45, Experiment1$X52, Experiment1$X53, Experiment1$X77, Experiment1$X79, Experiment1$X92, Experiment1$X97, Experiment1$X3,  Experiment1$X6,  Experiment1$X7,  Experiment1$X12, Experiment1$X15, Experiment1$X64, Experiment1$X67, Experiment1$X74, Experiment1$X86, Experiment1$X88, Experiment1$X90, Experiment1$X17, Experiment1$X20, Experiment1$X25, Experiment1$X36, Experiment1$X58)
DataRDMale=DataRD[MaCo, ]
DataRDFemale=DataRD[FeCo, ]

DataRDNew=data.frame(Experiment1$X27, Experiment1$X31, Experiment1$X34, Experiment1$X83, Experiment1$X94, Experiment1$X3, Experiment1$X6, Experiment1$X7, Experiment1$X12, Experiment1$X15, Experiment1$X64, Experiment1$X67, Experiment1$X74, Experiment1$X86, Experiment1$X88, Experiment1$X90, Experiment1$X17, Experiment1$X20, Experiment1$X25, Experiment1$X36, Experiment1$X58)
DataRDNewMale=DataRDNew[MaCo, ]
DataRDNewFemale=DataRDNew[FeCo, ]

SubjectCount.G = nrow(Experiment1)
AgeRange.G = paste0(range(Experiment1$Age)[1],"-",range(Experiment1$Age)[2])
FemaleCount.G = nrow(Experiment1[Experiment1$Sex=="Female",])
MaleCount.G = nrow(Experiment1[Experiment1$Sex=="Male",])

FemalePer.G = round(100*FemaleCount.G/SubjectCount.G,digits = 1)
MalePer.G = round(100*MaleCount.G/SubjectCount.G,digits = 1)

BDITest = t.test(BDI.II~Sex, data =Experiment1,var.equal=TRUE)
BAITest = t.test(BAI~Sex, data =Experiment1[Experiment1$SiteOfTesting!="Other",],var.equal=TRUE)
YOETest = t.test(Years.of.Education~Sex, data =Experiment1[Experiment1$SiteOfTesting!="Other",],var.equal=TRUE)

# Just to have a look at the differences
DataMeans=UMeans(Experiment1,"BAI","Sex","SiteOfTesting")
ggplot(DataMeans, aes(x=Sex,y=BAI,group=Sex, fill=Sex))+geom_bar(stat = "identity",color="black")+
  TypicalTheme+geom_errorbar(aes(ymin=BAI-SEM, ymax = BAI+SEM), width = 0.2)+
  facet_wrap(~SiteOfTesting)+scale_fill_manual(values = c("#Eb4C42","#0087BD"))

DataMeans=UMeans(Experiment1,"Years.of.Education","Sex","SiteOfTesting")
ggplot(DataMeans, aes(x=Sex,y=Years.of.Education,group=Sex, fill=Sex))+geom_bar(stat = "identity",color="black")+
  TypicalTheme+geom_errorbar(aes(ymin=Years.of.Education-SEM, ymax = Years.of.Education+SEM), width = 0.2)+
  facet_wrap(~SiteOfTesting)+scale_fill_manual(values = c("#Eb4C42","#0087BD"))

AllP = 0.2 #This is the number above which all p.values are present

AgeMean.M.G = round(mean(Experiment1$Age[Experiment1$Sex=="Male"]),digits = 1)
AgeMean.F.G = round(mean(Experiment1$Age[Experiment1$Sex=="Female"]),digits = 1)
AgeSD.M.G = round(sd(Experiment1$Age[Experiment1$Sex=="Male"]),digits = 1)
AgeSD.F.G = round(sd(Experiment1$Age[Experiment1$Sex=="Female"]),digits = 1)

EdMean.M.G = round(mean(Experiment1$Years.of.Education[Experiment1$Sex=="Male"]),digits = 1)
EdMean.F.G = round(mean(Experiment1$Years.of.Education[Experiment1$Sex=="Female"]),digits = 1)
EdSD.M.G = round(sd(Experiment1$Years.of.Education[Experiment1$Sex=="Male"]),digits = 1)
EdSD.F.G = round(sd(Experiment1$Years.of.Education[Experiment1$Sex=="Female"]),digits = 1)

MMSEMean.M.G = round(mean(Experiment1$MMSE[Experiment1$Sex=="Male"]),digits = 1)
MMSEMean.F.G = round(mean(Experiment1$MMSE[Experiment1$Sex=="Female"]),digits = 1)
MMSESD.M.G = round(sd(Experiment1$MMSE[Experiment1$Sex=="Male"]),digits = 1)
MMSESD.F.G = round(sd(Experiment1$MMSE[Experiment1$Sex=="Female"]),digits = 1)

WAISRMean.M.G = round(mean(Experiment1$WAIS.R[Experiment1$Sex=="Male"]),digits = 1)
WAISRMean.F.G = round(mean(Experiment1$WAIS.R[Experiment1$Sex=="Female"]),digits = 1)
WAISRSD.M.G = round(sd(Experiment1$WAIS.R[Experiment1$Sex=="Male"]),digits = 1)
WAISRSD.F.G = round(sd(Experiment1$WAIS.R[Experiment1$Sex=="Female"]),digits = 1)

BDIMean.M.G = round(mean(Experiment1$MMSE[Experiment1$Sex=="Male"]),digits = 1)
BDIMean.F.G = round(mean(Experiment1$MMSE[Experiment1$Sex=="Female"]),digits = 1)
BDISD.M.G = round(sd(Experiment1$MMSE[Experiment1$Sex=="Male"]),digits = 1)
BDISD.F.G = round(sd(Experiment1$MMSE[Experiment1$Sex=="Female"]),digits = 1)

BAIMean.M.G = round(mean(Experiment1$BAI[Experiment1$Sex=="Male"]),digits = 1)
BAIMean.F.G = round(mean(Experiment1$BAI[Experiment1$Sex=="Female"]),digits = 1)
BAISD.M.G = round(sd(Experiment1$BAI[Experiment1$Sex=="Male"]),digits = 1)
BAISD.F.G = round(sd(Experiment1$BAI[Experiment1$Sex=="Female"]),digits = 1)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                                   Experiment 2                                                #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
MDDCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"),])
HCCount = nrow(TPQData[(TPQData$Diagnosis=="HC"&TPQData$Session=="Test"),])
MDDMalesPer =round(100*nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Gender=="Male"),])/MDDCount,digits = 1)
MDDFemalesPer =round(100*nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Gender=="Female"),])/MDDCount,digits = 1)
HCMalesPer =round(100*nrow(TPQData[(TPQData$Diagnosis=="HC"&TPQData$Session=="Test"&TPQData$Gender=="Male"),])/HCCount,digits = 1)
HCFemalesPer =round(100*nrow(TPQData[(TPQData$Diagnosis=="HC"&TPQData$Session=="Test"&TPQData$Gender=="Female"),])/HCCount,digits = 1)

RetestPeriodRangeMDD = paste0(min(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="MDD"]),"-",max(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="MDD"]))
RetestPeriodRangeHC = paste0(min(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="HC"]),"-",max(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="HC"]))
RetestPeriodMeanMDD = round(mean(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
RetestPeriodMeanHC = round(mean(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)
RetestPeriodSDMDD = round(sd(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
RetestPeriodSDHC = round(sd(TPQData$RetestPeriod[TPQData$Session=="Test"&TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)

GADCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$GAD=="GAD"),])
GADPer =round(100*GADCount/MDDCount,digits = 1)

RespondersCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Responder"),])
RespondersPer =round(100*RespondersCount/MDDCount,digits = 1)

NonrespondersCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Non-responder"),])
NonrespondersPer =round(100*NonrespondersCount/MDDCount,digits = 1)

RespondersGADCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Responder"&TPQData$GAD=="GAD"),])
RespondersNoGADCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Responder"&TPQData$GAD=="NoGad"),])
NonrespondersGADCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Non-responder"&TPQData$GAD=="GAD"),])
NonrespondersNoGADCount = nrow(TPQData[(TPQData$Diagnosis=="MDD"&TPQData$Session=="Test"&TPQData$Response=="Non-responder"&TPQData$GAD=="NoGad"),])

AgeRangeMDD = paste0(min(TPQData$Age[TPQData$Diagnosis=="MDD"]), "-", max(TPQData$Age[TPQData$Diagnosis=="MDD"]))
AgeRangeHC = paste0(min(TPQData$Age[TPQData$Diagnosis=="HC"]), "-", max(TPQData$Age[TPQData$Diagnosis=="HC"]))

MeanAgeMDD=round(mean(TPQData$Age[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
SDAgeMDD=round(sd(TPQData$Age[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
MeanAgeHC=round(mean(TPQData$Age[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)
SDAgeHC=round(sd(TPQData$Age[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)


MeanMMSEMDD=round(mean(TPQData$MMSE[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
SDMMSEMDD=round(sd(TPQData$MMSE[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
MeanMMSEHC=round(mean(TPQData$MMSE[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)
SDMMSEHC=round(sd(TPQData$MMSE[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)

MeanWAISRMDD=round(mean(TPQData$WAIS.R[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
SDWAISRMDD=round(sd(TPQData$WAIS.R[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
MeanWAISRHC=round(mean(TPQData$WAIS.R[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)
SDWAISRHC=round(sd(TPQData$WAIS.R[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)

MeanEdMDD=round(mean(TPQData$Years.of.Education[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
SDEdMDD=round(sd(TPQData$Years.of.Education[TPQData$Diagnosis=="MDD"], na.rm = TRUE),digits = 1)
MeanEdHC=round(mean(TPQData$Years.of.Education[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)
SDEdHC=round(sd(TPQData$Years.of.Education[TPQData$Diagnosis=="HC"], na.rm = TRUE),digits = 1)

TPQData$SSRI=as.character(TPQData$SSRI)
TPQData$SSRI[TPQData$SSRI=="Paroxetine" | TPQData$SSRI=="Seroxat"]="1"
TPQData$SSRI[TPQData$SSRI=="Fluoxetine"]="2"
TPQData$SSRI[TPQData$SSRI=="Escitalopram" | TPQData$SSRI=="Cipralex"]="3"
TPQData$SSRI[TPQData$SSRI=="Sertraline"|TPQData$SSRI=="Citalopram"]="4"

#Drugs Descriptive
DrugCount =  sum(TPQData$SSRI=="4" | TPQData$SSRI=="3" | TPQData$SSRI=="2" | TPQData$SSRI=="1")
ParoxetinePer = round(100*sum(TPQData$SSRI=="1")/DrugCount,digits = 1)
FluoxetinePer = round(100*sum(TPQData$SSRI=="2")/DrugCount,digits = 1)
EscitalopramPer = round(100*sum(TPQData$SSRI=="3")/DrugCount,digits = 1)
SertralinePer = round(100*sum(TPQData$SSRI=="4")/DrugCount,digits = 1)


ParoxetineDoseMean = round(mean(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="1"])), na.rm = TRUE),digits = 1)
ParoxetineDoseSD = round(sd(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="1"])), na.rm = TRUE),digits = 1)
FluoxetineDoseMean = round(mean(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="2"])), na.rm = TRUE),digits = 1)
FluoxetineDoseSD = round(sd(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="2"])), na.rm = TRUE),digits = 1)
EscitalopramDoseMean = round(mean(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="3"])), na.rm = TRUE),digits = 1)
EscitalopramDoseSD = round(sd(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="3"])), na.rm = TRUE),digits = 1)
SertralineDoseMean = round(mean(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="4"])), na.rm = TRUE),digits = 1)
SertralineDoseSD = round(sd(as.numeric(as.character(TPQData$SSRIDose.mg.day[TPQData$SSRI=="4"])), na.rm = TRUE),digits = 1)
DrugDose=paste0(min(TPQData$Dose.SSRI.mg.day[complete.cases(TPQData$Dose.SSRI.mg.day)]),"-",max(TPQData$Dose.SSRI.mg.day[complete.cases(TPQData$Dose.SSRI.mg.day)]))



variables = c("MDDCount", "HCCount", "MDDMalesPer", "MDDFemalesPer", "HCMalesPer", "HCFemalesPer", "RetestPeriodRangeMDD", "RetestPeriodRangeHC", "RetestPeriodMeanMDD",
              "RetestPeriodMeanHC", "RetestPeriodSDMDD", "RetestPeriodSDHC", "GADCount", "GADPer", "RespondersCount", "RespondersPer", "NonrespondersCount",
              "NonrespondersPer", "RespondersGADCount", "RespondersNoGADCount", "NonrespondersGADCount", "NonrespondersNoGADCount", 
              "MeanMMSEMDD", "SDMMSEMDD", "MeanMMSEHC", "SDMMSEHC", "MeanWAISRMDD", "SDWAISRMDD", "MeanWAISRHC", "SDWAISRHC", "AgeRangeMDD", "AgeRangeHC",
              "MeanAgeMDD", "SDAgeMDD", "MeanAgeHC", "SDAgeHC", "MeanEdMDD", "SDEdMDD", "MeanEdHC", "SDEdHC", "DrugCount", "ParoxetinePer", "FluoxetinePer", "EscitalopramPer",
              "SertralinePer", "ParoxetineDoseMean", "ParoxetineDoseSD", "FluoxetineDoseMean", "FluoxetineDoseSD", "EscitalopramDoseMean", "EscitalopramDoseSD",
              "SertralineDoseMean", "SertralineDoseSD")

values = c(MDDCount, HCCount, MDDMalesPer, MDDFemalesPer, HCMalesPer, HCFemalesPer, RetestPeriodRangeMDD, RetestPeriodRangeHC, RetestPeriodMeanMDD,
           RetestPeriodMeanHC, RetestPeriodSDMDD, RetestPeriodSDHC, GADCount, GADPer, RespondersCount, RespondersPer, NonrespondersCount, 
           NonrespondersPer, RespondersGADCount, RespondersNoGADCount, NonrespondersGADCount, NonrespondersNoGADCount,
           MeanMMSEMDD, SDMMSEMDD, MeanMMSEHC, SDMMSEHC, MeanWAISRMDD, SDWAISRMDD, MeanWAISRHC, SDWAISRHC, AgeRangeMDD, AgeRangeHC,
           MeanAgeMDD, SDAgeMDD, MeanAgeHC, SDAgeHC, MeanEdMDD, SDEdMDD, MeanEdHC, SDEdHC, DrugCount, ParoxetinePer, FluoxetinePer, EscitalopramPer,
           SertralinePer, ParoxetineDoseMean, ParoxetineDoseSD, FluoxetineDoseMean, FluoxetineDoseSD, EscitalopramDoseMean, EscitalopramDoseSD,
           SertralineDoseMean, SertralineDoseSD)

write.table(data.frame(Variable=variables, Value = values),"/home/abdul-rahman/Dropbox/TPQ_Paper/Experiment2VArialbes.xls", sep="\t")
