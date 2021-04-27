# Title     : TODO
# Objective : TODO
# Created by: asawalma
# Created on: 12/7/20

####################################################################################################
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

  for (i in 2:length(as.character(Model$formula[[3]]))){
    local_factor= as.character(Model$formula[[3]])[i]
    data_frame_essential[[local_factor]]=Data[[local_factor]]
  }

  #Create prediction matrix
  PredMatrix=table(data_frame_essential$predicted.outcome,data_frame_essential$actual.outcome)
  rownames(PredMatrix)=c(paste0(BaseLevel,"(Predicted)"),paste0(FirstLevel,"(Predicted)"))
  colnames(PredMatrix)=c(paste0(BaseLevel,"(Actual)"),paste0(FirstLevel,"(Actual)"))

  LogValues=data.frame(Beta, SE, Lower.CI, odds_ratio, Upper.CI, Zvalues, Pvalues)
  DerivedValues=data.frame(chisq=modelChi, df=chidf, p.value.chi=chisq.prob, r2.hl=R2.hl, r2.cs=R2.cs, r2.n=R2.n, sensitivity=Sensitivity, specificity=Specificity)
  DataList=list(data_frame = data_frame_essential,
                log.values=LogValues, derived.values=DerivedValues,
                prob.formula = paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,")",sep = ""),
                PredMatrix = PredMatrix)

  print (Drawing)

  return (DataList)
}
TypicalTheme=theme_bw(base_size = 16,base_family = "Times")+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                                                                  plot.subtitle = element_text(hjust = 0.5,face = "italic"))

##############################################################################################################
# Linear Discriminant Analysis
# (the multiclass version of logistic regression)
# Author: Abdulrahman Sawalma


# This method is more widely used when we have mroe than two response classes
# it is more stable, especially if the distribution of the predictors is normal
# Reading is from the book: https://faculty.marshall.usc.edu/gareth-james/ISL/ISLR%20Seventh%20Printing.pdf
# Also, check this for ROC curve: https://towardsdatascience.com/understanding-the-roc-and-auc-curves-a05b68550b69
# using Bayes theorim (see the book) we can compute the probability of Y=k | X=x
# using the prior probability (pi_k) and f(x). pi_k is the overall probability that a randomly
# chosen observation comes from kth class. f(x) is the probability of X=x|Y=y
# This prability (pi_k) and function (f(x)) will be used to estimate the
# posterior probabiltiy (probability for each observation)

# Based on the Bayes classifier, we assign an observation to the class for which
# the p_k(x) is highest in the discriminant function
# LDA uses estimates for the mean for each class, prior probability for each class
# and sd for all classes to approximate the discriminant function. All what remains
# is to plug in the x value into the equation 4.17
# The boundary beyond which a variable is considered to belong to one of two classes
# is given by the equation (mu_1+mu_2)/2 ... for multivariate version, things are different

# In the multivariate case, the equation is the same, but the sd and mean is replaced with
# a covariance matrix and vector of means, respectively.
# The quadratic Discriminant Analysis assumes that each of the classes has different
# covariance matrix
# Shortly, if the classes shar common variance, LDA is better. If this assumption is far
# off, QDA is better


# Let's start with the examples
library(ISLR)

# get to know the data
?Smarket
summary(Smarket)

#find the nrows, ncolumns of Smarket
dim(Smarket)

#pairwise correlations of all variables
cor(Smarket[,-9])

#produce a testing data
train = Smarket$Year<2005
Smarket.2005 = Smarket[!train,]
Direction.2005 = Smarket$Direction[!train]

#First, let's do logistic regression
glm.fits=glm(Direction ~ Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=Smarket,
             family=binomial, subset = train)
summary(glm.fits)

# Predict theprobabilities
glm.probs = predict(glm.fits, Smarket.2005, type = "response")
glm.probs[1:10]

#Notice that the given values are the probability of (Y=1 | X).
# We can know what does "1" mean by looking at the contrasts(Smarket$Direction)
# where we see that Up was given a value of 1, and Down of 0

# calculate accuracy of prediction
glm.pred=rep("Down",dim(Smarket.2005)[1])
glm.pred[glm.probs >.5]="Up"

# For inspection
table(glm.pred,Direction.2005)

print (paste0("Percentage of correct = ", 100 * mean(Direction.2005 == glm.pred),"%"))



# Now the LDA
library(MASS)
lda.fit = lda (Direction ~ Lag1+Lag2, data = Smarket, subset = train)
lda.fit

#######################################################################################
# Cross validation
#######################################################################################

library(ISLR)
library(boot)

set.seed(17)
glm.fit=glm(mpg~horsepower,data=Auto)

Data = read.csv("/home/asawalma/Dropbox/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020.csv")

Data = Data[(complete.cases(Data$RD1) & Data$Session == "Test" & complete.cases(Data$RD2) & complete.cases(Data$RD3) &
  complete.cases(Data$RD1) & (Data$Diagnosis == "HC" | Data$Diagnosis == "MDD")),]
Data$Diagnosis = factor(Data$Diagnosis)

#Save data to open in python for comparison
write.csv(Data, "/home/asawalma/git/data_analysis/cv_trial.csv", row.names = FALSE)
glm.fit=glm(Diagnosis~RD1+RD2+RD3+RD4, data = Data, family = binomial())

summary(glm.fit)
cv.glm(data = Data ,glm.fit ,K=10)$delta[2]

LogisticFunction(glm.fit)


#######################################################################################
# Cross validation - Multinomial Logistic regression
#######################################################################################
require(foreign)
library(nnet)
require(ggplot2)
library(dplyr)
library(ISLR)
library(boot)
library(data.table)
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
Multinomial_LR = function(Model){
  coefficients  = summary(Model)$coefficients
  standard_errors = summary(Model)$standard.errors
  z_scores = coefficients/standard_errors
  p_values = pnorm(abs(z_scores),lower.tail = FALSE) * 2


  # Rename the first variable to become the same as the grouping variable
  grouping.var = all.vars(summary(Model)$call)[1]

  DataList = list()
  for (lvl in rownames(coefficients)){

    # exp_int is a 3d array, you can access it using exp_int[row, column, matrix]
    exp_int = confint(multinom.fit)

    #get matrix number
    matrix_no = match(lvl,dimnames(exp_int)[3][[1]])
    NewData = data.frame(Group = rep(lvl, ncol(coefficients)),
                         Predictor = colnames(coefficients),
                         coef = round(as.numeric(coefficients[lvl,]),digits = 4),
                         std.err = round(as.numeric(standard_errors[lvl,]),digits = 4),
                         z = round(as.numeric(z_scores[lvl,]),digits = 4),
                         p.value = round(as.numeric(p_values[lvl,]),digits = 4),
                         "0.025" = exp_int[,1,matrix_no],
                         "0.975" = exp_int[,2,matrix_no]
    )
    colnames(NewData)[1] = grouping.var
    DataList = append(DataList, list(NewData))
  }
  return (DataList)
}

# read the data
TPQData = read.csv("/home/asawalma/Dropbox/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020.csv")

#remove NAs
#TPQData = TPQData[complete.cases(TPQData$Initial.ID),]

# only include those with test and retest condition
#Test_ids = TPQData$Initial.ID[complete.cases(TPQData$Initial.ID) & TPQData$Session=="Test"]
#Retest_ids = TPQData$Initial.ID[complete.cases(TPQData$Initial.ID) & TPQData$Session=="Retest"]
#inclusion_condition = sapply(TPQData$Initial.ID[complete.cases(TPQData$Initial.ID)],
#                           function(x) (x %in% Test_ids & x %in% Retest_ids))

#TPQData = TPQData[inclusion_condition,]


#Calculate Date of Retest
#TPQData$DateRetest=TPQData$TestingDay
#TPQData=Distribute(TPQData,"DateRetest","Session","Retest","Test")
#TPQData$RetestPeriod=as.numeric(TPQData$DateRetest)-as.numeric(TPQData$TestingDay)

#Distribute response and GAD
#TPQData$Response = sapply(TPQData$Initial.ID,function(x) TPQData$Response[TPQData$Initial.ID==x & TPQData$Session=="Retest"])
#TPQData$GAD = sapply(TPQData$Initial.ID,function(x) TPQData$GAD[TPQData$Initial.ID==x & TPQData$Session=="Test"])

#Remove those with incomplete TPQ, diagnosis, response or retest period
#TPQData=TPQData[(TPQData$Session=="Retest"|(TPQData$RetestPeriod>27 & TPQData$RetestPeriod<76)),]
#TPQData=TPQData[complete.cases(TPQData$GAD=="NoGad"|TPQData$GAD=="GAD"|TPQData$Session=="Retest"|TPQData$Diagnosis=="HC"),]
#TPQData=TPQData[(TPQData$Response=="Non-responder"|TPQData$Response=="Responder"|TPQData$Response=="NoResponse"|TPQData$Session=="Retest"|TPQData$Diagnosis=="HC"),]

#Remove those who were on medications
#TPQData=TPQData[((as.character(TPQData$SSRI)=="None" & TPQData$Session=="Test" & TPQData$Diagnosis=="MDD") | TPQData$Session=="Retest" | TPQData$Diagnosis=="HC"),]

#Remove empty TPQs (sometimes they appear as all zeros)
TPQData$NS[(TPQData$NS+TPQData$HA+TPQData$RD)==0]=NA
TPQData=TPQData[(complete.cases(TPQData$NS)),]
#TPQData=MatchDelete(TPQData,"Session")

#####################################################################
#TPQData = TPQData[(complete.cases(TPQData$HA1) & TPQData$Session == "Test" & complete.cases(TPQData$HA2) & complete.cases(TPQData$HA3) &
#  complete.cases(TPQData$HA4) & (TPQData$Diagnosis == "HC" | TPQData$Diagnosis == "MDD" | TPQData$Diagnosis == "PTSD")),]
TPQData = TPQData[(complete.cases(TPQData$HA1) & TPQData$Session == "Test" & complete.cases(TPQData$HA2) & complete.cases(TPQData$HA3) &
  complete.cases(TPQData$HA4) & complete.cases(TPQData$Diagnosis)),]
TPQData$Diagnosis = factor(TPQData$Diagnosis)

included_vars = c("Diagnosis","HA1","HA2","HA3","HA4","NS1","NS2","NS3","NS4","RD1","RD2","RD3","RD4")
train = sample_frac(TPQData[,included_vars], 0.7)
sample_id = as.numeric(rownames(TPQData))
test = TPQData[-sample_id,included_vars]

#Save data to open in python for comparison
write.csv(train, "/home/asawalma/git/data_analysis/train.csv", row.names = FALSE)
write.csv(test, "/home/asawalma/git/data_analysis/test.csv", row.names = FALSE)
write.csv(TPQData, "/home/asawalma/git/data_analysis/all.csv", row.names = FALSE)


multinom.fit = multinom(Diagnosis~HA1+HA2+HA3+HA4+NS1+NS2+NS3+NS4+RD1+RD2+RD3+RD4, data = train)
summary(multinom.fit)
model_summary = Multinomial_LR(multinom.fit)
model_summary
y_pred = data.frame(multinom.fit$fitted.values)
y_pred$prediction = sapply(transpose(y_pred),function (x) match(max(x),x))
y_pred$prediction = sapply(y_pred$prediction, function(x) colnames(y_pred)[x])
y_pred$actual = train$Diagnosis


table(y_pred$actual,y_pred$prediction)

Accuracy = sum(y_pred$actual == y_pred$prediction)/length(y_pred$actual)
print (paste(round(100*Accuracy,1),"%"))


library(pROC)
pred = y_pred$MDD[(y_pred$actual=="GAD" | y_pred$actual=="HC")]
actual = factor(y_pred$actual[(y_pred$actual=="GAD" | y_pred$actual=="HC")])
roc(actual ~ pred, plot = TRUE, print.auc = TRUE)