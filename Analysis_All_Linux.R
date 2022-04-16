#The Ultimate File for Analysis
library(ggplot2)
library(extrafont)
library(reshape2)
library(car)
library(ppcor)
library(mlogit)
library(stringr)
library(dplyr)
library(ggthemes)
##############################################################################################################
#################################                                             ################################
#################################                   Functions                 ################################
#################################                                             ################################
##############################################################################################################
windowsFonts(Times=windowsFont("Amiri"))
TypicalTheme=theme_bw(base_size = 16,base_family = "Amiri")+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                                                                  plot.subtitle = element_text(hjust = 0.5,face = "italic"))

MinimalTheme=theme_minimal(base_size = 18,base_family = "Amiri")+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                                                                       plot.subtitle = element_text(hjust = 0.5,face = "italic"))

CPalette=c("#00A550","#Eb4C42","#0087BD",NA,
           "#C33E3B","#4EA3DF","#6cBE58","#808CA3","#B9B0AB",NA,
           "#B768A2","#FFD800","#E25822",NA,
           "#FF3800","#1B4D3E","#003153",NA,
           "#008000","#CC0000","#08457E",NA,
           "#682860","#FBEC5D","#FF6347",NA,NA,
           "#EFCC00","#A2A2D0","#FF0800","#D40000","#E62020",NA,
           "#E86100","#FF7F50",NA,"#FF355E","#FC5A8D","#E30B5D", NA,
           "#2F4F4F", "#CC6666", "#9999CC", "#66CC99")
FFF=as.factor(c(1:44))
ggplot(mapping=aes(x=FFF,y=rep(5,44),group=FFF,fill=FFF))+geom_bar(stat="identity",position = position_dodge())+scale_fill_manual(values = CPalette)+TypicalTheme+
  theme(legend.position = "none")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, Listing=FALSE) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  #I added the Listing variable so that if I provide a list in the function it can plot it (like the one used in "Correlate" function)
  if (isTRUE(Listing)){
    plots <- c(..., plotlist)
  }
  
  numPlots = length(plots)
  
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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
Attach.V <- function(data.frame1,data.frame2){
  #I will create a data frame with a length equal to the total length of the two data frames and put the Group variables in it
  Length1=nrow(data.frame1)
  Length2=nrow(data.frame2)
  FinalData=data.frame(G.Group=rep(NA,Length1+Length2))
  
  for (name in colnames(data.frame1)){
    if (name %in% colnames(data.frame2)){
      FinalData[[name]]=c(as.character(data.frame1[[name]]),as.character(data.frame2[[name]]))
      FinalData[[name]]=ReturnClass(data.frame1[[name]],FinalData[[name]])
    }
  }
  FinalData$G.Group=NULL
  return(FinalData)
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
Correlate <- function(Data,Var1,Var2,name1=NA_character_,name2=NA_character_,Method="lm",CorMethod=NA_character_,Factor=NA_character_,
                      Colors=NA_character_,Density=TRUE,FactorRelevel=NA_character_){
  #First we define the names as the variables' characters if they are NA_Character_s
  #if the factor is still Na, then we do the normal way of correlation
  #elseif there is a factor for segregating the factor, we make plots as many as the factor levels.
  #The following two lines to determine the levels of the factor, if the factor is not NA, then Levels=levels of the factor
  #otherwise, it is a one character vector
  if (is.na(Colors)){
    Colors=c("#006A4E","#FF4040","#1974D2","#E9D66B","#734F96","#8C92AC","#7B3F00")
  }
  
  Levels=Var1
  if (!is.na(Factor)){
    if (!is.na(FactorRelevel)){
      Levels=FactorRelevel
    }else{
      Levels=levels(as.factor(Data[[Factor]]))
    }
  }
  if (is.na(name1)){
    name1=Var1
  }
  if (is.na(name2)){
    name2=Var2
  }
  Text1=NULL; Text2=NULL; Text3=NULL
  
  FinalPlot=list()
  i=0
  for (item in Levels){
    if (is.na(Factor)){
      Variable1=Data[[Var1]]
      Variable2=Data[[Var2]]
    }else {
      Variable1=Data[[Var1]][Data[[Factor]]==item]
      Variable2=Data[[Var2]][Data[[Factor]]==item]
    }
    
    Shapiro1=shapiro.test(Variable1)
    Shapiro2=shapiro.test(Variable2)
    
    
    P1=round(Shapiro1$p.value,digits = 5); P2=round(Shapiro2$p.value,digits = 5)
    Est1 = round(Shapiro1$statistic,digits = 3); Est2 = round(Shapiro2$statistic,digits = 3)
    EstType=ifelse((P1>0.05 & P2>0.05),"pearson","spearman")
    if (!is.na(CorMethod)) {
      EstType=CorMethod
    }
    EstName=ifelse(EstType=="pearson","Pearson's r","Spearman's rho")
    if (EstType == "kendall"){
      EstName="Kendall's tau"
    }
    
    Length=length(Variable1)
    Text1=c(Text1,rep("",Length-1),paste("Shapiro-Wilk's W=",Est1,"\n p-value = ", P1,"\n\n\n\n"))# The many \n characters to make the text appear on a well-determined place no matter what is the height of the y-axis
    Text2=c(Text2,rep("",Length-1),paste("Shapiro-Wilk's W=",Est2,"\n p-value = ", P2,"\n\n\n\n"))
    
    Text1Position.x=(max(Data[[Var1]][complete.cases(Data[[Var1]])])+min(Data[[Var1]][complete.cases(Data[[Var1]])]))/2
    Text2Position.x=(max(Data[[Var2]][complete.cases(Data[[Var2]])])+min(Data[[Var2]][complete.cases(Data[[Var2]])]))/2
    
    Correlation=cor.test(Variable1,Variable2,method = EstType)
    print ("p-value is"); print (Correlation$p.value)
    CP=round(Correlation$p.value,digits = 3)
    CE=round(Correlation$estimate,digits = 3)
    Text3=c(Text3,rep("",Length-1),paste(EstName, "=",CE,"\n","p.value =",CP))
    Text3Position.x=(max(Data[[Var1]][complete.cases(Data[[Var1]])])+min(Data[[Var1]][complete.cases(Data[[Var1]])]))/2
    Text3Position.y=max(Data[[Var2]][complete.cases(Data[[Var2]])])*0.9
    
    Ncol=ifelse(is.na(Factor),2,1)
    
  }
  if (is.na(Factor)){
    g1=ggplot(mapping=aes(x=Data[[Var1]],y=Data[[Var2]]))+geom_point(size=2,shape=21,color="black",fill=Colors[1],alpha=0.8)+
      theme_bw(base_size = 15,base_family = "Amiri")+ggtitle(paste(name1, "Vs", name2))+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none")+
      scale_x_continuous(name = name1)+scale_y_continuous(name=name2)+scale_fill_manual(values = Colors)+
      geom_smooth(mapping = aes(color=Colors[1]),method=Method, size=1)+scale_color_manual(values = Colors)+
      geom_text(mapping = aes(x=Text3Position.x,y=Text3Position.y,label=Text3),family="Amiri",size=5)
    
    g2=ggplot(mapping = aes(x=Data[[Var1]]))+geom_density(alpha=0.8,fill=Colors[1])+theme_bw(base_size = 13,base_family = "Amiri")+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none",strip.text.x = element_blank())+
      scale_x_continuous(name = name1)+scale_y_continuous(name=NULL)+ggtitle("Density Plots")+
      geom_text(mapping=aes(x=Text1Position.x,y=0,label=Text1),family="Amiri",size=4,inherit.aes = FALSE)
    
    
    g3=ggplot(mapping = aes(x=Data[[Var2]]))+geom_density(alpha=0.8,fill=Colors[1])+theme_bw(base_size = 13,base_family = "Amiri")+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none",strip.text.x = element_blank())+
      scale_x_continuous(name = name2)+scale_y_continuous(name=NULL)+
      geom_text(mapping=aes(x=Text2Position.x,y=0,label=Text2),family="Amiri",size=4.5,inherit.aes = FALSE)
    
    if (isTRUE(Density)){
      LayOut <- matrix(c(1,1,1,1,2,3), nrow = 2, byrow = FALSE)
      print (multiplot(g1,g2,g3,layout = LayOut,cols=2))
    }else{
      print (g1)
    }
    
  }else{
    DataFrame=data.frame(VVar1=Data[[Var1]],VVar2=Data[[Var2]],VFactor=Data[[Factor]])
    if (!is.na(FactorRelevel)){
      DataFrame$VFactor=factor(DataFrame$VFactor,levels=FactorRelevel)
    }
    g1=ggplot(DataFrame,mapping=aes(x=VVar1,y=VVar2,fill=VFactor))+geom_point(size=2,shape=21,color="black",alpha=0.8)+theme_bw(base_size = 15,base_family = "Amiri")+
      ggtitle(paste("Scatter Plot for",name1, "and", name2))+facet_wrap(~VFactor,ncol = length(Levels))+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none")+
      scale_x_continuous(name = name1)+scale_y_continuous(name=name2)+scale_fill_manual(values = Colors)+
      geom_smooth(mapping = aes(color=VFactor),method=Method, size=1)+scale_color_manual(values = Colors)+
      geom_text(mapping = aes(x=Text3Position.x,y=Text3Position.y,label=Text3),family="Amiri",size=5)
    
    g2=ggplot(DataFrame,mapping = aes(x=VVar1,fill=VFactor))+geom_density(alpha=0.8)+theme_bw(base_size = 13,base_family = "Amiri")+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none",strip.text.x = element_blank())+
      facet_wrap(~VFactor,ncol = length(Levels))+scale_x_continuous(name = NULL)+scale_y_continuous(name=name1)+scale_fill_manual(values = Colors)+
      geom_text(mapping=aes(x=Text1Position.x,y=0,label=Text1),family="Amiri",size=4.5,inherit.aes = FALSE)
    
    g3=ggplot(DataFrame,mapping = aes(x=VVar2,fill=VFactor))+geom_density(alpha=0.8)+theme_bw(base_size = 13,base_family = "Amiri")+
      theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none",strip.text.x = element_blank())+
      facet_wrap(~VFactor,ncol = length(Levels))+scale_x_continuous(name = NULL)+scale_y_continuous(name=name2)+scale_fill_manual(values = Colors)+
      geom_text(mapping=aes(x=Text2Position.x,y=0,label=Text2),family="Amiri",size=4.5,inherit.aes = FALSE)
    
    if (isTRUE(Density)){
      FinalPlot=c(list(g1),list(g2), list(g3))
      LayOut <- matrix(c(1,1,1,1,2,3), nrow = 6, byrow = TRUE)
      print (multiplot(FinalPlot,cols=length(Levels),layout = LayOut,Listing=TRUE))
    }else{
      print (g1)
    }
  }
}
BinCorrelate <- function(Data,VariableList){
  #This function compares the correlations between non-parametric variables
  ExtendedVariableList={}
  
  ExcludedNumbers={}#This is a list of numbers that we do not want to include in the correlation...remember that each item will be
  #correlated with the others twice, we do not want that to happen, so, whenever an Item is compared to the rest, it should be
  #removed from future comparisons
  
  for (i in 1:length(VariableList)){
    ExcludedNumbers=append(ExcludedNumbers,i)
    for (k in (1:length(VariableList))[!((1:length(VariableList)) %in% ExcludedNumbers)]){
      ExtendedVariableList=append(ExtendedVariableList,c(VariableList[i],VariableList[k]))
    }
  }
  MatrixList=matrix(ExtendedVariableList,ncol = 2,byrow = TRUE)
  
  PValueList={}
  DfList={}
  PearsonChiList={}
  FactorNameList={}
  for (Row in 1:nrow(MatrixList)){
    ChiSquare=chisq.test(Data[[MatrixList[Row,1]]],Data[[MatrixList[Row,2]]],correct = FALSE)
    FactorNameList=append(FactorNameList,paste(MatrixList[Row,1]," - ",MatrixList[Row,2]))
    PValueList=append(PValueList,as.character(round(ChiSquare$p.value,6)))
    DfList=append(DfList,as.character(round(ChiSquare$parameter,1)))
    PearsonChiList=append(PearsonChiList,as.character(round(ChiSquare$statistic,1)))
  }
  DataFrame=data.frame(Variables=FactorNameList,Chi.Square=PearsonChiList,DF=DfList,p.value=PValueList)
  return(DataFrame)
}
Multimelt<-function(dataframe,...,Names=NULL, new.factor="Factor", common.var="Subject"){
  Variables=c(list(...))
  MeltedData=data.frame()
  DeletedVariables=NA_character_
  i=0
  for (item in Variables){
    i=i+1
    if (length(Names)<i){Names[i]=paste("var",i,sep = "")}
    
    if (item[1] == Variables[[1]][1]){
      #for the first set of variables, do a simple melt
      FinalData=melt(data=dataframe, measure.vars = item, value.name = Names[i])
    }else{
      #for the rest, attach horizontally
      FinalData=Attach.H(FinalData,melt(data=dataframe, measure.vars = item, value.name = Names[i]),common.var,common.var)
    }
    
    DeletedVariables=c(DeletedVariables,item)
  }
  #Finally, delete the remaining variables
  FinalData=FinalData[,!(colnames(FinalData) %in% DeletedVariables)]
  FinalData[[new.factor]]=FinalData$variable
  FinalData$variable=NULL
  return(FinalData)
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
    print (paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,")",sep = ""))
  }
  
  
  if (isTRUE(PrintList)){
    print (OriginalData)
  }
  PredMatrix=table(OriginalData[,length(OriginalData)-1],OriginalData[,(length(OriginalData)-2)])
  #colnames(PredMatrix)=c("0(Predicted)","1(Predicted)"); rownames(PredMatrix)=c("0(Actual)","1(Actual)")
  print (paste("Sensitivity =",round(PredMatrix[4]/(PredMatrix[3]+PredMatrix[4]),digits = 3)*100,"%"))
  print (paste("Specificity =",round(PredMatrix[1]/(PredMatrix[1]+PredMatrix[2]),digits = 3)*100,"%"))
  print (paste("Neg.Pred.Value =",round(PredMatrix[1]/(PredMatrix[1]+PredMatrix[3]),digits = 3)*100,"%"))
  print (paste("Pos.Pred.Value =",round(PredMatrix[4]/(PredMatrix[2]+PredMatrix[4]),digits = 3)*100,"%"))
  return (PredMatrix)
  
}
LogisticInfoFinder <- function(Model,iterations=10,Stepwise=TRUE,DeleterMean=100){
  #The deleter function is to make the function delete random variables to check better solutions
  #it depends on the rnorm function (rnorm (number, mean=DeleterMean)), so if the Deletermean is high (100) all values will be above 0,
  # which will make all the items in the TheSample list included without deletion
  #We start with formular Seperating: seperate the formula according to the "+" signs it contains
  TheData=Model$data
  DV=as.character(Model$formula[2][[1]])
  DV2=paste0(c(DV,"~"),collapse = "")
  MainFormula=Model$formula
  MainFormula[3]=1
  
  nothing=glm(as.formula(paste0(DV2,"1")),data=TheData,family = binomial())
  ModelNew=glm(Model$formula,data=TheData,family = binomial())
  StepForm=Model$formula
  if (isTRUE(Stepwise)){
    backward=step(ModelNew)
    StepForm=formula(backward)
  }
  StepVar=NULL
  
  
  #forward=step(nothing,scope = list(lower=formula(nothing),upper=formula(ModelNew)),direction = "forward")
  #bothways=step(nothing,scope = list(lower=formula(nothing),upper=formula(Model)),direction = "both",trace = 0)
  
  
  #ForForm=formula(forward); ForVar=NULL
  #BothForm=formula(both); BothVar=NULL
  
  #FormulaCharacter is the formula in string format (the second and third elemts of the formula pasted together)
  StepFormulaCharacter=strsplit(paste0(as.character(StepForm[3][[1]][2])," + ",as.character(StepForm[3][[1]][3])),split = "")
  StepFormulaCharacter=append(c("+",""),append(StepFormulaCharacter[[1]],c("","+")))
  StepCharacterNumbers=which(StepFormulaCharacter=="+")
  
  #now we will make the sensitivity variable (a variable to include the sequence of sensitivities produced by the loop below)
  #and Specificity variable, and a list of sampled randomized variables
  #The SensSpecVariable is the sum of sensitivities and specificities. I made it to see which value is greater, this means higher spec and sens
  SensVar=NULL; SpecVar= NULL; SensSpecVariable=NULL; SampleList=list(); PosPredValueVar=NULL; NegPredValueVar=NULL
  
  for (i in 2:length(StepCharacterNumbers)){
    StepVar=append(StepVar,paste0(StepFormulaCharacter[(StepCharacterNumbers[(i-1)]+2):(StepCharacterNumbers[i]-2)],collapse = ""))
  }
  for (i in 1:iterations){
    FinalFormula=NULL
    TheSample= sample(1:length(StepVar))
    TheSample=TheSample[rnorm(length(TheSample),mean = DeleterMean)>0]
    
    for (item in TheSample){
      FinalFormula=paste0(c(FinalFormula,StepVar[item]),collapse = "+")
    }
    FinalFormula=paste0(DV2,FinalFormula,collapse = "")
    FinalFormula=as.formula(FinalFormula)
    Matrix=LogisticInfo(glm(FinalFormula,data = TheData,family = binomial()),PrintAny = FALSE)
    Sens=round(Matrix[4]/(Matrix[3]+Matrix[4]),digits = 3)
    Spec=round(Matrix[1]/(Matrix[1]+Matrix[2]),digits = 3)
    NegPredValue=round(Matrix[1]/(Matrix[3]+Matrix[1]),digits = 3)
    PosPredValue=round(Matrix[4]/(Matrix[2]+Matrix[4]),digits = 3)
    #Now we append the values of spec, sens, and samples to their appropriate variables
    SensVar=append(SensVar,Sens); SpecVar=append(SpecVar,Spec); SampleList=append(SampleList,list(TheSample))
    NegPredValueVar=append(NegPredValueVar,NegPredValue); PosPredValueVar=append(PosPredValueVar,PosPredValue)
  }
  #I will put more wieght on the spec and sens
  SensSpecVariable=2*SensVar+2*SpecVar+NegPredValueVar+PosPredValueVar
  
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  #I will reproduce the list with the higher sens and spec
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  FinalFormula=NULL
  BestSampleNumber=which(SensSpecVariable==max(SensSpecVariable[complete.cases(SensSpecVariable)]))[1]
  TheSample=SampleList[[BestSampleNumber]]
  for (item in TheSample){
    FinalFormula=paste0(c(FinalFormula,StepVar[item]),collapse = "+")
  }
  FinalFormula=paste0(DV2,FinalFormula,collapse = "")
  FinalFormula=as.formula(FinalFormula)
  print ("#####################################################################################################################################")
  print ("The best sample in terms of specificity and sensitivity is:")
  print (TheSample)
  LogisticInfo(glm(FinalFormula,data = TheData,family = binomial()),PrintList = TRUE)
  
}

LogisticFunction = function(Data, DV, IVs,control = "HC" , Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE){
  
  #The Model is a logistic model (DV~IV1+IV2..., family = binomial())
  Data = Data[colnames(Data) %in% DV| colnames(Data) %in% IVs]
  Data = na.omit(Data)
  
  # create a training and a test data set
  train_inds = sample(1:nrow(Data), size = perc*nrow(Data))
  train = Data[train_inds,]
  test = Data[-train_inds,]
  
  Formula = as.formula(paste0(DV," ~ ",paste0(IVs,collapse = "+")))
  Model = glm(Formula, data = train, family = binomial())
  Factor = as.character(Formula[[2]])
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
  
  PredictedVar = factor(test[[DV]])
  PredictedVar = relevel(PredictedVar, ref = control)
  
  
  BaseLevel=levels(PredictedVar)[1]
  FirstLevel=levels(PredictedVar)[2]
  
  
  test$predicted.probability = predict(Model,newdata=test,type='response')
  test$predicted.outcome <-ifelse(test$predicted.probability<Threshold, BaseLevel, FirstLevel)
  
  TruePositive=nrow(test[(PredictedVar==FirstLevel&test$predicted.outcome==FirstLevel), ])
  FalsePositive=nrow(test[(PredictedVar==BaseLevel&test$predicted.outcome==FirstLevel), ])
  TrueNegative=nrow(test[(PredictedVar==BaseLevel&test$predicted.outcome==BaseLevel), ])
  FalseNegative=nrow(test[(PredictedVar==FirstLevel&test$predicted.outcome==BaseLevel), ])
  
  Sensitivity=round(TruePositive/(TruePositive+FalseNegative), digits = 3)
  Specificity=round(TrueNegative/(TrueNegative+FalsePositive), digits = 3)
  PPV=round(TruePositive/(TruePositive+FalsePositive), digits=3)
  NPV=round(TrueNegative/(TrueNegative+FalseNegative), digits=3)
  
  
  Xs=round(Coefficients[1,1],digits = 4)
  for (i in 2:length(Coefficients[,1])){
    VAr_name = names(Coefficients[,1])[i]
    L_coefficient = paste0("+",round(Coefficients[i,1],digits = 4))
    factor_value = L_coefficient
    factor_value = gsub("\\+-"," - ",L_coefficient)
    factor_value = gsub("\\+"," \\+ ",factor_value)
    Xs = paste0(Xs,factor_value,"*",VAr_name)
  }
  #The y.value finder, which is the maximum count on the plot divided by 2, this is for plotting reasons
  TList=NULL
  for (i in 0:20){
    interval=c((i-0.5)*0.05, (i+0.5)*0.05)
    Num=sum((predicted.probability>interval[1]&predicted.probability<interval[2]))
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
  
  Title = paste("Logistic Regression Function for", as.character(Formula[[2]]))
  Subtitle= paste("As Predicted by", paste0(IVs,collapse = "+"))
  if (plt_type == "histogram"){
    Drawing=ggplot(test, aes(x=predicted.probability, fill=PredictedVar))+geom_histogram(binwidth=0.05, color="black")+
      scale_fill_manual(name="Group", values=c("#08457E", "#FBEC5D"))+TypicalTheme+geom_vline(xintercept = Threshold)+
      scale_x_continuous("Predicted Probability",limits = c(-0.1, 1.1), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+
      annotate(geom = "text", size=5, family="Amiri", x = (Threshold/2), y = y.value*2.1, label = paste("Predicted to be", BaseLevel))+
      annotate(geom = "text", size=5, family="Amiri", x = ((1+Threshold)/2), y = y.value*2.1, label = paste("Predicted to be", FirstLevel))+
      annotate(geom = "text", size=5, family="Amiri", x = -0.05, y = y.value, label = GText_0, angle = 90)+
      annotate(geom = "text", size=5, family="Amiri", x = 1, y = y.value, label = GText_1, angle = 90)+
      ggtitle(Title, subtitle = Subtitle)+scale_y_continuous("Number of Cases")
  } else if (plt_type == "glm"){
    if (length(IVs)==1){
      test$Values = test$predicted.probability
      PredictedLevels = levels(test[[DV]])
      test$PredictedFactor = as.character(test[[DV]])
      test$PredictedFactor[test$PredictedFactor==PredictedLevels[1]]=0
      test$PredictedFactor[test$PredictedFactor==PredictedLevels[2]]=1
      test$PredictedFactor = as.numeric(test$PredictedFactor)
      Drawing =ggplot(data = test,mapping= aes(x = Values, y=PredictedFactor))+
        geom_jitter(shape = 21, color = "#7c0a02", fill="white",size = 4,width = 0.1, height = 0)+
        stat_smooth(method = "glm", method.args = list(family = "binomial"),se=F,color="#7c0a02",size=4)+
        scale_y_continuous(name = DV,breaks=c(0,1), labels= c(levels(test[[DV]])[1],levels(test[[DV]])[2]))+
        scale_x_continuous(name = paste0(IVs,collapse = "+"))+
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
  data_frame_essential = data.frame(predicted.probability= test$predicted.probability,
                                    predicted.outcome = test$predicted.outcome,
                                    actual.outcome = test[[DV]])
  
  for (i in 2:length(as.character(Model$formula[[3]]))){
    local_factor= as.character(Model$formula[[3]])[i]
    data_frame_essential[[local_factor]]=test[[local_factor]]
  }
  
  #Create prediction matrix
  PredMatrix=table(data_frame_essential$predicted.outcome,data_frame_essential$actual.outcome)
  #in case all cases are predicted to be 0 or 1, this will make PredMatrixa 1-row matrix. So, I will go through items one by one
  for (i in 1:length(rownames(PredMatrix))){
    rownames(PredMatrix)[i]=gsub(rownames(PredMatrix)[i],paste0(rownames(PredMatrix)[i],"(Predicted)"),rownames(PredMatrix)[i])
  }
  colnames(PredMatrix)=c("0(Actual)","1(Actual)")
  
  roc.info = "No ROC"
  if (nlevels(test[[DV]]) == 2){
    cases_probs = test$predicted.probability
    test_values = factor(test[[DV]])
    test_values = relevel(test_values, control)
    roc.info = roc(test_values, cases_probs, plot = plot.ROC, levels = levels(test_values),
                   percent = TRUE,legacy.axes = TRUE, print.auc = TRUE)
  }
  LogValues=data.frame(Beta, SE, Lower.CI, odds_ratio, Upper.CI, Zvalues, Pvalues)
  DerivedValues=data.frame(chisq=modelChi, df=chidf, p.value.chi=chisq.prob, r2.hl=R2.hl, r2.cs=R2.cs, r2.n=R2.n, sensitivity=Sensitivity, specificity=Specificity)
  DataList=list(data_frame = data_frame_essential,
                log.values=LogValues, derived.values=DerivedValues,
                prob.formula = paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,")",sep = ""),
                PredMatrix = PredMatrix,
                ROC = roc.info)
  
  print (Drawing)
  
  return (DataList)
  
}




Classes = c("Session"="factor","Subject"="factor","Original.ID"="factor","Age"="numeric","Sex"="factor","Years.of.Education"="numeric",
            "MMSE"="numeric","WAIS.R"="numeric",BDI.II="numeric","BAI"="numeric","Exclusion.Reason"="character","Retest.Period"="numeric",
            "Response"="factor","Include"="numeric","Comorbid.GAD"="character")

Original.Data=read.csv("/media/abdul-rahman/OS/Users/jsawa/Desktop/PNI/Analysis/Analysis.InfoSheet.csv",colClasses = Classes)
Quarters = read.csv("/home/abdul-rahman/Desktop/Quarters_HC_Ready_Old_1.csv")
NewQuarters = Attach.H(Original.Data,Quarters,SubjectID1 = "Original.ID", SubjectID2 = "Original.ID")

NewQuarters = NewQuarters[NewQuarters$Session=="Test",]
cor.test((NewQuarters$B1RewCor+NewQuarters$B2PunCor+NewQuarters$B3PunCor+NewQuarters$B4PunCor),NewQuarters$BAI)
cor.test((NewQuarters$B1RewCor+NewQuarters$B2PunCor+NewQuarters$B3PunCor+NewQuarters$B4PunCor),NewQuarters$BAI)
cor.test((NewQuarters$B1RewCor+NewQuarters$B2RewCor+NewQuarters$B3RewCor+NewQuarters$B4RewCor),NewQuarters$BDI.II)


PTSDData = read.csv("/home/abdul-rahman/Desktop/Fish_PTSD_Ready.csv")
TNPData = read.csv("/home/abdul-rahman/Desktop/Fish_TNP_Ready.csv")
HCData = read.csv("/home/abdul-rahman/Desktop/Fish_HC_Ready.csv")
PTSDData$Diagnosis = "PTSD"
TNPData$Diagnosis = "TNP"
HCData$Diagnosis = "HC"
AllData = Attach.V(PTSDData,Attach.V(TNPData,HCData))
AllData$Learning = AllData$Face1ErrorPh0+AllData$Face3ErrorPh0+AllData$Face1ErrorPh1+AllData$Face2ErrorPh1+AllData$Face3ErrorPh1+
  AllData$Face4ErrorPh1+AllData$Face1ErrorPh2+AllData$Face2ErrorPh2+AllData$Face3ErrorPh2+AllData$Face4ErrorPh2+AllData$Face1ErrorPh2Prime+AllData$Face3ErrorPh2Prime
AllData$Retention = AllData$Face1ErrorPh3+AllData$Face2ErrorPh3+AllData$Face3ErrorPh3+AllData$Face4ErrorPh3+AllData$Face1ErrorPh3Prime+AllData$Face3ErrorPh3Prime
AllData$Generalization = AllData$Face2ErrorPh3Prime+AllData$Face4ErrorPh3Prime

DataMeans_1=UMeans(AllData,"Learning","Diagnosis","Rule")
DataMeans_1$SEM_Learning = DataMeans_1$SEM
DataMeans_1=DataMeans_1[-2]
DataMeans_2=UMeans(AllData,"Retention","Diagnosis","Rule")
DataMeans_2$SEM_Retention = DataMeans_2$SEM
DataMeans_2=DataMeans_2[-2]
DataMeans_3=UMeans(AllData,"Generalization","Diagnosis","Rule")
DataMeans_3$SEM_Generalization = DataMeans_3$SEM
DataMeans_3=DataMeans_3[-2]

DataMeans = DataMeans_1
DataMeans$Retention = DataMeans_2$Retention
DataMeans$SEM_Retention = DataMeans_2$SEM_Retention
DataMeans$Generalization = DataMeans_3$Generalization
DataMeans$SEM_Generalization = DataMeans_3$SEM_Generalization


MeltedData = melt(data=DataMeans, measure.vars = c("Generalization","Retention","Learning"), value.name = "Errors",variable.name = "MeasureType")
MeltedData$SEM=ifelse(MeltedData$MeasureType=="Learning",MeltedData$SEM_Learning,ifelse(MeltedData$MeasureType=="Retention",MeltedData$SEM_Retention,MeltedData$SEM_Generalization))

ggplot(MeltedData,aes(x=Diagnosis,y=Errors,group=Diagnosis, fill=Diagnosis))+geom_bar(stat="identity")+theme_bw(base_size = 18,base_family = "Amiri")+
  geom_errorbar(aes(ymax=Errors+SEM,ymin=Errors-SEM),size=1.2,width=0.35)+theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Errors in Learning, Retention and Generalization According to Diagnosis")+#theme(legend.position = "none")+
  scale_y_continuous(name = "Erros")+scale_x_discrete(name="Diagnosis")+facet_wrap(~Rule*MeasureType)

two_way_fit <- summary(aov(AllData$Generalization~AllData$Rule*AllData$Diagnosis))
two_way_fit
AllData$yellowBoyGeneralization = ifelse(AllData$yellowBoy=="2",AllData$Face2ErrorPh3Prime,ifelse(AllData$yellowBoy=="4",AllData$Face4ErrorPh3Prime,0))
AllData$yellowWomanGeneralization = ifelse(AllData$yellowWoman=="2",AllData$Face2ErrorPh3Prime,ifelse(AllData$yellowWoman=="4",AllData$Face4ErrorPh3Prime,0))
AllData$brownGirlGeneralization = ifelse(AllData$brownGirl=="2",AllData$Face2ErrorPh3Prime,ifelse(AllData$brownGirl=="4",AllData$Face4ErrorPh3Prime,0))
AllData$brownManGeneralization = ifelse(AllData$brownMan=="2",AllData$Face2ErrorPh3Prime,ifelse(AllData$brownGirl=="4",AllData$Face4ErrorPh3Prime,0))

AllData$MaleGeneralization = AllData$brownManGeneralization + AllData$yellowBoyGeneralization
AllData$FemaleGeneralization = AllData$brownGirlGeneralization + AllData$yellowWomanGeneralization

#I need the gender to do the gender differences