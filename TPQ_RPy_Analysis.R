
# Data Preparation --------------------------------------------------------

library(reshape2)
library(ggplot2)
library(labelled)
library(corrplot)

library(reticulate)
library(readxl)

# Prepare projections
library(readxl)
library(boot)
library(e1071)
library(ggplot2)
library(sf)
library(pROC)
library(randomForest)
library(caret)
library(plotROC)
library(psy)

par(pty = "s") # to remove the side panels when plotting the ROC curve

if (dir.exists("/home/abdulrahman/anaconda3/envs/mne/bin/")){
  print("Home PC was found")
  use_python ("/home/abdulrahman/anaconda3/envs/mne/bin/python3", required = TRUE)
}else if(dir.exists("/home/asawalma/anaconda3/envs/mne/bin/")){
  print("Office PC was found")
  use_python("/home/asawalma/anaconda3/envs/mne/bin/python3", required = TRUE)
} else {
  print("Some Windows PC was found")
  use_python("C:/Users/jsawa/Anaconda3/envs/mne/python.exe", required = TRUE)
}

# load libraries
np = import("numpy")
pd = import("pandas")
os = import("os")


# define directories of main files and python functions file
path =        '/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine'
scores_path = "/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx"

if (dir.exists("/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive")) {
  path = gsub("/abdulrahman/abdulrahman.sawalma@gmail.com",
              "/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive",
              path)
  scores_path = gsub("/abdulrahman/abdulrahman.sawalma@gmail.com",
                     "/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive",
                     scores_path)
}
if (dir.exists("G:/My Drive")) {
  path = gsub("/home/abdulrahman/abdulrahman.sawalma@gmail.com",
              "G:/My Drive",
              path)
  scores_path = gsub("/home/abdulrahman/abdulrahman.sawalma@gmail.com",
                     "G:/My Drive",
                     scores_path)
}


# load original scores
TPQ_scores = as.data.frame(read_xlsx(scores_path, na = "NA"))
gen_info = TPQ_scores[c("Diagnosis", "Final ID", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session")]

# Load question list for each dimension, it will be used in data analysis below
TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))

npz_tpq_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz')


prepare_icasso <- function(npz_tpq_path){
  # prepare the data from the npz files sent to me by Juergen.
  # It creates the following five data frames/groups of data frames:
  # data_tpq, data_reco, scores, projections and sources
  # All DFs will have column names indicating what they represent
  # :param npz_tpq_path: the path to the tpq npz (which is the one that contains the basic info such as IDs)
  # :return: five data frames (see description)
  
  # load the tpq numpy array
  npz_tpq = np$load(npz_tpq_path, allow_pickle=TRUE)
  
  # load the results numpy array
  results_path = gsub("icasso_ICA-tpq",
                      "icasso-results_ICA-tpq",
                      npz_tpq_path)
  
  npz_results = np$load(results_path, allow_pickle=TRUE)
  
  
  # find the indices of original scores that correspond to the IDs in npz_tpq["IDs"]
  inc_ids = npz_tpq["IDs"] %in% gen_info$`Final ID` # IDs in numpy array that are found in the original scores file
  matched_inds = match(npz_tpq["IDs"],gen_info$`Final ID`)
  matched_inds = matched_inds[inc_ids]
  sub_info = data.frame(
    ID = factor(npz_tpq["IDs"][inc_ids],levels = npz_tpq["IDs"]),
    Diagnosis = factor(gen_info$Diagnosis[matched_inds]),
    Trauma = factor(gen_info$Trauma[matched_inds]),
    GAD = factor(gen_info$GAD[matched_inds]),
    Response = factor(gen_info$Response[matched_inds]),
    PTSD = factor(gen_info$PTSD[matched_inds]),
    MDD = factor(gen_info$MDD[matched_inds]),
    Session = factor(gen_info$Session[matched_inds])
  )
  
  data_tpq = data.frame(npz_tpq["data_tpq"])[inc_ids,]
  data_tpq[,c("ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session")] = sub_info
  
  colnames(data_tpq) = c(paste0("Q",1:100), c("ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD","MDD", "Session"))
  
  # rearrange the columns, just for aesthetics
  data_tpq = data_tpq[,c(c("ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD","MDD", "Session"),paste0("Q",1:100))]
  
  data_reco_df = npz_results["data_reco"]
  data_reco = list()
  for (i in 1:length(data_reco_df[,1,1])){
    new_df = cbind(sub_info,data.frame(data_reco_df[i,inc_ids,]))
    colnames(new_df) = c(c("ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD","MDD", "Session"),paste0("Q",1:100))
    data_reco[[i]] = new_df
  }
  
  scores = npz_results["scores"]
  projections = data.frame(npz_results['projection'])[inc_ids,]
  colnames(projections)  = paste0("IC",1:ncol(projections))
  projections = cbind(sub_info,projections)
  
  sources = data.frame(npz_results["sources"])
  colnames(sources) = paste0("Q",1:100)
  sources$IC = paste0("IC", 1:nrow(sources))
  
  df_list = list(tpq= data_tpq, reco = data_reco, scores = scores, projections = projections, sources = sources)
  return (df_list)
}

# Create a reconstructed data using all 15 ICs
npz_tpq_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz')
results_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD/icasso-results_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz')
npz_tpq = np$load(npz_tpq_path, allow_pickle=TRUE)
npz_results = np$load(results_path, allow_pickle=TRUE)
sources = npz_results["sources"]
sources = sources[1:15,]

icasso = npz_tpq['icasso'][[1]]
unmixing = icasso$get_centrotype_unmixing()# of shape: [n_cluster, n_chan]
unmixing = unmixing[1:15,]
All_reco=np$transpose(icasso$ica2data(sources, unmixing, idx_keep=0:14))
All_reco = All_reco[npz_tpq["IDs"] %in% gen_info$`Final ID`,]


# define paths for variuos decompositions
all_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz')
all_fast_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD_FastICA_MNE/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEfastica.npz')
all_fastsk_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD_FastICA_sklearn/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_SklearnFastICA.npz')
hc_path  = paste0(path, '/ICASSO_fixed/decomp_tpq/HC/icasso_ICA-tpq_HC_nsamp1202_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
mdd_path = paste0(path, '/ICASSO_fixed/decomp_tpq/MDD/icasso_ICA-tpq_MDD_nsamp455_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
gad_path = paste0(path, '/ICASSO_fixed/decomp_tpq/GAD/icasso_ICA-tpq_GAD_nsamp30_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
ptsd_path = paste0(path, '/ICASSO_fixed/decomp_tpq/PTSD/icasso_ICA-tpq_PTSD_nsamp57_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
tnp_path = paste0(path, '/ICASSO_fixed/decomp_tpq/TNP/icasso_ICA-tpq_TNP_nsamp78_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
disorders_path = paste0(path, '/ICASSO_fixed/decomp_tpq/MDD,PTSD,TNP,GAD/icasso_ICA-tpq_MDD,PTSD,TNP,GAD_nsamp620_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')

# Load all data into lists of variables
ICASSO_all = prepare_icasso(all_path)
ICASSO_all_fast = prepare_icasso(all_fast_path)
ICASSO_all_fastsk = prepare_icasso(all_fastsk_path)
ICASSO_hc = prepare_icasso(hc_path)
ICASSO_mdd = prepare_icasso(mdd_path)
ICASSO_gad = prepare_icasso(gad_path)
ICASSO_ptsd = prepare_icasso(ptsd_path)
ICASSO_tnp = prepare_icasso(tnp_path)
ICASSO_disorders = prepare_icasso(disorders_path)

df_cor <- function(df_1, df_2, Factor, max_components_1, max_components_2, estimate_threshold, stringent_alpha = TRUE,Projection = FALSE){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # A Function to find the correlation between two sources
  # Returns a matrix of spearman's rho values
  # -> stringent_alpha: if TRUE, the threshold will be
  #    divided by number of comparisons
  # -> Projection: set Projection to TRUE if you are
  #    including projections data frames
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  if (is.na(max_components_1)){
    max_components_1 = nrow(df_1)
    max_components_2 = nrow(df_2)
  }
  
  
  if (Projection) {
    df_1 = df_1[, 1:(max_components_1 + 8)]
    df_2 = df_2[, 1:(max_components_2 + 8)]
    
    # match the subjects within both data frames
    df_1 = df_1[!is.na(match(df_1$ID,df_2$ID)),]
    df_2 = df_2[!is.na(match(df_2$ID,df_1$ID)),]
    
    #arrange IDs to match in both data frames. Usually no need, but for extra assurance
    df_1 = df_1[match(df_1$ID,df_2$ID),]
    
    # flip them so that it has the same shape as the sources
    df_1 = t(df_1[,9:(max_components_1+8)])
    df_2 = t(df_2[,9:(max_components_2+8)])
    
  } else {
    df_1 = df_1[1:max_components_1, 1:100]
    df_2 = df_2[1:max_components_2, 1:100]
  }
  
  Cor_mat = matrix(0,nrow = max_components_1,ncol = max_components_2)
  dimnames(Cor_mat) = list(paste0("IC",1:max_components_1), paste0("IC",1:max_components_2))
  names(dimnames(Cor_mat)) = c("All",Factor)
  
  for (i in 1:max_components_1){
    threshold = 0.05
    if (stringent_alpha) {
      threshold = threshold / (max_components_2 * max_components_1)
    }
    
    estimates = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]))$estimate)
    p_values = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]))$p.value)
    estimates[p_values>threshold | abs(estimates)<estimate_threshold] = 0
    Cor_mat[i,] = estimates
  }
  return(Cor_mat)
}


simplify_mat <- function(sources_Cor){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # A function to simplify the correlation matrix into
  #   a 3-column matrix:
  #     1st col -> ICs of df1, 2 -> ICs of df2 and 3 -> spearman rho values
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  final_matrix = matrix(rep(0, 3 * nrow(sources_Cor)), nrow = 3)
  for (i in 1:nrow(sources_Cor)) {
    IC = paste0("IC", i)
    IC_all = IC
    
    #find the IC in the second df that has a correlation of >0
    # If you have two, choose the one with highest cor
    max_cor_IC = abs(sources_Cor[IC, ]) > 0 &
      abs(sources_Cor[IC, ]) == max(abs(sources_Cor[IC, ]))
    
    IC_df2 = colnames(sources_Cor)[max_cor_IC]
    Cor_value = round(sources_Cor[IC, max_cor_IC], digits = 2)
    if (sum(max_cor_IC) == 0) {
      final_matrix[, i] = c(NA, NA, NA)
    } else{
      final_matrix[, i] = c(IC_all, IC_df2, Cor_value)
    }
    dimnames(final_matrix) = list(c("All","DF2","Cor"), paste0("IC",1:10))
  }
  
  return(final_matrix)
}


shuffle_ICs <- function(melted_df,ICs_mat, Factor,Projection = FALSE){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Change the order of sources/ projections to reflect the IC of
  # greatest correlation with the sources of the All-decomposition
  # e.g. if IC3 of MDD group correlates with IC1 of All-decomposition
  # It will be put next to IC1 in the specified data frame
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  if (Projection){
    if (length(ICs_mat)>0){
      for (i in 1:ncol(ICs_mat)){
        melted_df$IC_group[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[1,i]] = ICs_mat[2,i]
        melted_df$Projections_group[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[1,i]] = melted_df$Projection[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[2,i]]
        melted_df$Cor[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[1,i]] = ICs_mat[3,i]
      }
    }
  }else {
    if (length(ICs_mat)>0){
      for (i in 1:ncol(ICs_mat)){
        melted_df$IC_group[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[1,i]] = ICs_mat[2,i]
        melted_df$Sources_group[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[1,i]] = melted_df$Source[melted_df$Group==Factor&melted_df$IC_all== ICs_mat[2,i]]
      }
    }
  }
  return(melted_df)
}

prep_proj <- function(df_projection){
  # Create all_na which is an empty data frame to fill in with the data of
  # each projections data frame. E.g. mdd_projections will have all subjects
  # but only the MDD subjects have values from mdd_projections, the rest is NA
  
  all_na = as.data.frame(matrix(NA, nrow = nrow(ref_projections), ncol = 10))
  all_na[ref_projections$ID %in%df_projection$ID, ] = df_projection[9:18]
  return(all_na)
}


test_roc = function(roc1, roc2, name1 = "ROC1", name2 = "ROC2", control = "HC",Subtitle = ""){
  d.len1 = length(roc1$sensitivities)
  d.len1 = length(roc2$sensitivities)
  perc1 = c(rep("",(d.len1-1)),paste0("AUC = ",round(as.numeric(roc1$auc),2),"%"))
  perc2 = c(rep("",(d.len1-1)),paste0("AUC = ",round(as.numeric(roc2$auc),2),"%"))
  
  case = levels(roc1$original.response)[levels(roc1$original.response)!=control]
  
  roc1_df = data.frame(D=factor(roc1$original.response, ordered = TRUE, levels = c(control, case)), M1 = roc1$original.predictor)
  roc2_df = data.frame(D=factor(roc2$original.response, ordered = TRUE, levels = c(control, case)), M1 = roc2$original.predictor)
  
  color1 = "#808CA3"
  color2 = "#2F4F4F"
  
  roc_test = roc.test(roc1,roc2)
  roc_p = round(roc_test$p.value,4)
  roc_stat = round(roc_test$statistic,3)
  roc_est = round(roc_test$parameter,3)
  roc_test_text = paste0("DeLong's test for two ROC curves\nD(",roc_est,") = ", roc_stat,", p value = ",roc_p)
  g_roc = ggplot()+
    geom_roc(data = roc1_df, mapping = aes(d = D, m = M1), color = color1, n.cuts = 0)+style_roc()+
    geom_roc(data = roc2_df, mapping = aes(d = D, m = M1), color = color2, n.cuts = 0)+
    theme_bw(base_size = 18,base_family = "Amiri")+
    theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5,face = "italic"))+
    geom_text(aes(x= 0.8,y=0.45,label = perc1), color = color1)+
    geom_text(aes(x= 0.8,y=0.4,label = perc2), color = color2)+
    geom_text(aes(x= 0.5,y=0.45,label = c(rep("",(d.len1-1)),name1)), color = color1)+
    geom_text(aes(x= 0.5,y=0.4,label = c(rep("",(d.len2-1)),name2)), color = color2)+
    geom_text(aes(x= 0.65,y=0.25,label = roc_test_text), color = "black")+
    #scale_x_continuous(name = "100 - Specificity (%)", labels = c(0,25,50,75,100))+
    #scale_y_continuous(name = "Sensitivity (%)", breaks = c(0.0,0.25,0.5,0.75,1.0), labels = c(0,25,50,75,100))+
    ggtitle(paste0("ROC Curves for ",name1," VS ",name2), subtitle = Subtitle)
  
  
  
  return(list(roc_plot = g_roc, roc_test = roc_test))
}

LogisticFunction = function(Data, DV, IVs,control = "HC" , Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE, k.fold = 10){
  #The Model is a logistic model (DV~IV1+IV2..., family = binomial())
  Data = Data[colnames(Data) %in% DV| colnames(Data) %in% IVs]
  Data = na.omit(Data)
  
  # create a training and a test data set
  train_inds = sample(1:nrow(Data), size = perc*nrow(Data))
  train = Data[train_inds,]
  test = Data[-train_inds,]
  
  
  Formula = as.formula(paste0(DV," ~ ",paste0(IVs,collapse = "+")))
  Model = glm(Formula, data = train, family = binomial())
  
  # calculate accuracy and prediction error
  folds = createFolds(train[[DV]], k = k.fold)
  accuracy = ""
  pred_error = ""
  cv = lapply(folds, function(fold) {
    train_fold = train[-fold,]
    test_fold = train[fold,]
    
    model_fold = glm(Formula, data = train_fold, family = binomial())
    preds_fold = predict(model_fold, test_fold,type='response')
    outcome = ifelse(test_fold[[DV]] == control,0,1)
    accuracy_fold = mean(outcome == ifelse(preds_fold>0.5,1,0))
    pred_error = mean((outcome- preds_fold)^2)
    
    # average ROC 
    roc.info_fold = "No ROC"
    if (nlevels(test_fold[[DV]]) == 2){
      test_values = factor(test_fold[[DV]])
      test_values = relevel(test_values, control)
      test_values_num = ifelse(test_values==control, 0,1)
      roc.info_fold = roc(test_values, preds_fold, plot = FALSE, levels = levels(test_values),
                          percent = TRUE,legacy.axes = TRUE, print.auc = TRUE, auc = TRUE)
    }
    return(list(length(test_values),accuracy_fold,roc.info_fold$auc,pred_error,roc.info_fold))
  })
  cv_values = as.data.frame(t(sapply(cv, function(x) as.numeric(x[1:4]))))
  roc_list = sapply(cv, function(x) x[5])
  accuracy = paste0(100*round(mean(cv_values[,2]),digits = 6),"%")
  pred_error = paste0(100*round(mean(cv_values[,4]),digits = 3),"%")
  AUC = paste0(round(mean(cv_values[,3]),digits = 3),"%")
  
  Factor = as.character(Formula[[2]])
  Summary=summary(Model)
  n.value=Model$df.null+1
  Coefficients = Summary$coefficients
  
  #The improvement we got (change in deviance) after including predictors (baseline deviance - model deviance)
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
    Num=sum((test$predicted.probability>interval[1]&test$predicted.probability<interval[2]))
    TList=append(TList, Num)
  }
  y.value=max(TList)/2
  
  GText_0=paste("Specificity =", (100*Specificity), "%  || ", "NPV =", (100*NPV), "%  || ","Accuracy =", accuracy)
  GText_1=paste("Sensitivity =", (100*Sensitivity), "%  || ", "PPV =", (100*PPV), "%  || ","Pred.Error =", pred_error)
  
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
                   percent = TRUE,legacy.axes = TRUE, print.auc = TRUE, auc = TRUE)
  }
  LogValues=data.frame(Beta, SE, Lower.CI, odds_ratio, Upper.CI, Zvalues, Pvalues)
  DerivedValues=data.frame(chisq=modelChi, df=chidf, p.value.chi=chisq.prob, r2.hl=R2.hl,
                           r2.cs=R2.cs, r2.n=R2.n, sensitivity=Sensitivity, specificity=Specificity,
                           cv_accuracy = accuracy, cv_MSE = pred_error, cv_AUC = AUC)
  DataList=list(data_frame = data_frame_essential,
                log.values=LogValues, derived.values=DerivedValues,
                prob.formula = paste("Probability of Y occuring = ","1/(1 + e ^ -(",Xs,")",sep = ""),
                PredMatrix = PredMatrix,
                ROC = roc.info,
                roc_cv_list = roc_list,
                Model = Model)
  
  print (Drawing)
  
  return (DataList)
  
}

show_svm = function(data_frame, DV, IVs, control = "HC", res = 75, k_type = "linear", perc = 0.8,
                    plot.ROC=TRUE, k.fold = 0, weighted = TRUE){
  
  # Function to calculate SVM. It will plot the SVM function if there were only two IVs
  # It will also plot an ROC curve, if specified, and if there are exactly two levels of the DV
  # Also, when k.fold > 1, there will be no plotting of any type
  # data_frame: the dataframe to be included,
  # DV: dependent variable
  # IVs: independent variables to predict the model from
  # control: the level name of controls in the DV, usually set to "HC"
  # res: resolution of the plotted polygons
  # k_type= kernel type, can be linear, polynomial, radial or sigmoid
  # perc: percentage of training data set
  # plot.ROC: boolean to plot ROC curve if needed
  # k.fold: number of folds to test the data on
  
  data_frame = data_frame[colnames(data_frame) %in% DV| colnames(data_frame) %in% IVs]
  data_frame = na.omit(data_frame)
  data_frame[,IVs] = scale(data_frame[,IVs])
  
  # create a test and train data sets
  
  train_inds = sample(1:nrow(data_frame), size = perc*nrow(data_frame))
  train = data_frame[train_inds,]
  test = data_frame[-train_inds,]
  
  make_grid = function(x1,x2,n = 75){  
    ## This function only creates a range of dots
    # These dots will be colored according to the predicted value based on our data
    x1 = seq(from = min(x1)-0.5, to = max(x1)+0.5, length = n)
    x2 = seq(from = min(x2)-0.5, to = max(x2)+0.5, length = n)
    
    new_df = expand.grid(X1 = x1, X2 = x2)
    colnames(new_df) = colnames(data_frame)[1:2]
    
    return(new_df)
  }
  
  convert_to_sf = function(new_df){
    # The next part converts the grid data frame to a group of polygons of class sf
    pre_raster_df = new_df
    pre_raster_df = cbind(pre_raster_df, cat = rep(1L, nrow(pre_raster_df)),stringsAsFactors = FALSE)
    
    # formula coordinates will be based on the colnames of the data input (e.g. ~X1+X2)
    cor_form = as.formula(paste0("~",colnames(new_df)[1],"+", colnames(new_df)[2]))
    coordinates(pre_raster_df) = cor_form
    gridded(pre_raster_df) <- TRUE
    
    # create raster df
    raster_df <- raster(pre_raster_df)
    
    # create spatial polygons
    sp_df = rasterToPolygons(raster_df, dissolve = TRUE)
    
    # convert polygons of class sf
    sf_polygons = st_as_sf(sp_df)
    
    return(sf_polygons)
  }
  
  if (weighted){
    classes = levels(train[[DV]])
    c.weights = sapply(classes, function(x) length(train[[DV]])/sum(train[[DV]] == x))
  }
  Formula = as.formula(paste0(DV," ~ ",paste0(IVs,collapse = "+")))
  svm_model = svm(Formula, data = train, kernel = k_type, cost = 1, scale = TRUE, probability = TRUE, cross = k.fold, class.weights = c.weights)
  
  if (length(IVs) == 2){
    grid = make_grid(train[[IVs[1]]],train[[IVs[2]]], n = res)
    preds = predict(svm_model, grid)
    predicted_df = data.frame(X1 = grid[,1], X2 = grid[,2], Y=preds)
    sf_polygons = convert_to_sf(predicted_df)
    
    Colors = c("#C33E3B","#4EA3DF","#6cBE58","#808CA3","#B9B0AB","#2F4F4F", "#CC6666", "#9999CC", "#66CC99","#682860","#FBEC5D","#FF6347","#FF3800","#1B4D3E","#E30B5D")
    # plot the model
    g_svm = ggplot()+
      geom_sf(data = sf_polygons, alpha = 0.25, mapping= aes(x=NULL,y=NULL,group = Y,fill = Y))+
      scale_fill_gradientn(breaks = 1:length(IVs),colors = Colors[1:length(IVs)])+
      geom_point(data = train,mapping = aes(x = .data[[IVs[1]]], y =.data[[IVs[2]]],color = .data[[DV]]),size = 3)+
      scale_color_manual(values = Colors)+
      ggtitle("Data points of the variables X0 and X1", subtitle = "Original Data with Predicted Values")+
      MinimalTheme
    
    print(g_svm)
  }
  
  
  preds = predict(svm_model, test, probability = TRUE)
  roc.info = "No ROC"
  if (nlevels(test[[DV]]) == 2){
    cases_probs = 1-as.data.frame(attr(preds, "probabilities"))[[control]]
    test_values = factor(test[[DV]])
    test_values = relevel(test_values, control)
    roc.info = roc(test_values, cases_probs, plot = plot.ROC, levels = levels(test_values),
                   percent = TRUE,legacy.axes = TRUE, print.auc = TRUE, auc = TRUE)
  }
  Accuracy = round(mean(preds== test[[DV]]),digits = 4)
  res_table = table(predicted = preds, actual = test[[DV]])
  if (k.fold>1){
    folds = createFolds(train[[DV]], k = k.fold)
    
    # in cv we are going to apply a created function to our 'folds'
    cv = lapply(folds, function(x) { # start of function
      training_fold = train[-x, ]
      test_fold = train[x, ] # here we describe the test fold individually
      # now apply (train) the classifer on the training_fold
      svm_fold = svm(formula = Formula,
                     data = training_fold,
                     type = 'C-classification',
                     kernel = k_type,
                     probability = TRUE, cost = 1, scale = TRUE,class.weights = c.weights)
      
      preds_fold = predict(svm_fold, test_fold, probability = TRUE)
      roc.info_fold = "No ROC"
      if (nlevels(test[[DV]]) == 2){
        cases_probs = 1-as.data.frame(attr(preds_fold, "probabilities"))[[control]]
        test_values = factor(test_fold[[DV]])
        test_values = relevel(test_values, control)
        test_values_num = ifelse(test_values==control, 0,1)
        roc.info_fold = roc(test_values, cases_probs, plot = FALSE, levels = levels(test_values),
                            percent = TRUE,legacy.axes = TRUE, print.auc = TRUE)
      }
      Accuracy_fold = round(mean(preds_fold== test_fold[[DV]]),digits = 4)
      pred_error_fold = round(mean((test_values_num-cases_probs)^2),digits = 4)
      
      return(c(Accuracy_fold,roc.info_fold$auc,pred_error_fold))
    })
    cv = t(data.frame(cv))
    colnames(cv) = c("Accuracy","AUC", "Pred.Error")
  }else{
    cv = "Data were not cross validated"
  }
  return(list(svm = svm_model,accuracy = svm_model$tot.accuracy, c_matrix = res_table, ROC = roc.info, cv = cv))
}


show_knn = function(data_frame, DV, IVs, perc, control = "HC", plot.ROC = TRUE, k.fold = 10, k=10){
  nor <-function(x) {(x -min(x))/(max(x)-min(x))}
  data_frame = data_frame[colnames(data_frame) %in% DV| colnames(data_frame) %in% IVs]
  data_frame = na.omit(data_frame)
  
  # create a test and train data sets
  
  train_inds = sample(1:nrow(data_frame), size = perc*nrow(data_frame))
  data_frame[,IVs] = as.data.frame(lapply(data_frame[,IVs], nor))
  train = data_frame[train_inds,]
  test = data_frame[-train_inds,]
  knn_model = caret::knn3(train[,IVs], train[[DV]], k = k)
  predicted = as.data.frame(predict(knn_model, test[,IVs]))
  case_name = "NA"
  if (nlevels(test[[DV]]) == 2){
    case_name = colnames(preds_fold)[colnames(preds_fold)!=control]
    
    cases_probs = 1-predicted[colnames(predicted) == control][[1]]
    test_values = factor(test[[DV]])
    test_values = relevel(test_values, control)
    roc.info = roc(test_values, cases_probs, plot = plot.ROC, levels = levels(test_values),
                   percent = TRUE,legacy.axes = TRUE, print.auc = TRUE, auc = TRUE)
  }
  if (k.fold>1){
    folds = createFolds(train[[DV]], k = k.fold)
    
    # in cv we are going to apply a created function to our 'folds'
    cv = lapply(folds, function(x) { # start of function
      training_fold = train[-x, ]
      test_fold = train[x, ] # here we describe the test fold individually
      # now apply (train) the classifer on the training_fold
      knn_fold = caret::knn3(training_fold[,IVs], training_fold[[DV]], k = k)
      
      preds_fold = predict(knn_fold, test_fold[,IVs], probability = TRUE)
      cases_probs = 1-preds_fold[,colnames(preds_fold) == control]
      controls = ifelse(preds_fold[,colnames(preds_fold) == control]>0.5,0,1)
      roc.info_fold = "No ROC"
      if (nlevels(test[[DV]]) == 2){
        
        test_values = factor(test_fold[[DV]])
        test_values = relevel(test_values, control)
        test_values_num = ifelse(test_values==control, 0,1)
        roc.info_fold = roc(test_values, cases_probs, plot = FALSE, levels = levels(test_values),
                            percent = TRUE,legacy.axes = TRUE, print.auc = TRUE)
      }
      preds_category_fold = ifelse(preds_fold[colnames(preds_fold) == control]>0.5,control,case_name)
      Accuracy_fold = round(mean(preds_category_fold== test_fold[[DV]]),digits = 4)
      pred_error_fold = round(mean((test_values_num-cases_probs)^2),digits = 4)
      
      
      return(c(Accuracy_fold,roc.info_fold$auc,pred_error_fold))
    })
    cv = t(data.frame(cv))
    colnames(cv) = c("Accuracy","AUC", "Pred.Error")
  }else{
    cv = "Data were not cross validated"
  }
  
  preds_categories = ifelse(predicted[colnames(predicted) == control]>0.5,control,case_name)
  Accuracy = mean(preds_categories == test_category)
  
  res_table = table(predicted = preds_categories, actual = test[[DV]])
  
  
  return(list(knn = knn_model,accuracy = Accuracy, c_matrix = res_table, ROC = roc.info, cv = cv))
}


# Correlations Between Sources of Different Data Decompositions ---------------

library(reshape2)
library(ggplot2)

ref_sources = ICASSO_all$sources
fast_sources = ICASSO_all_fast$sources
fastsk_sources = ICASSO_all_fastsk$sources
hc_sources = ICASSO_hc$sources
gad_sources = ICASSO_gad$sources
mdd_sources = ICASSO_mdd$sources
ptsd_sources = ICASSO_ptsd$sources
tnp_sources = ICASSO_tnp$sources


# Create the matrix of correlations' estimates
est_threshold = 0.3
AllFAST_sources_Cor = df_cor(ref_sources, fast_sources, "All_Fast", 10,10,est_threshold,Projection = FALSE)
AllFASTSK_sources_Cor = df_cor(ref_sources, fastsk_sources, "All_Fast_SK", 10,10,est_threshold,Projection = FALSE)
HC_sources_Cor = df_cor(ref_sources, hc_sources, "HC", 10,10,est_threshold,Projection = FALSE)
GAD_sources_Cor = df_cor(ref_sources, gad_sources, "GAD", 10,10,est_threshold,Projection = FALSE)
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, "MDD", 10,10,est_threshold,Projection = FALSE)
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, "PTSD", 10,10,est_threshold,Projection = FALSE)
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, "TNP", 10,10,est_threshold,Projection = FALSE)


# Find the ICs in the each sources df that has the highest correlation
AllFAST_sources_ICs = simplify_mat(AllFAST_sources_Cor)
AllFASTSK_sources_ICs = simplify_mat(AllFASTSK_sources_Cor)
HC_sources_ICs = simplify_mat(HC_sources_Cor)
GAD_sources_ICs = simplify_mat(GAD_sources_Cor)
MDD_sources_ICs = simplify_mat(MDD_sources_Cor)
PTSD_sources_ICs = simplify_mat(PTSD_sources_Cor)
TNP_sources_ICs = simplify_mat(TNP_sources_Cor)


all_sources = data.frame(
  IC_all = rep(paste0("IC", 1:10), 100),
  All = melt(ref_sources[1:10, ])$value,
  All_fast = melt(fast_sources[1:10, ])$value,
  All_fastsk = melt(fastsk_sources[1:10, ])$value,
  HC = melt(hc_sources[1:10, ])$value,
  GAD = melt(gad_sources[1:10, ])$value,
  MDD = melt(mdd_sources[1:10, ])$value,
  PTSD = melt(ptsd_sources[1:10, ])$value,
  TNP = melt(tnp_sources[1:10, ])$value
)

# Create the final data frames that is ready for plotting
all_sources$Question = rep(paste0("Q",1:100),each = 10)
all_sources_melted = melt(all_sources, id.vars = c("Question","IC_all"), value.name = "Source", variable.name = "Group")

all_sources_melted$IC_group = NA


all_sources_melted = shuffle_ICs(all_sources_melted, AllFAST_sources_ICs[,!is.na(AllFAST_sources_ICs[2,])],"All_fast")
all_sources_melted = shuffle_ICs(all_sources_melted, AllFASTSK_sources_ICs[,!is.na(AllFASTSK_sources_ICs[2,])],"All_fastsk")
all_sources_melted = shuffle_ICs(all_sources_melted, HC_sources_ICs[,!is.na(HC_sources_ICs[2,])],"HC")
all_sources_melted = shuffle_ICs(all_sources_melted, GAD_sources_ICs[,!is.na(GAD_sources_ICs[2,])],"GAD")
all_sources_melted = shuffle_ICs(all_sources_melted, PTSD_sources_ICs[,!is.na(PTSD_sources_ICs[2,])],"PTSD")
all_sources_melted = shuffle_ICs(all_sources_melted, TNP_sources_ICs[,!is.na(TNP_sources_ICs[2,])],"TNP")
all_sources_melted = shuffle_ICs(all_sources_melted, MDD_sources_ICs[,!is.na(MDD_sources_ICs[2,])],"MDD")
all_sources_melted$Sources_group[all_sources_melted$Group == "All"] = all_sources_melted$Source[all_sources_melted$Group == "All"]
all_sources_melted$IC_group_plot = as.vector(sapply(1:8, function(x)
  c(all_sources_melted$IC_group[all_sources_melted$Question == "Q1"][((x - 1) * 10 + 1):(x * 10)], rep(NA, 990))))
# Add correlation values only for one question (10 values, 1 per IC),for plotting purposes
all_sources_melted$Cor_plot = c(
  rep(NA, 1000),
  c(AllFAST_sources_ICs[3, ], rep(NA, 990)),
  c(AllFASTSK_sources_ICs[3, ], rep(NA, 990)),
  c(HC_sources_ICs[3, ], rep(NA, 990)),
  c(GAD_sources_ICs[3, ], rep(NA, 990)),
  c(MDD_sources_ICs[3, ], rep(NA, 990)),
  c(PTSD_sources_ICs[3, ], rep(NA, 990)),
  c(TNP_sources_ICs[3, ], rep(NA, 990))
)

all_sources_melted$IC_all = factor(all_sources_melted$IC_all, levels = paste0("IC",1:10))
all_sources_melted$Sources_all = rep(all_sources_melted$Sources[all_sources_melted$Group=="All"],nlevels(all_sources_melted$Group))

g1 = ggplot(data = all_sources_melted, aes(x=Question, y= Sources_group, group =Group))+
  geom_bar(stat = "identity")+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x= 30,y=5,label = IC_group_plot))+
  geom_text(aes(x= 70,y=5,label = Cor_plot))

ggplot(all_sources_melted[!(all_sources_melted$Group == "All_fast" | all_sources_melted$Group == "All_fastsk"| all_sources_melted$IC_all == "IC9"|all_sources_melted$IC_all == "IC10"),],
       aes (x= Sources_group, y= Sources_all, group = Group))+
  geom_point( size =1)+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x=-8,y=-5,label = IC_group_plot))+
  geom_text(aes(x= -2,y=-5,label = Cor_plot))


# Addition: To make it like the projections plot, I will scale the values and sort them by group
# and add a point of minimal value
all_sources_melted$Sources_zero = all_sources_melted$Sources_group

for (IC in paste0("IC",1:10)){
  for (Group in levels(all_sources_melted$Group)){
    subgroup_cond = all_sources_melted$Group ==Group & all_sources_melted$IC_all == IC
    pro_val = all_sources_melted$Sources_group[subgroup_cond]
    all_sources_melted$Sources_group_scaled[subgroup_cond] = pro_val/mean(abs(pro_val), na.rm = TRUE)
    min_value = min(abs(all_sources_melted$Sources_group_scaled[subgroup_cond]))
    all_sources_melted$Sources_zero[subgroup_cond][abs(all_sources_melted$Sources_group_scaled)[subgroup_cond] != min_value] = NA
  }
}

all_sources_melted$Sources_zero[!is.na(all_sources_melted$Sources_zero)] = all_sources_melted$Sources_group_scaled[!is.na(all_sources_melted$Sources_zero)]
all_sources_melted$Sources_all_scaled = rep(all_sources_melted$Sources_group_scaled[all_sources_melted$Group == "All"],8)
all_sources_melted$Sources_sorted_all = rep(all_sources_melted$Sources_sorted[all_sources_melted$Group == "All"],8)
all_sources_melted$Sources_sorted_all[is.na(all_sources_melted$Sources_sorted)] = NA
all_sources_melted$Sources_all_scaled[is.na(all_sources_melted$Sources_group_scaled)] = NA

for (IC in paste0("IC",1:10)){
  for (Group in levels(all_sources_melted$Group)){
    subgroup_cond = all_sources_melted$Group ==Group & all_sources_melted$IC_all == IC 
    all_sources_melted$Sources_sorted[subgroup_cond] = match(all_sources_melted$Sources_group_scaled[subgroup_cond], sort(all_sources_melted$Sources_group_scaled[subgroup_cond]))
  }
}

sources_melted_all = all_sources_melted[(all_sources_melted$Group =="All" | all_sources_melted$Group =="All_fast" | all_sources_melted$Group =="All_fastsk"),]
sources_melted_all$IC_Fast = sources_melted_all$IC_group_plot[sources_melted_all$Group == "All_fast"]
sources_melted_all$IC_Fastsk = sources_melted_all$IC_group_plot[sources_melted_all$Group == "All_fastsk"]
sources_melted_all$Cor_Fast = sources_melted_all$Cor_plot[sources_melted_all$Group == "All_fast"]
sources_melted_all$Cor_Fastsk = sources_melted_all$Cor_plot[sources_melted_all$Group == "All_fastsk"]

all_sources_melted_dia = melt(all_sources_melted, measure.vars = c("Sources_all_scaled","Sources_group_scaled"), variable.name = "Sources_type", value.name = "Sources_scaled")
all_sources_melted_dia = all_sources_melted_dia[(all_sources_melted_dia$Group !="All" & all_sources_melted_dia$Group !="All_fast"& all_sources_melted_dia$Group !="All_fastsk"),]
all_sources_melted_dia$Sources_type = as.character(all_sources_melted_dia$Sources_type)
all_sources_melted_dia$Sources_type[all_sources_melted_dia$Sources_type == "Sources_all_scaled"] = "Combined"
all_sources_melted_dia$Sources_type[all_sources_melted_dia$Sources_type == "Sources_group_scaled"] = "Separate"

ggplot(data = all_sources_melted_dia, mapping = aes(x = Sources_sorted_all,y=Sources_type))+
  geom_raster(aes(fill = Sources_scaled))+
  scale_x_continuous("Sorted Questions")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradientn(colors = c("#02422f","#1B4D3E","#FFFFFF","#85015d","#69026b"), name = "Scaled\nScores    ")+
  facet_grid(IC_all~Group, scale = "free")+
  geom_text(aes(x= 5,y=2.25,label = IC_group_plot))+
  geom_text(aes(x= 5,y=1.75,label = Cor_plot))+
  ggtitle("Sources Per Diagnosis Per IC", subtitle = "- Scaled by dividing on the mean of absolute values
          - Arranged by values of the \'Combined\' decompoistion
          - Each IC in the \'Separate\' decomposition represents the most correlating IC with the \'Combined\' decomposition
          - The name of the most correlating IC from the separate group and its value of correlation are on the left")+
  MinimalTheme


ggplot(data = sources_melted_all, mapping = aes(x = Sources_sorted_all,y= Group))+
  geom_raster(aes(fill = Sources_group_scaled))+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradientn(colors = c("#02422f","#1B4D3E","#FFFFFF","#85015d","#69026b"), name = "Scaled\nScores")+
  facet_grid(IC_all~., scale = "free")+
  geom_text(aes(x= 5,y=2.25,label = IC_Fast))+
  geom_text(aes(x= 5,y=1.75,label = Cor_Fast))+
  geom_text(aes(x= 5,y=3.25,label = IC_Fastsk))+
  geom_text(aes(x= 5,y=2.75,label = Cor_Fastsk))+
  ggtitle("Sources Per Decomposition Type")+
  MinimalTheme

# Correlations Between Projections of Different Data Decompositions -------

ref_projections = ICASSO_all$projections
fast_projections = ICASSO_all_fast$projections
fastsk_projections = ICASSO_all_fastsk$projections
hc_projections = ICASSO_hc$projections
gad_projections = ICASSO_gad$projections
mdd_projections = ICASSO_mdd$projections
ptsd_projections = ICASSO_ptsd$projections
tnp_projections = ICASSO_tnp$projections


# Create the matrix of correlations' estimates
est_threshold = 0.3
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, "All_Fast", 10,10,est_threshold,Projection = TRUE)
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)
HC_projections_Cor = df_cor(ref_projections, hc_projections, "HC", 10,10,est_threshold,Projection = TRUE)
GAD_projections_Cor = df_cor(ref_projections, gad_projections, "GAD", 10,10,est_threshold,Projection = TRUE)
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, "MDD", 10,10,est_threshold,Projection = TRUE)
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, "PTSD", 10,10,est_threshold,Projection = TRUE)
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, "TNP", 10,10,est_threshold,Projection = TRUE)


# Find the ICs in the each projections df that has the highest correlation
AllFAST_projections_ICs = simplify_mat(AllFAST_projections_Cor)
AllFASTSK_projections_ICs = simplify_mat(AllFASTSK_projections_Cor)
HC_projections_ICs = simplify_mat(HC_projections_Cor)
GAD_projections_ICs = simplify_mat(GAD_projections_Cor)
MDD_projections_ICs = simplify_mat(MDD_projections_Cor)
PTSD_projections_ICs = simplify_mat(PTSD_projections_Cor)
TNP_projections_ICs = simplify_mat(TNP_projections_Cor)




all_projections = data.frame(
  ID = rep(ref_projections$ID,10),
  IC_all = rep(paste0("IC", 1:10), each = nrow(ref_projections)),
  All = melt(ref_projections[ ,9:18])$value,
  All_fast = melt(fast_projections[ ,9:18])$value,
  All_fastsk = melt(fastsk_projections[ ,9:18])$value,
  HC = melt(prep_proj(hc_projections))$value,
  GAD = melt(prep_proj(gad_projections))$value,
  MDD = melt(prep_proj(mdd_projections))$value,
  PTSD = melt(prep_proj(ptsd_projections))$value,
  TNP = melt(prep_proj(tnp_projections))$value
)

# Create the final data frames that is ready for plotting
all_projections_melted = melt(all_projections, id.vars = c("ID","IC_all"), value.name = "Projection", variable.name = "Group")

all_projections_melted$IC_group = NA

all_projections_melted = shuffle_ICs(all_projections_melted, AllFAST_projections_ICs[,!is.na(AllFAST_projections_ICs[2,])],"All_fast", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, AllFASTSK_projections_ICs[,!is.na(AllFASTSK_projections_ICs[2,])],"All_fastsk", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, HC_projections_ICs[,!is.na(HC_projections_ICs[2,])],"HC", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, GAD_projections_ICs[,!is.na(GAD_projections_ICs[2,])],"GAD", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, PTSD_projections_ICs[,!is.na(PTSD_projections_ICs[2,])],"PTSD", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, TNP_projections_ICs[,!is.na(TNP_projections_ICs[2,])],"TNP", Projection = TRUE)
all_projections_melted = shuffle_ICs(all_projections_melted, MDD_projections_ICs[,!is.na(MDD_projections_ICs[2,])],"MDD", Projection = TRUE)
all_projections_melted$Projections_group[all_projections_melted$Group == "All"] = all_projections_melted$Projection[all_projections_melted$Group == "All"]
all_projections_melted$IC_group_plot = c(sapply(1:80, function(x)
  c(all_projections_melted$IC_group[all_projections_melted$ID == "C_32"][x], rep(NA, 1818))))
# Add correlation values only for one question (10 values, 1 per IC),for plotting purposes
all_projections_melted$Cor_plot = c(
  rep(NA, 18190),
  c(sapply(1:10, function(x) c(AllFAST_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(AllFASTSK_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(HC_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(GAD_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(MDD_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(PTSD_projections_ICs[3, ][x], rep(NA, 1818)))),
  c(sapply(1:10, function(x) c(TNP_projections_ICs[3, ][x], rep(NA, 1818))))
)

all_projections_melted$IC_all = factor(all_projections_melted$IC_all, levels = paste0("IC",1:10))
all_projections_melted$projections_all = rep(all_projections_melted$projections[all_projections_melted$Group=="All"],nlevels(all_projections_melted$Group))


gp = ggplot(data = all_projections_melted, aes(x=ID, y= Projections_group, group =Group))+
  geom_line(stat = "identity")+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x= 1000,y=2,label = IC_group_plot))+
  geom_text(aes(x= 1500,y=2,label = Cor_plot))



# Correlate Projections for All groups Combined ---------------------------

all_projections_grouped = data.frame(
  ID = rep(ref_projections$ID,10),
  IC_all = rep(paste0("IC", 1:10), each = nrow(ref_projections)),
  All = melt(ref_projections[ ,9:18])$value,
  Groups = melt(rbind(hc_projections[1:18],mdd_projections[1:18], ptsd_projections[1:18], tnp_projections[1:18],gad_projections[1:18]))$value
)

ap_grouped_melted = melt(all_projections_grouped, id.vars = c("ID","IC_all"), value.name = "Projection", variable.name = "Group")
ap_grouped_melted$Group = as.character(ap_grouped_melted$Group)

ap_grouped_melted$Group[ap_grouped_melted$Group != "All"] = as.character(rep(ref_projections$Diagnosis, 10))

ap_grouped_melted$IC_group = NA

ap_grouped_melted = shuffle_ICs(ap_grouped_melted, HC_projections_ICs[,!is.na(HC_projections_ICs[2,])],"HC", Projection = TRUE)
ap_grouped_melted = shuffle_ICs(ap_grouped_melted, GAD_projections_ICs[,!is.na(GAD_projections_ICs[2,])],"GAD", Projection = TRUE)
ap_grouped_melted = shuffle_ICs(ap_grouped_melted, PTSD_projections_ICs[,!is.na(PTSD_projections_ICs[2,])],"PTSD", Projection = TRUE)
ap_grouped_melted = shuffle_ICs(ap_grouped_melted, TNP_projections_ICs[,!is.na(TNP_projections_ICs[2,])],"TNP", Projection = TRUE)
ap_grouped_melted = shuffle_ICs(ap_grouped_melted, MDD_projections_ICs[,!is.na(MDD_projections_ICs[2,])],"MDD", Projection = TRUE)
ap_grouped_melted$Projections_group[ap_grouped_melted$Group == "All"] = ap_grouped_melted$Projection[ap_grouped_melted$Group == "All"]

for (h_lvl in levels(factor(ap_grouped_melted[ap_grouped_melted$Group!="All","Group"]))){
  for (l_lvl in levels(factor(ap_grouped_melted[ap_grouped_melted$Group!="All","IC_all"]))){
    IC_n = ap_grouped_melted$IC_group[ap_grouped_melted$Group == h_lvl & ap_grouped_melted$IC_all == l_lvl]
    if (is.na(IC_n[1])){IC_n[1] = "NA"}
    ap_grouped_melted$IC_group_plot[ap_grouped_melted$Group == h_lvl & ap_grouped_melted$IC_all == l_lvl] = c(IC_n[1],rep(NA,length(IC_n)-1))
  } 
}
ap_grouped_melted$Cor_plot = ap_grouped_melted$Cor
ap_grouped_melted$Cor_plot[is.na(ap_grouped_melted$IC_group_plot)]=NA
# ap_grouped_melted$IC_group_plot[ap_grouped_melted$IC_group_plot=="NA"]=NA
# Add correlation values only for one question (10 values, 1 per IC),for plotting purposes


ap_grouped_melted$IC_all = factor(ap_grouped_melted$IC_all, levels = paste0("IC",1:10))
ap_grouped_melted$ID = factor(ap_grouped_melted$ID, levels = ref_projections$ID)
ap_grouped_melted$projections_all = rep(ap_grouped_melted$projections[ap_grouped_melted$Group=="All"],nlevels(ap_grouped_melted$Group))
ap_grouped_melted$Gen_group = rep(c("Combind","Separate"),each = 18190)

for (IC in paste0("IC",1:10)){
  IC_combined = ap_grouped_melted$IC_group_plot[ap_grouped_melted$Group !="All" & ap_grouped_melted$IC_all == IC]
  IC_Cor = paste0(IC_combined[!is.na(IC_combined)],collapse = "-")
  ap_grouped_melted$Gen_IC_group[ap_grouped_melted$Group !="All" & ap_grouped_melted$IC_all == IC] = c(IC_Cor,rep(NA,1818))
  
  Cor_combined = ap_grouped_melted$Cor_plot[ap_grouped_melted$Group !="All" & ap_grouped_melted$IC_all == IC]
  Cor_val = paste0(Cor_combined[!is.na(Cor_combined)],collapse = "-")
  ap_grouped_melted$Gen_Cor_group[ap_grouped_melted$Group !="All" & ap_grouped_melted$IC_all == IC] = c(Cor_val,rep(NA,1818))
}
ggplot(data = ap_grouped_melted, aes(x=ID, y= Projections_group, group =Gen_group, color = Group))+
  geom_line(stat = "identity")+
  facet_grid(Gen_group~IC_all,scales = "free")+
  geom_text(aes(x= 1000,y=0.5,label = Gen_IC_group))+
  geom_text(aes(x= 1000,y=0.25,label = Gen_Cor_group))



# Check the differences in Projections types  --------------------------------------------
ref_projections = ICASSO_all$projections
fast_projections = ICASSO_all_fast$projections
fastsk_projections = ICASSO_all_fastsk$projections

df_cor(fastsk_projections, fastsk_projections, "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)

AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, "All_Fast", 10,10,est_threshold,Projection = TRUE)
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)

AllFAST_projections_ICs = simplify_mat(AllFAST_projections_Cor)
AllFASTSK_projections_ICs = simplify_mat(AllFASTSK_projections_Cor)

# Create a list of the separate diagnoses data frames
all_projections = rbind(hc_projections[,1:18],mdd_projections[,1:18], ptsd_projections[,1:18], tnp_projections[,1:18], gad_projections[,1:18])
ref_projections_melted = melt(ref_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
dia_projections_melted = melt(dia_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")

# Plot projections in a matrix --------------------------------------------

ref_projections = ICASSO_all$projections
fast_projections = ICASSO_all_fast$projections
fastsk_projections = ICASSO_all_fastsk$projections
hc_projections = ICASSO_hc$projections
gad_projections = ICASSO_gad$projections
mdd_projections = ICASSO_mdd$projections
ptsd_projections = ICASSO_ptsd$projections
tnp_projections = ICASSO_tnp$projections


# Create the matrix of correlations' estimates
est_threshold = 0.3
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, "All_Fast", 10,10,est_threshold,Projection = TRUE)
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)
HC_projections_Cor = df_cor(ref_projections, hc_projections, "HC", 10,10,est_threshold,Projection = TRUE)
GAD_projections_Cor = df_cor(ref_projections, gad_projections, "GAD", 10,10,est_threshold,Projection = TRUE)
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, "MDD", 10,10,est_threshold,Projection = TRUE)
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, "PTSD", 10,10,est_threshold,Projection = TRUE)
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, "TNP", 10,10,est_threshold,Projection = TRUE)


# Find the ICs in the each projections df that has the highest correlation
AllFAST_projections_ICs = simplify_mat(AllFAST_projections_Cor)
AllFASTSK_projections_ICs = simplify_mat(AllFASTSK_projections_Cor)
HC_projections_ICs = simplify_mat(HC_projections_Cor)
GAD_projections_ICs = simplify_mat(GAD_projections_Cor)
MDD_projections_ICs = simplify_mat(MDD_projections_Cor)
PTSD_projections_ICs = simplify_mat(PTSD_projections_Cor)
TNP_projections_ICs = simplify_mat(TNP_projections_Cor)



# Create a list of the separate diagnoses data frames
dia_projections = rbind(hc_projections[,1:18],mdd_projections[,1:18], ptsd_projections[,1:18], tnp_projections[,1:18], gad_projections[,1:18])
fast_projections_melted = melt(fast_projections[,1:18],id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
fastsk_projections_melted = melt(fastsk_projections[,1:18],id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
ref_projections_melted = melt(ref_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
dia_projections_melted = melt(dia_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")

projections_melted = rbind(ref_projections_melted,fast_projections_melted,fastsk_projections_melted,dia_projections_melted)
projections_melted$Grouping = rep(c("Combined","All_Fast","All_Fast_Sickit","Separate"),each = nrow(ref_projections_melted))
projections_melted$Grouping = factor(projections_melted$Grouping, levels = c("Combined","All_Fast","All_Fast_Sickit","Separate"))

# Do a diagnosis-specific sorting
for (Factor in levels(projections_melted$Diagnosis)){
  for (IC in paste0("IC",1:10)){
    for (Group in levels(projections_melted$Grouping)){
      projection_grouping = projections_melted$Projection[projections_melted$IC == IC &
                                                            projections_melted$Diagnosis == Factor &
                                                            projections_melted$Grouping == Group]
      projections_melted$Sort[projections_melted$IC == IC &
                                projections_melted$Diagnosis == Factor &
                                projections_melted$Grouping == Group] = match(projection_grouping, sort(projection_grouping))
    }
  }
}

# Do a sorting for each of the All-subject decomposition
for (IC in paste0("IC",1:10)){
  for (Group in levels(projections_melted$Grouping)){
    projection_grouping = projections_melted$Projection[projections_melted$IC == IC &
                                                          projections_melted$Grouping == Group]
    projections_melted$Sort_main[projections_melted$IC == IC &
                              projections_melted$Grouping == Group] = match(projection_grouping, sort(projection_grouping))
  }
}


# Some sort values are similar in two subjects (if they have exactly the same value). Fix this
for (k in 1:2){
  for (IC in paste0("IC",1:10)){
    for (Group in levels(projections_melted$Grouping)){
      for (Factor in levels(projections_melted$Diagnosis)){
        Sorting_Cond = projections_melted$Diagnosis== Factor & projections_melted$Grouping == Group & projections_melted$IC == IC
        for (i in 1:length(projections_melted$Sort[Sorting_Cond])){
          if (sum (projections_melted$Sort[Sorting_Cond] == projections_melted$Sort[Sorting_Cond][i]) >1){
            projections_melted$Sort[Sorting_Cond][i] = projections_melted$Sort[Sorting_Cond][i]+1
          }
        }
      }
    }
  }
}

for (k in 1:2){
  for (IC in paste0("IC",1:10)){
    for (Group in levels(projections_melted$Grouping)){
      Sorting_Cond = projections_melted$Grouping == Group & projections_melted$IC == IC
      for (i in 1:length(projections_melted$Sort_main[Sorting_Cond])){
        if (sum (projections_melted$Sort_main[Sorting_Cond] == projections_melted$Sort_main[Sorting_Cond][i]) >1){
          projections_melted$Sort_main[Sorting_Cond][i] = projections_melted$Sort_main[Sorting_Cond][i]+1
          print(c(IC,Group,projections_melted$Sort_main[Sorting_Cond][i]))
        }
      }
    }
  }
}

# add a variable to include only the close-to-zero values
projections_melted$Projection_zero = projections_melted$Projection
projections_melted$Projection_zero_main = projections_melted$Projection

for (Group in levels(factor(projections_melted$Grouping))){
  for (IC in levels(projections_melted$IC)){
    for (Diagnosis in levels(projections_melted$Diagnosis)){
      subgroup_cond = projections_melted$Grouping==Group & projections_melted$IC==IC & projections_melted$Diagnosis==Diagnosis
      pro_val = projections_melted$Projection[subgroup_cond]
      projections_melted$Projection_scaled[subgroup_cond] = pro_val/mean(abs(pro_val), na.rm = TRUE)
      min_value = min(abs(projections_melted$Projection_scaled[subgroup_cond]))
      projections_melted$Projection_zero[subgroup_cond][abs(projections_melted$Projection_scaled)[subgroup_cond] != min_value] = NA
    }
  }
}

projections_melted$Projection_zero[!is.na(projections_melted$Projection_zero)] = projections_melted$Projection_scaled[!is.na(projections_melted$Projection_zero)]

for (Group in levels(factor(projections_melted$Grouping))){
  for (IC in levels(projections_melted$IC)){
    subgroup_cond = projections_melted$Grouping==Group & projections_melted$IC==IC
    pro_val = projections_melted$Projection[subgroup_cond]
    projections_melted$Projection_scaled_main[subgroup_cond] = pro_val/mean(abs(pro_val), na.rm = TRUE)
    min_value = min(abs(projections_melted$Projection_scaled_main[subgroup_cond]))
    projections_melted$Projection_zero_main[subgroup_cond][abs(projections_melted$Projection_scaled_main)[subgroup_cond] != min_value] = NA
  }
}

# Change the IC designation to match that of the highest correlation for each IC from the all-subjects decompoisition

for (ICn in 1:10){
  for (Fn in 1:nlevels(projections_melted$Diagnosis)){
    Factor = levels(projections_melted$Diagnosis)[Fn]
    Pro_ICs = list(GAD_projections_ICs, HC_projections_ICs, MDD_projections_ICs, PTSD_projections_ICs, TNP_projections_ICs)[[Fn]]
    projections_melted$IC_group[projections_melted$Diagnosis == Factor & projections_melted$Grouping == "Separate"& projections_melted$IC == paste0("IC",ICn)] = Pro_ICs[2,ICn]
    projections_melted$Cor_arranged[projections_melted$Diagnosis == Factor & projections_melted$Grouping == "Separate"& projections_melted$IC == paste0("IC",ICn)] = Pro_ICs[3,ICn]
  }
  projections_melted$IC_group[projections_melted$Grouping == "All_Fast"& projections_melted$IC == paste0("IC",ICn)] = AllFAST_projections_ICs[2,ICn]
  projections_melted$Cor_arranged[projections_melted$Grouping == "All_Fast"& projections_melted$IC == paste0("IC",ICn)] = AllFAST_projections_ICs[3,ICn]
  
  projections_melted$IC_group[projections_melted$Grouping == "All_Fast_Sickit"& projections_melted$IC == paste0("IC",ICn)] = AllFASTSK_projections_ICs[2,ICn]
  projections_melted$Cor_arranged[projections_melted$Grouping == "All_Fast_Sickit"& projections_melted$IC == paste0("IC",ICn)] = AllFASTSK_projections_ICs[3,ICn]
}



for (IC in paste0("IC", 1:10)) {
  for (Factor in levels(projections_melted$Diagnosis)) {
    cond_all_ic = (projections_melted$IC == IC & projections_melted$Diagnosis == Factor & projections_melted$Grouping == "Separate")
    new_ic = projections_melted$IC_group[cond_all_ic][1]
    cond_new_ic = (projections_melted$IC == new_ic & projections_melted$Diagnosis == Factor & projections_melted$Grouping == "Separate")
    cond_new_ic_all = (projections_melted$IC == IC & projections_melted$Diagnosis == Factor & projections_melted$Grouping == "Combined")
    
    projections_melted$Projections_group[cond_all_ic] = projections_melted$Projection[cond_new_ic]
    projections_melted$Projections_all[cond_all_ic] = projections_melted$Projection[cond_new_ic_all] # Projection value from the Combined grouping that has most correlation with this projection in the "Separate" Grouping
    projections_melted$Projections_scaled_group[cond_all_ic] = projections_melted$Projection_scaled[cond_new_ic]
    projections_melted$Sort_group[cond_all_ic] = projections_melted$Sort[cond_new_ic]
    projections_melted$Projection_zero_group[cond_all_ic] = projections_melted$Projection_zero[cond_new_ic]
  }
  cond_ic_fast = (projections_melted$IC == IC & projections_melted$Grouping == "All_Fast")
  cond_ic_fastsk = (projections_melted$IC == IC & projections_melted$Grouping == "All_Fast_Sickit")
  new_ic_fast = projections_melted$IC_group[cond_ic_fast][1]
  new_ic_fastsk = projections_melted$IC_group[cond_ic_fastsk][1]
  cond_new_ic_fast = (projections_melted$IC == new_ic_fast & projections_melted$Grouping == "All_Fast")
  cond_new_ic_fastsk = (projections_melted$IC == new_ic_fastsk & projections_melted$Grouping == "All_Fast_Sickit")
  cond_new_ic_fast_all = (projections_melted$IC == IC & projections_melted$Grouping == "Combined")
  cond_new_ic_fastsk_all = (projections_melted$IC == IC & projections_melted$Grouping == "Combined")
  
  
  projections_melted$Projections_group[cond_ic_fast] = projections_melted$Projection[cond_new_ic_fast]
  projections_melted$Projections_group[cond_ic_fastsk] = projections_melted$Projection[cond_new_ic_fastsk]
  
  projections_melted$Projections_all[cond_ic_fast] = projections_melted$Projection[cond_new_ic_fast_all]
  projections_melted$Projections_all[cond_ic_fastsk] = projections_melted$Projection[cond_new_ic_fastsk_all]
  
  projections_melted$Projections_scaled_group_main[cond_ic_fast] = projections_melted$Projection_scaled_main[cond_new_ic_fast]
  projections_melted$Projections_scaled_group_main[cond_ic_fastsk] = projections_melted$Projection_scaled_main[cond_new_ic_fastsk]

  projections_melted$Sort_main_group[cond_ic_fast] = projections_melted$Sort_main[cond_new_ic_fast]
  projections_melted$Sort_main_group[cond_ic_fastsk] = projections_melted$Sort_main[cond_new_ic_fastsk]
  
  projections_melted$Sort_main_diagnosis[cond_ic_fast] = projections_melted$Sort[cond_new_ic_fast]
  projections_melted$Sort_main_diagnosis[cond_ic_fastsk] = projections_melted$Sort[cond_new_ic_fastsk]

  projections_melted$Projection_zero_group_main[cond_ic_fast] = projections_melted$Projection_zero_main[cond_new_ic_fast]
  projections_melted$Projection_zero_group_main[cond_ic_fastsk] = projections_melted$Projection_zero_main[cond_new_ic_fastsk]
}

var_label(projections_melted$Projections_all) <- "Projections of this 'IC' in the 'Combined' Grouping"
var_label(projections_melted$Projections_group) <- "Projections of the 'IC_group' in this Grouping"

projections_melted$Projections_all[projections_melted$Grouping == "Combined"] = projections_melted$Projection[projections_melted$Grouping == "Combined"]
projections_melted$Projections_scaled_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection_scaled[projections_melted$Grouping == "Combined"]
projections_melted$Projections_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection[projections_melted$Grouping == "Combined"]
projections_melted$Sort_group[projections_melted$Grouping == "Combined"] = projections_melted$Sort[projections_melted$Grouping == "Combined"]
projections_melted$Sort_main_group[projections_melted$Grouping == "Combined"] = projections_melted$Sort_main[projections_melted$Grouping == "Combined"]
projections_melted$Sort_main_diagnosis[projections_melted$Grouping == "Combined"] = projections_melted$Sort[projections_melted$Grouping == "Combined"]
projections_melted$Sort_main_diagnosis[projections_melted$Grouping == "Separate"] = projections_melted$Sort_group[projections_melted$Grouping == "Separate"]

#projections_melted$Projections_scaled_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA
#projections_melted$Projections_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA
#projections_melted$Sort_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA

projections_melted$Projections_scaled_group_main[projections_melted$Grouping == "Combined"] = projections_melted$Projection_scaled_main[projections_melted$Grouping == "Combined"]
projections_melted$Projections_scaled_group_main[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "All_Fast"])] = NA
projections_melted$Projections_scaled_group_main[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "All_Fast"])] = NA
projections_melted$Sort_main_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "All_Fast"])] = NA
projections_melted$Sort_main_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "All_Fast_Sickit"])] = NA

projections_melted$Projection_zero_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection_zero[projections_melted$Grouping == "Combined"]


# remove all IC designations, except for the first one, for plotting purposes
# find the indices to keep
Fast_inds = 18190 + seq(from = 1, to=1819*10, by = 1819)
FastSK_inds = 18190*2 + seq(from = 1, to=1819*10, by = 1819)
HC_indices = 18190*3 + seq(from = 1, to=1819*10, by = 1819)
MDD_indices = HC_indices + nrow(projections_melted[projections_melted$Diagnosis == "HC",])/40
PTSD_indices = MDD_indices + nrow(projections_melted[projections_melted$Diagnosis == "MDD",])/40
TNP_indices = PTSD_indices + nrow(projections_melted[projections_melted$Diagnosis == "PTSD",])/40
GAD_indices = TNP_indices + nrow(projections_melted[projections_melted$Diagnosis == "TNP",])/40
Inc_inds = c(Fast_inds, FastSK_inds, HC_indices, MDD_indices, PTSD_indices, TNP_indices, GAD_indices)
projections_melted$IC_group_plot = projections_melted$IC_group
projections_melted$Cor_arranged_plot = projections_melted$Cor_arranged
projections_melted$IC_group_plot[-Inc_inds] = NA
projections_melted$Cor_arranged_plot[-Inc_inds] = NA

# make sure that IC_group is not empty when Grouping == "Combined"
projections_melted$IC_group[projections_melted$Grouping == "Combined"] = as.character(projections_melted$IC[projections_melted$Grouping == "Combined"])
projections_melted$IC_group = factor(projections_melted$IC_group, levels = paste0("IC",1:10))

# Finally make a df for diagnoses alone, and one for the three methods of IC decomposition (Infomax, fast, and fast_sickit)
projections_melted_dia = projections_melted[(projections_melted$Grouping == "Combined" | projections_melted$Grouping == "Separate"),]
projections_melted_all = projections_melted[!projections_melted$Grouping == "Separate",]

# All, non arranged
ggplot(data = projections_melted, mapping = aes(x = Sort,y= IC))+
  geom_raster(aes(fill = Projection_scaled))+
  scale_fill_gradient2(high = "#D40000",mid = "white", low = "#08457E", midpoint = 0)+
  geom_point(data = projections_melted[!is.na(projections_melted$Projection_zero),], aes(color = Projection_zero))+
  scale_color_gradient2(high = "black", mid = "black", low = "black")+
  facet_grid(Grouping~Diagnosis, scale = "free")+
  MinimalTheme

# Significant correlation, arranged
ggplot(data = projections_melted_dia, mapping = aes(x = Sort_group,y= Grouping))+
  geom_raster(aes(fill = Projections_group))+
  scale_fill_gradient2(high = "#a10000",mid = "white", low = "#002138", midpoint = 0, name = "Projection")+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  geom_point(data = projections_melted_dia[!is.na(projections_melted_dia$Projection_zero_group),], aes(color = Projection_zero_group))+
  scale_color_gradient2(high = "black", mid = "black", low = "black", name = "Zero", labels = c())+
  facet_grid(IC~Diagnosis, scale = "free")+
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "- Scaled by dividing on the mean of absolute values
          - Each IC in the \'Separate\' decomposition represents the most correlating IC with the \'Combined\' decomposition
          - The Dot represents the zero value for that IC")+
  MinimalTheme

projections_melted_dia$Sort_group_comb = rep(projections_melted_dia$Sort_group[projections_melted_dia$Grouping == "Combined"],2)
projections_melted_dia$Diagnosis = factor(projections_melted_dia$Diagnosis, levels = c("HC","GAD","MDD","PTSD","TNP"))
ggplot(data = projections_melted_dia, mapping = aes(x = Sort_group_comb,y= Grouping))+
  geom_raster(aes(fill = Projections_scaled_group))+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradient2(high = "#a10000",mid = "white", low = "#002138", midpoint = 0, name = "Scaled\nProjection")+
  facet_grid(IC~Diagnosis, scale = "free")+
  geom_text(aes(x= 3,y=2.25,label = IC_group_plot))+
  geom_text(aes(x= 3,y=1.75,label = Cor_arranged_plot))+
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "- Scaled by dividing on the mean of absolute values
          - Arranged by values of the \'Combined\' decompoistion
          - Each IC in the \'Separate\' decomposition represents the most correlating IC with the \'Combined\' decomposition
          - The name of the most correlating IC from the separate group and its value of correlation are on the left")+
  MinimalTheme


# Todo: make an all sorting variable, also make projections_group not NA for All
projections_melted_all$IC_group_Fast_Sickit = NA
projections_melted_all$IC_group_Fast = NA
projections_melted_all$IC_group_Fast_Sickit[projections_melted_all$Grouping != "All_Fast"] = projections_melted_all$IC_group_plot[projections_melted_all$Grouping != "All_Fast"]
projections_melted_all$IC_group_Fast[projections_melted_all$Grouping != "All_Fast_Sickit"] = projections_melted_all$IC_group_plot[projections_melted_all$Grouping != "All_Fast_Sickit"]

projections_melted_all$Cor_arranged_Fast_Sickit = NA
projections_melted_all$Cor_arranged_Fast = NA
projections_melted_all$Cor_arranged_Fast_Sickit[projections_melted_all$Grouping != "All_Fast"] = projections_melted_all$Cor_arranged_plot[projections_melted_all$Grouping != "All_Fast"]
projections_melted_all$Cor_arranged_Fast[projections_melted_all$Grouping != "All_Fast_Sickit"] = projections_melted_all$Cor_arranged_plot[projections_melted_all$Grouping != "All_Fast_Sickit"]

projections_melted_all$Sort_main_group_comb = rep(projections_melted_all$Sort_main_group[projections_melted_all$Grouping == "Combined"],3)

ggplot(data = projections_melted_all, mapping = aes(x = Sort_main_group_comb,y= Grouping))+
  geom_raster(aes(fill = Projections_scaled_group_main))+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradient2(high = "#a10000",mid = "white", low = "#002138", midpoint = 0, name = "Scaled\nProjection")+
  facet_grid(IC~., scale = "free")+
  geom_text(aes(x= 3,y=2.25,label = IC_group_Fast))+
  geom_text(aes(x= 3,y=3.25,label = IC_group_Fast_Sickit))+
  geom_text(aes(x= 3,y=1.75,label = Cor_arranged_Fast))+
  geom_text(aes(x= 3,y=2.75,label = Cor_arranged_Fast_Sickit))+
  ggtitle("Projection Values Per Decomposition Type")+
  MinimalTheme



# Plot the correlation figure
projections_melted_cor = projections_melted
projections_melted_cor$Projection_comb[projections_melted_cor$Grouping == "Separate"] = projections_melted_cor$Projection[projections_melted_cor$Grouping == "Combined"]
projections_melted_cor$Projections_group_comb[projections_melted_cor$Grouping == "Separate"] = projections_melted_cor$Projections_group[projections_melted_cor$Grouping == "Combined"]
projections_melted_cor = projections_melted_cor[projections_melted_cor$Grouping == "Separate",]

# remove the variables that are no longer needed
projections_melted_cor$Projection_zero = NULL;projections_melted_cor$Projection_zero_group = NULL; projections_melted_cor$Sort_main = NULL; projections_melted_cor$Projection_scaled = NULL
projections_melted_cor$Projections_scaled_group = NULL; projections_melted_cor$Sort = NULL ; projections_melted_cor$Sort_group = NULL

Correlate(projections_melted_cor[projections_melted_cor$Diagnosis == "HC",],Var1 = "Projections_group",Var2 = "Projections_group_comb", Factor = "IC")

ggplot(data = projections_melted_cor, mapping = aes(x = Projections_group,y= Projections_group_comb))+
  geom_point()+
  facet_grid(IC~Diagnosis, scale = "free")+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "Scaled by dividing on the mean of absolute values")+
  MinimalTheme+
  geom_smooth(method="lm", size=1, se = FALSE)+
  geom_text(aes(x= -0.1,y=0.25,label = IC_group))+
  geom_text(aes(x= -0.1,y=-0.25,label = Cor_arranged_plot))+
  facet_grid(IC~Diagnosis, scale = "free")
  

# Additional part: Present the data the way Juergen asked for. Use same data, just make it into a square
# First, add a column to represent rows of the figure, basically a number close to the square root of subject count (43)
# create a new sorting factor, which is basically the IDs arranged in the order of diagnosis

projections_melted$ID = factor(projections_melted$ID, levels = projections_melted$ID[projections_melted$Grouping == "Combined" & projections_melted$IC == "IC1"])

for (IC in paste0("IC",1:10)){
#  for (Diagnosis in levels(projections_melted$Diagnosis)){
    for (Grouping in levels(projections_melted$Grouping)){
      cond = projections_melted$Grouping == Grouping & projections_melted$IC == IC
      max_num = max(as.numeric(projections_melted$ID[cond]))
      mat_seq = c(seq(from = 1, to = max_num, by = 45),max_num+1)
      for (num in 2:length(mat_seq)){
        projections_melted$row_seq[cond & as.numeric(projections_melted$ID)>mat_seq[num-1]-1 & as.numeric(projections_melted$ID)<(mat_seq[num])] = num-1
        
        #create a matched ID, that arranges the IDs by projection value
        projections_melted$Sort_mod[cond&projections_melted$Diagnosis =="HC"] = projections_melted$Sort_main_diagnosis[cond&projections_melted$Diagnosis =="HC"]
        projections_melted$Sort_mod[cond&projections_melted$Diagnosis =="MDD"] = projections_melted$Sort_main_diagnosis[cond&projections_melted$Diagnosis =="MDD"]+1202
        projections_melted$Sort_mod[cond&projections_melted$Diagnosis =="PTSD"] = projections_melted$Sort_main_diagnosis[cond&projections_melted$Diagnosis =="PTSD"]+1654
        projections_melted$Sort_mod[cond&projections_melted$Diagnosis =="TNP"] = projections_melted$Sort_main_diagnosis[cond&projections_melted$Diagnosis =="TNP"]+1711
        projections_melted$Sort_mod[cond&projections_melted$Diagnosis =="GAD"] = projections_melted$Sort_main_diagnosis[cond&projections_melted$Diagnosis =="GAD"]+1789
        projections_melted$row_seq_arr[cond & projections_melted$Sort_mod[cond]>mat_seq[num-1]-1 & projections_melted$Sort_mod[cond]<(mat_seq[num])] = num-1
      }
    }
#  }
}


# now the columns
projections_melted$col_seq = as.numeric(projections_melted$ID) - (45*(projections_melted$row_seq-1))
projections_melted$col_seq_arr = projections_melted$Sort_mod - (45*(projections_melted$row_seq_arr-1))

# I will create a separating line for each of the diagnoses
for (i in 1:5){
  lvl = c("HC","MDD","PTSD","TNP","GAD")[i]
  projections_melted$row_seq[projections_melted$Diagnosis == lvl] = projections_melted$row_seq[projections_melted$Diagnosis == lvl]+i*3
  projections_melted$row_seq_arr[projections_melted$Diagnosis == lvl] = projections_melted$row_seq_arr[projections_melted$Diagnosis == lvl]+i*3
}

plot_matrix <- function(Grouping, Group, scale_values=FALSE){
  pro_data = projections_melted[(projections_melted$Grouping == Grouping),]
  pro_data$Projections_group[is.na(pro_data$Projections_group)] = 0
  
  # make all IC projections of 0 mean and sd of 1
  if (scale_values){
    for (IC in paste0("IC",1:10)){
      pro_data$Projections_group[pro_data$IC == IC] = (pro_data$Projections_group[pro_data$IC == IC]-
                                                         mean(pro_data$Projections_group[pro_data$IC == IC]))/
        sd(pro_data$Projections_group[pro_data$IC == IC])
      
      
      if(sum(is.na(pro_data$Projections_group[pro_data$IC == IC])) < (182)){
        estimate = ifelse(cor.test(pro_data$Projections_group[pro_data$IC == IC], pro_data$Projections_all[pro_data$IC == IC])$estimate>0,1,-1)
        pro_data$Projections_group[pro_data$IC == IC] = pro_data$Projections_group[pro_data$IC == IC]*estimate
        
        if (estimate==-1){
          for (Diagnosis in c("HC","MDD","PTSD","TNP","GAD")){
            pro_data$Sort_main_diagnosis[pro_data$IC == IC & pro_data$Diagnosis ==Diagnosis] = 
              min(pro_data$Sort_main_diagnosis[pro_data$IC == IC& pro_data$Diagnosis ==Diagnosis],na.rm = TRUE)+
              max(pro_data$Sort_main_diagnosis[pro_data$IC == IC& pro_data$Diagnosis ==Diagnosis],na.rm = TRUE)-
              pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis ==Diagnosis]
          }
          
          max_num = max(as.numeric(pro_data$ID[pro_data$IC == IC]))
          mat_seq = c(seq(from = 1, to = max_num, by = 45),max_num+1)
          for (num in 2:length(mat_seq)){
            #create a matched ID, that arranges the IDs by projection value
            pro_data$Sort_mod[pro_data$IC == IC&pro_data$Diagnosis =="HC"] = pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis =="HC"]
            pro_data$Sort_mod[pro_data$IC == IC&pro_data$Diagnosis =="MDD"] = pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis =="MDD"]+1202
            pro_data$Sort_mod[pro_data$IC == IC&pro_data$Diagnosis =="PTSD"] = pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis =="PTSD"]+1654
            pro_data$Sort_mod[pro_data$IC == IC&pro_data$Diagnosis =="TNP"] = pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis =="TNP"]+1711
            pro_data$Sort_mod[pro_data$IC == IC&pro_data$Diagnosis =="GAD"] = pro_data$Sort_main_diagnosis[pro_data$IC == IC&pro_data$Diagnosis =="GAD"]+1789
            pro_data$row_seq_arr[pro_data$IC == IC & pro_data$Sort_mod[pro_data$IC == IC]>mat_seq[num-1]-1 & pro_data$Sort_mod[pro_data$IC == IC]<(mat_seq[num])] = num-1
          }
          pro_data$col_seq_arr[pro_data$IC == IC] = pro_data$Sort_mod[pro_data$IC == IC] - (45*(pro_data$row_seq_arr[pro_data$IC == IC]-1))
          for (i in 1:5){
            lvl = c("HC","MDD","PTSD","TNP","GAD")[i]
            pro_data$row_seq_arr[pro_data$IC == IC &pro_data$Diagnosis == lvl] = pro_data$row_seq_arr[pro_data$IC == IC &pro_data$Diagnosis == lvl]+i*3
          }
        }
      }
    }
  }
  pro_data$Projections_group[is.na(pro_data$Projections_group)] = 0
  
  
  sources_data = sources_melted_all[(sources_melted_all$Group == Group),]
  # make all IC sources of 0 mean and sd of 1
  if (scale_values){
    for (IC in paste0("IC",1:10)){
      sources_data$Sources_group[sources_data$IC_all == IC] = (sources_data$Sources_group[sources_data$IC_all == IC]-
                                                                 mean(sources_data$Sources_group[sources_data$IC_all == IC]))/
        sd(sources_data$Sources_group[sources_data$IC_all == IC])
      if(sum(is.na(sources_data$Sources_group[sources_data$IC_all == IC])) < (10)){
        estimate = ifelse(cor.test(sources_data$Sources_group[sources_data$IC_all == IC], sources_data$Sources_all[sources_data$IC_all == IC])$estimate>0,1,-1)
        sources_data$Sources_group[sources_data$IC_all == IC] = sources_data$Sources_group[sources_data$IC_all == IC]*estimate
      }
    }
  }
  sources_data$Question = factor(sources_data$Question, levels = paste0("Q",1:100))
  
  #make sure that IC_group is not empty (which is the case when sources_data$Group == "All")
  sources_data$IC_group[sources_data$Group == Group] = as.character(sources_data$IC_all[sources_data$Group == Group])
  sources_data$IC_group = factor(sources_data$IC_group, levels = paste0("IC",1:10))
  
  # make sure there are no NAs
  sources_data$Sources_group[is.na(sources_data$Sources_group)] = 0
  # for aesthetic purposes, I will remove those above 4
  
  if (scale_values){
    sources_data$Sources_group[sources_data$Sources_group>4] = 4
    sources_data$Sources_group[sources_data$Sources_group< (-4)] = -4
    pro_data$Projections_group[pro_data$Projections_group>4] = 4
    pro_data$Projections_group[pro_data$Projections_group< (-4)] = -4
  }
  
  # create labels from the Diagnosis column, only the first of each group will be included
  #First, create the True/False condition
  dia_cond = unlist(sapply(c("HC","MDD","PTSD","TNP","GAD"), function(x) c(TRUE,rep(FALSE,sum(pro_data$Diagnosis == x)/10-1))))
  
  pro_labels = pro_data$Diagnosis
  pro_labels_arr = pro_data$Diagnosis
  
  pro_labels[!dia_cond] = NA
  pro_labels_arr[pro_data$Sort_main_diagnosis != 1]= NA
  
  
  g1_1 = ggplot(pro_data[pro_data$IC %in% paste0("IC",1:5),], aes(x= col_seq, y=row_seq))+
    geom_raster(aes(fill = Projections_group))+
    scale_fill_gradientn(colors = c("#002138", "#FFFFFF","#a10000"), name = "Proj")+
    geom_text(aes(x=0,label = pro_labels[pro_data$IC %in% paste0("IC",1:5)]),size=3,hjust = 1,vjust=0)+
    scale_y_discrete("")+
    scale_x_continuous("", limits = c(-4.5,45), breaks =NULL)+
    facet_grid(IC~.,scales = "fixed")+
    TypicalTheme#+theme(legend.position = "none")
  
  
  g1_2 = ggplot(pro_data[pro_data$IC %in% paste0("IC",6:10),], aes(x= col_seq, y=row_seq))+
    geom_raster(aes(fill = Projections_group))+
    scale_fill_gradientn(colors = c("#002138", "#FFFFFF","#a10000"), name = "Proj")+
    geom_text(aes(x=0,label = pro_labels[pro_data$IC %in% paste0("IC",6:10)]),size=3,hjust = 1,vjust=0)+
    scale_y_discrete("")+
    scale_x_continuous("", limits = c(-4.5,45), breaks =NULL)+
    facet_grid(IC~.,scales = "fixed")+
    TypicalTheme#+theme(legend.position = "none")
  
  g1_1_arr = ggplot(pro_data[pro_data$IC %in% paste0("IC",1:5),], aes(x= col_seq_arr, y=row_seq_arr))+
    geom_raster(aes(fill = Projections_group))+
    scale_fill_gradientn(colors = c("#002138", "#FFFFFF","#a10000"),  name = "Proj")+
    geom_text(aes(x=0,label = pro_labels_arr[pro_data$IC %in% paste0("IC",1:5)]),size=3,hjust = 1,vjust=0)+
    facet_grid(IC~.,scales = "fixed")+
    scale_y_discrete("")+
    scale_x_continuous("", limits = c(-4.5,45), breaks = NULL)+
    TypicalTheme#+ theme(legend.position = "none")
  
  g1_2_arr = ggplot(pro_data[pro_data$IC %in% paste0("IC",6:10),], aes(x= col_seq_arr, y=row_seq_arr))+
    geom_raster(aes(fill = Projections_group))+
    scale_fill_gradientn(colors = c("#002138", "#FFFFFF","#a10000"),name = "Proj")+
    geom_text(aes(x=0,label = pro_labels_arr[pro_data$IC %in% paste0("IC",6:10)]),size=3,hjust = 1,vjust=0)+
    facet_grid(IC~.,scales = "fixed")+
    scale_y_discrete("")+
    scale_x_continuous("", limits = c(-4.5,45), breaks = NULL)+
    TypicalTheme#+theme(legend.position = "none")
  
  g2_1 = ggplot(data = sources_data[sources_data$IC_all %in% paste0("IC",1:5),], mapping = aes(x = Question,y= Sources_group))+
    geom_bar(stat = "identity")+
    scale_x_discrete("Question")+
    scale_y_continuous("")+
    scale_fill_gradientn(colors = c("#02422f","#1B4D3E","#FFFFFF","#85015d","#69026b"), name = "Scores")+
    MinimalTheme+
    facet_grid(IC_all~.,scales = "fixed")+
    theme(axis.text.x = element_text(angle = 90,size=5))
  
  g2_2 = ggplot(data = sources_data[sources_data$IC_all %in% paste0("IC",6:10),], mapping = aes(x = Question,y= Sources_group))+
    geom_bar(stat = "identity")+
    scale_x_discrete("Question")+
    scale_y_continuous("")+
    scale_fill_gradientn(colors = c("#02422f","#1B4D3E","#FFFFFF","#85015d","#69026b"), name = "Scores")+
    MinimalTheme+
    facet_grid(IC_all~.,scales = "fixed")+
    theme(axis.text.x = element_text(angle = 90,size=5))
  
  multiplot(g1_1,g2_1, g1_2, g2_2, cols = 4, layout = matrix(c(1,1,2,2,2,3,3,4,4,4), nrow = 1))
  multiplot(g1_1_arr,g2_1, g1_2_arr, g2_2, cols = 4, layout = matrix(c(1,1,2,2,2,3,3,4,4,4), nrow = 1))
}
plot_matrix("Combined","All", scale_values = FALSE)
plot_matrix("All_Fast","All_fast", scale_values = FALSE)
plot_matrix("All_Fast_Sickit","All_fastsk", scale_values =FALSE)
plot_matrix("Separate","All", scale_values = FALSE)

ggplot(all_sources_melted_dia, aes(x=Question, y = Sources_group))+
  geom_bar(stat = "identity")+
  facet_grid(Group~IC_group)+
  MinimalTheme

### Find the Differences between ICs
#First, make everything range from -1 to 1
for(IC in paste0("IC",1:10)){
  print(IC)
  Grouping = "All_Fast_Sickit"
  fast_projections = projections_melted$Projections_group[projections_melted$Grouping == "All_Fast" & projections_melted$IC == IC]
  fastsk_projections = projections_melted$Projections_group[projections_melted$Grouping == "All_Fast_Sickit" & projections_melted$IC == IC]
  
  all_projections = projections_melted$Projections_all[projections_melted$Grouping == "All_Fast" & projections_melted$IC == IC]
  print(mean(abs((fast_projections-all_projections)/fast_projections)))
  if (IC!="IC2"){
    # if the correlation is above 0, define the estimate as 1, otherwise as -1,
    # this is done to reverse the values of negative correlations
    estimate_fast = ifelse(cor.test(fast_projections, all_projections)$estimate>0,1,-1)
    estimate_fastsk = ifelse(cor.test(fastsk_projections, all_projections)$estimate>0,1,-1)
    
    # convert everything to have a mean of 0 and sd of 1
    all_projections = (all_projections-mean(all_projections))/sd(all_projections)
    
    fast_projections = estimate_fast*(fast_projections-mean(fast_projections))/sd(fast_projections)
    fastsk_projections = estimate_fastsk*(fastsk_projections-mean(fastsk_projections))/sd(fastsk_projections)
    
    #(fast_projections-all_projections)/all_projections
    print(mean(abs((fast_projections-all_projections)/fast_projections)))
  }
  
}

dd = data.frame(fast = fast_projections, fastsk = fastsk_projections, all = all_projections, fast_dif = (all_projections-fast_projections))
ggplot(dd, aes(x= all,y = fast))+
  geom_point()



# Now for the sources
IC = "IC1"
Grouping = "All_fast"
group_sources = sources_melted_all$Sources_group[sources_melted_all$Group == Grouping & sources_melted_all$IC_all == IC]
all_sources = sources_melted_all$Sources_all[sources_melted_all$Group == Grouping & sources_melted_all$IC_all == IC]

# if the correlation is above 0, define the estimate as 1, otherwise as -1,
# this is done to reverse the values of negative correlations
estimate = ifelse(cor.test(group_sources, all_sources)$estimate>0,1,-1)

# convert everything to have a mean of 0 and sd of 1
all_sources = (all_sources-mean(all_sources))/sd(all_sources)
group_sources = estimate*(group_sources-mean(group_sources))/sd(group_sources)
mean(abs(all_sources))
mean(abs(all_sources-group_sources))
dd = data.frame(group = group_sources, all = all_sources)
ggplot(dd, aes(x= all,y = group))+
  scale_fill_gradientn(colors = c("#002138", "#FFFFFF","#a10000"), name = "Proj")

# Combine sources and projections analyses --------------------------------


ref_sources = ICASSO_all$sources
fast_sources = ICASSO_all_fast$sources
fastsk_sources = ICASSO_all_fastsk$sources
hc_sources = ICASSO_hc$sources
gad_sources = ICASSO_gad$sources
mdd_sources = ICASSO_mdd$sources
ptsd_sources = ICASSO_ptsd$sources
tnp_sources = ICASSO_tnp$sources


# Create the matrix of correlations' estimates
est_threshold = 0.3
HC_sources_Cor = df_cor(ref_sources, hc_sources, "HC", 10,10,est_threshold,Projection = FALSE)
GAD_sources_Cor = df_cor(ref_sources, gad_sources, "GAD", 10,10,est_threshold,Projection = FALSE)
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, "MDD", 10,10,est_threshold,Projection = FALSE)
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, "PTSD", 10,10,est_threshold,Projection = FALSE)
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, "TNP", 10,10,est_threshold,Projection = FALSE)


HC_projections_Cor = df_cor(ref_projections, hc_projections, "HC", 10,10,est_threshold,Projection = TRUE)
GAD_projections_Cor = df_cor(ref_projections, gad_projections, "GAD", 10,10,est_threshold,Projection = TRUE)
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, "MDD", 10,10,est_threshold,Projection = TRUE)
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, "PTSD", 10,10,est_threshold,Projection = TRUE)
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, "TNP", 10,10,est_threshold,Projection = TRUE)

HC_sources_ICs = simplify_mat(HC_sources_Cor)
GAD_sources_ICs = simplify_mat(GAD_sources_Cor)
MDD_sources_ICs = simplify_mat(MDD_sources_Cor)
PTSD_sources_ICs = simplify_mat(PTSD_sources_Cor)
TNP_sources_ICs = simplify_mat(TNP_sources_Cor)


HC_projections_ICs = simplify_mat(HC_projections_Cor)
GAD_projections_ICs = simplify_mat(GAD_projections_Cor)
MDD_projections_ICs = simplify_mat(MDD_projections_Cor)
PTSD_projections_ICs = simplify_mat(PTSD_projections_Cor)
TNP_projections_ICs = simplify_mat(TNP_projections_Cor)

All_ICs = data.frame(IC_all = rep(paste0("IC",1:10),10),
                     Diagnosis = rep(c("HC","GAD","MDD","PTSD", "TNP"), each = 10),
                     Type = rep(c("Sources","Projections"),each = 50),
                     IC=c(HC_sources_ICs[2,],GAD_sources_ICs[2,],MDD_sources_ICs[2,],PTSD_sources_ICs[2,],TNP_sources_ICs[2,],
                          HC_projections_ICs[2,],GAD_projections_ICs[2,],MDD_projections_ICs[2,],PTSD_projections_ICs[2,],TNP_projections_ICs[2,]),
                     Cor=as.numeric(c(HC_sources_ICs[3,],GAD_sources_ICs[3,],MDD_sources_ICs[3,],PTSD_sources_ICs[3,],TNP_sources_ICs[3,],
                                      HC_projections_ICs[3,],GAD_projections_ICs[3,],MDD_projections_ICs[3,],PTSD_projections_ICs[3,],TNP_projections_ICs[3,])))
All_ICs$IC = as.character(All_ICs$IC)
All_ICs$IC[rep(All_ICs$IC[All_ICs$Type == "Sources"] != All_ICs$IC[All_ICs$Type == "Projections"],2)]=NA
All_ICs$IC[rep(is.na(All_ICs$IC[All_ICs$Type == "Sources"]),2)]=NA
All_ICs$IC[rep(is.na(All_ICs$IC[All_ICs$Type == "Projections"]),2)]=NA
All_ICs$IC = factor(All_ICs$IC,levels = paste0("IC",1:10))

All_ICs$Cor[is.na(All_ICs$IC)]=NA
All_ICs$IC_all = factor(All_ICs$IC_all,levels = paste0("IC",1:10))

All_ICs_wide = reshape(All_ICs[1:50,1:4], idvar = c("Diagnosis","Type"), timevar = "IC_all", direction = "wide")
colnames(All_ICs_wide) = c("Diagnosis", "Type", paste0("IC",1:10))

ggplot(data = All_ICs[All_ICs$Type == "Sources",], mapping = aes(x = 1,y= Type))+
  geom_text(aes(y=1,label = IC))+
  geom_text(aes(y=0,label = Cor))+
  facet_grid(IC_all~Diagnosis)+
  TypicalTheme+
  scale_y_discrete("IC")+
  scale_x_discrete("Diagnosis")+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(), axis.ticks.y=element_blank())+
  ggtitle("ICs with Most Correlations from the Separate-Group Decomposition", subtitle = "All Correlations are with a p-value of <0.0005 and Spearman rho >0.3")


# Choose Questions --------------------------------------------------------


ref_sources = ICASSO_all$sources
fast_sources = ICASSO_all_fast$sources
fastsk_sources = ICASSO_all_fastsk$sources
hc_sources = ICASSO_hc$sources
gad_sources = ICASSO_gad$sources
mdd_sources = ICASSO_mdd$sources
ptsd_sources = ICASSO_ptsd$sources
tnp_sources = ICASSO_tnp$sources

# Create the matrix of correlations' estimates
est_threshold = 0.3
HC_sources_Cor = df_cor(ref_sources, hc_sources, "HC", 10,10,est_threshold,Projection = FALSE)
GAD_sources_Cor = df_cor(ref_sources, gad_sources, "GAD", 10,10,est_threshold,Projection = FALSE)
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, "MDD", 10,10,est_threshold,Projection = FALSE)
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, "PTSD", 10,10,est_threshold,Projection = FALSE)
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, "TNP", 10,10,est_threshold,Projection = FALSE)


HC_sources_ICs = simplify_mat(HC_sources_Cor)
GAD_sources_ICs = simplify_mat(GAD_sources_Cor)
MDD_sources_ICs = simplify_mat(MDD_sources_Cor)
PTSD_sources_ICs = simplify_mat(PTSD_sources_Cor)
TNP_sources_ICs = simplify_mat(TNP_sources_Cor)


All_ICs = data.frame(IC_all = rep(paste0("IC",1:10),10),
                     Diagnosis = rep(c("HC","GAD","MDD","PTSD", "TNP"), each = 10),
                     Type = rep(c("Sources","Projections"),each = 50),
                     IC=c(HC_sources_ICs[2,],GAD_sources_ICs[2,],MDD_sources_ICs[2,],PTSD_sources_ICs[2,],TNP_sources_ICs[2,],
                          HC_projections_ICs[2,],GAD_projections_ICs[2,],MDD_projections_ICs[2,],PTSD_projections_ICs[2,],TNP_projections_ICs[2,]),
                     Cor=as.numeric(c(HC_sources_ICs[3,],GAD_sources_ICs[3,],MDD_sources_ICs[3,],PTSD_sources_ICs[3,],TNP_sources_ICs[3,],
                                      HC_projections_ICs[3,],GAD_projections_ICs[3,],MDD_projections_ICs[3,],PTSD_projections_ICs[3,],TNP_projections_ICs[3,])))
All_ICs$IC = as.character(All_ICs$IC)
All_ICs$IC[rep(All_ICs$IC[All_ICs$Type == "Sources"] != All_ICs$IC[All_ICs$Type == "Projections"],2)]=NA
All_ICs$IC[rep(is.na(All_ICs$IC[All_ICs$Type == "Sources"]),2)]=NA
All_ICs$IC[rep(is.na(All_ICs$IC[All_ICs$Type == "Projections"]),2)]=NA
All_ICs$IC = factor(All_ICs$IC,levels = paste0("IC",1:10))

All_ICs$Cor[is.na(All_ICs$IC)]=NA
All_ICs$IC_all = factor(All_ICs$IC_all,levels = paste0("IC",1:10))

All_ICs_wide = reshape(All_ICs[1:50,1:4], idvar = c("Diagnosis","Type"), timevar = "IC_all", direction = "wide")
colnames(All_ICs_wide) = c("Diagnosis", "Type", paste0("IC",1:10))


# Prepare the data frame

all_sources = data.frame(
  IC_all = rep(paste0("IC", 1:10), 100),
  All = melt(ref_sources[1:10, ])$value,
  All_fast = melt(fast_sources[1:10, ])$value,
  All_fastsk = melt(fastsk_sources[1:10, ])$value,
  HC = melt(hc_sources[1:10, ])$value,
  GAD = melt(gad_sources[1:10, ])$value,
  MDD = melt(mdd_sources[1:10, ])$value,
  PTSD = melt(ptsd_sources[1:10, ])$value,
  TNP = melt(tnp_sources[1:10, ])$value
)

# Create the final data frames that is ready for plotting
all_sources$Question = rep(paste0("Q",1:100),each = 10)
all_sources_melted = melt(all_sources, id.vars = c("Question","IC_all"), value.name = "Source", variable.name = "Group")

all_sources_melted$IC_group = NA


all_sources_melted = shuffle_ICs(all_sources_melted, AllFAST_sources_ICs[,!is.na(AllFAST_sources_ICs[2,])],"All_fast")
all_sources_melted = shuffle_ICs(all_sources_melted, AllFASTSK_sources_ICs[,!is.na(AllFASTSK_sources_ICs[2,])],"All_fastsk")
all_sources_melted = shuffle_ICs(all_sources_melted, HC_sources_ICs[,!is.na(HC_sources_ICs[2,])],"HC")
all_sources_melted = shuffle_ICs(all_sources_melted, GAD_sources_ICs[,!is.na(GAD_sources_ICs[2,])],"GAD")
all_sources_melted = shuffle_ICs(all_sources_melted, PTSD_sources_ICs[,!is.na(PTSD_sources_ICs[2,])],"PTSD")
all_sources_melted = shuffle_ICs(all_sources_melted, TNP_sources_ICs[,!is.na(TNP_sources_ICs[2,])],"TNP")
all_sources_melted = shuffle_ICs(all_sources_melted, MDD_sources_ICs[,!is.na(MDD_sources_ICs[2,])],"MDD")
all_sources_melted$Sources_group[all_sources_melted$Group == "All"] = all_sources_melted$Source[all_sources_melted$Group == "All"]
all_sources_melted$IC_group_plot = as.vector(sapply(1:8, function(x)
  c(all_sources_melted$IC_group[all_sources_melted$Question == "Q1"][((x - 1) * 10 + 1):(x * 10)], rep(NA, 990))))

# Add correlation values only for one question (10 values, 1 per IC),for plotting purposes
all_sources_melted$Cor_plot = c(
  rep(NA, 1000),
  c(AllFAST_sources_ICs[3, ], rep(NA, 990)),
  c(AllFASTSK_sources_ICs[3, ], rep(NA, 990)),
  c(HC_sources_ICs[3, ], rep(NA, 990)),
  c(GAD_sources_ICs[3, ], rep(NA, 990)),
  c(MDD_sources_ICs[3, ], rep(NA, 990)),
  c(PTSD_sources_ICs[3, ], rep(NA, 990)),
  c(TNP_sources_ICs[3, ], rep(NA, 990))
)

all_sources_melted$IC_all = factor(all_sources_melted$IC_all, levels = paste0("IC",1:10))
all_sources_melted$Sources_all = rep(all_sources_melted$Sources[all_sources_melted$Group=="All"],nlevels(all_sources_melted$Group))
all_sources_melted$IC_group[all_sources_melted$Group == "All"]=as.character(all_sources_melted$IC_all[all_sources_melted$Group == "All"])
all_sources_melted$IC_group = factor(all_sources_melted$IC_group, levels = paste0("IC",1:10))

all_sources_melted$Sources_high = all_sources_melted$Sources_group
for (IC in all_sources_melted$IC_group){
  SS = all_sources_melted$Sources_group[all_sources_melted$IC_group == IC]
  qq = quantile(abs(SS),0.85, na.rm = TRUE)
  all_sources_melted$Sources_high[all_sources_melted$IC_group == IC & abs(all_sources_melted$Sources_group) < qq] = NA
}
df_plot = all_sources_melted[(all_sources_melted$Group!="All" & all_sources_melted$Group!="All_fast" & all_sources_melted$Group!="All_fastsk"),]
df_plot$IC_group = factor(df_plot$IC_group, levels = paste0("IC",1:10))

ggplot(df_plot, aes(x=Question, y= Sources_group))+
  geom_bar(stat = "identity")+
  facet_grid(Group~IC_all, scales =  "free")+
  geom_text(aes(x=50,y=5, label = IC_group_plot))
      






# Logistic Regression for Diagnosis ----------------------------------
# read projections data
library(readxl)
library(boot)
library(e1071)
library(ggplot2)
library(sf)
library(pROC)
library(randomForest)
library(caret)
par(pty = "s") # to remove the side panels when plotting the ROC curve

Projections = ICASSO_all$projections
Projections = Projections[(Projections$Diagnosis == "MDD" | Projections$Diagnosis == "HC"),]
Projections$Diagnosis = factor(Projections$Diagnosis, levels = c("HC","MDD"))


corrplot::corrplot(cor(Projections[,paste0("IC",1:10)]))

# Add "Sources_reconstructed
Sources_reco = data.frame(All_reco)
colnames(Sources_reco) = paste0("Q",1:100)
Sources_reco$Diagnosis = ICASSO_all$projections$Diagnosis
Sources_reco$ID = ICASSO_all$projections$ID
Sources_reco = Sources_reco[as.character(Sources_reco$ID) %in% as.character(Projections$ID),]
Sources_reco$Diagnosis = factor(Sources_reco$Diagnosis, levels=c("HC","MDD"))
corrplot::corrplot(cor(Sources_reco[,paste0("Q",1:100)]))

# read TPQ data
ex_dir = '/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx'

TPQ = as.data.frame(read_excel(ex_dir))

# convert all Questions to numeric
TPQ[,paste0("Q",1:100)] = apply(TPQ[,paste0("Q",1:100)],2, as.numeric)
TPQ[,c(paste0("NS",1:4),paste0("HA",1:4),paste0("RD",1:4),c("NS","HA","RD"))] = 
  apply(TPQ[,c(paste0("NS",1:4),paste0("HA",1:4),paste0("RD",1:4),c("NS","HA","RD"))],2, as.numeric)


# recalculate main dimensions, because they are not calculated by default
TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))
for (item in names(TPQQuestions)){
  # convert Question values (0s and 1s) to numeric
  TPQ[,paste0("Q",TPQQuestions[[item]])] = apply(TPQ[,paste0("Q",TPQQuestions[[item]])],2,as.numeric)
  
  # convret Question answers (T and F) to numeric
  TPQ[,paste0("QO",TPQQuestions[[item]])][TPQ[,paste0("QO",TPQQuestions[[item]])] == "T"]=1
  TPQ[,paste0("QO",TPQQuestions[[item]])][TPQ[,paste0("QO",TPQQuestions[[item]])] == "F"]=0
  TPQ[,paste0("QO",TPQQuestions[[item]])]= apply(TPQ[,paste0("QO",TPQQuestions[[item]])], 2, as.numeric)
  
  # calculate Cloninger's subscales by just summing True values (because the values given to Jurgen were not flipped,
  # so, I will calculate a Data-driven equivalent for these subscales)
  TPQ[[paste0("O_",item)]]=apply(TPQ[,paste0("QO",TPQQuestions[[item]])], 1,sum, na.rm = TRUE)
}

# Sum questions to reproduce scales and subscales
for (item in names(TPQQuestions)){
  TPQ[[item]] = apply(TPQ[,paste0("Q",TPQQuestions[[item]])],1,sum, na.rm = TRUE)
}
# Make sure that TPQ includes the same subjects as the "Projections" data set, and no NAs are included
TPQ = TPQ[as.character(TPQ$`Final ID`) %in% as.character(Projections$ID),]
TPQ$Diagnosis = factor(TPQ$Diagnosis, levels = c("HC","MDD"))

# check correlation between all dimensions
Main_dim = c("NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")
subscales = c("NS1","NS2","NS3","NS4","HA1","HA2","HA3","HA4","RD1","RD2","RD3","RD4")
subscales_or = paste0("O_",subscales)
scales = c("NS","HA","RD")
scales_or = paste0("O_",scales)

included_questions =paste0("Q",1:100)
included_questions = included_questions[included_questions!="Q61"&included_questions!="Q71"]


# prepare a data frame with only questions and scales of original data, plus the data driven subscales
Sources = ICASSO_all$sources
melted_s = melt(Sources)
melted_s$IC = factor(melted_s$IC, levels = paste0("IC",1:15))
Sources_w = tidyr::spread(melted_s, IC, value)

Sources_w[,2:16]= ifelse(abs(Sources_w[,2:16])>0.75,1,0)
q_num = as.numeric(substr(Sources_w$variable,2,nchar(as.character(Sources_w$variable))))
for(i in 1:12){
  subscale = names(TPQQuestions)[i]
  Sources_w$subscale[q_num %in% TPQQuestions[[i]]] = subscale
}
Sources_l = melt(Sources_w, id.vars = c("variable","subscale"), variable.name = "IC", value.name = "IC_loading")
ggplot(Sources_l, aes(x=subscale,y=IC_loading))+
  geom_bar(stat = "identity")+
  facet_grid(IC~.)


TPQ_source = TPQ[,c(paste0("QO",1:100),subscales_or)]

for (item in paste0("IC",1:15)){
  cor(Sources_w[[item]])
}


corrplot::corrplot(cor(TPQ[,Main_dim]))
corrplot::corrplot(cor(cbind(TPQ[,subscales], Projections[,paste0("IC",1:10)])))

corrplot::corrplot(cor(cbind(TPQ[,subscales_or], Projections[,paste0("IC",1:10)])))

Combined_data = cbind(TPQ[,subscales], Projections[,paste0("IC",1:10)])
Combined_data_or = cbind(TPQ[,subscales_or], Projections[,paste0("IC",1:10)])

corrplot::corrplot(cor(Combined_data))
corrplot::corrplot(cor(Combined_data_or))

cor_mat = matrix(rep(c(0),120), nrow = 10)
p_mat = matrix(rep(c(0),120), nrow = 10)
rownames(cor_mat) = paste0("IC",1:10)
colnames(cor_mat)= subscales
rownames(p_mat) = paste0("IC",1:10)
colnames(p_mat)= subscales

for (i in 1:10){
  IC =paste0("IC",1:10)[i]
  IC_p = sapply(subscales,function(item) cor.test(Combined_data[[item]],Combined_data[[IC]], method = "spearman")$p.value)
  est =  sapply(subscales,function(item) cor.test(Combined_data[[item]],Combined_data[[IC]], method = "spearman")$estimate)
  cor_mat[i,] = est
  p_mat[i,] = IC_p
}
corrplot::corrplot(
  cor_mat,
  method = "number",
  p.mat = p_mat,
  sig.level = 0.005,
  insig = "blank",
  title = "\n\nCorrelations Between Cloninger's and Data-Driven Decomposition of TPQ\nSpearman Rho Values for Significant Results Only"
)


for (i in 1:10){
  IC =paste0("IC",1:10)[i]
  IC_p = sapply(subscales_or,function(item) cor.test(Combined_data_or[[item]],Combined_data_or[[IC]], method = "spearman")$p.value)
  est =  sapply(subscales_or,function(item) cor.test(Combined_data_or[[item]],Combined_data_or[[IC]], method = "spearman")$estimate)
  cor_mat[i,] = est
  p_mat[i,] = IC_p
}
corrplot::corrplot(
  cor_mat,
  method = "number",
  p.mat = p_mat,
  sig.level = 0.0005,
  insig = "blank",
  title = "\n\nCorrelations Between Cloninger's and Data-Driven Decomposition of TPQ\nSpearman Rho Values for Significant Results Only"
)

# Now we do the subject groups


# create a projections glm model
Form = as.formula(paste0("Diagnosis ~ ",paste0("IC",1:10,collapse =  "+")))
ModelProjections = glm(Form, data = Projections, family = binomial())
#ModelProjections = step(ModelProjections,direction = "back")
summary(ModelProjections)

Projections_glm = LogisticFunction(Data = Projections, DV = "Diagnosis", IVs =paste0("IC",1:10), control = "HC",
                 Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)
ProjectionsPE = cv.glm(data = Projections ,ModelProjections ,K=10)$delta[2]
print(paste("IC Model accuracy based on 10 fold cv = ", 100*(1-round(ProjectionsPE,3)),"%"))

# create Sources-reconstructed glm model
Form = as.formula(paste0("Diagnosis ~ ",paste0("Q",1:100,collapse =  "+")))
ModelSources = glm(Form, data = Sources_reco, family = binomial())
#ModelSources = step(ModelSources,direction = "back")
#ModelSources = glm(Diagnosis~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q11 + Q12 + Q13 + Q14, data = Sources_reco, family = binomial())
summary(ModelSources)
LogisticFunction(ModelSources, plt_type = "histogram", Threshold = 0.5)
SourcesPE = cv.glm(data = Sources_reco ,ModelSources ,K=10)$delta[2]
print(paste("IC Model accuracy based on 10 fold cv = ", 100*(1-round(SourcesPE,3)),"%"))

LogisticFunction(ModelProjections, plt_type = "histogram", Threshold = 0.5)
ProjectionsPE = cv.glm(data = Projections ,ModelProjections ,K=10)$delta[2]
print(paste("IC Model accuracy based on 10 fold cv = ", 100*(1-round(ProjectionsPE,3)),"%"))


# Create a TPQ glm model
Form = as.formula(paste0("Diagnosis ~ ",paste0(subscales,collapse =  "+")))
ModelTPQ = glm(Form, data = TPQ, family = binomial())
#ModelTPQ = step(ModelTPQ,direction = "back")
summary(ModelTPQ)

TPQ_subscales_glm = LogisticFunction(Data = TPQ, DV = "Diagnosis", IVs =subscales, control = "HC",
                                     Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)

TPQ_scales_glm = LogisticFunction(Data = TPQ, DV = "Diagnosis", IVs =c("NS","HA","RD"), control = "HC",
                                     Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)

PETPE = cv.glm(data = TPQ ,ModelTPQ ,K=10)$delta[2]
print(paste("TPQ Model accuracy based on 10 fold cv = ", 100*(1-round(PETPE,3)),"%"))

# compare cloninger to the new method
test_roc(Projections_glm$ROC,TPQ_subscales_glm$ROC, name1 = "IC-based", name2 = "Traditional Personality Subscales",
         Subtitle = "For the Ability of Logistic Regression Models to Predict the Presence of Depression")
test_roc(Projections_glm$ROC,TPQ_scales_glm$ROC, name1 = "IC-based", name2 = "Traditional Personality Scales",
         Subtitle = "For the Ability of Logistic Regression Models to Predict the Presence of Depression")

# Create a TPQ questions glm model
TPQ_NA = na.exclude(TPQ[,c("Diagnosis",included_questions)])
Form = as.formula(paste0("Diagnosis ~ ",paste0(included_questions,collapse =  "+")))
ModelQuestions = glm(Form, data = TPQ_NA, family = binomial())
ModelQuestions = step(ModelQuestions,direction = "back")
#BackQuestions = c("Q2", "Q3", "Q5", "Q7", "Q8", "Q9", "Q13", "Q18", "Q19", "Q22", "Q26", "Q38", "Q39", "Q40", "Q41", "Q43", "Q45", "Q46", "Q49", "Q51", "Q54", "Q55", "Q56", "Q57", "Q60", "Q62", "Q64", "Q66", "Q69", "Q75", "Q77", "Q83", "Q85", "Q86", "Q88", "Q92", "Q93", "Q95")
#Form = as.formula(paste0("Diagnosis ~ ",paste0(BackQuestions,collapse =  "+")))
#ModelQuestions = glm(Form, data = TPQ_NA, family = binomial())
summary(ModelQuestions)
LogisticFunction(Data = TPQ, DV = "Diagnosis", IVs =included_questions, control = "HC",
                 Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)
# I will use the caret package to confirm the results
library(caret)
require(mlbench)

Form = as.formula(paste0("Diagnosis ~ ",paste0("IC",1:10,collapse =  "+")))
fitControl_data = trainControl(method = "cv", number = 10, savePredictions = T)
mod_fitcv_data = train(Form, data = Projections, method = "glm", family = "binomial", trControl = fitControl_data)
summary(mod_fitcv_data)
confusionMatrix(table((mod_fitcv_data$pred)$pred, (mod_fitcv_data$pred)$obs))

#Now for Sources

Form = Diagnosis~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q11 + Q12 + Q13 + Q14
fitControl_data = trainControl(method = "cv", number = 10, savePredictions = T)
mod_fitcv_data = train(Form, data = Sources_reco, method = "glm", family = "binomial", trControl = fitControl_data)
summary(mod_fitcv_data)
confusionMatrix(table((mod_fitcv_data$pred)$pred, (mod_fitcv_data$pred)$obs))


# Now for Cloninger's TPQ
Form = as.formula(paste0("Diagnosis ~ ",paste0(subscales,collapse =  "+")))
fitControl_TPQ = trainControl(method = "cv", number = 10, savePredictions = T)
mod_fitcv_TPQ = train(Form, data = TPQ, method = "glm", family = "binomial", trControl = fitControl_TPQ)
summary(mod_fitcv_TPQ)
confusionMatrix(table((mod_fitcv_TPQ$pred)$pred, (mod_fitcv_TPQ$pred)$obs))



# Now for Questions
BackQuestions = c("Q2", "Q3", "Q5", "Q7", "Q8", "Q9", "Q13", "Q18", "Q19", "Q22", "Q26", "Q38", "Q39", "Q40", "Q41", "Q43", "Q45", "Q46", "Q49", "Q51", "Q54", "Q55", "Q56", "Q57", "Q60", "Q62", "Q64", "Q66", "Q69", "Q75", "Q77", "Q83", "Q85", "Q86", "Q88", "Q92", "Q93", "Q95")
Form = as.formula(paste0("Diagnosis ~ ",paste0(BackQuestions,collapse =  "+")))
fitControl_Questions = trainControl(method = "cv", number = 10, savePredictions = T)
mod_fitcv_Questions = train(Form, data = TPQ_NA, method = "glm", family = "binomial", trControl = fitControl_Questions)
summary(mod_fitcv_Questions)
confusionMatrix(table((mod_fitcv_Questions$pred)$pred, (mod_fitcv_Questions$pred)$obs))

# Apparently, the two "Accuracy" concepts are different. I will use this one as it seems better


# SVM -------------------------------------------------------------


ref_projections = ICASSO_all$projections
fast_projections = ICASSO_all_fast$projections
fastsk_projections = ICASSO_all_fastsk$projections
hc_projections = ICASSO_hc$projections
gad_projections = ICASSO_gad$projections
mdd_projections = ICASSO_mdd$projections
ptsd_projections = ICASSO_ptsd$projections
tnp_projections = ICASSO_tnp$projections


# Create the matrix of correlations' estimates
est_threshold = 0.3
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, "All_Fast", 10,10,est_threshold,Projection = TRUE)
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)
HC_projections_Cor = df_cor(ref_projections, hc_projections, "HC", 10,10,est_threshold,Projection = TRUE)
GAD_projections_Cor = df_cor(ref_projections, gad_projections, "GAD", 10,10,est_threshold,Projection = TRUE)
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, "MDD", 10,10,est_threshold,Projection = TRUE)
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, "PTSD", 10,10,est_threshold,Projection = TRUE)
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, "TNP", 10,10,est_threshold,Projection = TRUE)


# Find the ICs in the each projections df that has the highest correlation
AllFAST_projections_ICs = simplify_mat(AllFAST_projections_Cor)
AllFASTSK_projections_ICs = simplify_mat(AllFASTSK_projections_Cor)
HC_projections_ICs = simplify_mat(HC_projections_Cor)
GAD_projections_ICs = simplify_mat(GAD_projections_Cor)
MDD_projections_ICs = simplify_mat(MDD_projections_Cor)
PTSD_projections_ICs = simplify_mat(PTSD_projections_Cor)
TNP_projections_ICs = simplify_mat(TNP_projections_Cor)



# Create a list of the separate diagnoses data frames
dia_projections = rbind(hc_projections[,1:18],mdd_projections[,1:18], ptsd_projections[,1:18], tnp_projections[,1:18], gad_projections[,1:18])
fast_projections = melt(fast_projections[,1:18],id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
fastsk_projections = melt(fastsk_projections[,1:18],id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
ref_projections = melt(ref_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")
dia_projections = melt(dia_projections,id.vars = c("ID","Diagnosis"), measure.vars = paste0("IC",1:10), value.name = "Projection", variable.name = "IC")

projections = rbind(ref_projections,fast_projections,fastsk_projections,dia_projections)
projections$Group = rep(c("All","All_Fast","All_Fast_Sickit","Separate"),each = nrow(ref_projections))
projections$Group = factor(projections$Group, levels = c("All","All_Fast","All_Fast_Sickit","Separate"))


for (ICn in 1:10){
  for (Fn in 1:nlevels(projections$Diagnosis)){
    Factor = levels(projections$Diagnosis)[Fn]
    Pro_ICs = list(GAD_projections_ICs, HC_projections_ICs, MDD_projections_ICs, PTSD_projections_ICs, TNP_projections_ICs)[[Fn]]
    projections$IC_group[projections$Diagnosis == Factor & projections$Group == "Separate"& projections$IC == paste0("IC",ICn)] = Pro_ICs[2,ICn]
  }
  projections$IC_group[projections$Group == "All_Fast"& projections$IC == paste0("IC",ICn)] = AllFAST_projections_ICs[2,ICn]
  projections$IC_group[projections$Group == "All_Fast_Sickit"& projections$IC == paste0("IC",ICn)] = AllFASTSK_projections_ICs[2,ICn]
}


for (Factor in levels(projections$Diagnosis)) {
  for (IC in paste0("IC", 1:10)) {
    for (Group in c("All_Fast", "All_Fast_Sickit", "Separate")) {
      cond_all_ic = (projections$IC == IC & projections$Diagnosis == Factor & projections$Group == Group)
      new_ic = projections$IC_group[cond_all_ic][1]
      cond_new_ic = (projections$IC == new_ic & projections$Diagnosis == Factor & projections$Group == Group)
      
      projections$Projections_group[cond_all_ic] = projections$Projection[cond_new_ic]
      projections$Projections_scaled_group[cond_all_ic] = projections$Projection_scaled[cond_new_ic]
      projections$Sort_group[cond_all_ic] = projections$Sort[cond_new_ic]
      projections$Projection_zero_group[cond_all_ic] = projections$Projection_zero[cond_new_ic]
    }
  }
}

### ### ### ### ### ### ### ### 
# Prepare sources
### ### ### ### ### ### ### ### 
library(data.table)
all_sources = transpose(ICASSO_all$sources, make.names = "IC")
fast_sources = transpose(ICASSO_all_fast$sources, make.names = "IC")
fastsk_sources = transpose(ICASSO_all_fastsk$sources, make.names = "IC")
hc_sources = transpose(ICASSO_hc$sources, make.names = "IC")
gad_sources = transpose(ICASSO_gad$sources, make.names = "IC")
mdd_sources = transpose(ICASSO_mdd$sources, make.names = "IC")
ptsd_sources = transpose(ICASSO_ptsd$sources, make.names = "IC")
tnp_sources = transpose(ICASSO_tnp$sources, make.names = "IC")

TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))

q_num = as.numeric(substr(sources_melted_all$Question,2,nchar(sources_melted_all$Question)))
for(i in 1:12){
  subscale = names(TPQQuestions)[i]
  sources_melted_all$subscale[q_num %in% TPQQuestions[[i]]] = subscale
}

sources_melted_all$q_num = q_num


# Experiments ------------------------------------------------------------
find_scale = function(scale_list,Q_num){
  Question = NA
  for (item in names(scale_list)){
    if (Q_num %in% scale_list[[item]]){
      Question = item
    }
  }
  return(Question)
}

subscales = sapply(1:100, function(x) find_scale(TPQQuestions[1:12],x))
scales = sapply(1:100, function(x) find_scale(TPQQuestions[13:15],x))
sc_df = data.frame(scales = as.factor(scales), subscales = as.factor(subscales), Questions = factor(paste0("Q",1:100),levels = paste0("Q",1:100)))

all_sources[,c("scale","subscale","Question")] = sc_df
fast_sources[,c("scale","subscale","Question")] = sc_df
fastsk_sources[,c("scale","subscale","Question")] = sc_df
hc_sources[,c("scale","subscale","Question")] = sc_df
gad_sources[,c("scale","subscale","Question")] = sc_df
mdd_sources[,c("scale","subscale","Question")] = sc_df
ptsd_sources[,c("scale","subscale","Question")] = sc_df
tnp_sources[,c("scale","subscale","Question")] = sc_df






### ### ###
### ### ###
#   SVM   #
### ### ###
### ### ###
library(e1071)
library(ggplot2)
library(sf)
library(pROC)
library(randomForest)
library(caret)
par(pty = "s") # to remove the side panels when plotting the ROC curve

pro_all = ICASSO_all$projections[,c(9:18,2)]
pro_fast = ICASSO_all_fast$projections[,c(9:18,2)]
pro_fastsk = ICASSO_all_fastsk$projections[,c(9:18,2)]

#pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"), colnames(pro_fast) %in% c("Diagnosis","IC1","IC2")]
pro_all = pro_all[pro_all$Diagnosis %in% c("HC","MDD"),]
pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"),]
pro_fastsk = pro_fastsk[pro_fastsk$Diagnosis %in% c("HC","MDD"),]

pro_all$Diagnosis = factor(pro_all$Diagnosis)
pro_fast$Diagnosis = factor(pro_fast$Diagnosis)
pro_fastsk$Diagnosis = factor(pro_fastsk$Diagnosis)





svm_info_lin = show_svm(data_frame = pro_fast, k_type = "linear",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_rad = show_svm(data_frame = pro_fast, k_type = "radial",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_pol = show_svm(data_frame = pro_fast, k_type = "polynomial",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_sig = show_svm(data_frame = pro_fast, k_type = "sigmoid",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)


svm_lin_auc = mean(svm_info_lin$cv[,2])
svm_rad_auc = mean(svm_info_rad$cv[,2])
svm_pol_auc = mean(svm_info_pol$cv[,2])
svm_sig_auc = mean(svm_info_sig$cv[,2])



# Can be deleted later
# This is a comparison with python's SVM
# Seems the answers make sense
train_test_split = import("sklearn")$model_selection$train_test_split
py_svm= import("sklearn")$svm
np = import("numpy")
pd = import("pandas")


show_py_svm = function(data_frame, perc = 0.8, svm_method = "linear"){
  # create a test and train data sets
  train_inds = sample(1:nrow(data_frame), size = perc*nrow(data_frame))
  train = data_frame[train_inds,]
  test = data_frame[-train_inds,]
  
  svm_model = py_svm$SVC(kernel=svm_method, C=1, decision_function_shape='ovo')$fit(train[1:10], train[[11]])
  res_table = table(svm_model$predict(test[1:10]),test$Diagnosis)
  Accuracy = svm_model$score(test[1:10], test[[11]])
  
  return (list(res_table,Accuracy))
}
show_py_svm(data_frame = pro_fast, perc = 0.8, svm_method = "linear") #"rbf","poly","sigmoid"
# End delete
# now we do cloninger
tpq = ICASSO_all$tpq

or_tpq = as.data.frame(readxl::read_xlsx("/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx", na = "NA"))[c("Final ID", "Diagnosis", paste0("QO",1:100))]
TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))



or_tpq = or_tpq[or_tpq$`Final ID` %in% tpq$ID,]
or_tpq[,paste0("QO",1:100)] = ifelse(or_tpq[,paste0("QO",1:100)] == "T", 1, ifelse(or_tpq[,paste0("QO",1:100)] == "F",0,NA))
or_tpq = or_tpq[or_tpq$Diagnosis %in% c("MDD","HC"),]
or_tpq$Diagnosis = factor(or_tpq$Diagnosis)

for (item in names(TPQQuestions)){
  question_labels = paste0("QO",TPQQuestions[[item]])
  or_tpq[[item]] = apply(or_tpq[,question_labels], 1, sum, na.rm = TRUE)
}


svm_cloninger_subscales_lin = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "linear", k.fold = 10)
svm_cloninger_subscales_rad = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "radial", k.fold = 10)
svm_cloninger_subscales_pol = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "polynomial", k.fold = 10)

svm_cloninger_scales_lin = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "linear", k.fold = 10)
svm_cloninger_scales_rad = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "radial", k.fold = 10)
svm_cloninger_scales_pol = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "polynomial", k.fold = 10)

# plot ROC
test_roc(svm_info_lin$ROC,svm_cloninger_subscales_lin$ROC, name1 = "IC-based", name2 = "Traditional Personality Subscales",
         Subtitle = "For the Ability of SVM Models to Predict the Presence of Depression")
# plot ROC
test_roc(svm_info_lin$ROC,svm_cloninger_scales_lin$ROC, name1 = "IC-based", name2 = "Traditional Personality Scales",
         Subtitle = "For the Ability of SVM Models to Predict the Presence of Depression")

# to plot average roc curve with a data frame: https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html

# CLoninger appears to have very low predictive value for MDD

# SVM for scales and subscales


all_source_svm = mean(sapply(1:50, function(x) show_svm(all_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
fast_source_svm = mean(sapply(1:50, function(x) show_svm(fast_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
fastsk_source_svm = mean(sapply(1:50, function(x) show_svm(fastsk_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
hc_source_svm = mean(sapply(1:50, function(x) show_svm(hc_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
gad_source_svm = mean(sapply(1:50, function(x) show_svm(gad_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
mdd_source_svm = mean(sapply(1:50, function(x) show_svm(mdd_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
ptsd_source_svm = mean(sapply(1:50, function(x) show_svm(ptsd_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
tnp_source_svm = mean(sapply(1:50, function(x) show_svm(tnp_sources, DV = "subscale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))

all_source_svm = mean(sapply(1:50, function(x) show_svm(all_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
fast_source_svm = mean(sapply(1:50, function(x) show_svm(fast_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
fastsk_source_svm = mean(sapply(1:50, function(x) show_svm(fastsk_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
hc_source_svm = mean(sapply(1:50, function(x) show_svm(hc_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
gad_source_svm = mean(sapply(1:50, function(x) show_svm(gad_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
mdd_source_svm = mean(sapply(1:50, function(x) show_svm(mdd_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
ptsd_source_svm = mean(sapply(1:50, function(x) show_svm(ptsd_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))
tnp_source_svm = mean(sapply(1:50, function(x) show_svm(tnp_sources, DV = "scale",IVs = paste0("IC",1:10), k_type = "radial")[[2]]))


# KNN -------------------------------------------------------------
library(class)
df = pro_all
DV = "Diagnosis"
IVs = paste0("IC",1:10)

## ## ## ## ## ## ## ## ## ## 
# projections
## ## ## ## ## ## ## ## ## ## 
pro_all = ICASSO_all$projections[,c(9:18,2)]
pro_fast = ICASSO_all_fast$projections[,c(9:18,2)]
pro_fastsk = ICASSO_all_fastsk$projections[,c(9:18,2)]

pro_all = pro_all[pro_all$Diagnosis %in% c("HC","MDD"),]
pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"),]
pro_fastsk = pro_fastsk[pro_fastsk$Diagnosis %in% c("HC","MDD"),]

pro_all$Diagnosis = factor(pro_all$Diagnosis)
pro_fast$Diagnosis = factor(pro_fast$Diagnosis)
pro_fastsk$Diagnosis = factor(pro_fastsk$Diagnosis)




knn_all = show_knn(data_frame = pro_all, DV = "Diagnosis", IVs = paste0("IC",1:10),perc = 0.8,control = "HC", plot.ROC = TRUE, k.fold = 10, k=11)
knn_fast = show_knn(data_frame = pro_fast, DV = "Diagnosis", IVs = paste0("IC",1:10),perc = 0.8,control = "HC", plot.ROC = TRUE, k.fold = 10, k=11)
knn_fastsk = show_knn(data_frame = pro_fastsk, DV = "Diagnosis", IVs = paste0("IC",1:10),perc = 0.8,control = "HC", plot.ROC = TRUE, k.fold = 10, k=11)


## ## ## ## ## ## ## ## ## ## 
# Cloninger
## ## ## ## ## ## ## ## ## ## 
tpq = ICASSO_all$tpq

or_tpq = as.data.frame(readxl::read_xlsx("/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx", na = "NA"))[c("Final ID", "Diagnosis", paste0("QO",1:100))]
TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))

or_tpq = or_tpq[or_tpq$`Final ID` %in% tpq$ID,]
or_tpq[,paste0("QO",1:100)] = ifelse(or_tpq[,paste0("QO",1:100)] == "T", 1, ifelse(or_tpq[,paste0("QO",1:100)] == "F",0,NA))
or_tpq = or_tpq[or_tpq$Diagnosis %in% c("MDD","HC"),]
or_tpq$Diagnosis = factor(or_tpq$Diagnosis)

for (item in names(TPQQuestions)){
  question_labels = paste0("QO",TPQQuestions[[item]])
  or_tpq[[item]] = apply(or_tpq[,question_labels], 1, sum, na.rm = TRUE)
}



knn_tpq_subscale = show_knn(data_frame = or_tpq, DV = "Diagnosis", IVs =  c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),
         perc = 0.8,control = "HC", plot.ROC = TRUE, k.fold = 10, k=13)

knn_tpq_scale = show_knn(data_frame = or_tpq, DV = "Diagnosis", IVs =  c("HA","NS","RD"),
                        perc = 0.8,control = "HC", plot.ROC = TRUE, k.fold = 10, k=13)


test_roc(knn_all$ROC,knn_tpq_subscale$ROC, name1 = "IC-based", name2 = "Traditional Personality Subscales",
         Subtitle = "For the Ability of  KNN Model to Predict the Presence of Depression")

test_roc(knn_all$ROC,knn_tpq_scale$ROC, name1 = "IC-based", name2 = "Traditional Personality Scales",
         Subtitle = "For the Ability of KNN Model to Predict the Presence of Depression")

#Cloninger is a little bit worse
## ## ## ## ## ## ## ## ## ## 
# sources
## ## ## ## ## ## ## ## ## ## 
library(data.table)
all_sources = transpose(ICASSO_all$sources, make.names = "IC")
fast_sources = transpose(ICASSO_all_fast$sources, make.names = "IC")
fastsk_sources = transpose(ICASSO_all_fastsk$sources, make.names = "IC")
hc_sources = transpose(ICASSO_hc$sources, make.names = "IC")
gad_sources = transpose(ICASSO_gad$sources, make.names = "IC")
mdd_sources = transpose(ICASSO_mdd$sources, make.names = "IC")
ptsd_sources = transpose(ICASSO_ptsd$sources, make.names = "IC")
tnp_sources = transpose(ICASSO_tnp$sources, make.names = "IC")

TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))

###
find_scale = function(scale_list,Q_num){
  Question = NA
  for (item in names(scale_list)){
    if (Q_num %in% scale_list[[item]]){
      Question = item
    }
  }
  return(Question)
}

subscales = sapply(1:100, function(x) find_scale(TPQQuestions[1:12],x))
scales = sapply(1:100, function(x) find_scale(TPQQuestions[13:15],x))
sc_df = data.frame(scales = as.factor(scales), subscales = as.factor(subscales), Questions = factor(paste0("Q",1:100),levels = paste0("Q",1:100)))

all_sources[,c("scale","subscale","Question")] = sc_df
fast_sources[,c("scale","subscale","Question")] = sc_df
fastsk_sources[,c("scale","subscale","Question")] = sc_df

all_dc = mean(sapply(1:50, function(x) show_decision_tree(data_frame = all_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
fast_dc = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fast_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
fastsk_dc = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fastsk_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))

all_dc_subscale = mean(sapply(1:50, function(x) show_decision_tree(data_frame = all_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
fast_dc_subscale = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fast_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
fastsk_dc_subscale = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fastsk_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))

show_decision_tree(data_frame = all_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = TRUE, max_depth = 5, use_control = TRUE)
fast_dc = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fast_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
fastsk_dc = mean(sapply(1:50, function(x) show_decision_tree(data_frame = fastsk_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))

all_dc
fast_dc
fastsk_dc

all_dc_subscale
fast_dc_subscale
fastsk_dc_subscale



# Machine Learning -------------------------------------------------------------


Projections = ICASSO_all$projections
Projections = Projections[(Projections$Diagnosis == "MDD" | Projections$Diagnosis == "HC"),]
Projections$Diagnosis = factor(Projections$Diagnosis, levels = c("HC","MDD"))

# prepare TPQ Data
ex_dir = '/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx'

TPQ = as.data.frame(read_excel(ex_dir))

# convert all Questions to numeric
TPQ[,paste0("Q",1:100)] = apply(TPQ[,paste0("Q",1:100)],2, as.numeric)
TPQ[,c(paste0("NS",1:4),paste0("HA",1:4),paste0("RD",1:4),c("NS","HA","RD"))] = 
  apply(TPQ[,c(paste0("NS",1:4),paste0("HA",1:4),paste0("RD",1:4),c("NS","HA","RD"))],2, as.numeric)


# recalculate main dimensions, because they are not calculated by default
TPQQuestions = list(NS1=c(2, 4, 9, 11, 40, 43, 85, 93, 96),
                    NS2=c(30, 46, 48, 50, 55, 56, 81, 99),
                    NS3=c(32, 66, 70, 72, 76, 78, 87),
                    NS4=c(13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA1=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98),
                    HA2=c(18, 19, 23, 26, 29, 47, 51),
                    HA3=c(33, 37, 38, 42, 44, 89, 100),
                    HA4=c(49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD1=c(27, 31, 34, 83, 94),
                    RD2=c(39, 41, 45, 52, 53, 77, 79, 92, 97),
                    RD3=c(3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90),
                    RD4=c(17, 20, 25, 36, 58),
                    NS=c(2, 4, 9, 11, 40, 43, 85, 93, 96, 30, 46, 48, 50, 55, 56, 81, 99, 32, 66, 70, 72, 76, 78, 87, 13, 16, 21, 22, 24, 28, 35, 60, 62, 65),
                    HA=c(1, 5, 8, 10, 14, 82, 84, 91, 95, 98, 18, 19, 23, 26, 29, 47, 51, 33, 37, 38, 42, 44, 89, 100, 49, 54, 57, 59, 63, 68, 69, 73, 75, 80),
                    RD=c(27, 31, 34, 83, 94, 39, 41, 45, 52, 53, 77, 79, 92, 97, 3, 6, 7, 12, 15, 64, 67, 74, 86, 88, 90, 17, 20, 25, 36, 58))
for (item in names(TPQQuestions)){
  # convert Question values (0s and 1s) to numeric
  TPQ[,paste0("Q",TPQQuestions[[item]])] = apply(TPQ[,paste0("Q",TPQQuestions[[item]])],2,as.numeric)
  
  # convret Question answers (T and F) to numeric
  TPQ[,paste0("QO",TPQQuestions[[item]])][TPQ[,paste0("QO",TPQQuestions[[item]])] == "T"]=1
  TPQ[,paste0("QO",TPQQuestions[[item]])][TPQ[,paste0("QO",TPQQuestions[[item]])] == "F"]=0
  TPQ[,paste0("QO",TPQQuestions[[item]])]= apply(TPQ[,paste0("QO",TPQQuestions[[item]])], 2, as.numeric)
  
  # calculate Cloninger's subscales by just summing True values (because the values given to Jurgen were not flipped,
  # so, I will calculate a Data-driven equivalent for these subscales)
  TPQ[[paste0("O_",item)]]=apply(TPQ[,paste0("QO",TPQQuestions[[item]])], 1,sum, na.rm = TRUE)
}

# Sum questions to reproduce scales and subscales
for (item in names(TPQQuestions)){
  TPQ[[item]] = apply(TPQ[,paste0("Q",TPQQuestions[[item]])],1,sum, na.rm = TRUE)
}
# Make sure that TPQ includes the same subjects as the "Projections" data set, and no NAs are included
TPQ = TPQ[as.character(TPQ$`Final ID`) %in% as.character(Projections$ID),]
TPQ$Diagnosis = factor(TPQ$Diagnosis, levels = c("HC","MDD"))


# Logistic Regression -------------------------------------------------------------
Projections_glm = LogisticFunction(Data = Projections, DV = "Diagnosis", IVs =paste0("IC",1:10), control = "HC",
                                   Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)

TPQ_subscales_glm = LogisticFunction(Data = TPQ, DV = "Diagnosis", IVs =subscales, control = "HC",
                                     Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)

TPQ_scales_glm = LogisticFunction(Data = TPQ, DV = "Diagnosis", IVs =c("NS","HA","RD"), control = "HC",
                                  Threshold = 0.5, plt_type = "histogram", perc= 0.8, plot.ROC = TRUE)

# compare cloninger to the new method
test_roc(Projections_glm$ROC,TPQ_subscales_glm$ROC, name1 = "IC-based", name2 = "Traditional Personality Subscales",
         Subtitle = "For the Ability of Logistic Regression Models to Predict the Presence of Depression")
test_roc(Projections_glm$ROC,TPQ_scales_glm$ROC, name1 = "IC-based", name2 = "Traditional Personality Scales",
         Subtitle = "For the Ability of Logistic Regression Models to Predict the Presence of Depression")


# SVM -------------------------------------------------------------

# # # # # # # # # # # # 
# IC-based SVM
# # # # # # # # # # # # 
pro_all = ICASSO_all$projections[,c(9:18,2)]
pro_fast = ICASSO_all_fast$projections[,c(9:18,2)]
pro_fastsk = ICASSO_all_fastsk$projections[,c(9:18,2)]

#pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"), colnames(pro_fast) %in% c("Diagnosis","IC1","IC2")]
pro_all = pro_all[pro_all$Diagnosis %in% c("HC","MDD"),]
pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"),]
pro_fastsk = pro_fastsk[pro_fastsk$Diagnosis %in% c("HC","MDD"),]

pro_all$Diagnosis = factor(pro_all$Diagnosis)
pro_fast$Diagnosis = factor(pro_fast$Diagnosis)
pro_fastsk$Diagnosis = factor(pro_fastsk$Diagnosis)

svm_info_lin = show_svm(data_frame = pro_all, k_type = "linear",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_rad = show_svm(data_frame = pro_all, k_type = "radial",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_pol = show_svm(data_frame = pro_all, k_type = "polynomial",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)
svm_info_sig = show_svm(data_frame = pro_all, k_type = "sigmoid",DV = "Diagnosis",IVs = paste0("IC",1:10), k.fold = 10)


# # # # # # # # # # # # 
# Traditional SVM
# # # # # # # # # # # # 
svm_cloninger_subscales_lin = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "linear", k.fold = 10)
svm_cloninger_subscales_rad = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "radial", k.fold = 10)
svm_cloninger_subscales_pol = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),k_type = "polynomial", k.fold = 10)

svm_cloninger_scales_lin = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "linear", k.fold = 10)
svm_cloninger_scales_rad = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "radial", k.fold = 10)
svm_cloninger_scales_pol = show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"),k_type = "polynomial", k.fold = 10)

# plot ROC
test_roc(svm_info_lin$ROC,svm_cloninger_subscales_lin$ROC, name1 = "IC-based", name2 = "Traditional Personality Subscales",
         Subtitle = "For the Ability of SVM Models to Predict the Presence of Depression")
# plot ROC
test_roc(svm_info_lin$ROC,svm_cloninger_scales_lin$ROC, name1 = "IC-based", name2 = "Traditional Personality Scales",
         Subtitle = "For the Ability of SVM Models to Predict the Presence of Depression")

# Basic TPQ Analysis --------------------------------------------------------

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#             Crobnach Alpha                  #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

return_alpha = function(item) {
  # returns 3 items, the name of the subscale,
  # c.alpha of males and c.alpha of females
  
  questions = paste0("Q", TPQQuestions[[item]])
  df = TPQ_scores[questions]
  M_alpha = cronbach(df[TPQ_scores$Gender == "Male", ])$alpha
  F_alpha = cronbach(df[TPQ_scores$Gender == "Female", ])$alpha
  return(c( M_alpha, F_alpha))
}
c_alpha_values = sapply(names(TPQQuestions), function(x) return_alpha(x))
c_alpha_df = as.data.frame(t(c_alpha_values))
colnames(c_alpha_df) = c("Males","Females")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#             Effect of Treatment on Temperaments                  #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ##
scales = c("NS", "HA", "RD")
subscales = c(paste0("NS",1:4),paste0("HA",1:4),paste0("RD",1:4))
all_scales = c("NS",paste0("NS",1:4),"HA",paste0("HA",1:4),"RD",paste0("RD",1:4))

scale_names = c("NS: Novelty Seeking",
                "NS1: Exploratory excitability vs stoic rigidity",
                "NS2: Impulsiveness vs reflection",
                "NS3: Extravagance vs reserve",
                "NS4: Disorderliness vs regimentation",
                "HA: Harm avoidance",
                "HA1: Anticipatory worry vs uninhibited optimism",
                "HA2: Fear of uncertainty vs confidence",
                "HA3: Shyness with strangers vs gregariousness",
                "HA4: Fatigability and asthenia vs vigor",
                "RD: Reward Dependence",
                "RD1: Semtimentality vs insensitiveness",
                "RD2: Persistence vs irresoluteness",
                "RD3: Attachment vs detachment",
                "RD4: Dependence vs independence")

scale_names = c("NS: Novelty Seeking","NS1: Exploratory excitability","NS2: Impulsiveness","NS3: Extravagance","NS4: Disorderliness",
                "HA: Harm avoidance", "HA1: Anticipatory worry", "HA2: Fear of uncertainty", "HA3: Shyness with strangers", "HA4: Fatigability and asthenia",
                "RD: Reward Dependence", "RD1: Semtimentality", "RD2: Persistence", "RD3: Attachment", "RD4: Dependence")

# Create the test-retest data set
TPQ_scores_TR = TPQ_scores[complete.cases(TPQ_scores$`Initial ID`),]
TPQ_scores_TR = TPQ_scores_TR[complete.cases(TPQ_scores_TR$Response),]
TPQ_scores_TR = MatchDelete(TPQ_scores_TR,"Session","Initial ID")

# Create array for the Test-retest data set(TPQ_scores_TR)
Means_array = sapply(all_scales, function(x)
  SMeans(TPQ_scores_TR, DV = x, IVs = c("Session", "Response"),
    GroupBy = "Response"), simplify = FALSE)

# Bind all data frames into one
Means_df = Means_array[[1]]
colnames(Means_df) = c("Value","SEM","Count","Session","Response","Groups")
for (i in 2:length(Means_array)){
  new_df = Means_array[[i]]
  colnames(new_df) = c("Value","SEM","Count","Session","Response","Groups")
  Means_df = rbind(Means_df,new_df)
}

# Add scales and subscales columns
Means_df$subscale = factor(rep(c("Scale",paste("Sub_",1:4)), each = 6))
Means_df$scale = factor(rep(scales, each = 30))

# Add scale and subscale names
Means_df$scale_names = ""
Means_df$scale_names[seq(from = 1, to = 90, by = 6)] = scale_names
# make session into a factor with levels = "Test"&"Retest"
Means_df$Session = factor(Means_df$Session, levels = c("Test","Retest"))

# plot all
ggplot(Means_df, aes(x=Session, y=Value, group = Groups, color = Groups))+
  geom_line(position = position_dodge(0.05),size=1)+
  geom_errorbar(aes(ymax = Value+SEM,ymin = Value-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  geom_text(mapping = aes(x=1.5,y=0,label = scale_names), color = "black")+
  scale_color_manual(values = c("#6cBE58","#C33E3B","#4EA3DF"))+
  ggtitle("Harm Avoidance Comparison According to Response")+
  facet_grid(scale~subscale,scale = "free")


# I will see the reconstructed data
IC_length = length(ICASSO_all$reco)
IC_sums = sapply(1:IC_length, function(x) apply(ICASSO_all$reco[[x]][,9:108],1,sum))
TPQ_reco = ICASSO_all$reco[[1]][,1:8]

TPQ_reco[paste0("IC",1:IC_length)] = IC_sums

sapply(TPQ_reco$ID, function(id) TPQ_scores$`Initial ID`[TPQ_scores$`Final ID` == id])
TPQ_reco$`Initial ID` = 

# Create the test-retest data set
ICA_scores_TR = TPQ_reco[complete.cases(TPQ_reco$`Initial ID`),]
ICA_scores_TR = ICA_scores_TR[complete.cases(ICA_scores_TR$Response),]
ICA_scores_TR = MatchDelete(ICA_scores_TR,"Session","Initial ID")

# Create array for the Test-retest data set(TPQ_scores_TR)
Means_array = sapply(all_scales, function(x)
  SMeans(TPQ_scores_TR, DV = x, IVs = c("Session", "Response"),
         GroupBy = "Response"), simplify = FALSE)

# Bind all data frames into one
Means_df = Means_array[[1]]
colnames(Means_df) = c("Value","SEM","Count","Session","Response","Groups")
for (i in 2:length(Means_array)){
  new_df = Means_array[[i]]
  colnames(new_df) = c("Value","SEM","Count","Session","Response","Groups")
  Means_df = rbind(Means_df,new_df)
}

# Add scales and subscales columns
Means_df$subscale = factor(rep(c("Scale",paste("Sub_",1:4)), each = 6))
Means_df$scale = factor(rep(scales, each = 30))

# Add scale and subscale names
Means_df$scale_names = ""
Means_df$scale_names[seq(from = 1, to = 90, by = 6)] = scale_names
# make session into a factor with levels = "Test"&"Retest"
Means_df$Session = factor(Means_df$Session, levels = c("Test","Retest"))

# plot all
ggplot(Means_df, aes(x=Session, y=Value, group = Groups, color = Groups))+
  geom_line(position = position_dodge(0.05),size=1)+
  geom_errorbar(aes(ymax = Value+SEM,ymin = Value-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  geom_text(mapping = aes(x=1.5,y=0,label = scale_names), color = "black")+
  scale_color_manual(values = c("#6cBE58","#C33E3B","#4EA3DF"))+
  ggtitle("Harm Avoidance Comparison According to Response")+
  facet_grid(scale~subscale,scale = "free")
