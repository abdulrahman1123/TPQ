# Data Preparation --------------------------------------------------------
library(abind)
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
library(psych)
library(doParallel)
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

# load original scores
TPQ_scores = as.data.frame(read_xlsx(scores_path, na = "NA"))
for (item in names(TPQQuestions)){
  TPQ_scores[[item]] = apply(TPQ_scores[,paste0("Q",TPQQuestions[[item]])],1,sum, na.rm = TRUE)
}

# find which subjects have more than 10 unanswerd questions
na_values = apply(is.na(TPQ_scores[,paste0("QO",1:100)]),1,sum) > 10
TPQ_scores[na_values,names(TPQQuestions)] = NA

TPQ_scores[,paste0("QO",1:100)] = ifelse(TPQ_scores[,paste0("QO",1:100)] =="T",TRUE,FALSE)

gen_info = TPQ_scores[c("Diagnosis", "Final ID", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session")]

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
All_reco=np$transpose(icasso$ica2data(sources, unmixing, idx_keep=0:9))
All_reco = All_reco[npz_tpq["IDs"] %in% gen_info$`Final ID`,]


# define paths for variuos decompositions
all_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz')
all_fast_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD_FastICA_MNE/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEfastica.npz')
all_fastsk_path = paste0(path, '/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD_FastICA_sklearn/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_SklearnFastICA.npz')
hc_path  = paste0(path, '/ICASSO_fixed/decomp_tpq/HC/icasso_ICA-tpq_HC_nsamp1202_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
hc_retest_path  = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC/retest/icasso_ICA-tpq_HC_nsamp90_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
hc_test_path  = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC/test/icasso_ICA-tpq_HC_nsamp90_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
hc_test_retest_path ='' # fill later
# TODO for MDD as well
mdd_retest_path  = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/retest/icasso_ICA-tpq_MDD_nsamp91_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
mdd_test_path  = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/test/icasso_ICA-tpq_MDD_nsamp91_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
mdd_path = paste0(path, '/ICASSO_fixed/decomp_tpq/MDD/icasso_ICA-tpq_MDD_nsamp455_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
gad_path = paste0(path, '/ICASSO_fixed/decomp_tpq/GAD/icasso_ICA-tpq_GAD_nsamp30_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
ptsd_path = paste0(path, '/ICASSO_fixed/decomp_tpq/PTSD/icasso_ICA-tpq_PTSD_nsamp57_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
tnp_path = paste0(path, '/ICASSO_fixed/decomp_tpq/TNP/icasso_ICA-tpq_TNP_nsamp78_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
disorders_path = paste0(path, '/ICASSO_fixed/decomp_tpq/MDD,PTSD,TNP,GAD/icasso_ICA-tpq_MDD,PTSD,TNP,GAD_nsamp620_n_comp13_n_iter100_dist0.30_MNEinfomax.npz')
noHC_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD,PTSD,TNP,GAD/icasso_ICA-tpq_MDD,PTSD,TNP,GAD_nsamp620_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
noMDD_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC,PTSD,TNP,GAD/icasso_ICA-tpq_HC,PTSD,TNP,GAD_nsamp1367_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
noTNP_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC,MDD,PTSD,GAD/icasso_ICA-tpq_HC,MDD,PTSD,GAD_nsamp1744_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
noPTSD_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC,MDD,TNP,GAD/icasso_ICA-tpq_HC,MDD,TNP,GAD_nsamp1765_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
noGAD_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/HC,MDD,TNP,PTSD/icasso_ICA-tpq_HC,MDD,TNP,PTSD_nsamp1792_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
test_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/All/test/icasso_ICA-tpq_All_nsamp181_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
retest_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/All/retest/icasso_ICA-tpq_All_nsamp181_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
resp_test_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/responders/test/icasso_ICA-tpq_MDD_nsamp51_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
resp_retest_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/responders/retest/icasso_ICA-tpq_MDD_nsamp51_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
nonresp_test_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/nonresponder/test/icasso_ICA-tpq_MDD_nsamp27_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'
nonresp_retest_path = '/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/ICASSO/decomp_tpq/MDD/nonresponder/retest/icasso_ICA-tpq_MDD_nsamp27_n_comp13_n_iter100_dist0.30_MNEinfomax.npz'

# Load all data into lists of variables
ICASSO_all = prepare_icasso(all_path)
ICASSO_all_fast = prepare_icasso(all_fast_path)
ICASSO_all_fastsk = prepare_icasso(all_fastsk_path)
ICASSO_hc = prepare_icasso(hc_path)
ICASSO_hc_test = prepare_icasso(hc_test_path)
ICASSO_hc_retest = prepare_icasso(hc_retest_path)
ICASSO_mdd = prepare_icasso(mdd_path)
ICASSO_mdd_test = prepare_icasso(mdd_test_path)
ICASSO_mdd_retest = prepare_icasso(mdd_retest_path)
ICASSO_gad = prepare_icasso(gad_path)
ICASSO_ptsd = prepare_icasso(ptsd_path)
ICASSO_tnp = prepare_icasso(tnp_path)
ICASSO_disorders = prepare_icasso(disorders_path)
ICASSO_noHC = prepare_icasso(noHC_path)
ICASSO_noMDD = prepare_icasso(noMDD_path)
ICASSO_noTNP = prepare_icasso(noTNP_path)
ICASSO_noPTSD = prepare_icasso(noPTSD_path)
ICASSO_noGAD = prepare_icasso(noGAD_path)
ICASSO_test = prepare_icasso(test_path)
ICASSO_retest = prepare_icasso(retest_path)
ICASSO_resp_test = prepare_icasso(resp_test_path)
ICASSO_resp_retest = prepare_icasso(resp_retest_path)
ICASSO_nonresp_test = prepare_icasso(nonresp_test_path)
ICASSO_nonresp_retest = prepare_icasso(nonresp_retest_path)

#df_2 = pro_hc_retest
#df_1 = pro_hc_test
df_cor <- function(df_1, df_2, Factor1, Factor2, max_components_1, max_components_2, estimate_threshold,
                   stringent_alpha = TRUE,Projection = FALSE, plot = TRUE, cor_method = "pearson"){
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
    if (Projection){
      max_components_1 = ncol(df_1)-8
      max_components_2 = ncol(df_2)-8
    }
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
  names(dimnames(Cor_mat)) = c(Factor1,Factor2)
  
  p_mat = matrix(0,nrow = max_components_1,ncol = max_components_2)
  dimnames(p_mat) = list(paste0("IC",1:max_components_1), paste0("IC",1:max_components_2))
  names(dimnames(p_mat)) = c(Factor1,Factor2)
  
  for (i in 1:max_components_1){
    threshold = 0.05
    if (stringent_alpha) {
      threshold = threshold / (max_components_2 * max_components_1)
    }
    
    estimates = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]), method = cor_method)$estimate)
    p_values = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]), method = cor_method)$p.value)
    #estimates[p_values>threshold | abs(estimates)<estimate_threshold] = 0
    p_mat[i,] = p_values
    Cor_mat[i,] = estimates
  }
  #Cor_mat[p_values>threshold |abs(Cor_mat)<estimate_threshold] =0
  if(plot){
    
    corrplot::corrplot(Cor_mat, method = "number", p.mat=p_mat, sig.level = threshold,
                       insig = "blank", title = paste0("\n\nCorrelations Between ",Factor1," and ",Factor2))
  }
  return(list(correlation = Cor_mat, pvalue = p_mat))
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



find_threshold_old = function(or_data_1, or_data_2, n_shuffle, cor_method = "spearman"){
  
  # A function to find the correlation threshold of two data frames
  # both data frames should have the same number of columns
  shuffled_data_1 = array(data = 0, dim = c(n_shuffle,nrow(or_data_1),ncol(or_data_1)))
  shuffled_data_2 = array(data = 0, dim = c(n_shuffle,nrow(or_data_2),ncol(or_data_2)))
  
  
  for (perm in 1:n_shuffle){
    # shuffle the two data frames
    data_copy_1 = sample(or_data_1)
    data_copy_2 = sample(or_data_2)
    
    # add the shuffled data frames to the shuffled_data_1 and shuffled_data_2
    shuffled_data_1[perm,,] = as.matrix(data_copy_1)
    shuffled_data_2[perm,,] = as.matrix(data_copy_2)
  }
  
  
  # find the cross-correlation for each of the rows of the two data frames
  # I should have a final data frame with a sahpe of (n_shuffle, or_data_1.shape[0], or_data_2.shape[0])
  cor_matrix = array(c(0), dim = c(n_shuffle, nrow(or_data_1), nrow(or_data_2)))
  
  
  for (m in 1:ncol(shuffled_data_1)){
    for (n in 1:ncol(shuffled_data_2)){
      cor_r = sapply(1:n_shuffle, function(shuf) cor.test(shuffled_data_1[shuf, m, ],shuffled_data_2[shuf, n, ],method = cor_method)$estimate)
      cor_matrix[, m, n] = abs(cor_r)
    }
  }
  
  return (cor_matrix)
  
}


find_threshold = function(or_data_1, or_data_2, method = "spearman", n_shuffle = 1000, threshold_choice = 0.999,
                          Comparison = "Sources", G1, G2, plot=TRUE, simulate = FALSE){
  # find the threshold of correlation for two data frames
  # This is done by shuffling the data a number of times
  # The correlation values that result represent all the random correlations possible
  # so, we are sure that anything within the range of these values
  # is not acceptable. Any value above the threshold can be considered as good
  
  # remove non-numeric columns
  or_data_1 = or_data_1[,sapply(or_data_1,is.numeric)]
  or_data_2 = or_data_2[,sapply(or_data_2,is.numeric)]
  
  
  # shuffle the data frames
  local_shuffle = function(or_data_1, or_data_2, method, n_shuffle){
    # A function to create data shuffles from original data
    # Create two lists from the original data frames
    data_1 = abind(c(or_data_1))
    data_2 = abind(c(or_data_2))
    
    # normalize  
    data_1 = (data_1-mean(data_1))/sd(data_1)
    data_2 = (data_2-mean(data_2))/sd(data_2)
    
    # create a vector containing the max values of each correlation
    threshold_vector = c()
    for (i in 1:n_shuffle){
      shuffled_d1 = data.frame(matrix(sample(data_1), ncol = ncol(or_data_1)))
      shuffled_d2 = data.frame(matrix(sample(data_2), ncol = ncol(or_data_2)))
      
      colnames(shuffled_d1) = colnames(or_data_1)
      colnames(shuffled_d2) = colnames(or_data_2)
      
      cor_matrix = abs(corr.test(shuffled_d1,shuffled_d2 , method = method)$r)
      local_threshold = quantile(cor_matrix,0.99)#max(cor_matrix)
      threshold_vector = append(threshold_vector,local_threshold)
    }
    return (threshold_vector)
  }
  
  
  # the simulate function
  local_simulate = function(or_data_1, or_data_2, method = "spearman", n_shuffle = 100){
    # find sd and mean for each of the data frames
    sd_1 = sd(abind(c(or_data_1)))
    sd_2 = sd(abind(c(or_data_2)))
    mean_1 = mean(abind(c(or_data_1)))
    mean_2 = mean(abind(c(or_data_2)))
    
    # create a vector for thresholds
    threshold_vector = c()
    for (i in 1:n_shuffle){
      # create the simulated data
      sim_data_1 = data.frame(sapply(colnames(or_data_1)[1], function(icol)
        rnorm(nrow(or_data_1), mean = mean_1, sd = sd_1)))
      
      sim_data_2 = data.frame(sapply(colnames(or_data_2)[1], function(icol)
        rnorm(nrow(or_data_2), mean = mean_2,sd = sd_2)))
      
      
      cor_matrix = abs(corr.test(sim_data_1,sim_data_2 , method = method)$r)
      local_threshold = max(cor_matrix)#quantile(cor_matrix,0.95)
      threshold_vector = append(threshold_vector,local_threshold)
    }
    return(threshold_vector)
  }
  
  
  # run the local simulator/shuffler function. Change the cores argument to the number of cores in the machine
  if (simulate){
    registerDoParallel(cores=8)
    threshold_list = foreach(i = 1:8) %dopar% local_simulate (
      or_data_1 = or_data_1, or_data_2 = or_data_2, method = method, n_shuffle = ceiling(n_shuffle / 8))
  } else{
    registerDoParallel(cores=8)
    threshold_list = foreach(i = 1:8) %dopar% local_shuffle (
      or_data_1 = or_data_1, or_data_2 = or_data_2, method = method, n_shuffle = ceiling(n_shuffle / 8)
    )
  }
  
  # The result is 8 threshold lists (because of parallelization). Now we combine them
  thresholds = unlist(threshold_list)
  cor_threshold = quantile(thresholds,threshold_choice)
  standard_thresholds = sapply(c(0.1,0.5,0.9,0.95,0.99), function(x) quantile(thresholds,x))
  # norm_data_1 = (or_data_1-mean(data_1))/sd(data_1)
  # norm_data_2 = (or_data_2-mean(data_2))/sd(data_2)
  
  cor_matrix = abs(corr.test(or_data_1,or_data_2 , method = method)$r)
  
  
  cor_list_final_1 = cor_matrix
  cor_list_final_1[cor_list_final_1<cor_threshold] = 0
  
  # cor_list_2 is the list passing twice as the threshold
  cor_list_final_2 = cor_matrix
  cor_list_final_2[cor_list_final_2<(2*cor_threshold)] = 0
  
  if (plot){
    C_palette = colorRampPalette(c("#FFFFFF","#0087BD","#C33E3B"))(200)
    # anything above 0 will have red color "#C33E3B", anything above two thresholds will have blue color "#0087BD"
    
    first_position = round(200*min(cor_matrix)/max(cor_matrix),0)
    threshold_position = round(200*(cor_threshold-min(cor_matrix))/(max(cor_matrix)-min(cor_matrix)),0)
    if (threshold_position<0){
      threshold_position = 0}
    threshold_two_position = round(200*(2*(cor_threshold)-min(cor_matrix))/(max(cor_matrix)-min(cor_matrix)),0)
    if (threshold_two_position<0){
      threshold_two_position =0
    }else if (threshold_two_position>200){
      threshold_two_position = 200
    }
    C_palette[threshold_position:threshold_two_position] = "#C33E5B"
    C_palette[1:threshold_position] = "#FFFFFF"
    C_palette[threshold_two_position:200] = "#0085BD"
    # This happens if the 2*threshold is higher than the max
    if (threshold_two_position == 200){
      C_palette[threshold_two_position:200] = C_palette[199]
    }
    C_palette[199] = gsub("5","4",C_palette[200])
    C_palette[200] = gsub("4","8",C_palette[199])
    
    corrplot(cor_matrix, col = C_palette,
             title = paste0("\n\nCorrelation of ",Comparison," - ",G1," with ",G2," -",
                            if (simulate){"Simulation"}else{"Shuffling"}," Threshold = ",
                            round(cor_threshold,3), " (", threshold_choice*100,"th centile) - ",
                            n_shuffle," Iterations", "\n Red = Passing the Threshold ... Blue = passing 2*Threshold"),
             method = c("n"),cl.lim = c(min(cor_matrix),max(cor_matrix)))
  }
  return(list(cor_list_threshold = cor_list_final_1,cor_list_twice = cor_list_final_2, cor_list = cor_matrix, threshold = cor_threshold, standard_thresholds = standard_thresholds))
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
AllFAST_sources_Cor = df_cor(ref_sources, fast_sources, Factor1 = "All",Factor2 = "All_Fast", 10,10,est_threshold,Projection = FALSE)$correlation
AllFASTSK_sources_Cor = df_cor(ref_sources, fastsk_sources, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = FALSE)$correlation
HC_sources_Cor = df_cor(ref_sources, hc_sources, Factor1 = "All",Factor2 = "HC", 10,10,est_threshold,Projection = FALSE)$correlation
GAD_sources_Cor = df_cor(ref_sources, gad_sources, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = FALSE)$correlation
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = FALSE)$correlation
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = FALSE)$correlation
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = FALSE)$correlation


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
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, Factor1 = "All",Factor2 =  "All_Fast", 10,10,est_threshold,Projection = TRUE)$correlation
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)$correlation
HC_projections_Cor = df_cor(ref_projections, hc_projections, Factor1 = "All",Factor2 =  "HC", 10,10,est_threshold,Projection = TRUE)$correlation
GAD_projections_Cor = df_cor(ref_projections, gad_projections, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = TRUE)$correlation
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = TRUE)$correlation
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = TRUE)$correlation
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = TRUE)$correlation


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

df_cor(fastsk_projections, fastsk_projections, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)$correlation

AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, Factor1 = "All",Factor2 =  "All_Fast", 10,10,est_threshold,Projection = TRUE)$correlation
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)$correlation

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
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, Factor1 = "All",Factor2 =  "All_Fast", 10,10,est_threshold,Projection = TRUE)$correlation
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)$correlation
HC_projections_Cor = df_cor(ref_projections, hc_projections, Factor1 = "All",Factor2 =  "HC", 10,10,est_threshold,Projection = TRUE)$correlation
GAD_projections_Cor = df_cor(ref_projections, gad_projections, Factor1 = "All",Factor2 = "GAD", 10,10,est_threshold,Projection = TRUE)$correlation
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, Factor1 = "All",Factor2 = "MDD", 10,10,est_threshold,Projection = TRUE)$correlation
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = TRUE)$correlation
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = TRUE)$correlation


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
HC_sources_Cor = df_cor(ref_sources, hc_sources, Factor1 = "All",Factor2 =  "HC", 10,10,est_threshold,Projection = FALSE)$correlation
GAD_sources_Cor = df_cor(ref_sources, gad_sources, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = FALSE)$correlation
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = FALSE)$correlation
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = FALSE)$correlation
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = FALSE)$correlation


HC_projections_Cor = df_cor(ref_projections, hc_projections, Factor1 = "All",Factor2 =  "HC", 10,10,est_threshold,Projection = TRUE)$correlation
GAD_projections_Cor = df_cor(ref_projections, gad_projections, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = TRUE)$correlation
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = TRUE)$correlation
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = TRUE)$correlation
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = TRUE)$correlation

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
  ggtitle("ICs with Most Correlations from the Separate-Group Decomposition for Sources", subtitle = "All Correlations are with a p-value of <0.0005 and Spearman rho >0.3")


ggplot(data = All_ICs[All_ICs$Type == "Projections",], mapping = aes(x = 1,y= Type))+
  geom_text(aes(y=1,label = IC))+
  geom_text(aes(y=0,label = Cor))+
  facet_grid(IC_all~Diagnosis)+
  TypicalTheme+
  scale_y_discrete("IC")+
  scale_x_discrete("Diagnosis")+
  theme(axis.text.y=element_blank(),axis.text.x=element_blank(), axis.ticks.y=element_blank())+
  ggtitle("ICs with Most Correlations from the Separate-Group Decomposition for Projections", subtitle = "All Correlations are with a p-value of <0.0005 and Spearman rho >0.3")


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
HC_sources_Cor = df_cor(ref_sources, hc_sources, Factor1 = "All",Factor2 =  "HC", 10,10,est_threshold,Projection = FALSE)$correlation
GAD_sources_Cor = df_cor(ref_sources, gad_sources, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = FALSE)$correlation
MDD_sources_Cor = df_cor(ref_sources, mdd_sources, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = FALSE)$correlation
PTSD_sources_Cor = df_cor(ref_sources, ptsd_sources, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = FALSE)$correlation
TNP_sources_Cor = df_cor(ref_sources, tnp_sources, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = FALSE)$correlation


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
df_plot$Group = factor(df_plot$Group, levels = c("GAD","HC","MDD","PTSD","TNP"))

ggplot(df_plot, aes(x=Question, y= Sources_group))+
  geom_bar(stat = "identity")+
  facet_grid(IC_all~Group, scales =  "free")+
  geom_text(aes(x=50,y=5, label = IC_group_plot))+
  geom_text(aes(x=50,y=-5, label = Cor_plot))
      


# Questions with highest loading 
all_sources = all_sources_melted[(all_sources_melted$Group=="All"),]
all_sources 

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
AllFAST_projections_Cor = df_cor(ref_projections, fast_projections, Factor1 = "All",Factor2 = "All_Fast", 10,10,est_threshold,Projection = TRUE)$correlation
AllFASTSK_projections_Cor = df_cor(ref_projections, fastsk_projections, Factor1 = "All",Factor2 =  "All_Fast_SK", 10,10,est_threshold,Projection = TRUE)$correlation
HC_projections_Cor = df_cor(ref_projections, hc_projections, Factor1 = "All",Factor2 = "HC", 10,10,est_threshold,Projection = TRUE)$correlation
GAD_projections_Cor = df_cor(ref_projections, gad_projections, Factor1 = "All",Factor2 =  "GAD", 10,10,est_threshold,Projection = TRUE)$correlation
MDD_projections_Cor = df_cor(ref_projections, mdd_projections, Factor1 = "All",Factor2 =  "MDD", 10,10,est_threshold,Projection = TRUE)$correlation
PTSD_projections_Cor = df_cor(ref_projections, ptsd_projections, Factor1 = "All",Factor2 =  "PTSD", 10,10,est_threshold,Projection = TRUE)$correlation
TNP_projections_Cor = df_cor(ref_projections, tnp_projections, Factor1 = "All",Factor2 =  "TNP", 10,10,est_threshold,Projection = TRUE)$correlation


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

# To find c.alpha for ICs, we need to assign questions to each of the ICs
# Find which IC each question loads most on
sources_all = ICASSO_all$sources[ICASSO_all$scores>0.85,1:100]
n_ICs = sum(ICASSO_all$scores>0.85)
max_source = matrix(rep(apply(abs(sources_all),2,max),n_ICs),nrow = n_ICs, byrow = TRUE)
IC_questions = as.data.frame(ifelse(max_source == abs(sources_all),1,0))
colnames(IC_questions) = paste0("Q",1:100)

# Create a list for each of the ICs and the questions that belong to it
ICs = ICASSO_all$sources$IC[ICASSO_all$scores>0.85]
TPQ_IC_Questions = sapply(ICs, function(x) IC="", simplify = FALSE)
for (IC in ICs){
  Questions = colnames(IC_questions)[IC_questions[ICs == IC,] == 1]
  if (length(Questions)>2){
    TPQ_IC_Questions[[IC]] = Questions
  } else{
    TPQ_IC_Questions[[IC]] = NULL
  }
}


#plot the results to show the questions chosen for each of the ICs
low_q_count = apply(IC_questions,1,sum)<3
IC_Questions_plot = IC_questions
IC_Questions_plot$IC = factor(ICs, levels = ICs)
IC_Questions_plot = reshape2::melt(IC_Questions_plot, id.vars = "IC", value.name = "Value", variable.name="Question")

IC_q_g = ggplot(IC_Questions_plot,aes(y=Question,x=IC,fill = Value))+geom_raster()+
  scale_fill_gradientn(colors = c("#FFFFFF","#02422f"))+
  TypicalTheme+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+
  ggtitle("Chosen Questions for ICs based on Sources",subtitle = "ICs with low Question count are crossed with a red line")+
  geom_vline(xintercept=which(low_q_count), color = "red")

#Find factors
factanal <- fa(r=TPQ_scores[,paste0("Q",c(1:60,62:70,72:100))],
               nfactors = 12, 
               # covar = FALSE
               #SMC = TRUE,
               fm="minres",
               max.iter=50,
               rotate="varimax")

loadings = as.data.frame(t(matrix(factanal$loadings, byrow = FALSE, nrow = 98)))
loadings = sapply(1:ncol(loadings), function(x) ifelse(abs(loadings[,x])==max(abs(loadings[,x])),1,0))
colnames(loadings) = paste0("Q",c(1:60,62:70,72:100))
# Assess correspondance between these factors and Cloninger's
Factors = paste0("F",1:12)
Factors_list = sapply(Factors, function(x) Factor="", simplify = FALSE)
for (Factor in Factors){
  Questions = colnames(loadings)[loadings[Factors == Factor,] == 1]
  if (length(Questions)>2){
    Factors_list[[Factor]] = Questions
  } else{
    Factors_list[[Factor]] = NULL
  }
}



#plot the results to show the questions chosen for each of the ICs
low_q_count_factor = apply(loadings,1,sum)<3
F_Questions_plot = as.data.frame(loadings)
F_Questions_plot$Factor = factor(Factors, levels = Factors)
F_Questions_plot = reshape2::melt(F_Questions_plot, id.vars = "Factor", value.name = "Value", variable.name="Question")

factor_q_g = ggplot(F_Questions_plot,aes(y=Question,x=Factor,fill = Value))+geom_raster()+
  scale_fill_gradientn(colors = c("#FFFFFF","#02422f"))+
  TypicalTheme+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+
  ggtitle("Chosen Questions for Factors based on Q. Loading",subtitle = "Factors with low Question count are crossed with a red line")+
  geom_vline(xintercept=which(low_q_count_factor), color = "red")



#plot TPQ Question count
tpq_q_cloninger = as.data.frame(sapply(names(TPQQuestions), function(item) ifelse(1:100 %in% TPQQuestions[[item]],1,0)))
tpq_q_cloninger = t(tpq_q_cloninger)
colnames(tpq_q_cloninger) = paste0("Q",1:100)
low_q_count_factor = apply(tpq_q_cloninger,1,sum)<3

Scale_Questions_plot = as.data.frame(tpq_q_cloninger)
Scale_Questions_plot$Scale = factor(names(TPQQuestions), levels = names(TPQQuestions))
Scale_Questions_plot = reshape2::melt(Scale_Questions_plot, id.vars = "Scale", value.name = "Value", variable.name="Question")

TPQ_q_g = ggplot(Scale_Questions_plot,aes(y=Question,x=Scale,fill = Value))+geom_raster()+
  scale_fill_gradientn(colors = c("#FFFFFF","#02422f"))+
  TypicalTheme+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+
  ggtitle("Chosen Questions for Subscales\nAccording to Cloninger")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#       Effect of Diagnosis and Treatment on Temperaments          #
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
TPQ_TR_g = ggplot(Means_df, aes(x=Session, y=Value, group = Groups, color = Groups))+
  geom_line(position = position_dodge(0.05),size=1)+
  geom_errorbar(aes(ymax = Value+SEM,ymin = Value-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  geom_text(mapping = aes(x=1.5,y=0,label = scale_names), color = "black")+
  scale_color_manual(values = c("#6cBE58","#C33E3B","#4EA3DF"))+
  ggtitle("TPQ Scales", subtitle = "Comparing MDD Groups by Response")+
  facet_grid(scale~subscale,scale = "free")

# Plot everything for all groups
TPQ_scores_all = TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`) & TPQ_scores$Session=="Test"),]

# Create array for the Test-retest data set(TPQ_scores_TR)
Means_array_all = sapply(all_scales, function(x)
  SMeans(TPQ_scores_all, DV = x, IVs = c("Diagnosis"),
         GroupBy = "Diagnosis"), simplify = FALSE)

# Bind all data frames into one
Means_df_all= Means_array_all[[1]]
colnames(Means_df_all) = c("Value","SEM","Count","Diagnosis","Groups")
for (i in 2:length(Means_array_all)){
  new_df = Means_array_all[[i]]
  colnames(new_df) = c("Value","SEM","Count","Diagnosis","Groups")
  Means_df_all = rbind(Means_df_all,new_df)
}

# Add scales and subscales columns
Means_df_all$subscale = factor(rep(c("Scale",paste("Sub_",1:4)), each = 5))
Means_df_all$scale = factor(rep(scales, each = 25))
Means_df_all$scale_names = ""
Means_df_all$scale_names[seq(from = 1, to = 75, by = 5)] = scale_names

# plot the result

TPQ_all_g = ggplot(Means_df_all, aes(x=Diagnosis, y=Value, group = Groups, fill = Groups))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax = Value+SEM,ymin = Value-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  geom_text(mapping = aes(x=3,y=20,label = scale_names), color = "black")+
  scale_fill_manual(values = c("#4EA3DF","#6cBE58","#C33E3B","#808CA3","#B9B0AB"))+
  ggtitle("TPQ Scales", subtitle = "Comparing Groups by Diagnosis")+
  facet_grid(scale~subscale,scale = "free")








## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#            Effect of Diagnosis and Treatment on ICs              #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ##
#projections = ICASSO_all$projections
#projections$`Initial ID` = TPQ_scores$`Initial ID`[match(projections$ID,TPQ_scores$`Final ID`)]
#projections = projections[complete.cases(projections$`Initial ID`),]
#projections_melted = melt(projections, measure.vars = paste0("IC",1:15), variable.name = "IC", value.name = "Projection")
#All_means_pro = SMeans(projections_melted, DV = "Projection", IVs = c("Diagnosis","IC"), GroupBy = "Diagnosis")
#All_means_pro$IC = factor(All_means$IC, levels=paste0("IC",1:15))

# sourced_TPQ
sourced_TPQ = TPQ_scores[,c("Initial ID","Session","Diagnosis","Response", paste0("Q",1:100))]
for (x in 1:15){
  mult_matrix = matrix(rep(sources[x,],nrow(sourced_TPQ)), nrow = nrow(sourced_TPQ), byrow = TRUE)
  sourced_TPQ[paste0("S",x)] =  apply(sourced_TPQ[,5:104]*mult_matrix,1,sum, na.rm= TRUE)
}
sourced_melted = melt(sourced_TPQ, measure.vars = paste0("S",1:15), variable.name = "Source", value.name = "Value")
All_means_s = SMeans(sourced_melted, DV = "Value", IVs = c("Diagnosis","Source"), GroupBy = "Diagnosis")
All_means_s$Source = factor(All_means_s$Source, levels=paste0("S",1:15))


IC_Qs = as.data.frame(sapply(TPQ_IC_Questions, function(item) rowSums(TPQ_scores[,item], na.rm = TRUE)))
IC_Qs[,c("Initial ID","Session","Diagnosis","Response")] = TPQ_scores[,c("Initial ID","Session","Diagnosis","Response")]
IC_Qs = IC_Qs[complete.cases(IC_Qs$`Initial ID`),]
IC_melted = melt(IC_Qs, id.vars = c("Initial ID","Session","Diagnosis","Response"), variable.name = "IC", value.name = "Score")
All_means = SMeans(IC_melted, DV = "Score", IVs = c("Diagnosis","IC"), GroupBy = "Diagnosis")
All_means$IC = factor(All_means$IC, levels=names(TPQ_IC_Questions))

IC_TR=IC_Qs[complete.cases(IC_Qs$Response),]
IC_TR = MatchDelete(IC_TR,"Session",common="Initial ID")
IC_melted_TR = melt(IC_TR, measure.vars = names(TPQ_IC_Questions), variable.name = "IC", value.name = "Score")
All_TR_means = SMeans(IC_melted_TR, DV = "Score", IVs = c("Session","Response","IC"), GroupBy = "Response")
All_TR_means$IC = factor(All_TR_means$IC, levels=names(TPQ_IC_Questions))

IC_all_s = ggplot(All_means_s, aes(x= Diagnosis, y= Value, group = Diagnosis, fill = Diagnosis))+
  geom_bar(stat = "identity")+geom_errorbar(aes(ymax = Value+SEM,ymin = Value-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  scale_fill_manual(values = c("#4EA3DF","#6cBE58","#C33E3B","#808CA3","#B9B0AB"))+
  ggtitle("ICs", subtitle = "Comparing TPQ Scores of Diagnoses Groups")+
  facet_wrap(~Source,scale = "free")



IC_all_g = ggplot(All_means, aes(x= Diagnosis, y= Score, group = Diagnosis, fill = Diagnosis))+
  geom_bar(stat = "identity")+geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  scale_fill_manual(values = c("#4EA3DF","#6cBE58","#C33E3B","#808CA3","#B9B0AB"))+
  ggtitle("ICs", subtitle = "Comparing TPQ Scores of Diagnoses Groups")+
  facet_wrap(~IC,scale = "free")

IC_TR_g = ggplot(All_TR_means, aes(x= Session, y= Score, group = Response, color = Response))+
  geom_line(stat = "identity", size =2,position = position_dodge(0.05))+
  geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1.2,position = position_dodge(0.05))+
  TypicalTheme+
  scale_color_manual(values = c("#00A550","#Eb4C42","#0087BD"))+
  ggtitle("ICs", subtitle = "Comparing MDD Groups by Response")+
  facet_wrap(~IC,scale = "free")


IC_all_pro = ggplot(All_means_pro, aes(x= Diagnosis, y= Projection, group = Diagnosis, fill = Diagnosis))+
  geom_bar(stat = "identity")+geom_errorbar(aes(ymax = Projection+SEM,ymin = Projection-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  scale_fill_manual(values = c("#4EA3DF","#6cBE58","#C33E3B","#808CA3","#B9B0AB"))+
  ggtitle("ICs", subtitle = "Comparing TPQ Scores of Diagnoses Groups")+
  facet_wrap(~IC,scale = "free")

IC_TR_g_pro = ggplot(All_TR_means_pro, aes(x= Session, y= Projection, group = Response, color = Response))+
  geom_line(stat = "identity", size =2,position = position_dodge(0.05))+
  geom_errorbar(aes(ymax = Projection+SEM,ymin = Projection-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1.2,position = position_dodge(0.05))+
  TypicalTheme+
  scale_color_manual(values = c("#00A550","#Eb4C42","#0087BD"))+
  ggtitle("ICs", subtitle = "Comparing MDD Groups by Response")+
  facet_wrap(~IC,scale = "free")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#          Effect of Diagnosis and Treatment on Factors            #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ##

factored_TPQ = TPQ_scores[,c("Initial ID","Session","Diagnosis","Response", paste0("Q",1:100))]
for (x in 1:12){
  mult_matrix = matrix(rep(factanal$loadings[x,],nrow(factored_TPQ)), nrow = nrow(factored_TPQ), byrow = TRUE)
  factored_TPQ[paste0("F",x)] =  apply(factored_TPQ[,5:104]*mult_matrix,1,sum, na.rm= TRUE)
}
factored_melted = melt(factored_TPQ, measure.vars = paste0("F",1:12), variable.name = "Factor", value.name = "Value")
All_means_s = SMeans(sourced_melted, DV = "Value", IVs = c("Diagnosis","Source"), GroupBy = "Diagnosis")
All_means_s$Source = factor(All_means_s$Source, levels=paste0("S",1:15))



Factor_Qs = as.data.frame(sapply(Factors_list, function(item) rowSums(TPQ_scores[,item], na.rm = TRUE)))
Factor_Qs[,c("Initial ID","Session","Diagnosis","Response")] = TPQ_scores[,c("Initial ID","Session","Diagnosis","Response")]
Factor_Qs = Factor_Qs[complete.cases(Factor_Qs$`Initial ID`),]
Factor_melted = melt(Factor_Qs, measure.vars = paste0("F",1:11), variable.name = "Factor", value.name = "Score")
All_means = SMeans(Factor_melted, DV = "Score", IVs = c("Diagnosis","Factor"), GroupBy = "Diagnosis")
All_means$Factor = factor(All_means$Factor, levels=paste0("F",1:11))



Factors_TR=Factor_Qs[complete.cases(Factor_Qs$Response),]
Factors_TR = MatchDelete(Factors_TR,"Session",common="Initial ID")
Factors_melted_TR = melt(Factors_TR, measure.vars = paste0("F",1:11), variable.name = "Factor", value.name = "Score")
All_TR_means = SMeans(Factors_melted_TR, DV = "Score", IVs = c("Session","Response","Factor"), GroupBy = "Response")
All_TR_means$Factor = factor(All_TR_means$Factor, levels=paste0("F",1:11))



factor_all_g = ggplot(All_means, aes(x= Diagnosis, y= Score, group = Groups, fill = Groups))+
  geom_bar(stat = "identity")+geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1,position = position_dodge(0.05))+
  TypicalTheme+
  scale_fill_manual(values = c("#4EA3DF","#6cBE58","#C33E3B","#808CA3","#B9B0AB"))+
  ggtitle("Factors", subtitle = "Comparing Groups by Diagnosis")+
  facet_wrap(~Factor,scale = "fixed")

factor_TR_g = ggplot(All_TR_means, aes(x= Session, y= Score, group = Groups, color = Groups))+
  geom_line(stat = "identity", size =2,position = position_dodge(0.05))+
  geom_errorbar(aes(ymax = Score+SEM,ymin = Score-SEM), position = position_dodge(0.05), width = 0.15,color="black")+
  geom_point(shape = 21, color = "black",fill = "white",size = 1.2,position = position_dodge(0.05))+
  TypicalTheme+
  scale_color_manual(values = c("#00A550","#Eb4C42","#0087BD"))+
  ggtitle("Factors", subtitle = "Comparing MDD Groups by Response")+
  facet_wrap(~Factor,scale = "free")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#           plot the chosen questions with the results             #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ##
multiplot(TPQ_q_g,TPQ_TR_g, layout = matrix(c(1,2,2,2), nrow = 1))
multiplot(IC_q_g,IC_TR_g, layout = matrix(c(1,1,2,2,2,2), nrow = 1))
multiplot(factor_q_g,factor_TR_g, layout = matrix(c(1,1,2,2,2,2), nrow = 1))

multiplot(TPQ_q_g,TPQ_all_g, layout = matrix(c(1,2,2,2), nrow = 1))
multiplot(IC_q_g,IC_all_g, layout = matrix(c(1,1,2,2,2,2), nrow = 1))
multiplot(factor_q_g,factor_all_g, layout = matrix(c(1,1,2,2,2,2), nrow = 1))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ## 
#       Correspondance between Factors, ICs and temperaments       #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## ##  ## ##  ## ##


# Now we find the correlations between factors and ICs, and factors with subtemperaments
Factors_Temps = data.frame(Factor = paste0("F",1:11), Temp_coocurrance_ratio = NA, Temp = NA,
                           IC_coocurrance_ratio = NA, IC = NA)
TPQQuestions_2 = sapply(TPQQuestions[1:12], function(x) paste0("Q", x))
for (i in 1:length(Factors_list)){
  item = Factors_list[[i]]
  sim_count = sapply(TPQQuestions_2, function(TPQ_Q_group) sum(item %in% TPQ_Q_group))
  max_ratio =max(sim_count)/length(item)
  max_temp = names(TPQQuestions_2)[which.max(sim_count)]
  Factors_Temps[i,c("Temp_coocurrance_ratio", "Temp")] = c(round(max_ratio,3),max_temp)
}

for (i in 1:length(Factors_list)){
  item = Factors_list[[i]]
  sim_count = sapply(TPQ_IC_Questions, function(TPQ_IC_group) sum(item %in% TPQ_IC_group))
  max_ratio =max(sim_count)/length(item)
  max_IC = names(TPQ_IC_Questions)[which.max(sim_count)]
  Factors_Temps[i,c("IC_coocurrance_ratio", "IC")] = c(round(max_ratio,3),max_IC)
}

# Find correlations
Temp_sums = sapply(TPQQuestions_2, function(item) rowSums(TPQ_scores[,item], na.rm = TRUE))
IC_sums = sapply(TPQ_IC_Questions, function(item) rowSums(TPQ_scores[,item], na.rm = TRUE))
Factor_sums = sapply(Factors_list, function(item) rowSums(TPQ_scores[,item], na.rm = TRUE))
All_sums = data.frame(cbind(Temp_sums,IC_sums))#All_sums = data.frame(cbind(Temp_sums,IC_sums,Factor_sums))
# find the threshold of inclusion, uncomment if you want to calculate again

IC_cooc = find_threshold(All_sums[1:12], All_sums[13:22], Comparison = "Correlations",G1="TPQ Questions",
               G2="IC-related Questions", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)

apply(IC_cooc$cor_list,2,max)

cor_mat=cor(All_sums)
cor_mat = ifelse(abs(cor_mat)>2*0.2,cor_mat,0)
corrplot(cor_mat, sig = 0.0005,insig = "p-value", hclust.method = "average", title = "\n\n\nCorrelation Between ICs & Factors Scores With TPQ Scores")
corrplot(cor_mat[(length(TPQQuestions_2)+1):ncol(cor_mat),1:length(TPQQuestions_2)], sig = 0.0005,insig = "p-value", hclust.method = "average", title = "\n\n\nCorrelation Between ICs Scores With TPQ Scores")
corrplot(cor_mat[1:length(TPQQuestions_2),(length(TPQQuestions_2)+1):ncol(cor_mat)], sig = 0.0005,insig = "p-value", hclust.method = "average", title = "\n\n\nCorrelation Between Factors Scores With TPQ Scores\n\n")

# Show this as percentage of co-occurance matrices

loadings_cor = matrix(nrow = 12, ncol = 15)
dimnames(loadings_cor) = list(paste0("F",1:12), names(TPQQuestions))
names(dimnames(loadings_cor)) = c("Factor Loading","Scales")

for (i in 1:nrow(loadings)){
  cors = sapply(1:15, function(x) mean(which(loadings[i,]==1)%in% which(tpq_q_cloninger[x,]==1) ))
  loadings_cor[i,] = cors
}

loadings_cor = t(ifelse(loadings_cor<0.7, 0,round(loadings_cor,2)))

IC_questions_cor = matrix(nrow = nrow(IC_questions), ncol =15)
dimnames(IC_questions_cor) = list(paste0("IC",1:n_ICs), names(TPQQuestions)[1:15])
names(dimnames(IC_questions_cor)) = c("ICs","Scales")


for (i in 1:nrow(IC_questions)){
  cors = sapply(1:15, function(x) mean(which(IC_questions[i,]==1) %in% which(tpq_q_cloninger[x,]==1)))
  IC_questions_cor[i,] = cors
}

IC_questions_cor = t(ifelse(IC_questions_cor<0.7, 0,round(IC_questions_cor,2)))

#plot the final result
corrplot(cbind(IC_questions_cor[,!low_q_count[1:15]],loadings_cor[,!low_q_count_factor[1:12]]), title = "\n\n\nPercentage of Questions of ICs and Factors Included in TPQ Scales")
corrplot(IC_questions_cor[,!low_q_count[1:15]], title = "\n\n\nPercentage of Questions of ICs Included in TPQ Scales")
corrplot(loadings_cor[,!low_q_count_factor[1:12]], title = "\n\n\nPercentage of Questions of Factors Included in TPQ Scales")


# Find the correlations between factor loadings and TPQ scores

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#             Crobnach Alpha                  #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

return_alpha = function(questions) {
  # returns 3 items, the name of the subscale,
  # c.alpha of males and c.alpha of females
  
  df = TPQ_scores[questions]
  M_alpha = cronbach(df[TPQ_scores$Gender == "Male", ])$alpha
  F_alpha = cronbach(df[TPQ_scores$Gender == "Female", ])$alpha
  All_alpha = cronbach(df)$alpha
  return(c( M_alpha, F_alpha, All_alpha))
}
TPQ_alpha_values = sapply(TPQQuestions, function(x) return_alpha(paste0("Q", x)))
TPQ_alpha_cloninger = as.data.frame(t(TPQ_alpha_values))
colnames(TPQ_alpha_cloninger) = c("Males","Females", "All")


# Find c.alpha for ICs
IC_alpha_values = sapply(ICs,function(x) return_alpha(TPQ_IC_Questions[[x]]))
IC_alpha_values = as.data.frame(t(IC_alpha_values))
colnames(IC_alpha_values) = c("Males","Females", "All")

# Find c.alpha for ICs
F_alpha_values = sapply(paste0("F",1:11),function(x) return_alpha(Factors_list[[x]]))
F_alpha_values = as.data.frame(t(F_alpha_values))
colnames(F_alpha_values) = c("Males","Females", "All")


# Find c_alpha for factored data
cronbach(factored_TPQ[,paste0("F",1:12)])

# compare both cronbachs
apply(IC_alpha_values, 2, mean, na.rm = TRUE)
apply(F_alpha_values, 2, mean, na.rm = TRUE)
apply(TPQ_alpha_cloninger[1:12,], 2, mean, na.rm = TRUE)






# Comare different ICASSOs (Test vs Retest and leave-one-out) ----------------------------------------------
# Import HC data
hc_pro_test = ICASSO_hc_test$projections[,1:(8+sum(ICASSO_hc_test$scores>0.89))]
hc_pro_retest = ICASSO_hc_retest$projections[,1:(8+sum(ICASSO_hc_retest$scores>0.89))]
hc_pro_test$Initial_ID = TPQ_scores$`Initial ID`[match(hc_pro_test$ID,TPQ_scores$`Final ID`)]
hc_pro_retest$Initial_ID = TPQ_scores$`Initial ID`[match(hc_pro_retest$ID,TPQ_scores$`Final ID`)]
hc_pro_test = hc_pro_test[match(hc_pro_retest$Initial_ID,hc_pro_test$Initial_ID),]

hc_src_test = data.frame(t(ICASSO_hc_test$sources[ICASSO_hc_test$scores>0.89,-101]))
hc_src_retest = data.frame(t(ICASSO_hc_retest$sources[ICASSO_hc_retest$scores>0.89,-101]))
colnames(hc_src_test) = paste0("IC",1:ncol(hc_src_test))
colnames(hc_src_retest) = paste0("IC",1:ncol(hc_src_retest))

# Import MDD data
mdd_pro_test = ICASSO_mdd_test$projections[,1:(8+sum(ICASSO_mdd_test$scores>0.89))]
mdd_pro_retest = ICASSO_mdd_retest$projections[,1:(8+sum(ICASSO_mdd_retest$scores>0.89))]
mdd_pro_test$Initial_ID = TPQ_scores$`Initial ID`[match(mdd_pro_test$ID,TPQ_scores$`Final ID`)]
mdd_pro_retest$Initial_ID = TPQ_scores$`Initial ID`[match(mdd_pro_retest$ID,TPQ_scores$`Final ID`)]
mdd_pro_test = mdd_pro_test[match(mdd_pro_retest$Initial_ID,mdd_pro_test$Initial_ID),]

mdd_src_test = data.frame(t(ICASSO_mdd_test$sources[ICASSO_mdd_test$scores>0.89,-101]))
mdd_src_retest = data.frame(t(ICASSO_mdd_retest$sources[ICASSO_mdd_retest$scores>0.89,-101]))
colnames(mdd_src_test) = paste0("IC",1:ncol(mdd_src_test))
colnames(mdd_src_retest) = paste0("IC",1:ncol(mdd_src_retest))


# Projection correlations
cor_hc_pro = find_threshold(rbind(hc_pro_test,hc_pro_test), rbind(hc_pro_retest,hc_pro_retest),Comparison = "Projections",G1="HC at Test",
                            G2="HC at Retest", method = "spearman", n_shuffle = 1000, threshold_choice = 0.95)


cor_mdd_pro = find_threshold(mdd_pro_test, mdd_pro_retest,Comparison = "Projections",G1="MDD at Test",
                            G2="MDD at Retest", method = "spearman", n_shuffle = 400, threshold_choice = 0.95)

# Sources correlations
cor_hc_src = find_threshold(hc_src_test, hc_src_retest,Comparison = "Sources",G1="HC at Test",
                            G2="HC at Retest", method = "spearman", n_shuffle = 400, threshold_choice = 0.95)

cor_mdd_src = find_threshold(mdd_src_test, mdd_src_retest,Comparison = "Sources",G1="MDD at Test",
                             G2="MDD at Retest", method = "spearman", n_shuffle = 400, threshold_choice = 0.95)

# sources correlation between MDD and HC - validity check
cor_test_src = find_threshold(mdd_src_test, hc_src_test,Comparison = "Sources",G1="MDD at Test",
                              G2="HC at Test", method = "spearman", n_shuffle = 400, threshold_choice = 0.95)
cor_retest_src = find_threshold(mdd_src_retest, hc_src_retest,Comparison = "Sources",G1="MDD at Retest",
                              G2="HC at Retest", method = "spearman", n_shuffle = 400, threshold_choice = 0.95)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Leave-one out correlations
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# import data
all_src = data.frame(t(ICASSO_all$sources[ICASSO_all$scores>0.849,-101]))
colnames(all_src) = paste0("IC",1:ncol(all_src))

hc_src = data.frame(t(ICASSO_hc$sources[ICASSO_hc$scores>0.849,-101]))
no_hc_src = data.frame(t(ICASSO_noHC$sources[ICASSO_noHC$scores>0.849,-101]))
colnames(hc_src) = paste0("IC",1:ncol(hc_src))
colnames(no_hc_src) = paste0("IC",1:ncol(no_hc_src))

mdd_src = data.frame(t(ICASSO_mdd$sources[ICASSO_mdd$scores>0.849,-101]))
no_mdd_src = data.frame(t(ICASSO_noMDD$sources[ICASSO_noMDD$scores>0.849,-101]))
colnames(mdd_src) = paste0("IC",1:ncol(mdd_src))
colnames(no_mdd_src) = paste0("IC",1:ncol(no_mdd_src))

gad_src = data.frame(t(ICASSO_gad$sources[ICASSO_gad$scores>0.849,-101]))
no_gad_src = data.frame(t(ICASSO_noGAD$sources[ICASSO_noGAD$scores>0.849,-101]))
colnames(gad_src) = paste0("IC",1:ncol(gad_src))
colnames(no_gad_src) = paste0("IC",1:ncol(no_gad_src))

ptsd_src = data.frame(t(ICASSO_ptsd$sources[ICASSO_ptsd$scores>0.849,-101]))
no_ptsd_src = data.frame(t(ICASSO_noPTSD$sources[ICASSO_noPTSD$scores>0.849,-101]))
colnames(ptsd_src) = paste0("IC",1:ncol(ptsd_src))
colnames(no_ptsd_src) = paste0("IC",1:ncol(no_ptsd_src))

tnp_src = data.frame(t(ICASSO_tnp$sources[ICASSO_tnp$scores>0.849,-101]))
no_tnp_src = data.frame(t(ICASSO_noTNP$sources[ICASSO_noTNP$scores>0.849,-101]))
colnames(tnp_src) = paste0("IC",1:ncol(tnp_src))
colnames(no_tnp_src) = paste0("IC",1:ncol(no_tnp_src))

# for HC
cor_hc_src = find_threshold(all_src, no_hc_src,Comparison = "Sources",G1="Decomposition with HC Alone",
                              G2="Decomposition without HC", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)

# for MDD
cor_mdd_src = find_threshold(all_src, no_mdd_src,Comparison = "Sources",G1="Decomposition with MDD Alone",
                              G2="Decomposition without MDD", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)

# for GAD
cor_gad_src = find_threshold(all_src, no_gad_src,Comparison = "Sources",G1="Decomposition with GAD Alone",
                             G2="Decomposition without GAD", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)

# for TNP
cor_tnp_src = find_threshold(all_src, no_tnp_src,Comparison = "Sources",G1="Decomposition with TNP Alone",
                             G2="Decomposition without TNP", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)

# for PTSD
cor_ptsd_src = find_threshold(all_src, no_ptsd_src,Comparison = "Sources",G1="Decomposition with PTSD Alone",
                             G2="Decomposition without PTSD", method = "spearman", n_shuffle = 1000, threshold_choice = 0.999)





# Reconstructed Data Correlation (Test vs Retest) --------------------------------------------
hc_test_reco = ICASSO_hc_test$reco[ICASSO_hc_test$scores>0.849]
hc_retest_reco = ICASSO_hc_retest$reco[ICASSO_hc_retest$scores>0.849]

# add initial ID to hc_retest_reco
for (i in 1:length(hc_retest_reco)){
  hc_retest_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(hc_retest_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}


# add initial ID to hc_test_reco
for (i in 1:length(hc_test_reco)){
  hc_test_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(hc_test_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}

# rearrange_inds is the indices needed to shuffle the test data frames to match those of retest
rearrange_retest = match(hc_retest_reco[[1]]$Initial_ID,hc_test_reco[[1]]$Initial_ID)

# rearrange the IDs at test to match those of retest
for (i in 1:length(hc_test_reco)){
  hc_test_reco[[i]] = hc_test_reco[[i]][rearrange_retest,]
}


# Create the matrix of correlation. It should consist, for now, from vectorized data frames
test_dfs = sapply(1:length(hc_test_reco), function(x) melt(hc_test_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_hc_test = data.frame(Subject = test_dfs[[1]]$ID, IC1 = test_dfs[[1]]$value)
for (i in 2:length(test_dfs)){
  final_hc_test[[paste0("IC",i)]] = test_dfs[[i]]$value
}

retest_dfs = sapply(1:length(hc_retest_reco), function(x) melt(hc_retest_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_hc_retest = data.frame(Subject = retest_dfs[[1]]$ID, IC1 = retest_dfs[[1]]$value)
for (i in 2:length(retest_dfs)){
  final_hc_retest[[paste0("IC",i)]] = retest_dfs[[i]]$value
}


cor_matrix = find_threshold(final_hc_test, final_hc_retest,Comparison = "Reconstructed Data",G1="HC at Test",
                            G2="HC at Retest", method = "spearman", n_shuffle = 50, threshold_choice = 0.999)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# For MDD
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

mdd_test_reco = ICASSO_mdd_test$reco[ICASSO_mdd_test$scores>0.849]
mdd_retest_reco = ICASSO_mdd_retest$reco[ICASSO_mdd_retest$scores>0.849]

# add initial ID to mdd_retest_reco
for (i in 1:length(mdd_retest_reco)){
  mdd_retest_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(mdd_retest_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}


# add initial ID to mdd_test_reco
for (i in 1:length(mdd_test_reco)){
  mdd_test_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(mdd_test_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}

# rearrange_inds is the indices needed to shuffle the test data frames to match those of retest
rearrange_retest = match(mdd_retest_reco[[1]]$Initial_ID,mdd_test_reco[[1]]$Initial_ID)

# rearrange the IDs at test to match those of retest
for (i in 1:length(mdd_test_reco)){
  mdd_test_reco[[i]] = mdd_test_reco[[i]][rearrange_retest,]
}


# Create the matrix of correlation. It should consist, for now, from vectorized data frames
test_dfs = sapply(1:length(mdd_test_reco), function(x) melt(mdd_test_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_mdd_test = data.frame(Subject = test_dfs[[1]]$ID, IC1 = test_dfs[[1]]$value)
for (i in 2:length(test_dfs)){
  final_mdd_test[[paste0("IC",i)]] = test_dfs[[i]]$value
}

retest_dfs = sapply(1:length(mdd_retest_reco), function(x) melt(mdd_retest_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_mdd_retest = data.frame(Subject = retest_dfs[[1]]$ID, IC1 = retest_dfs[[1]]$value)
for (i in 2:length(retest_dfs)){
  final_mdd_retest[[paste0("IC",i)]] = retest_dfs[[i]]$value
}

cor_matrix_mdd = find_threshold(final_mdd_test, final_mdd_retest,Comparison = "Reconstructed Data",G1="MDD at Test",
                            G2="MDD at Retest", method = "spearman", n_shuffle = 50, threshold_choice = 0.999)


#old_threshold = find_threshold_old(t(final_mdd_test[,2:length(final_mdd_test)]), t(final_mdd_retest[,2:length(final_mdd_retest)]),
#                   150, cor_method = "spearman")

# Correlate a list of HC at test and MDD at test - Confirmatory analysis
conf_mdd = final_mdd_test[-seq(91, 9100, by=91),]# remove the last subject of each dataset
conf_hc = final_hc_test
cor_mat_test = find_threshold(conf_mdd, conf_hc,Comparison = "Reconstructed Data",G1="MDD",
                              G2="HC", method = "spearman", n_shuffle = 50, threshold_choice = 0.999)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# For MDD-responders
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
resp_test_reco = ICASSO_resp_test$reco[ICASSO_resp_test$scores>0.849]
resp_retest_reco = ICASSO_resp_retest$reco[ICASSO_resp_retest$scores>0.849]

# add initial ID to mdd_retest_reco
for (i in 1:length(resp_retest_reco)){
  resp_retest_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(resp_retest_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}


# add initial ID to resp_test_reco
for (i in 1:length(resp_test_reco)){
  resp_test_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(resp_test_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}

# rearrange_inds is the indices needed to shuffle the test data frames to match those of retest
rearrange_retest = match(resp_retest_reco[[1]]$Initial_ID,resp_test_reco[[1]]$Initial_ID)

# rearrange the IDs at test to match those of retest
for (i in 1:length(resp_test_reco)){
  resp_test_reco[[i]] = resp_test_reco[[i]][rearrange_retest,]
}


# Create the matrix of correlation. It should consist, for now, from vectorized data frames
test_dfs = sapply(1:length(resp_test_reco), function(x) melt(resp_test_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_resp_test = data.frame(Subject = test_dfs[[1]]$ID, IC1 = test_dfs[[1]]$value)
for (i in 2:length(test_dfs)){
  final_resp_test[[paste0("IC",i)]] = test_dfs[[i]]$value
}

retest_dfs = sapply(1:length(resp_retest_reco), function(x) melt(resp_retest_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_resp_retest = data.frame(Subject = retest_dfs[[1]]$ID, IC1 = retest_dfs[[1]]$value)
for (i in 2:length(retest_dfs)){
  final_resp_retest[[paste0("IC",i)]] = retest_dfs[[i]]$value
}

cor_matrix_resp = find_threshold(final_resp_test, final_resp_retest,Comparison = "Reconstructed Data",G1="Responders at Test",
                                 G2="Responders at Retest", method = "spearman", n_shuffle = 50, threshold_choice = 0.999)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# For MDD-nonresponders
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nonresp_test_reco = ICASSO_nonresp_test$reco[ICASSO_nonresp_test$scores>0.849]
nonresp_retest_reco = ICASSO_nonresp_retest$reco[ICASSO_nonresp_retest$scores>0.849]

# add initial ID to mdd_retest_reco
for (i in 1:length(nonresp_retest_reco)){
  nonresp_retest_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(nonresp_retest_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}


# add initial ID to nonresp_test_reco
for (i in 1:length(nonresp_test_reco)){
  nonresp_test_reco[[i]]$Initial_ID = TPQ_scores$`Initial ID`[match(nonresp_test_reco[[i]]$ID,TPQ_scores$`Final ID`)]
}

# rearrange_inds is the indices needed to shuffle the test data frames to match those of retest
rearrange_retest = match(nonresp_retest_reco[[1]]$Initial_ID,nonresp_test_reco[[1]]$Initial_ID)

# rearrange the IDs at test to match those of retest
for (i in 1:length(nonresp_test_reco)){
  nonresp_test_reco[[i]] = nonresp_test_reco[[i]][rearrange_retest,]
}


# Create the matrix of correlation. It should consist, for now, from vectorized data frames
test_dfs = sapply(1:length(nonresp_test_reco), function(x) melt(nonresp_test_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_nonresp_test = data.frame(Subject = test_dfs[[1]]$ID, IC1 = test_dfs[[1]]$value)
for (i in 2:length(test_dfs)){
  final_nonresp_test[[paste0("IC",i)]] = test_dfs[[i]]$value
}

retest_dfs = sapply(1:length(nonresp_retest_reco), function(x) melt(nonresp_retest_reco[[x]][,c(1,9:108)], id.vars = "ID"), simplify = FALSE)
final_nonresp_retest = data.frame(Subject = retest_dfs[[1]]$ID, IC1 = retest_dfs[[1]]$value)
for (i in 2:length(retest_dfs)){
  final_nonresp_retest[[paste0("IC",i)]] = retest_dfs[[i]]$value
}

cor_matrix_nonresp = find_threshold(final_nonresp_test, final_nonresp_retest,Comparison = "Reconstructed Data",G1="NonResponders at Test",
                                 G2="NonResponders at Retest", method = "spearman", n_shuffle = 50, threshold_choice = 0.999)


# IC Names--------------------------------------------
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#                     IC Names                      #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##

# Test the Treshold Method ------------------------------------------------------
or_data_1 = ICASSO_hc$sources[ICASSO_hc$scores > 0.85,1:100]
or_data_2 = ICASSO_noHC$sources[ICASSO_noHC$scores > 0.85,1:100]
n_shuffle=1250

# registerDoParallel(cores=8)
# cor_matrices = foreach(i=1:8) %dopar% find_threshold(or_data_1,or_data_2, n_shuffle)
# cor_matrix = abind(cor_matrices, along = 1)

cor_threshold = sapply(1:nrow(or_data_2), function (n) sapply(1:nrow(or_data_1), function (m) quantile(cor_matrix[,m,n],0.999)))

# delete later
# find the shuffle threshold and whether or not there is a cieling for this
threshold_array = rep(list(matrix(nrow = nrow(or_data_1), ncol = nrow(or_data_2))),2)
shuffles = c(15,50,100,500,1000,2500,5000,7500)
for (i in 1:length(shuffles)){
  n_shuffle = shuffles[i]
  registerDoParallel(cores=8)
  cor_matrices = foreach(i=1:8) %dopar% find_threshold(or_data_1,or_data_2, n_shuffle)
  cor_matrix = abind(cor_matrices, along = 1)
  
  cor_threshold = sapply(1:nrow(or_data_2), function (n)
    sapply(1:nrow(or_data_1), function (m)
      quantile(cor_matrix[, m, n], 0.999)))
  threshold_array[[i]] = cor_threshold
}

data_threshold = data.frame(correlation=c(abind(threshold_array)))
data_threshold$Data1 = factor(rep(paste0("C",1:nrow(or_data_1)),nrow(or_data_2)*length(thresholds)), levels = paste0("C",1:nrow(or_data_1)))
data_threshold$Data2 = factor(rep(rep(paste0("C",1:nrow(or_data_2)),each = nrow(or_data_1)),length(thresholds)),paste0("C",1:nrow(or_data_2)))
data_threshold$n_shuflle = rep(shuffles, each = nrow(or_data_1)*nrow(or_data_2))
for(d1 in levels(factor(data_threshold$Data1))){
  for (d2 in levels(factor(data_threshold$Data2))){
    data_threshold$max[data_threshold$Data1 == d1&data_threshold$Data2 == d2] =
      max(data_threshold$correlation[data_threshold$Data1 == d1&data_threshold$Data2 == d2])
  }
}


ggplot(data_threshold,aes(x=n_shuflle,y= correlation))+geom_line(size = 1.2)+
  geom_hline(aes(yintercept = max), color = "red")+
  facet_grid(Data1~Data2)+theme()
# end delete
# Improve upon the threshold method -------------------------------------------
# A function to find the correlation threshold of two data frames
# both data frames should have the same number of columns
sources_hc_test = data.frame(t(ICASSO_hc_test$sources[ICASSO_hc$scores > 0.85,1:100]))
sources_hc_retest = data.frame(t(ICASSO_hc_retest$sources[ICASSO_hc$scores > 0.85,1:100]))

sources_mdd_test = data.frame(t(ICASSO_mdd_test$sources[ICASSO_mdd_test$scores > 0.85,1:100]))
sources_mdd_retest = data.frame(t(ICASSO_mdd_retest$sources[ICASSO_mdd_retest$scores > 0.85,1:100]))

sources_nonresp_test = data.frame(t(ICASSO_nonresp_test$sources[ICASSO_nonresp_test$scores > 0.85,1:100]))
sources_nonresp_retest = data.frame(t(ICASSO_nonresp_retest$sources[ICASSO_nonresp_retest$scores > 0.85,1:100]))


#sources_hc_test = read.csv("/home/asawalma/Desktop/sources_hc_retest.xlsx")
#sources_hc_retest = read.csv("/home/asawalma/Desktop/sources_hc_retest.xlsx")


sh_th = find_threshold(or_data_1 =sources_hc_test,or_data_2 = sources_hc_retest,
               Comparison = "Sources",G1="HC at test", G2="HC at retest",
               method = "spearman", n_shuffle =1000, threshold_choice = 0.95, simulate = FALSE)



sim_th = find_threshold(or_data_1 =sources_hc_test,or_data_2 = sources_hc_retest,
                        Comparison = "Sources",G1="HC at test", G2="HC at retest",
                        method = "spearman", n_shuffle = 100, threshold_choice = 0.999, simulate = TRUE)



hc_mdd_th = find_threshold(or_data_1 =sources_mdd_test,or_data_2 = sources_hc_test,
                           Comparison = "Sources",G1="MDD at test", G2="HC at test",
                           method = "spearman", n_shuffle = 10000, threshold_choice = 0.999, simulate = TRUE)

mdd_th = find_threshold(or_data_1 =sources_mdd_test,or_data_2 = sources_mdd_retest,
                        Comparison = "Sources",G1="MDD at test", G2="MDD at retest",
                        method = "spearman", n_shuffle = 10000, threshold_choice = 0.999, simulate = TRUE)

resp_th = find_threshold(or_data_1 =sources_nonresp_test,or_data_2 = sources_nonresp_retest,
                        Comparison = "Sources",G1="NonResponders at test", G2="NonResponders at retest",
                        method = "spearman", n_shuffle = 10000, threshold_choice = 0.999, simulate = TRUE)








# Delete later
TPQ_responders = TPQ_scores[complete.cases(TPQ_scores$Response) & complete.cases(TPQ_scores$QO1)& TPQ_scores$Response == "Responder",
                            c("Initial ID","Final ID","Session","Diagnosis",paste0("QO",1:100))]
TPQ_responders = MatchDelete(TPQ_responders, "Session","Initial ID")

TPQ_nonresponders = TPQ_scores[complete.cases(TPQ_scores$Response) & complete.cases(TPQ_scores$QO1)& TPQ_scores$Response == "Non-responder",
                            c("Initial ID","Final ID","Session","Diagnosis",paste0("QO",1:100))]
TPQ_nonresponders = MatchDelete(TPQ_nonresponders, "Session","Initial ID")


write.csv(TPQ_nonresponders,"/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/TPQData_ICA_nonrespoder.xlsx")
write.csv(TPQ_responders,"/home/asawalma/git/Data/TPQ-Analysis/TPQ_data/TPQData_ICA_respoder.xlsx")


TPQ_scores_tr = TPQ_scores[c("Initial ID","Diagnosis","Response","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
TPQ_scores_tr = TPQ_scores_tr[complete.cases(TPQ_scores_tr$`Initial ID`) & complete.cases(TPQ_scores_tr$NS1),]
TPQ_scores_tr = MatchDelete(TPQ_scores_tr, "Session","Initial ID")
TPQ_scores_t = TPQ_scores_tr[TPQ_scores_tr$Session == "Test",]
TPQ_scores_r = TPQ_scores_tr[TPQ_scores_tr$Session == "Retest",]
TPQ_scores_t = TPQ_scores_t[match(TPQ_scores_r$`Initial ID`,TPQ_scores_t$`Initial ID`),]


cor_all = find_threshold(or_data_1 = TPQ_scores_t,or_data_2 = TPQ_scores_r,
                        Comparison = "Reconstructed Data",G1="HC at test", G2="HC at retest",
                        method = "spearman", n_shuffle = 1000, threshold_choice = 0.999, simulate = TRUE)

cor_hc = find_threshold(or_data_1 = TPQ_scores_t[TPQ_scores_t$Diagnosis == "HC",],
                        or_data_2 = TPQ_scores_r[TPQ_scores_t$Diagnosis == "HC",],
                        Comparison = "Reconstructed Data",G1="HC at test", G2="HC at retest",
                        method = "spearman", n_shuffle = 200, threshold_choice = 0.95, simulate = TRUE)

cor_mdd = find_threshold(or_data_1 = TPQ_scores_t[TPQ_scores_t$Diagnosis == "MDD",],
                         or_data_2 = TPQ_scores_r[TPQ_scores_r$Diagnosis == "MDD",],
                         Comparison = "Reconstructed Data",G1="MDD at test", G2="MDD at retest",
                         method = "spearman", n_shuffle = 200, threshold_choice = 0.95, simulate = TRUE)

cor_mdd_n_resp = find_threshold(or_data_1 = TPQ_scores_t[(TPQ_scores_t$Diagnosis == "MDD"&TPQ_scores_t$Response == "Non-responder"&complete.cases(TPQ_scores_t$Response)),],
                                or_data_2 = TPQ_scores_r[(TPQ_scores_r$Diagnosis == "MDD"&TPQ_scores_r$Response == "Non-responder"&complete.cases(TPQ_scores_r$Response)),],
                                Comparison = "Reconstructed Data",G1="MDD-NonResponders at test", G2="MDD-NonResponders at retest",
                                method = "spearman", n_shuffle = 200, threshold_choice = 0.95, simulate = TRUE)

cor_mdd_resp = find_threshold(or_data_1 = TPQ_scores_t[(TPQ_scores_t$Diagnosis == "MDD"&TPQ_scores_t$Response == "Responder"&complete.cases(TPQ_scores_t$Response)),],
                                or_data_2 = TPQ_scores_r[(TPQ_scores_r$Diagnosis == "MDD"&TPQ_scores_r$Response == "Responder"&complete.cases(TPQ_scores_r$Response)),],
                                Comparison = "Reconstructed Data",G1="MDD-Resp at test", G2="MDD-Resp at retest",
                                method = "spearman", n_shuffle = 200, threshold_choice = 0.95, simulate = TRUE)





# Statistics on the threshold method, sample size and numbre of iterations -------------------------------------------

# Difference between simulation and shuffling
# for data size of 100 - same 100 tested and retested
shuffles = c(50,100,500,1000,1500,2000,3000,4000,6000,8000,10000)#,12500, 15000)
TPQ_scores_r = TPQ_scores_r[1:100,]
TPQ_scores_t = TPQ_scores_t[1:100,]

sim_values_100 = sapply(shuffles, function(shuf)  find_threshold(TPQ_scores_t, TPQ_scores_r,Comparison = "Reconstructed Data",G1="HC at test",
                                                                 G2="HC at retest", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = TRUE)$standard_thresholds)
shuf_values_100 = sapply(shuffles, function(shuf)  find_threshold(TPQ_scores_t, TPQ_scores_r,Comparison = "Reconstructed Data",G1="HC at test",
                                                                  G2="HC at retest", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds)

sim_shuffle_100 = data.frame(simulations = c(sim_values_100), shuffles = c(shuf_values_100), perm = rep(shuffles, each = 5),
                             threshold_lvl = rep(c(10,50,90,95,99), length(shuffles)))
sim_dif_100 = melt(sim_shuffle_100,
                   id.vars = c("perm","threshold_lvl"),value.name = "threshold", variable.name = "type")


# for data size of 200 - same 200 tested and retested
sim_values_200 = sapply(shuffles, function(shuf)  find_threshold(rbind(TPQ_scores_t,TPQ_scores_t), rbind(TPQ_scores_r,TPQ_scores_r),Comparison = "Reconstructed Data",G1="HC at test",
                                                                 G2="HC at retest", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = TRUE)$standard_thresholds)
shuf_values_200 = sapply(shuffles, function(shuf)  find_threshold(rbind(TPQ_scores_t,TPQ_scores_t), rbind(TPQ_scores_r,TPQ_scores_r),Comparison = "Reconstructed Data",G1="HC at test",
                                                                  G2="HC at retest", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds)

sim_shuffle_200 = data.frame(simulations = c(sim_values_200), shuffles = c(shuf_values_200), perm = rep(shuffles, each = 5),
                             threshold_lvl = rep(c(10,50,90,95,99), length(shuffles)))
sim_dif_200 = melt(sim_shuffle_200,
                   id.vars = c("perm","threshold_lvl"),value.name = "threshold", variable.name = "type")




#for a size of 500 -different subjects
TPQ_scores_1 = TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`)&complete.cases(TPQ_scores$Q1)),][1:500,c("Initial ID","Diagnosis","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
TPQ_scores_2= TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`)&complete.cases(TPQ_scores$Q1)),][501:1000,c("Initial ID","Diagnosis","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
sim_values_500 = sapply(shuffles, function(shuf)
  find_threshold(TPQ_scores_1,
                 TPQ_scores_2,
                 Comparison = "Reconstructed Data",
                 G1 = "HC at test",
                 G2 = "HC at retest",
                 method = "spearman",
                 n_shuffle = shuf,
                 threshold_choice = 0.95,
                 simulate = TRUE)$standard_thresholds)
shuf_values_500 = sapply(shuffles, function(shuf)  find_threshold(TPQ_scores_1, TPQ_scores_2,Comparison = "Reconstructed Data",G1="random 500",
                                                                  G2="another random 500", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds)

sim_shuffle_500 = data.frame(simulations = c(sim_values_500), shuffles = c(shuf_values_500), perm = rep(shuffles, each = 5),
                             threshold_lvl = rep(c(10,50,90,95,99), length(shuffles)))
sim_dif_500 = melt(sim_shuffle_500,
                   id.vars = c("perm","threshold_lvl"),value.name = "threshold", variable.name = "type")




#for a size of 500 - exactly the same subjects
TPQ_scores_1 = TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`)&complete.cases(TPQ_scores$Q1)),][1:500,c("Initial ID","Diagnosis","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
sim_values_500_2 = sapply(shuffles, function(shuf)
  find_threshold(TPQ_scores_1,
                 TPQ_scores_1,
                 Comparison = "Reconstructed Data",
                 G1 = "HC at test",
                 G2 = "HC at test",
                 method = "spearman",
                 n_shuffle = shuf,
                 threshold_choice = 0.95,
                 simulate = TRUE)$standard_thresholds)
shuf_values_500_2 = sapply(shuffles, function(shuf)  find_threshold(TPQ_scores_1, TPQ_scores_1,Comparison = "Reconstructed Data",G1="random 500",
                                                                   G2="the same random 500", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds)

sim_shuffle_500_2 = data.frame(simulations = c(sim_values_500_2), shuffles = c(shuf_values_500_2), perm = rep(shuffles, each = 5),
                              threshold_lvl = rep(c(10,50,90,95,99), length(shuffles)))
sim_dif_500_2 = melt(sim_shuffle_500_2, id.vars = c("perm","threshold_lvl"),value.name = "threshold", variable.name = "type")


#for a size of 900 - different subjects
TPQ_scores_1 = TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`)&complete.cases(TPQ_scores$Q1)),][1:900,c("Initial ID","Diagnosis","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
TPQ_scores_2= TPQ_scores[(complete.cases(TPQ_scores$`Initial ID`)&complete.cases(TPQ_scores$Q1)),][901:1800,c("Initial ID","Diagnosis","Session","NS1","NS2","NS3","NS4","NS","HA1","HA2","HA3","HA4","HA","RD1","RD2","RD3","RD4","RD")]
sim_values_900 = sapply(shuffles, function(shuf)
  find_threshold(TPQ_scores_1,
                 TPQ_scores_2,
                 Comparison = "Reconstructed Data",
                 G1 = "HC at test",
                 G2 = "HC at retest",
                 method = "spearman",
                 n_shuffle = shuf,
                 threshold_choice = 0.95,
                 simulate = TRUE)$standard_thresholds)
shuf_values_900 = sapply(shuffles, function(shuf)  find_threshold(TPQ_scores_1, TPQ_scores_2,Comparison = "Reconstructed Data",G1="Random 900",
                                                                   G2="Another Random 900", method = "spearman", n_shuffle = shuf, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds)

sim_shuffle_900 = data.frame(simulations = c(sim_values_900), shuffles = c(shuf_values_900), perm = rep(shuffles, each = 5),
                              threshold_lvl = rep(c(10,50,90,95,99), length(shuffles)))
sim_dif_900 = melt(sim_shuffle_900,
                    id.vars = c("perm","threshold_lvl"),value.name = "threshold", variable.name = "type")




sim_dif_fin = rbind(sim_dif_100,sim_dif_200,sim_dif_500,sim_dif_500_2,sim_dif_900)
sim_dif_fin$size = rep(c("100","200", "500_dif", "500_same", "900"), each = nrow(sim_dif_100))
ggplot(sim_dif_fin, aes(x = perm, y= threshold, color = type))+geom_line()+
  geom_point()+scale_color_manual(values = c("#Eb4C42","#0087BD","#6cBE58","#00A550","#808CA3"))+
  TypicalTheme+ggtitle("Threshold values at different Number of Simulations/Permutations")+
  facet_grid(size~threshold_lvl)

# For threshold 0.95 only
ggplot(sim_dif_fin[sim_dif_fin$threshold_lvl == "95",], aes(x = perm, y= threshold, color =size ))+geom_line(size = 1.5)+
  geom_point(size = 2.5,fill = "white",shape = 21)+scale_color_manual(values = c("#Eb4C42","#0087BD","#6cBE58","#00A550","#808CA3"))+
  TypicalTheme+ggtitle("Threshold values at different Number of Subjects - Threshold = 0.95")+
  facet_grid(.~type)





find_threshold(TPQ_scores_t, TPQ_scores_r,Comparison = "Reconstructed Data",G1="HC at test",
               G2="HC at retest", method = "pearson", n_shuffle = 1000, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds

find_threshold(TPQ_scores_1, TPQ_scores_2,Comparison = "Reconstructed Data",G1="HC at test",
               G2="HC at retest", method = "pearson", n_shuffle = 1000, threshold_choice = 0.95, simulate = FALSE)$standard_thresholds

find_threshold(TPQ_scores_1, TPQ_scores_2,Comparison = "Reconstructed Data",G1="HC at test",
               G2="HC at retest", method = "pearson", n_shuffle = 1000, threshold_choice = 0.95, simulate = TRUE)$standard_thresholds



# create two random data sets and correlate them using the shuffle method
Data_1 = data.frame(matrix(sample(1000,1000,replace = TRUE),ncol = 10))
Data_2 = data.frame(matrix(sample(1000,1000,replace = TRUE),ncol = 10))
colnames(Data_1) = paste0("C1_",1:10)
colnames(Data_2) = paste0("C2_",1:10)
rownames(Data_1) = paste0("R1_",1:100)
rownames(Data_2) = paste0("R2_",1:100)

Data_1_t = data.frame(t(Data_1))
Data_2_t = data.frame(t(Data_2))


Data_1 = TPQ_scores_t[(TPQ_scores_t$Diagnosis == "MDD"&TPQ_scores_t$Response == "Responder"&complete.cases(TPQ_scores_t$Response)),]
Data_2 = TPQ_scores_r[(TPQ_scores_r$Diagnosis == "MDD"&TPQ_scores_r$Response == "Responder"&complete.cases(TPQ_scores_r$Response)),]

find_threshold(Data_1_t, Data_2_t,Comparison = "Reconstructed Data",G1="Random 200",
               G2="Random 200", method = "spearman", n_shuffle = 100, threshold_choice = 0.95, simulate = FALSE)
find_threshold(Data_1, Data_2,Comparison = "Reconstructed Data",G1="Random 200",
               G2="Random 200", method = "spearman", n_shuffle = 100, threshold_choice = 0.95, simulate = FALSE)

Cor_mat = matrix(0,nrow = 10,ncol = 10)
dimnames(Cor_mat) = list(paste0("C1_",1:10), paste0("C2_",1:10))
names(dimnames(Cor_mat)) = c("Data1","Data2")

p_mat = matrix(0,nrow = 10,ncol = 10)
dimnames(p_mat) = list(paste0("C1_",1:10), paste0("C2_",1:10))
names(dimnames(p_mat)) = c("Data1","Data2")


for (i in 1:10){
  threshold = 0.05
  if (stringent_alpha) {
    threshold = 0.05/100
  }
  
  estimates = sapply(1:10, function(x) cor.test(as.numeric(Data_1[,i]),as.numeric(Data_2[,x]), method = "spearman")$estimate)
  p_values = sapply(1:10, function(x) cor.test(as.numeric(Data_1[,i]),as.numeric(Data_2[,x]), method = "spearman")$p.value)
  #estimates[p_values>threshold | abs(estimates)<estimate_threshold] = 0
  p_mat[i,] = p_values
  Cor_mat[i,] = estimates
}

corrplot::corrplot(Cor_mat, method = "number", p.mat=p_mat, sig.level = threshold,
                   insig = "blank", title = "\n\nCorrelations Between Data_1 and Data_2")
find_threshold(Data_1, Data_2,Comparison = "Reconstructed Data",G1="Random 200",
               G2="Random 200", method = "spearman", n_shuffle = 1000, threshold_choice = 0.15, simulate = FALSE)

# Co-occurance of questions between ICASSO at Test and ICASSO at retest

sources_hc_test = ICASSO_hc_test$sources[ICASSO_hc_test$scores>0.85,1:100]
n_ICs_hc_test = sum(ICASSO_hc_test$scores>0.85)
max_source_hc_test = matrix(rep(apply(abs(sources_hc_test),2,max),n_ICs_hc_test),nrow = n_ICs_hc_test, byrow = TRUE)
thresholds = sapply(1:nrow(sources_hc_test), function(x) quantile(abs(sources_hc_test[x,]),0.9))
IC_questions_hc_test = data.frame(t(sapply(1:nrow(sources_hc_test), function(x) ifelse(sources_hc_test[x,]>thresholds[x],1,0))))
colnames(IC_questions_hc_test) = paste0("Q",1:100)

sources_hc_retest = ICASSO_hc_retest$sources[ICASSO_hc_retest$scores>0.85,1:100]
n_ICs_hc_retest = sum(ICASSO_hc_retest$scores>0.85)
max_source_hc_retest = matrix(rep(apply(abs(sources_hc_retest),2,max),n_ICs_hc_retest),nrow = n_ICs_hc_retest, byrow = TRUE)
thresholds = sapply(1:nrow(sources_hc_retest), function(x) quantile(abs(sources_hc_retest[x,]),0.9))
IC_questions_hc_retest = data.frame(t(sapply(1:nrow(sources_hc_retest), function(x) ifelse(sources_hc_retest[x,]>thresholds[x],1,0))))
colnames(IC_questions_hc_retest) = paste0("Q",1:100)

sources_mdd_test = ICASSO_mdd_test$sources[ICASSO_mdd_test$scores>0.85,1:100]
n_ICs_mdd_test = sum(ICASSO_mdd_test$scores>0.85)
max_source_mdd_test = matrix(rep(apply(abs(sources_mdd_test),2,max),n_ICs_mdd_test),nrow = n_ICs_mdd_test, byrow = TRUE)
thresholds = sapply(1:nrow(sources_mdd_test), function(x) quantile(abs(sources_mdd_test[x,]),0.9))
IC_questions_mdd_test = data.frame(t(sapply(1:nrow(sources_mdd_test), function(x) ifelse(sources_mdd_test[x,]>thresholds[x],1,0))))
colnames(IC_questions_mdd_test) = paste0("Q",1:100)

sources_mdd_retest = ICASSO_mdd_retest$sources[ICASSO_mdd_retest$scores>0.85,1:100]
n_ICs_mdd_retest = sum(ICASSO_mdd_retest$scores>0.85)
max_source_mdd_retest = matrix(rep(apply(abs(sources_mdd_retest),2,max),n_ICs_mdd_retest),nrow = n_ICs_mdd_retest, byrow = TRUE)
thresholds = sapply(1:nrow(sources_mdd_retest), function(x) quantile(abs(sources_mdd_retest[x,]),0.9))
IC_questions_mdd_retest = data.frame(t(sapply(1:nrow(sources_mdd_retest), function(x) ifelse(sources_mdd_retest[x,]>thresholds[x],1,0))))
colnames(IC_questions_mdd_retest) = paste0("Q",1:100)

sapply(1:n_ICs_hc_test, function(IC) colnames(IC_questions_hc_test)[IC_questions_hc_test[IC,] == 1])
sapply(1:n_ICs_hc_retest, function(IC) colnames(IC_questions_hc_retest)[IC_questions_hc_retest[IC,] == 1])
sapply(1:n_ICs_mdd_test, function(IC) colnames(IC_questions_mdd_test)[IC_questions_mdd_test[IC,] == 1])
sapply(1:n_ICs_mdd_retest, function(IC) colnames(IC_questions_mdd_retest)[IC_questions_mdd_retest[IC,] == 1])



IC_questions_cooc_hcTPQ = matrix(nrow = nrow(IC_questions_hc_test), ncol =nrow(tpq_q_cloninger))
dimnames(IC_questions_cooc_hcTPQ) = list(paste0("IC",1:nrow(IC_questions_hc_test)), names(TPQQuestions))

for (i in 1:nrow(IC_questions_hc_test)){
  cooc = sapply(1:nrow(tpq_q_cloninger), function(x) mean(which(IC_questions_hc_test[i,]==1) %in% which(tpq_q_cloninger[x,]==1)))
  IC_questions_cooc_hcTPQ[i,] = cooc
}
corrplot(IC_questions_cooc_hcTPQ, title = "\n\n\nPercentage of Questions of ICs in HC at Test and at Retest",
         cl.lim = c(0,1), col = col <- colorRampPalette(c( "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                                           "#4393C3", "#2166AC", "#053061"))(200))




IC_questions_cooc_hc = matrix(nrow = nrow(IC_questions_hc_test), ncol =nrow(IC_questions_hc_retest))
dimnames(IC_questions_cooc_hc) = list(paste0("IC",1:nrow(IC_questions_hc_test)), paste0("IC",1:nrow(IC_questions_hc_retest)))

for (i in 1:nrow(IC_questions_hc_test)){
  cooc = sapply(1:nrow(IC_questions_hc_retest), function(x) mean(which(IC_questions_hc_test[i,]==1) %in% which(IC_questions_hc_retest[x,]==1)))
  IC_questions_cooc_hc[i,] = cooc
}
corrplot(IC_questions_cooc_hc, title = "\n\n\nPercentage of Questions of ICs in HC at Test and at Retest",
         cl.lim = c(0,1), col = col <- colorRampPalette(c( "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                                           "#4393C3", "#2166AC", "#053061"))(200))





IC_questions_cooc_mdd = matrix(nrow = nrow(IC_questions_mdd_test), ncol =nrow(IC_questions_mdd_retest))
dimnames(IC_questions_cooc_mdd) = list(paste0("IC",1:nrow(IC_questions_mdd_test)), paste0("IC",1:nrow(IC_questions_mdd_retest)))

for (i in 1:nrow(IC_questions_mdd_test)){
  cooc = sapply(1:nrow(IC_questions_mdd_retest), function(x) mean(which(IC_questions_mdd_test[i,]==1) %in% which(IC_questions_mdd_retest[x,]==1)))
  IC_questions_cooc_mdd[i,] = cooc
}


corrplot(IC_questions_cooc_mdd, title = "\n\n\nPercentage of Questions of ICs in MDD at Test and at Retest",
         cl.lim = c(0,1), col = col <- colorRampPalette(c( "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                                           "#4393C3", "#2166AC", "#053061"))(200))




IC_questions_cooc_hcmdd = matrix(nrow = nrow(IC_questions_hc_test), ncol =nrow(IC_questions_mdd_test))
dimnames(IC_questions_cooc_hcmdd) = list(paste0("IC",1:nrow(IC_questions_hc_test)), paste0("IC",1:nrow(IC_questions_mdd_test)))

for (i in 1:nrow(IC_questions_hc_test)){
  cooc = sapply(1:nrow(IC_questions_mdd_test), function(x) mean(which(IC_questions_hc_test[i,]==1) %in% which(IC_questions_mdd_test[x,]==1)))
  IC_questions_cooc_hcmdd[i,] = cooc
}


corrplot(IC_questions_cooc_hcmdd, title = "Percentage of Questions of ICs in HC at Test and MDD at Test",
         cl.lim = c(0,1), col = col <- colorRampPalette(c( "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                                           "#4393C3", "#2166AC", "#053061"))(200))


View(TPQ_scores)
matched_inds = match(npz_tpq["IDs"],TPQ_scores$`Final ID`)
new_tpq = TPQ_scores[matched_inds,]
new_tpq = new_tpq[complete.cases(new_tpq$`Final ID`),]
new_tpq[paste0("QO",1:100)] = ifelse(new_tpq[paste0("QO",1:100)] == "TRUE","T","F")

write.table(new_tpq,'/home/asawalma/git/tpq_analysis/data/Tpq Scores.xlsx',row.names = FALSE)





# Simple CFA Model
