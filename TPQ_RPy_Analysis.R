##############################################################################
###                            Prepare Data                                ###
##############################################################################
library(reticulate)
library(ggplot2)


# define directories of main files and python functions file
path = '/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine'
py_function_path = '/home/abdulrahman/git/TPQ/TPQ_functions.py'

if (!dir.exists(path)) {
  path = gsub("/abdulrahman/abdulrahman.sawalma@gmail.com",
              "/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive",
              path)
  py_function_path = gsub('home/abdulrahman',
                          'home/asawalma',
                          py_function_path)
}


# import python functions
source_python(py_function_path)

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


prepare_icasso_list <- function(df_list){
  # prepare data by changing the names of the list items and changing the type of numeric variables to numeric
  # 0. Change the names of list variables so that they can be easily called
  names(df_list) = c("tpq", "reco", "scores", "projections", "sources")
  
  # 1. write a function to prepare single data sets
  prep_single_df <- function(df, python_object = FALSE){
    ##############################################
    #   A Function to clear NAs and change the   #
    # types of variables to the appropriate ones #
    ##############################################
    # 0. Make sure you have an R object
    if (python_object) {
      df = py_to_r(df)
    }
    
    # 1. remove NAs
    df = df[df$Diagnosis != "NA", ]
    for (i in 1:8) {
      df[, i][df[, i] == "NA"] = NA
    }
    
    # 2. change variable types
    df[, 2:8] = as.data.frame(apply(df[, 2:8], 2, as.factor))
    df[, 1] = factor(df[, 1], levels = df[, 1])
    df[, 9:ncol(df)] = as.data.frame(apply(df[, 9:ncol(df)],2,as.numeric))
    
    return(df)
  }
  
  
  # 2. Prepare single data frames
  # 2.A. Prepare original data (tpq)
  df_list$tpq = prep_single_df(df_list$tpq)
  
  
  # 2.B. Prepare reconstructed data
  for (i in 1:length(df_list$reco)) {
    df_list$reco[[i]] = prep_single_df(df_list$reco[[i]], python_object = TRUE)
  }
  
  
  # 2.C. Prepare Projections
  df_list$projections = prep_single_df(df_list$projections)
  
  # 2.D. prepare sources
  df_list$sources$IC = paste0("IC", 1:nrow(df_list$sources))
  
  return(df_list)
}


ICASSO_all = prepare_icasso_list(ICASSO_all)
ICASSO_all_fast = prepare_icasso_list(ICASSO_all_fast)
ICASSO_all_fastsk = prepare_icasso_list(ICASSO_all_fastsk)
ICASSO_hc = prepare_icasso_list(ICASSO_hc)
ICASSO_mdd = prepare_icasso_list(ICASSO_mdd)
ICASSO_gad = prepare_icasso_list(ICASSO_gad)
ICASSO_ptsd = prepare_icasso_list(ICASSO_ptsd)
ICASSO_tnp = prepare_icasso_list(ICASSO_tnp)
ICASSO_disorders = prepare_icasso_list(ICASSO_disorders)


##############################################################################
###                   Find Correlations Between ICs of                     ###
###                     Different Data Decompositions                      ###
##############################################################################


ref_sources = ICASSO_all$sources
fast_sources = ICASSO_all_fast$sources
fastsk_sources = ICASSO_all_fastsk$sources
hc_sources = ICASSO_hc$sources
gad_sources = ICASSO_gad$sources
mdd_sources = ICASSO_mdd$sources
ptsd_sources = ICASSO_ptsd$sources
tnp_sources = ICASSO_tnp$sources
source_cor <- function(df_1, df_2, Factor, max_components_1, max_components_2, estimate_threshold, stringent_alpha = TRUE){
  ############################################################
  # A Function to find the correlation between two sources
  # Returns a matrix of spearman's rho values
  ############################################################
  if (is.na(max_components_1)){
    max_components_1 = nrow(df_1)
    max_components_2 = nrow(df_2)
  }
  df_1 = df_1[1:max_components_1,1:100]
  df_2 = df_2[1:max_components_2,1:100]
  
  Sources_cor = matrix(0,nrow = max_components_1,ncol = max_components_2)
  dimnames(Sources_cor) = list(paste0("IC",1:max_components_1), paste0("IC",1:max_components_2))
  names(dimnames(Sources_cor)) = c("All",Factor)
  
  for (i in 1:max_components_1){
    threshold = 0.05
    if (stringent_alpha) {
      threshold = threshold / (max_components_2 * max_components_1)
    }
    
    estimates = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]))$estimate)
    p_values = sapply(1:max_components_2, function(x) cor.test(as.numeric(df_1[i,]),as.numeric(df_2[x,]))$p.value)
    estimates[p_values>threshold | abs(estimates)<estimate_threshold] = 0
    Sources_cor[i,] = estimates
  }
  return(Sources_cor)
}

simplify_mat <- function(sources_Cor){
  ##########################################################################
  # A function to simplify the correlation matrix into
  #   a 3-column matrix:
  #     1st col -> ICs of df1, 2 -> ICs of df2 and 3 -> spearman rho values
  ##########################################################################
  final_matrix = matrix(rep(0, 3 * nrow(sources_Cor)), nrow = 3)
  for (i in 1:nrow(sources_Cor)) {
    IC = paste0("IC", i)
    IC_all = IC
    
    #find the IC in the second df that has a correlation of >0
    # If you have two, choose the one with highest cor
    max_cor_IC = abs(sources_Cor[IC, ]) > 0 &
      abs(sources_Cor[IC, ]) == max(abs(sources_Cor[IC, ]))
    
    IC_df2 = colnames(sources_Cor)[max_cor_IC]
    Cor_value = round(sources_Cor[IC, max_cor_IC], digits = 3)
    if (sum(max_cor_IC) == 0) {
      final_matrix[, i] = c(NA, NA, NA)
    } else{
      final_matrix[, i] = c(IC_all, IC_df2, Cor_value)
    }
    dimnames(final_matrix) = list(c("All","DF2","Cor"), paste0("IC",1:10))
  }
  
  return(final_matrix)
}


est_threshold = 0.6
AllFAST_sources_Cor = source_cor(ref_sources, fast_sources, "All_Fast", 10,10,est_threshold)
AllFASTSK_sources_Cor = source_cor(ref_sources, fastsk_sources, "All_Fast_SK", 10,10,est_threshold)
HC_sources_Cor = source_cor(ref_sources, hc_sources, "HC", 10,10,est_threshold)
GAD_sources_Cor = source_cor(ref_sources, gad_sources, "GAD", 10,10,est_threshold)
MDD_sources_Cor = source_cor(ref_sources, mdd_sources, "MDD", 10,10,est_threshold)
PTSD_sources_Cor = source_cor(ref_sources, ptsd_sources, "PTSD", 10,10,est_threshold)
TNP_sources_Cor = source_cor(ref_sources, tnp_sources, "TNP", 10,10,est_threshold)


# now find the ICs in the each sources df that 
AllFAST_sources_ICs = simplify_mat(AllFAST_sources_Cor)
AllFASTSK_sources_ICs = simplify_mat(AllFASTSK_sources_Cor)
HC_sources_ICs = simplify_mat(HC_sources_Cor)
GAD_sources_ICs = simplify_mat(GAD_sources_Cor)
MDD_sources_ICs = simplify_mat(MDD_sources_Cor)
PTSD_sources_ICs = simplify_mat(PTSD_sources_Cor)
TNP_sources_ICs = simplify_mat(TNP_sources_Cor)


all_sources = data.frame(
  Question = melt(hc_sources[1:10, ])$variable,
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

all_sources$Question = rep(paste0("Q",1:100),each = 10)
all_sources_melted = melt(all_sources, id.vars = c("Question","IC_all"), value.name = "Source", variable.name = "Group")

all_sources_melted$IC_group = NA
shuffle_ICs <- function(Sources_IC, Factor){
  if (length(Sources_IC)>0){
    for (i in 1:ncol(Sources_IC)){
      all_sources_melted$IC_group[all_sources_melted$Group==Factor&all_sources_melted$IC_all== Sources_IC[1,i]] = Sources_IC[2,i]
      all_sources_melted$Sources_group[all_sources_melted$Group==Factor&all_sources_melted$IC_all== Sources_IC[1,i]] = all_sources_melted$Source[all_sources_melted$Group==Factor&all_sources_melted$IC_all== Sources_IC[2,i]]
    }
  }
  return(all_sources_melted)
}

all_sources_melted = shuffle_ICs(AllFAST_sources_ICs[,!is.na(AllFAST_sources_ICs[2,])],"All_fast")
all_sources_melted = shuffle_ICs(AllFASTSK_sources_ICs[,!is.na(AllFASTSK_sources_ICs[2,])],"All_fastsk")
all_sources_melted = shuffle_ICs(HC_sources_ICs[,!is.na(HC_sources_ICs[2,])],"HC")
all_sources_melted = shuffle_ICs(GAD_sources_ICs[,!is.na(GAD_sources_ICs[2,])],"GAD")
all_sources_melted = shuffle_ICs(PTSD_sources_ICs[,!is.na(PTSD_sources_ICs[2,])],"PTSD")
all_sources_melted = shuffle_ICs(TNP_sources_ICs[,!is.na(TNP_sources_ICs[2,])],"TNP")
all_sources_melted = shuffle_ICs(MDD_sources_ICs[,!is.na(MDD_sources_ICs[2,])],"MDD")
all_sources_melted$Sources_group[all_sources_melted$Group == "All"] = all_sources_melted$Source[all_sources_melted$Group == "All"]
all_sources_melted$IC_group_plot = as.vector(sapply(1:8, function(x)
  c(all_sources_melted$IC_group[all_sources_melted$Question == "Q1"][((x - 1) * 10 + 1):(x * 10)], rep(NA, 990))))
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

g1 = ggplot(data = all_sources_melted, aes(x=Question, y= Sources_group, group =Group))+
  geom_line()+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x= 50,y=6,label = IC_group_plot))+
  geom_text(aes(x= 50,y=10,label = Cor_plot))

g1


# delete later
#read valence data
library(readxl)
valenced = data.frame(read_excel(paste0(path,"/TPQ_DataAndAnalysis/FactorCoding_19.07.2021.xlsx")))
valenced[,1:7] = as.data.frame(apply(valenced[,1:7],2,as.factor))
valenced[,1] = factor(valenced[,1], levels = valenced[,1])

#add valence to all_sources_melted
all_sources_melted_new = all_sources_melted
for (CN in colnames(valenced)[3:7]){
  all_sources_melted[[CN]] = rep(valenced[[CN]],each = 10)
}
ggplot(data = all_sources_melted, aes(x=Question, y= Sources_group, group =Group, color = Active_Passive))+
  geom_bar(stat = "identity")+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x= 50,y=6,label = IC_group_plot))+
  geom_text(aes(x= 50,y=10,label = Cor_plot))
# end delete