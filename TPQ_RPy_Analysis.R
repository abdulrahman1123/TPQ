
# Data Preparation --------------------------------------------------------

library(reshape2)
library(ggplot2)

library(reticulate)
library(readxl)
if (dir.exists("/home/abdulrahman/anaconda3/envs/mne/bin/")){
  use_python ("/home/abdulrahman/anaconda3/envs/mne/bin/python3", required = TRUE)
}else if(dir.exists("/home/asawalma/anaconda3/envs/mne/bin/python3")){
  use_python("/home/asawalma/anaconda3/envs/mne/bin/python3", required = TRUE)
} else {
  use_python("C:/Users/jsawa/Anaconda3/envs/mne/python.exe", required = TRUE)
}

# load libraries
np = import("numpy")
pd = import("pandas")
os = import("os")


# define directories of main files and python functions file
path =        '/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine'
scores_path = "/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx"

if (dir.exists("/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive")) {
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

or_scores = as.data.frame(read_xlsx(scores_path, na = "NA"))[c("Diagnosis", "Final ID", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session")]

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
  inc_ids = npz_tpq["IDs"] %in% or_scores$`Final ID` # IDs in numpy array that are found in the original scores file
  matched_inds = match(npz_tpq["IDs"],or_scores$`Final ID`)
  matched_inds = matched_inds[inc_ids]
  sub_info = data.frame(
    ID = factor(npz_tpq["IDs"][inc_ids],levels = npz_tpq["IDs"]),
    Diagnosis = factor(or_scores$Diagnosis[matched_inds]),
    Trauma = factor(or_scores$Trauma[matched_inds]),
    GAD = factor(or_scores$GAD[matched_inds]),
    Response = factor(or_scores$Response[matched_inds]),
    PTSD = factor(or_scores$PTSD[matched_inds]),
    MDD = factor(or_scores$MDD[matched_inds]),
    Session = factor(or_scores$Session[matched_inds])
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
est_threshold = 0.7
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

# Some sort values are similar in two subjects (if they have exactly the same value). Fix this
for (k in 1:2){
  for (Factor in levels(projections_melted$Diagnosis)){
    for (IC in paste0("IC",1:10)){
      for (Group in levels(projections_melted$Grouping)){
        Sorting_Cond = projections_melted$Diagnosis== Factor & projections_melted$Grouping == Group & projections_melted$IC == IC
        for (i in 1:length(projections_melted$Sort[Sorting_Cond])){
          if (sum (projections_melted$Sort[Sorting_Cond] == projections_melted$Sort[Sorting_Cond][i]) >1){
            projections_melted$Sort[Sorting_Cond][i] = projections_melted$Sort[Sorting_Cond][i]+1
            print(projections_melted$Sort[Sorting_Cond][i])
          }
        }
      }
    }
  }
}




# add a variable to include only the close-to-zero values
projections_melted$Projection_zero = projections_melted$Projection

for (Diagnosis in levels(projections_melted$Diagnosis)){
  for (Group in levels(factor(projections_melted$Grouping))){
    for (IC in levels(projections_melted$IC)){
      subgroup_cond = projections_melted$Grouping==Group & projections_melted$IC==IC & projections_melted$Diagnosis==Diagnosis
      pro_val = projections_melted$Projection[subgroup_cond]
      projections_melted$Projection_scaled[subgroup_cond] = pro_val/mean(abs(pro_val), na.rm = TRUE)
      min_value = min(abs(projections_melted$Projection_scaled[subgroup_cond]))
      projections_melted$Projection_zero[subgroup_cond][abs(projections_melted$Projection_scaled)[subgroup_cond] != min_value] = NA
    }
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



for (Factor in levels(projections_melted$Diagnosis)) {
  for (IC in paste0("IC", 1:10)) {
    for (Group in c("All_Fast", "All_Fast_Sickit", "Separate")) {
      cond_all_ic = (projections_melted$IC == IC & projections_melted$Diagnosis == Factor & projections_melted$Grouping == Group)
      new_ic = projections_melted$IC_group[cond_all_ic][1]
      cond_new_ic = (projections_melted$IC == new_ic & projections_melted$Diagnosis == Factor & projections_melted$Grouping == Group)

      projections_melted$Projections_group[cond_all_ic] = projections_melted$Projection[cond_new_ic]
      projections_melted$Projections_scaled_group[cond_all_ic] = projections_melted$Projection_scaled[cond_new_ic]
      projections_melted$Sort_group[cond_all_ic] = projections_melted$Sort[cond_new_ic]
      projections_melted$Projection_zero_group[cond_all_ic] = projections_melted$Projection_zero[cond_new_ic]
    }
  }
}

projections_melted$Projections_scaled_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection_scaled[projections_melted$Grouping == "Combined"]
projections_melted$Projections_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection[projections_melted$Grouping == "Combined"]
projections_melted$Sort_group[projections_melted$Grouping == "Combined"] = projections_melted$Sort[projections_melted$Grouping == "Combined"]
projections_melted$Projections_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA
projections_melted$Projections_scaled_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA
projections_melted$Sort_group[projections_melted$Grouping == "Combined" & is.na(projections_melted$Projections_group[projections_melted$Grouping == "Separate"])] = NA
projections_melted$Projection_zero_group[projections_melted$Grouping == "Combined"] = projections_melted$Projection_zero[projections_melted$Grouping == "Combined"]
# remove all IC designations, except for the first one, for plotting purposes
# find the indices to keep
Fast_inds = 18190 + seq(from = 1, to=1819*10, by = 1819)
FastSK_inds = 18190*2 + seq(from = 1, to=1819*10, by = 1819)
HC_indices = 18190*3 + seq(from = 1, to=1819*10, by = 1819)
MDD_indices = HC_indices + nrow(projections_melted[projections_melted$Diagnosis == "HC",])/20
PTSD_indices = MDD_indices + nrow(projections_melted[projections_melted$Diagnosis == "MDD",])/20
TNP_indices = PTSD_indices + nrow(projections_melted[projections_melted$Diagnosis == "PTSD",])/20
GAD_indices = TNP_indices + nrow(projections_melted[projections_melted$Diagnosis == "TNP",])/20
Inc_inds = c(Fast_inds, FastSK_inds, HC_indices, MDD_indices, PTSD_indices, TNP_indices, GAD_indices)
projections_melted$IC_group[-Inc_inds] = NA
projections_melted$Cor_arranged[-Inc_inds] = NA


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
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "")+
  MinimalTheme

ggplot(data = projections_melted_dia, mapping = aes(x = Sort_group,y= Grouping))+
  geom_raster(aes(fill = Projections_scaled_group))+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradient2(high = "#a10000",mid = "white", low = "#002138", midpoint = 0, name = "Scaled\nProjection")+
  geom_point(data = projections_melted_dia[!is.na(projections_melted_dia$Projection_zero_group),], aes(color = Projection_zero_group))+
  scale_color_gradient2(high = "black", mid = "black", low = "black", name = "Zero", labels = c())+
  facet_grid(IC~Diagnosis, scale = "free")+
  geom_text(aes(x= 3,y=2.25,label = IC_group))+
  geom_text(aes(x= 3,y=1.75,label = Cor_arranged))+
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "Scaled by dividing on the mean of absolute values")+
  MinimalTheme


# Todo: make an all sorting variable
ggplot(data = projections_melted_all, mapping = aes(x = Sort_group,y= Grouping))+
  geom_raster(aes(fill = Projections_scaled_group))+
  scale_x_continuous("Sorted Subjects")+
  scale_y_discrete("Decomposition Type")+
  scale_fill_gradient2(high = "#a10000",mid = "white", low = "#002138", midpoint = 0, name = "Scaled\nProjection")+
  geom_point(data = projections_melted_all[!is.na(projections_melted_all$Projection_zero_group),], aes(color = Projection_zero_group))+
  scale_color_gradient2(high = "black", mid = "black", low = "black", name = "Zero", labels = c())+
  facet_grid(IC~Grouping, scale = "free")+
  geom_text(aes(x= 3,y=2.25,label = IC_group))+
  geom_text(aes(x= 3,y=1.75,label = Cor_arranged))+
  ggtitle("Projection Values Per Diagnosis Per IC", subtitle = "Scaled by dividing on the mean of absolute values")+
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
  geom_text(aes(x= -0.1,y=-0.25,label = Cor_arranged))+
  facet_grid(IC~Diagnosis, scale = "free")
  

# Logistic Regression for Diagnosis ----------------------------------
# read projections data
library(readxl)
library(boot)
Projections = ICASSO_all$projections
Projections = Projections[(Projections$Diagnosis == "MDD" | Projections$Diagnosis == "HC"),]
Projections$Diagnosis = factor(Projections$Diagnosis, levels = c("HC","MDD"))


# read TPQ data
ex_dir = "G:/My Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx"
TPQ = as.data.frame(read_excel(ex_dir))

# convert all Questions to numeric
TPQ[,paste0("Q",1:100)] = apply(TPQ[,paste0("Q",1:100)],2, as.numeric)


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
  
  # calculate Cloninger's subscales as suggested by Cloninger (some answers are flipped)
  TPQ[[item]]=apply(TPQ[,paste0("Q",TPQQuestions[[item]])], 1,sum, na.rm = TRUE)
  
  # calculate Cloninger's subscales by just summing True values (because the values given to Jurgen were not flipped,
  # so, I will calculate a Data-driven equivalent for these subscales)
  TPQ[[paste0("O_",item)]]=apply(TPQ[,paste0("QO",TPQQuestions[[item]])], 1,sum, na.rm = TRUE)
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


# create a projections glm model
Form = as.formula(paste0("Diagnosis ~ ",paste0("IC",1:10,collapse =  "+")))
ModelProjections = glm(Form, data = Projections, family = binomial())
#ModelProjections = step(ModelProjections,direction = "back")
summary(ModelProjections)
LogisticFunction(ModelProjections, plt_type = "histogram", Threshold = 0.5)
ProjectionsPE = cv.glm(data = Projections ,ModelProjections ,K=10)$delta[2]
print(paste("IC Model accuracy based on 10 fold cv = ", 100*(1-round(ProjectionsPE,3)),"%"))


# Create a TPQ glm model
Form = as.formula(paste0("Diagnosis ~ ",paste0(subscales,collapse =  "+")))
ModelTPQ = glm(Form, data = TPQ, family = binomial())
#ModelTPQ = step(ModelTPQ,direction = "back")
summary(ModelTPQ)
LogisticFunction(ModelTPQ, plt_type = "histogram", Threshold = 0.5)
PETPE = cv.glm(data = TPQ ,ModelTPQ ,K=10)$delta[2]
print(paste("TPQ Model accuracy based on 10 fold cv = ", 100*(1-round(PETPE,3)),"%"))

# Create a TPQ questions glm model
TPQ_NA = na.exclude(TPQ[,c("Diagnosis",included_questions)])
Form = as.formula(paste0("Diagnosis ~ ",paste0(included_questions,collapse =  "+")))
ModelQuestions = glm(Form, data = TPQ_NA, family = binomial())
ModelQuestions = step(ModelQuestions,direction = "back")
#BackQuestions = c("Q2", "Q3", "Q5", "Q7", "Q8", "Q9", "Q13", "Q18", "Q19", "Q22", "Q26", "Q38", "Q39", "Q40", "Q41", "Q43", "Q45", "Q46", "Q49", "Q51", "Q54", "Q55", "Q56", "Q57", "Q60", "Q62", "Q64", "Q66", "Q69", "Q75", "Q77", "Q83", "Q85", "Q86", "Q88", "Q92", "Q93", "Q95")
#Form = as.formula(paste0("Diagnosis ~ ",paste0(BackQuestions,collapse =  "+")))
#ModelQuestions = glm(Form, data = TPQ_NA, family = binomial())
summary(ModelQuestions)
LogisticFunction(ModelQuestions, plt_type = "histogram", Threshold = 0.5)

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


show_svm = function(data_frame, DV, IVs, res = 75, m_type = "linear", perc = 0.8){
  data_frame = data_frame[colnames(data_frame) %in% DV| colnames(data_frame) %in% IVs]
  data_frame = na.omit(data_frame)
  make_grid = function(x1,x2,n = 75){  
    ## This function only creates a range of dots
    # These dots will be colored according to the predicted value based on our data
    x1 = seq(from = min(x1)-0.5, to = max(x1)+0.5, length = n)
    x2 = seq(from = min(x2)-0.5, to = max(x2)+0.5, length = n)
    
    new_df = expand.grid(X1 = x1, X2 = x2)
    colnames(new_df) = colnames(x)[1:2]
    
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
  
  # create a test and train data sets
  train_inds = sample(1:nrow(data_frame), size = perc*nrow(data_frame))
  train = data_frame[train_inds,]
  test = data_frame[-train_inds,]
  
  
  Formula = as.formula(paste0(DV," ~ ",paste0(IVs,collapse = "+")))
  svm_model = svm(Formula, data = train, kernel = m_type, cost = 1, scale = TRUE)

  if (length(IVs) == 2){
    grid = make_grid(train[[IVs[1]]],train[[IVs[2]]], n = res)
    preds = predict(svm_model, grid)
    predicted_df = data.frame(X1 = grid[,1], X2 = grid[,2], Y=preds)
    sf_polygons = convert_to_sf(predicted_df)
    
    Colors = c("#C33E3B","#4EA3DF","#6cBE58","#808CA3","#B9B0AB","#2F4F4F", "#CC6666", "#9999CC", "#66CC99","#682860","#FBEC5D","#FF6347","#FF3800","#1B4D3E","#E30B5D")
    # plot the model
    g_svm = ggplot()+
      geom_sf(data = sf_polygons, alpha = 0.25, mapping= aes(x=NULL,y=NULL,group = sf_polygons$Y,fill = sf_polygons$Y))+
      scale_fill_gradientn(breaks = 1:length(IVs),colors = Colors[1:length(IVs)])+
      geom_point(data = train,mapping = aes(x = train[[IVs[1]]], y =train[[IVs[2]]],color = train[[DV]]),size = 3)+
      scale_color_manual(values = Colors)+
      ggtitle("Data points of the variables X0 and X1", subtitle = "Original Data with Predicted Values")+
      MinimalTheme
    
    print(g_svm)
  }
    
  
  preds = predict(svm_model, test)
  Accuracy = round(mean(preds== test[[DV]]),digits = 4)
  res_table = table(predicted = preds, actual = test[[DV]])
  return(list(svm_model, Accuracy, res_table))
}


show_svm(data_frame = pro_fast, m_type = "linear",DV = "Diagnosis",IVs = paste0("IC",1:10))
lin_info = mean(sapply(1:50, function(x) show_svm(data_frame = pro_all, m_type = "linear", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
rad_info = mean(sapply(1:50, function(x) show_svm(data_frame = pro_all, m_type = "radial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
pol_info = mean(sapply(1:50, function(x) show_svm(data_frame = pro_all, m_type = "polynomial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))

lin_fastmne = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast, m_type = "linear", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
rad_fastmne = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast, m_type = "radial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
pol_fastmne = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast, m_type = "polynomial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))

lin_fastsk = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast_sk, m_type = "linear", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
rad_fastsk = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast_sk, m_type = "radial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))
pol_fastsk = mean(sapply(1:50, function(x) show_svm(data_frame = pro_fast_sk, m_type = "polynomial", DV = "Diagnosis",IVs = paste0("IC",1:10))[[2]]))

lin_info
rad_info
pol_info

lin_fastmne
rad_fastmne
pol_fastmne

lin_fastsk
rad_fastsk
pol_fastsk

mean(accuracies)
# 10 ICs = 80.592 - linear
# 20 ICs = 81.467 - linear
# 10 ICs = 80.818 - radial
# 20 ICs = 82.43 - radial
# 10 ICs = 79.203 - polynomial
# 20 ICs = 80.563  - polynomial
# 10 ICs = 75.22  - sigmoid
# 20 ICs = 73.23 - sigmoid

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

or_tpq = as.data.frame(readxl::read_xlsx("G:/My Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx", na = "NA"))[c("Final ID", "Diagnosis", paste0("QO",1:100))]
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


tpq_radial = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),m_type = "radial")[[2]]))
tpq_linear = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),m_type = "linear")[[2]]))
tpq_poly = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)),m_type = "polynomial")[[2]]))

tpq_radial_scale = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("NS","HA","RD"),m_type = "radial")[[2]]))
tpq_linear_scale = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("NS","HA","RD"),m_type = "linear")[[2]]))
tpq_poly_scale = mean(sapply(1:50, function(x) show_svm(data = or_tpq,DV = "Diagnosis",IVs = c("NS","HA","RD"),m_type = "polynomial")[[2]]))

# CLoninger appears to have very low predictive value for MDD

# SVM for scales and subscales


all_source_svm = mean(sapply(1:50, function(x) show_svm(all_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
fast_source_svm = mean(sapply(1:50, function(x) show_svm(fast_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
fastsk_source_svm = mean(sapply(1:50, function(x) show_svm(fastsk_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
hc_source_svm = mean(sapply(1:50, function(x) show_svm(hc_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
gad_source_svm = mean(sapply(1:50, function(x) show_svm(gad_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
mdd_source_svm = mean(sapply(1:50, function(x) show_svm(mdd_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
ptsd_source_svm = mean(sapply(1:50, function(x) show_svm(ptsd_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
tnp_source_svm = mean(sapply(1:50, function(x) show_svm(tnp_sources, DV = "subscale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))

all_source_svm = mean(sapply(1:50, function(x) show_svm(all_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
fast_source_svm = mean(sapply(1:50, function(x) show_svm(fast_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
fastsk_source_svm = mean(sapply(1:50, function(x) show_svm(fastsk_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
hc_source_svm = mean(sapply(1:50, function(x) show_svm(hc_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
gad_source_svm = mean(sapply(1:50, function(x) show_svm(gad_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
mdd_source_svm = mean(sapply(1:50, function(x) show_svm(mdd_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
ptsd_source_svm = mean(sapply(1:50, function(x) show_svm(ptsd_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))
tnp_source_svm = mean(sapply(1:50, function(x) show_svm(tnp_sources, DV = "scale",IVs = paste0("IC",1:10), m_type = "radial")[[2]]))


# KNN -------------------------------------------------------------
library(class)
show_knn = function(df, DV, IVs, perc){
  nor <-function(x) {(x -min(x))/(max(x)-min(x))}
  
  df = na.omit(df[,c(IVs,DV)])
  data_norm <- as.data.frame(lapply(df[,IVs], nor))
  
  train_inds = sample(1:nrow(df), size = perc*nrow(df))
  df_train = data_norm[train_inds,] 
  df_test = data_norm[-train_inds,] 
  
  train_category = df[train_inds,colnames(df) == DV]
  
  test_category <- df[-train_inds,colnames(df) == DV]
  
  predicted = knn(df_train,df_test,cl=train_category,k=13)
  
  tab = table(predicted,test_category)
  Accuracy = mean(predicted == test_category)
  
  return(list(tab, Accuracy))
}

###########################
# projections
###########################

pro_all = ICASSO_all$projections[,c(9:18,2)]
pro_fast = ICASSO_all_fast$projections[,c(9:18,2)]
pro_fastsk = ICASSO_all_fastsk$projections[,c(9:18,2)]

pro_all = pro_all[pro_all$Diagnosis %in% c("HC","MDD"),]
pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"),]
pro_fastsk = pro_fastsk[pro_fastsk$Diagnosis %in% c("HC","MDD"),]

pro_all$Diagnosis = factor(pro_all$Diagnosis)
pro_fast$Diagnosis = factor(pro_fast$Diagnosis)
pro_fastsk$Diagnosis = factor(pro_fastsk$Diagnosis)


knn_all = mean(sapply(1:50, function(x) show_knn(df = pro_all, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
knn_fast = mean(sapply(1:50, function(x) show_knn(df = pro_fast, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
knn_fastsk = mean(sapply(1:50, function(x) show_knn(df = pro_fastsk, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))


###########################
# Cloninger
###########################
tpq = ICASSO_all$tpq

or_tpq = as.data.frame(readxl::read_xlsx("G:/My Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx", na = "NA"))[c("Final ID", "Diagnosis", paste0("QO",1:100))]
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

tpq_subscale = mean(sapply(1:50, function(x) show_knn(df = or_tpq,DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)), perc = 0.8)[[2]]))
tpq_scale = mean(sapply(1:50, function(x) show_knn(df = or_tpq,DV = "Diagnosis",IVs = c("HA","NS","RD"), perc = 0.8)[[2]]))

#Cloninger is a little bit worse
###########################
# sources
###########################
all_knn = mean(sapply(1:50, function(x) show_knn(df = all_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
fast_knn = mean(sapply(1:50, function(x) show_knn(df = fast_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
fastsk_knn = mean(sapply(1:50, function(x) show_knn(df = fastsk_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
hc_knn = mean(sapply(1:50, function(x) show_knn(df = hc_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
gad_knn = mean(sapply(1:50, function(x) show_knn(df = gad_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
mdd_knn = mean(sapply(1:50, function(x) show_knn(df = mdd_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
ptsd_knn = mean(sapply(1:50, function(x) show_knn(df = ptsd_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
tnp_knn = mean(sapply(1:50, function(x) show_knn(df = tnp_sources, DV = "scale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))

all_knn = mean(sapply(1:50, function(x) show_knn(df = all_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
fast_knn = mean(sapply(1:50, function(x) show_knn(df = fast_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
fastsk_knn = mean(sapply(1:50, function(x) show_knn(df = fastsk_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
hc_knn = mean(sapply(1:50, function(x) show_knn(df = hc_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
gad_knn = mean(sapply(1:50, function(x) show_knn(df = gad_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
mdd_knn = mean(sapply(1:50, function(x) show_knn(df = mdd_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
ptsd_knn = mean(sapply(1:50, function(x) show_knn(df = ptsd_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))
tnp_knn = mean(sapply(1:50, function(x) show_knn(df = tnp_sources, DV = "subscale", IVs = paste0("IC",1:10), perc = 0.8)[[2]]))


# decision trees -------------------------------------------------------------

library(dplyr)
library(rpart.plot)
library(rpart)



show_decision_tree = function(data_frame, IVs, DV, perc = 0.8, dc_method = "class", min_split = 4,
                              min_bucket = 2, max_depth = 3, c_p = 0, use_control = FALSE, plot = TRUE){
  #hyperparameters
  control = rpart.control(minsplit = min_split, # min number of observations before the algorith splits
                          minbucket = min_bucket, # min number of observations in the final node
                          maxdepth = max_depth, # max depth of a node (with the root being node 0)
                          cp = c_p)
  
  
  data_frame = data_frame[,c(IVs,DV)]
  data_frame = na.omit(data_frame)
  
  train_inds = sample(1:nrow(data_frame), size = perc*nrow(data_frame))
  train = data_frame[train_inds,]
  test = data_frame[-train_inds,]
  
  Formula = as.formula(paste0(DV," ~ ",paste0(IVs,collapse = "+")))
  if (use_control){
    fit = rpart(Formula, data = train, method = dc_method, control = control) # you can use class for classification, and anova for regression
  } else {
    fit = rpart(Formula, data = train, method = dc_method)
  }

  if (plot){
    rpart.plot(fit, extra = "auto") # extra is for extra information, refer to: https://cran.r-project.org/web/packages/rpart.plot/rpart.plot.pdf
  }
  
  # make prediction
  predict_unseen = predict(fit, test, type = "class")
  pred_table = table(predicted = predict_unseen, actual = test[[DV]])
  Accuracy = mean(predict_unseen == test[[DV]])
  
  return(list(pred_table, Accuracy))
  
}

##############################
# Projections
##############################

pro_all = ICASSO_all$projections[,c(9:18,2)]
pro_fast = ICASSO_all_fast$projections[,c(9:18,2)]
pro_fast_sk = ICASSO_all$projections[,c(9:18,2)]

pro_all = pro_all[pro_all$Diagnosis %in% c("HC","MDD"),]
pro_fast = pro_fast[pro_fast$Diagnosis %in% c("HC","MDD"),]
pro_fast_sk = pro_fast_sk[pro_fast_sk$Diagnosis %in% c("HC","MDD"),]

pro_all$Diagnosis = factor(pro_all$Diagnosis)
pro_fast$Diagnosis = factor(pro_fast$Diagnosis)
pro_fast_sk$Diagnosis = factor(pro_fast_sk$Diagnosis)

dc_all = mean(sapply(1:50, function(x) show_decision_tree(data_frame = pro_all, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
dc_fast = mean(sapply(1:50, function(x) show_decision_tree(data_frame = pro_fast, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))
dc_fast_sk = mean(sapply(1:50, function(x) show_decision_tree(data_frame = pro_all, DV = "Diagnosis", IVs = paste0("IC",1:10), perc = 0.8, plot = FALSE)[[2]]))


##############################
# Cloninger
##############################
tpq = ICASSO_all$tpq

or_tpq = as.data.frame(readxl::read_xlsx("G:/My Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx", na = "NA"))[c("Final ID", "Diagnosis", paste0("QO",1:100))]
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


dc_subscales = mean(sapply(1:50, function(x) show_decision_tree(data_frame = or_tpq, DV = "Diagnosis",IVs = c(paste0("HA",1:4),paste0("NS",1:4),paste0("RD",1:4)), perc = 0.8, plot = FALSE)[[2]]))
dc_scales = mean(sapply(1:50, function(x) show_decision_tree(data_frame = or_tpq, DV = "Diagnosis",IVs = c("HA", "NS","RD"), perc = 0.8, plot = FALSE)[[2]]))

##############################
# sources
##############################
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

# Experiments -------------------------------------------------------------

# delete later
#read valence data
library(readxl)
valenced = data.frame(read_excel(paste0(path,"/TPQ_DataAndAnalysis/FactorCoding_19.07.2021.xlsx")))
valenced[,1:7] = as.data.frame(apply(valenced[,1:7],2,as.factor))
valenced[,1] = factor(valenced[,1], levels = valenced[,1])

#add valence to sources
sources_new = sources
for (CN in colnames(valenced)[3:7]){
  all_sources_melted[[CN]] = rep(valenced[[CN]],each = 10)
}
ggplot(data = all_sources_melted, aes(x=Question, y= Sources_group, group =Group, color = Active_Passive))+
  geom_bar(stat = "identity")+
  facet_grid(Group~IC_all,scales = "free")+
  geom_text(aes(x= 50,y=6,label = IC_group_plot))+
  geom_text(aes(x= 50,y=10,label = Cor_plot))
# end delete




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

