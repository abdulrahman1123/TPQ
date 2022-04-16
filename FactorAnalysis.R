# Pythoon code for factor analysis ----------------------------------------
library(reticulate)
library(ggplot2)
library(psych)
library(readxl)
library(lavaan)

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

pd = import("pandas")
np = import ("numpy")
FactorAnalysis = import("sklearn.decomposition")$FactorAnalysis
cross_val_score = import("sklearn.model_selection")$cross_val_score

n_comp_list = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25)

name_xls = 'TPQData_ICA_25.11.2020.xlsx'
sheet_name = 'Sheet1'

stage = "/home/asawalma/git/Data/TPQ-Analysis/"
path_data = paste0(stage, 'TPQ_data/')
path_fig  = paste0(path_data, 'figures/')
path_results = paste0(path_data, 'results/')
fname_xls = paste0(path_data, name_xls)
fndata = gsub(".xlsx", '.npz',fname_xls)

sessions = c(1)
groups = c('HC')
kfolds = c(3) #c(3, 5, 10)

n_perm = 1        # apply CV after permuting the subjects

df = as.data.frame(readxl::read_excel(fname_xls))
df = df[(df$Session %in% sessions & df$Group %in% groups),4:103]
data = ifelse(df == "T",1,ifelse(df == "F",-1,0))
data[is.na(data)] = 0

# Info about data
n_samples = nrow(data)

# subtract mean value
dmean = sapply(apply(data,2,mean), function(x) rep(x,n_samples)) # create a means data frame with same shape as original data
data = data - dmean
# Python part
g_plots = list()
plot_result = FALSE
for (cv in kfolds){
  data_cv = t(data)

  rank_list = vector()
  n_samples = ncol(data_cv)  # data must be of shape [n_features, n_samples]
  
  fold_size = as.integer(n_samples / cv)
  
  for (iperm in 1:n_perm){
    data_perm = data_cv
  
    # Do factor analysis
    f_a = FactorAnalysis()
    fa_scores = vector()
    for (n in n_comp_list) {
      f_a$n_components = as.integer(n)
      fa_scores = append(fa_scores, np$mean(cross_val_score(f_a, t(data_cv), cv=as.integer(cv))))
      f_a$fit(t(data_cv))
      print(f_a$components_[1:6,1])
      if (plot_result){
        ss = f_a$components_
        ss = sapply(1:ncol(ss), function(x) ifelse(abs(ss[,x])==max(abs(ss[,x])),1,0))
        X = 1:nrow(ss)
        Y= apply(ss,1,sum)
        df = data.frame(X=X,Y=Y)
        gg=  ggplot(df, aes(x = X, y = Y))+ geom_line()+ geom_point()+
          geom_text(aes(y=Y+1),label = cumsum(Y))
        g_plots = append(g_plots,list(gg))
      }
    }
    
    n_components_fa = n_comp_list[fa_scores == max(fa_scores)]
    print(fa_scores)
    
    rank = n_components_fa
    
    line = paste0('FA rank estimation (cv=',cv,', n_samples=', nrow(data_perm),', fold size=',
                 fold_size,', perm.',iperm, '): rank = ',rank)
    rank_list = append(rank_list,rank)
  }
  rank_list_fa = rank_list
  fold_sz = as.integer(n_samples/cv)
  line1 = paste0(' FA rank estimation (cv=',cv ,', n_samples=',n_samples ,', fold size=',fold_sz ,', perm.',n_perm ,'):')
  line2 = paste0('   mean rank:  ',mean(rank_list_fa) ,' (+- ',sd(rank_list_fa) ,')')
  line3 = paste0('   range:  [',min(rank_list_fa) ,' - ',max(rank_list_fa) ,']')
  print (paste0(line1,line2,line3))
}
multiplot(plotlist=g_plots, Listing=TRUE, cols = 2)


# R part for comparison
g_plots = list()
for (cv in kfolds){
  data_cv = t(data)
  
  rank_list = vector()
  n_samples = ncol(data_cv)  # data must be of shape [n_features, n_samples]
  
  fold_size = as.integer(n_samples / cv)
  
  for (iperm in 1:n_perm){
    data_perm = data_cv
    
    # Do factor analysis
    f_a = FactorAnalysis()
    fa_scores = vector()
    for (n in n_comp_list) {
      f_a$n_components = as.integer(n)
      fa_scores = append(fa_scores, np$mean(cross_val_score(f_a, t(data_cv), cv=as.integer(cv))))
      f_a$fit(t(data_cv))
      ss = f_a$components_
      ss = sapply(1:ncol(ss), function(x) ifelse(abs(ss[,x])==max(abs(ss[,x])),1,0))
      X = 1:nrow(ss)
      Y= apply(ss,1,sum)
      df = data.frame(X=X,Y=Y)
      gg=  ggplot(df, aes(x = X, y = Y))+ geom_line()+ geom_point()+
        geom_text(aes(y=apply(ss,1,sum)+1),label = cumsum(Y))
      g_plots = append(g_plots,list(gg))
    }
    
    
    n_components_fa = n_comp_list[fa_scores == max(fa_scores)]
    print(fa_scores)
    
    rank = n_components_fa
    
    line = paste0('FA rank estimation (cv=',cv,', n_samples=', nrow(data_perm),', fold size=',
                  fold_size,', perm.',iperm, '): rank = ',rank)
    rank_list = append(rank_list,rank)
  }
  rank_list_fa = rank_list
  fold_sz = as.integer(n_samples/cv)
  line1 = paste0(' FA rank estimation (cv=',cv ,', n_samples=',n_samples ,', fold size=',fold_sz ,', perm.',n_perm ,'):')
  line2 = paste0('   mean rank:  ',mean(rank_list_fa) ,' (+- ',sd(rank_list_fa) ,')')
  line3 = paste0('   range:  [',min(rank_list_fa) ,' - ',max(rank_list_fa) ,']')
  print (paste0(line1,line2,line3))
}

multiplot(plotlist=g_plots, Listing=TRUE, cols = n_perm)


# Factor Analysis using R -------------------------------------------------


for (n in 2:20){
  factanal.none = fa(r=data+dmean, 
                      nfactors = n, 
                      # covar = FALSE
                      #SMC = TRUE,
                      fm="ml",
                      max.iter=50,
                      rotate="none")
  print(factanal.none$RMSEA)
  loadings = as.data.frame(t(matrix(factanal.none$loadings, byrow = FALSE, nrow = 100)))
  #head(loadings[,1:5])
  loadings = sapply(1:ncol(loadings), function(x) ifelse(abs(loadings[,x])==max(abs(loadings[,x])),1,0))
  X = 1:nrow(loadings)
  Y= apply(loadings,1,sum)
  print(cumsum(Y))
}
n=25
#factanal.none = factanal(data+dmean, factors=n, rotation = "none")


# Exploratory/Confirmatory factor analysis -------------------------------------------------
library(readxl)
library(lavaan)

# Exploratory Factor analysis
subscales_df = data.frame(read_xlsx('/home/asawalma/git/tpq_analysis/data/Subscales.xlsx'))
fit = factanal(mydata, 3, rotation = 'varimax')
print(fit, digits=3, cutoff=.3, sort=FALSE)
# plot factor 1 by factor 2
load = fit$loadings[,1:3]

# Confirmatory factor anlaysis
HS.model = ' NS =~ NS1+NS2+NS3+NS4
              HA =~ HA1+HA2+HA3+HA4
              RD =~ RD1+RD2+RD3+RD4'

fit = cfa(HS.model, data=subscales_df)
ss = summary(fit, fit.measures=TRUE)
sum_df = data.frame(df = ss$FIT['df'])
for (item in c('chisq','pvalue','cfi','rmsea','srmr')){
  sum_df[item] = round(ss$FIT[item],3)
}


# replaice the psych library function for parallel analysis -------------------------------------------------
library(psych)
library(readxl)

data = as.data.frame(read_xlsx('/home/asawalma/git/tpq_analysis/data/Tpq Scores.xlsx'))
data = data[data$Session == "Test",]
data = data[,113:212]
data[data=="T"] = 1
data[data=="F"] = -1
data[is.na(data) | data == "NA"] = 0

data = as.data.frame(apply(data,2,as.numeric))

fa.parallel(data, fm = 'ml', nfactors = 1)

# The function in psych is this
fm = "minres"
nfactors = 1
main = "Parallel Analysis Scree Plots"
n.iter = 20
quant = 0.95



nsub = dim(x)[1]
nvariables = dim(x)[2]

x = data

rx = cor(x, use = 'pairwise', method = 'pearson') 

valuesx = eigen(rx)$values
fa.valuesx = fa(rx, nfactors = 1, rotate = "none", 
                 fm = fm, warnings = FALSE)$values


sim_list = list()
for (i in 1:n.iter){
  sampledata = matrix(apply(x, 2, function(y) sample(y, nsub, replace = TRUE)), ncol = nvariables)
  colnames(sampledata) = colnames(x)
  C = cor(sampledata, use = 'pairwise', method = "pearson")

  samp = eigen(C)$values
  samp.fa = fa(C, fm = fm, nfactors = 1, SMC = FALSE)$values
  
  simdata = matrix(rnorm(nsub * nvariables), nrow = nsub, 
                   ncol = nvariables)
  sim.cor = cor(simdata)
  simulated_eigen = eigen(sim.cor)$values
  sim.fa = fa(sim.cor, fm = fm, nfactors = 1, 
                        rotate = "none", SMC = FALSE, warnings = FALSE)$values

  sim_list[[i]] = list(samp = samp, samp.fa = samp.fa, sim = simulated_eigen, sim.fa = sim.fa)
}



ylabel = "eigenvalues of principal components and factor analysis"

values = t(matrix(unlist(sim_list), ncol = n.iter))
values.sim.mean = colMeans(values, na.rm = TRUE)
values.ci = apply(values, 2, function(x) quantile(x, quant))

values.sim.se = apply(values, 2, sd, na.rm = TRUE)

ymin = min(valuesx, values.sim.mean)
ymax = max(valuesx, values.sim.mean)
sim.pcr =NA
sim.far = NA

plot(valuesx, type = "b", main = main, ylab = ylabel, 
     ylim = c(ymin, ymax), xlab = "Factor/Component Number", 
     pch = 4, col = "blue")
points(fa.valuesx, type = "b", pch = 2, col = "blue")


sim.pcr = values.sim.mean[1:nvariables]
sim.pcr.ci = values.ci[1:nvariables]
sim.se.pcr = values.sim.se[1:nvariables]
sim.far = values.sim.mean[(nvariables + 1):(2 * nvariables)]
sim.se.far = values.sim.se[(nvariables + 1):(2 * nvariables)]
sim.far.ci = values.ci[(nvariables + 1):(2 * nvariables)]
sim.pc = values.sim.mean[(2 * nvariables + 1):(3 * nvariables)]
sim.pc.ci = values.ci[(2 * nvariables + 1):(3 * nvariables)]
sim.se.pc = values.sim.se[(2 * nvariables + 1):(3 * nvariables)]
sim.fa = values.sim.mean[(3 * nvariables + 1):(4 * nvariables)]
sim.fa.ci = values.ci[(3 * nvariables + 1):(4 * nvariables)]
sim.se.fa = values.sim.se[(3 * nvariables + 1):(4 * nvariables)]
pc.test = which(!(valuesx > sim.pcr.ci))[1] - 1
fa.test = which(!(fa.valuesx > sim.far.ci))[1] - 1

points(sim.pc, type = "l", lty = "dotted", pch = 4,  col = "red")
points(sim.fa, type = "l", lty = "dotted", pch = 4, col = "red")
points(sim.pcr, type = "l", lty = "dashed", pch = 2, col = "red")
points(sim.far, type = "l", lty = "dashed", pch = 2, col = "red")

pc.test = which(!(valuesx > sim.pc.ci))[1] - 1
fa.test = which(!(fa.valuesx > sim.fa.ci))[1] - 1



points(sim.pcr, type = "l", lty = "dashed", pch = 4, col = "red")
points(sim.far, type = "l", lty = "dashed", pch = 4, col = "red")



legend("topright", c("  PC  Actual Data", 
                     "  PC  Simulated Data", " PC  Resampled Data", 
                     "  FA  Actual Data", "  FA  Simulated Data", 
                     " FA  Resampled Data"), col = c("blue", 
                                                     "red", "red", "blue", "red", "red"), pch = c(4, 
                                                                                                  NA, NA, 2, NA, NA), text.col = "green4", 
       lty = c("solid", "dotted", "dashed", "solid", 
               "dotted", "dashed"), merge = TRUE, bg = "gray90")

colnames(values) = paste0("Sim", 1:ncol(values))
abline(h = 1)
results = list(fa.values = fa.valuesx, pc.values = valuesx, 
                pc.sim = sim.pc, pc.simr = sim.pcr, fa.sim = sim.fa, 
                fa.simr = sim.far, nfact = fa.test, ncomp = pc.test)

colnames(values)[1:(2 * nvariables)] = c(paste0("C", 
                                                 1:nvariables), paste0("F", 1:nvariables))
colnames(values)[(2 * nvariables + 1):ncol(values)] = c(paste0("CSim", 
                                                                1:nvariables), paste0("Fsim", 1:nvariables))

results$nfact = fa.test

results$ncomp = pc.test
results$values = values
cat("Parallel analysis suggests that ")
cat("the number of factors = ", fa.test, " and the number of components = ", 
    pc.test, "\n")
class(results) = c("psych", "parallel")
final_result = invisible(results)
