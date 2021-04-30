import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = '/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine/ICASSO/icasso_tpq_reco/bootstrap/'

fname_reco_all = path + 'icasso_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.50_bootstrap_MNEinfomax_data_reco.npz'
fname_reco_hc = path + 'icasso_HC_nsamp1202_n_comp13_n_iter100_dist0.50_bootstrap_MNEinfomax_data_reco.npz'
fname_reco_mdd = path + 'icasso_MDD_nsamp455_n_comp06_n_iter100_dist0.45_bootstrap_MNEinfomax_data_reco.npz'
ids_all = path + "../ids_all_sess1_2.npy"
# load data
npz_all = np.load(fname_reco_all, allow_pickle=True)
npz_hc = np.load(fname_reco_hc, allow_pickle=True)
npz_mdd = np.load(fname_reco_mdd, allow_pickle=True)
npz_ids = np.load(ids_all, allow_pickle=True)

#rearrange groups and IDs according to the IDs in npz_ids
df = pd.read_excel("/home/abdulrahman/abdulrahman.sawalma@gmail.com/PhD/Data/TPQ_DataAndAnalysis/TPQData_ICA_25.11.2020.xlsx")
sub_inds = [np.where(df["ID"] == npz_ids[i])[0][0] for i in range(len(npz_ids))]
groups = [df["Group"][item] for item in sub_inds]
IDs = [df["ID"][item] for item in sub_inds]
session = [df["Session"][item] for item in sub_inds]
#create info_array which contains groups, IDs and session
info_array = np.array([[groups[i], IDs[i], session[i]] for i in range(len(IDs))])

#load subject_data, make subjects into a 2d array of (n_subjects,1) so that it is appended
# to a data of shape (n_subjects,100)

#params = npz['params'].item()              # info about subject groups and ICASSO parameters
#ica_params = npz['ica_params'].item()      # info about additional ICA parameters
#cluster_idx = npz['cluster_idx']           # indices of selected clusters
#data_reco = npz['data_reco']               # data array of reconstructed ICA components



average_all = np.average(npz_all['data_reco'],axis = 1)
average_hc = np.average(npz_hc['data_reco'],axis = 1)
average_mdd = np.average(npz_mdd['data_reco'],axis = 1)


pd.DataFrame(np.transpose(average_all)).to_csv(path+"npz_all_avg.csv", index=False)
pd.DataFrame(np.transpose(average_hc)).to_csv(path+"npz_hc_avg.csv", index=False)
pd.DataFrame(np.transpose(average_mdd)).to_csv(path+"npz_mdd_avg.csv", index=False)

# create data frames for each of the ICs and save it
for i in range(npz_all['data_reco'].shape[0]):
    IC_all = npz_all['data_reco'][i,:,:]

    # add info arrays to data frames
    IC_all = np.append(info_array,IC_all,axis = 1)

    # convert to data frame
    IC_all = pd.DataFrame(IC_all)

    #change column names
    colnames = ["Group","ID","Session"]
    colnames+=list(range(1,101))
    IC_all.columns = colnames

    #save to csv
    IC_all.to_csv(f"{path}IC_arrays/npz_all_{i}.csv", index=False)

for i in range(npz_mdd['data_reco'].shape[0]):
    IC_mdd = npz_mdd['data_reco'][i,:,:]
    pd.DataFrame(IC_mdd).to_csv(f"{path}IC_arrays/npz_mdd_{i}.csv", index=False)

for i in range(npz_hc['data_reco'].shape[0]):
    IC_hc = npz_hc['data_reco'][i,:,:]
    pd.DataFrame(IC_hc).to_csv(f"{path}IC_arrays/npz_hc_{i}.csv", index=False)


# Next issue: does Fast ICA differ between bootstraping and no bootstraping?
path_n_bootstrap = path + "../no_bootstrap/icasso_HC_nsamp1202_n_comp13_n_iter100_dist0.40_no-bootstrap_MNEfastica_data_reco.npz"
path_bootstrap = path + "icasso_HC_nsamp1202_n_comp13_n_iter100_dist0.50_bootstrap_MNEfastica_data_reco.npz"

npz_boot = np.load(path_bootstrap, allow_pickle=True)
npz_n_boot = np.load(path_n_bootstrap, allow_pickle=True)


data_boot = npz_boot['data_reco']
data_n_boot = npz_n_boot['data_reco']


for i in range(data_boot.shape[0]):
    IC_boot = data_boot[i,:,:]
    pd.DataFrame(IC_boot).to_csv(f"{path}../no_bootstrap/IC_arrays/npz_boot_{i}.csv", index=False)


for i in range(data_n_boot.shape[0]):
    IC_n_boot = data_n_boot[i,:,:]
    pd.DataFrame(IC_n_boot).to_csv(f"{path}../no_bootstrap/IC_arrays/npz_n_boot_{i}.csv", index=False)
