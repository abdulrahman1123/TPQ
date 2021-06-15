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


# The New ICASSO with fixed results
import numpy as np
import pandas as pd
import os

# create the main dir, which is different in different PCs
main_dir = "/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine"
if not os.path.exists(main_dir):
    main_dir = main_dir.replace("/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive","/home/abdulrahman/abdulrahman.sawalma@gmail.com")

def match_value(Column,IDs):
    # Search for matching IDs in the original data and return the matched information of the column "Column"
    final_arranged = []
    for item in IDs:
        if item in original_data["Final ID"].tolist():
            final_arranged.append(
                np.array(original_data.loc[original_data["Final ID"] == item, original_data.columns == Column])[
                    0, 0])
        else:
            final_arranged.append("NA")
    return (final_arranged)


original_data = pd.read_excel(main_dir + "/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx")


# read tpq data, which includes general information about the ICA and the IDs

def extract_data(tpq_npz_dir):
    results_npz_dir = tpq_npz_dir.replace("icasso_ICA-tpq","icasso-results_ICA-tpq")

    data_tpq = np.load(tpq_npz_dir, allow_pickle=True)
    data_results = np.load(results_npz_dir, allow_pickle=True)

    IDs = data_tpq["IDs"]

    # match the values of the new columns (Diagnosis, Trauma ... etc.) with the IDs within the data
    sub_info = np.array([[match_value(item,IDs)]*data_results["data_reco"].shape[0] for item in ['Diagnosis', 'Trauma', 'GAD', 'Response', 'PTSD','MDD', 'Session', 'Final ID']])
    sub_info = sub_info.transpose(1,2,0)

    # put the information in data frames after attaching the sub_info to the beginning
    data_reco_array = np.append(sub_info,data_results["data_reco"],axis=2)
    data_reco = {"IC"+str(i) : pd.DataFrame(data_reco_array[i]) for i in range(len(data_reco_array))}
    sources = pd.DataFrame(data_results["sources"])
    scores = pd.DataFrame(data_results["scores"])
    projection = pd.DataFrame(np.append(sub_info[0,:,:],data_results["projection"],axis=1))

    return(data_reco, sources,scores,projection)

tpq_npz_dir_GAD =main_dir + "/ICASSO_fixed/decomp_tpq/GAD/icasso_ICA-tpq_GAD_nsamp30_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"
tpq_npz_dir_hc =main_dir + "/ICASSO_fixed/decomp_tpq/HC/icasso_ICA-tpq_HC_nsamp1202_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"
tpq_npz_dir_all =main_dir + "/ICASSO_fixed/decomp_tpq/HC,MDD,PTSD,TNP,GAD_infomax/icasso_ICA-tpq_HC,MDD,PTSD,TNP,GAD_nsamp1822_n_comp13_n_iter100_dist0.40_MNEinfomax.npz"
tpq_npz_dir_mdd =main_dir + "/ICASSO_fixed/decomp_tpq/MDD/icasso_ICA-tpq_MDD_nsamp455_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"
tpq_npz_dir_disorders =main_dir + "/ICASSO_fixed/decomp_tpq/MDD,PTSD,TNP,GAD/icasso_ICA-tpq_MDD,PTSD,TNP,GAD_nsamp620_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"
tpq_npz_dir_ptsd =main_dir + "/ICASSO_fixed/decomp_tpq/PTSD/icasso_ICA-tpq_PTSD_nsamp57_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"
tpq_npz_dir_tnp =main_dir + "/ICASSO_fixed/decomp_tpq/TNP/icasso_ICA-tpq_TNP_nsamp78_n_comp13_n_iter100_dist0.30_MNEinfomax.npz"

data_reco_GAD, sources_GAD,scores_GAD,projection_GAD = extract_data(tpq_npz_dir = tpq_npz_dir_GAD)
data_reco_hc, sources_hc,scores_hc,projection_hc = extract_data(tpq_npz_dir = tpq_npz_dir_hc)
data_reco_all, sources_all,scores_all,projection_all = extract_data(tpq_npz_dir = tpq_npz_dir_all)
data_reco_mdd, sources_mdd,scores_mdd,projection_mdd = extract_data(tpq_npz_dir = tpq_npz_dir_mdd)
data_reco_disorders, sources_disorders,scores_disorders,projection_disorders = extract_data(tpq_npz_dir = tpq_npz_dir_disorders)
data_reco_ptsd, sources_ptsd,scores_ptsd,projection_ptsd = extract_data(tpq_npz_dir = tpq_npz_dir_ptsd)
data_reco_tnp, sources_tnp,scores_tnp,projection_tnp = extract_data(tpq_npz_dir = tpq_npz_dir_tnp)


tpq_npz_dir = tpq_npz_dir_all
results_npz_dir = tpq_npz_dir.replace("icasso_ICA-tpq","icasso-results_ICA-tpq")

npz_tpq = np.load(tpq_npz_dir, allow_pickle=True)
data_results = np.load(results_npz_dir, allow_pickle=True)

IDs = npz_tpq["IDs"]

# niterrupt
# Search for matching IDs in the original data and return the matched information of the column "Column"
Column = "Diagnosis"
final_arranged = []
for item in IDs:
    if item in original_data["Final ID"].tolist():
        final_arranged.append(
            np.array(original_data.loc[original_data["Final ID"] == item, original_data.columns == Column])[
                0, 0])
    else:
        final_arranged.append("NA")
#end einterrupt
# match the values of the new columns (Diagnosis, Trauma ... etc.) with the IDs within the data
sub_info = np.array([[match_value(item,IDs)]*data_results["data_reco"].shape[0] for item in ['Diagnosis', 'Trauma', 'GAD', 'Response', 'PTSD','MDD', 'Session', 'Final ID']])
sub_info = sub_info.transpose(1,2,0)

data_tpq = pd.DataFrame(npz_tpq["data_tpq"])
data_tpq[["Diagnosis","Trauma", "GAD", "Response", "PTSD", "MDD", "Session","ID"]] = sub_info[0]
# put the information in data frames after attaching the sub_info to the beginning
data_reco_array = np.append(sub_info,data_results["data_reco"],axis=2)
data_reco = {"IC"+str(i) : pd.DataFrame(data_reco_array[i]) for i in range(len(data_reco_array))}
sources = pd.DataFrame(data_results["sources"])
scores = pd.DataFrame(data_results["scores"])
projection = pd.DataFrame(np.append(sub_info[0,:,:],data_results["projection"],axis=1))

return(data_reco, sources,scores,projection)