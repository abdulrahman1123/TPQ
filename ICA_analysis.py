import pandas as pd
import numpy as np
import statistics

name_hc = "/home/abdul-rahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine/TPQ_Final_ABD_22.04.2020_HC_ica-logcosh_ncomp_13.npz"
name_mdd = "/home/abdul-rahman/abdulrahman.sawalma@gmail.com/PhD/Data/Palestine/TPQ_Final_ABD_22.04.2020_MDD_ica-logcosh_ncomp_13.npz"
dinfo_hc = np.load(name_hc, allow_pickle= True)['dinfo'].item()
dinfo_mdd = np.load(name_mdd, allow_pickle= True)['dinfo'].item()

mdd_ica = dinfo_mdd['ic_data']
hc_ica = dinfo_hc['ic_data']


ica_thresh_hc = dinfo_hc['data']
ica_thresh_mdd = dinfo_mdd['data']

Question_indices_hc = [list(np.where(abs(column)>0)[0]) for column in ica_thresh_hc]
Question_indices_mdd = [list(np.where(abs(column)>0)[0]) for column in ica_thresh_mdd]

#This is just a small code to extract the required code for R (So that I don't write them by hand and check too many times)
Question_labels_hc = ["TPQ_Q$IC_"+str(i)+"_HC = "+"".join([")+as.numeric(TPQ_Q$QO"+ str(int(item)+1) for item in (np.where(abs(ica_thresh_hc[i])>0)[0]).astype(str).tolist()])+")" for i in range(13)]
Question_labels_mdd = ["TPQ_Q$IC_"+str(i)+"_MDD = "+"".join([")+as.numeric(TPQ_Q$QO"+ str(int(item)+1) for item in (np.where(abs(ica_thresh_mdd[i])>0)[0]).astype(str).tolist()])+")" for i in range(13)]
print ("\n".join(Question_labels_hc).replace("= )+","= ").replace("= )","= NA"))
print ("\n".join(Question_labels_mdd).replace("= )+","= ").replace("= )","= NA"))


# The following is the same code as above, but I will be adding the values of each ICA
Question_labels_hc = ["TPQ_Q$Complete_IC_"+str(i)+"_HC = "+"".join([")+as.numeric(TPQ_Q$QO"+ str(int(item)+1) for item in (np.where(abs(mdd_ica[i])>0)[0]).astype(str).tolist()])+")" for i in range(13)]
Question_labels_mdd = ["TPQ_Q$Complete_IC_"+str(i)+"_MDD = "+"".join([")+as.numeric(TPQ_Q$QO"+ str(int(item)+1) for item in (np.where(abs(mdd_ica[i])>0)[0]).astype(str).tolist()])+")" for i in range(13)]
print ("\n".join(Question_labels_hc).replace("= )+","= "))
print ("\n".join(Question_labels_mdd).replace("= )+","= "))


Question_labels_hc = ["A"+str(i)"".join([")+as.numeric(TPQ_Q$QO"+ str(int(item)+1) for item in (np.where(abs(mdd_ica[i])>0)[0]).astype(str).tolist()])+")" for i in range(13)]

[print (component) for component in mdd_ica[i] for i in range(13)]

Question_labels_hc_complete = ["TPQ_Q$Complete_IC_"+str(k)+"_HC = " + "".join([")+ "+str(round(hc_ica[k][i],3))+" * as.numeric(TPQ_Q$QO"+str(1+i) for i in range(100)])+")" for k in range(13)]
Question_labels_mdd_complete = ["TPQ_Q$Complete_IC_"+str(k)+"_MDD = "+"".join([")+ "+str(round(mdd_ica[k][i],3))+" * as.numeric(TPQ_Q$QO"+str(1+i) for i in range(100)]) +")"for k in range(13)]
print ("\n".join(Question_labels_hc_complete).replace("= )+","= "))
print ("\n".join(Question_labels_mdd_complete).replace("= )+","= "))



statistics.mean(hc_ica[4])
statistics.mean(mdd_ica[4])