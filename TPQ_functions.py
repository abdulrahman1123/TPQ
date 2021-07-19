import numpy as np
import pandas as pd
import os

def prepare_icasso(npz_tpq_path):
    """
    prepare the data from the npz files sent to me by Juergen.
    It creates the following five data frames/groups of data frames:
    data_tpq, data_reco, scores, projections and sources
    All DFs will have column names indicating what they represent
    :param npz_tpq_path: the path to the tpq npz (which is the one that contains the basic info such as IDs)
    :return: five data frames (see description)
    """
    # Import original scores, there you can find information about diagnosis and other stuff
    scores_path = "/home/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive/PhD/Data/Palestine/TPQ_DataAndAnalysis/TPQ_Analysis_All_25.11.2020_modified.xlsx"
    if not os.path.exists(scores_path):
        scores_path = scores_path.replace("/abdulrahman/abdulrahman.sawalma@gmail.com","/asawalma/Insync/abdulrahman.sawalma@gmail.com/Google Drive")
    or_scores = pd.read_excel(scores_path)[
        ["Diagnosis", "Final ID", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session"]]
    or_scores.loc[pd.isna(or_scores["Trauma"]), or_scores.columns == "Trauma"] = "NA"
    or_scores.loc[pd.isna(or_scores["GAD"]), or_scores.columns == "GAD"] = "NA"
    or_scores.loc[pd.isna(or_scores["Response"]), or_scores.columns == "Response"] = "NA"
    or_scores.loc[pd.isna(or_scores["PTSD"]), or_scores.columns == "PTSD"] = "NA"
    or_scores.loc[pd.isna(or_scores["MDD"]), or_scores.columns == "MDD"] = "NA"

    npz_tpq = np.load(npz_tpq_path, allow_pickle=True)
    results_path = npz_tpq_path.replace("icasso_ICA-tpq", "icasso-results_ICA-tpq").replace("icasso_ICA-subjects",
                                                                                            "icasso-results_ICA-subjects")
    npz_results = np.load(results_path, allow_pickle=True)
    # extract IDs
    IDs = npz_tpq["IDs"]
    diagnoses = []
    trauma = []
    gad = []
    response = []
    ptsd = []
    mdd = []
    session = []
    for i in range(len(IDs)):
        if IDs[i] in list(or_scores["Final ID"]):
            diagnoses.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "Diagnosis"])[0])
            trauma.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "Trauma"])[0])
            gad.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "GAD"])[0])
            response.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "Response"])[0])
            ptsd.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "PTSD"])[0])
            mdd.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "MDD"])[0])
            session.append(list(or_scores.loc[or_scores["Final ID"] == IDs[i], "Session"])[0])
        else:
            diagnoses.append("NA")
            trauma.append("NA")
            gad.append("NA")
            response.append("NA")
            ptsd.append("NA")
            mdd.append("NA")
            session.append("NA")
    # create two data frames to be added to the data inside each dictionary item
    sub_info_2d = [[IDs[i], diagnoses[i], trauma[i], gad[i], response[i], ptsd[i], mdd[i], session[i]] for i in
                   range(len(IDs))]
    sub_info_3d = [sub_info_2d] * npz_results["data_reco"].shape[0]
    sub_info = pd.DataFrame(sub_info_2d,
                            columns=["ID", "Diagnosies", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session"])
    # read all data
    data_tpq = pd.DataFrame(npz_tpq["data_tpq"])
    data_tpq[["ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session"]] = sub_info
    # rename the columns
    data_tpq.columns = ["Q" + str(i) for i in range(1, 101)] + ["ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD",
                                                                "MDD", "Session"]
    # rearrange the columns, just for aesthetics
    new_colnames = ["ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session"] + ["Q" + str(i) for i in
                                                                                                 range(1, 101)]
    data_tpq = data_tpq[new_colnames]
    data_reco = npz_results["data_reco"]
    data_reco = [pd.DataFrame(item, columns=new_colnames) for item in
                 np.append(np.array(sub_info_3d), data_reco, axis=2)]
    scores = npz_results["scores"]
    projections = npz_results['projection']
    projections = pd.DataFrame(np.append(np.array(sub_info_2d), projections, axis=1))
    projections.columns = ["ID", "Diagnosis", "Trauma", "GAD", "Response", "PTSD", "MDD", "Session"] + ["IC" + str(i)
                                                                                                        for i in
                                                                                                        range(1,
                                                                                                              len(projections.columns) - 7)]
    sources = npz_results["sources"]
    sources = pd.DataFrame(sources, columns=["Q" + str(i) for i in range(1, 101)])
    return (data_tpq, data_reco, scores, projections, sources)



