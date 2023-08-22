import pandas as pd

datatypes = ["all", "ki", "kd"]
modeltypes = [
    "linearRegression",
    "Ridge",
    "Lasso",
    "ElasticNet",
    "SVR",
    "DecisionTree",
    "RandomForest",
]
additional_information = ["final"]
scoretypes = ["delta_G", "Affinity_Data_Value", "pKd_pKi_pIC50"]
signs = ["positive", "negative"]

for sign in signs:
    data_list = []
    model_list = []
    add_info_list = []
    score_list = []
    start = False
    percent1avg_list = []
    percent5avg_list = []
    percent10avg_list = []
    percent1cluster_list = []
    percent5cluster_list = []
    percent10cluster_list = []
    percent1can_list = []
    percent5can_list = []
    percent10can_list = []


    '''
    Summary of the forward screening power: =========================================================
    Average enrichment factor among top 1% = 0.60
    Average enrichment factor among top 5% = 0.89
    Average enrichment factor among top 10% = 0.60
    The best ligand is found among top 1% candidates for  0 cluster(s); success rate = 0.0%
    The best ligand is found among top 5% candidates for  0 cluster(s); success rate = 0.0%
    The best ligand is found among top 10% candidates for  0 cluster(s); success rate = 0.0%
    =================================================================================================
    '''

    for data in datatypes:
        for model in modeltypes:
            for add_info in additional_information:
                for score in scoretypes:
                    output = open(
                        f"results_forward_screening_{sign}/{data}{model}{add_info}{score}.out", "r"
                    )
                    lines = output.readlines()
                    for i, line in enumerate(lines):
                        if line[:39] == "Summary of the forward screening power:":
                            data_list.append(data)
                            model_list.append(model)
                            add_info_list.append(add_info)
                            score_list.append(score)
                            tmp = i
                            start = True
                        elif start == True:
                            if i == tmp + 1:
                                tmp_line = str(line).split(sep=" ")
                                percent1avg = tmp_line[7]
                                percent1avg_list.append(percent1avg)
                            if i == tmp + 2:
                                tmp_line = str(line).split(sep=" ")
                                percent5avg = tmp_line[7]
                                percent5avg_list.append(percent5avg)
                            if i == tmp + 3:
                                tmp_line = str(line).split(sep=" ")
                                percent10avg = tmp_line[7]
                                percent10avg_list.append(percent10avg)
                            if i == tmp + 4:
                                tmp_line = str(line).split(sep=" ")
                                percent1cluster = tmp_line[11]
                                percent1cluster_list.append(percent1cluster)
                                percent1can = tmp_line[16]
                                percent1can_list.append(percent1can)
                            if i == tmp + 5:
                                tmp_line = str(line).split(sep=" ")
                                percent5cluster = tmp_line[11]
                                percent5cluster_list.append(percent5cluster)
                                percent5can = tmp_line[16]
                                percent5can_list.append(percent5can)
                            if i == tmp + 6:
                                tmp_line = str(line).split(sep=" ")
                                percent10cluster = tmp_line[11]
                                percent10cluster_list.append(percent10cluster)
                                percent10can = tmp_line[16]
                                percent10can_list.append(percent10can)

                    print(data,model,add_info,score,"Done")


    data = {
        "Datatype": data_list,
        "Modeltype": model_list,
        "Featuretype": add_info_list,
        "Scoretype": score_list,
        "Top 1'%' average enrichment factor:": percent1avg_list,
        "Top 5'%' average enrichment factor:": percent5avg_list,
        "Top 10'%' average enrichment factor:": percent10avg_list,
        "1'%' cluster number:": percent1cluster_list,
        "1'%' success rate:": percent1can_list,
        "5'%' cluster number:": percent5cluster_list,
        "5'%' success rate:": percent5can_list,
        "10'%' cluster number:": percent10cluster_list,
        "10'%' success rate:": percent10can_list,
    }
    df = pd.DataFrame(data)
    df.to_csv(f"results_forward_screening_{sign}.csv")
