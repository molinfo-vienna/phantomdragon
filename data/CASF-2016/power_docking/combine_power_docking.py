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
signs = ["positive","negative"]

data_list = []
model_list = []
add_info_list = []
score_list = []
start = False
top1_number_list = []
top1_percent_list = []
top2_number_list = []
top2_percent_list = []
top3_number_list = []
top3_percent_list = []
spearman_02_list = []
spearman_03_list = []
spearman_04_list = []
spearman_05_list = []
spearman_06_list = []
spearman_07_list = []
spearman_08_list = []
spearman_09_list = []
spearman_010_list = []
conf_02_list = []
conf_03_list = []
conf_04_list = []
conf_05_list = []
conf_06_list = []
conf_07_list = []
conf_08_list = []
conf_09_list = []
conf_010_list = []

for sign in signs:
    for data in datatypes:
        for model in modeltypes:
            for add_info in additional_information:
                for score in scoretypes:
                    output = open(
                        f"results_power_docking_{sign}/{data}{model}{add_info}{score}.out", "r"
                    )
                    lines = output.readlines()
                    for i, line in enumerate(lines):
                        if line[:29] == "Summary of the docking power:":
                            data_list.append(data)
                            model_list.append(model)
                            add_info_list.append(add_info)
                            score_list.append(score)
                            tmp = i
                            start = True
                        elif start == True:
                            if i == tmp + 2:
                                tmp_line = str(line).split(sep=" ")
                                top1_number = tmp_line[6]
                                top1_number = top1_number[0:-1]
                                top1_percent = tmp_line[-1]
                                top1_percent = top1_percent[0:-1]
                                top1_number_list.append(top1_number)
                                top1_percent_list.append(top1_percent)
                            if i == tmp + 4:
                                tmp_line = str(line).split(sep=" ")
                                top2_number = tmp_line[6]
                                top2_number = top2_number[0:-1]
                                top2_percent = tmp_line[-1]
                                top2_percent = top2_percent[0:-1]
                                top2_number_list.append(top2_number)
                                top2_percent_list.append(top2_percent)
                            if i == tmp + 6:
                                tmp_line = str(line).split(sep=" ")
                                top3_number = tmp_line[6]
                                top3_number = top3_number[0:-1]
                                top3_percent = tmp_line[-1]
                                top3_percent = top3_percent[0:-1]
                                top3_number_list.append(top3_number)
                                top3_percent_list.append(top3_percent)
                            if i == tmp + 7:
                                tmp_line = str(line).split(sep=" ")
                                spearman_02 = tmp_line[-1]
                                spearman_02 = spearman_02[0:-1]
                                spearman_02_list.append(spearman_02)
                            if i == tmp + 8:
                                tmp_line = str(line).split()
                                conf_02 = tmp_line[-3:]
                                conf_02 = str(conf_02[0]) + str(conf_02[1]) + str(conf_02[2])
                                conf_02_list.append(conf_02)
                            if i == tmp + 9:
                                tmp_line = str(line).split(sep=" ")
                                spearman_03 = tmp_line[-1]
                                spearman_03 = spearman_03[0:-1]
                                spearman_03_list.append(spearman_03)
                            if i == tmp + 10:
                                tmp_line = str(line).split()
                                conf_03 = tmp_line[-3:]
                                conf_03 = str(conf_03[0]) + str(conf_03[1]) + str(conf_03[2])
                                conf_03_list.append(conf_03)
                            if i == tmp + 11:
                                tmp_line = str(line).split(sep=" ")
                                spearman_04 = tmp_line[-1]
                                spearman_04 = spearman_04[0:-1]
                                spearman_04_list.append(spearman_04)
                            if i == tmp + 12:
                                tmp_line = str(line).split()
                                conf_04 = tmp_line[-3:]
                                conf_04 = str(conf_04[0]) + str(conf_04[1]) + str(conf_04[2])
                                conf_04_list.append(conf_04)
                            if i == tmp + 13:
                                tmp_line = str(line).split(sep=" ")
                                spearman_05 = tmp_line[-1]
                                spearman_05 = spearman_05[0:-1]
                                spearman_05_list.append(spearman_05)
                            if i == tmp + 14:
                                tmp_line = str(line).split()
                                conf_05 = tmp_line[-3:]
                                conf_05 = str(conf_05[0]) + str(conf_05[1]) + str(conf_05[2])
                                conf_05_list.append(conf_05)
                            if i == tmp + 15:
                                tmp_line = str(line).split(sep=" ")
                                spearman_06 = tmp_line[-1]
                                spearman_06 = spearman_06[0:-1]
                                spearman_06_list.append(spearman_06)
                            if i == tmp + 16:
                                tmp_line = str(line).split()
                                conf_06 = tmp_line[-3:]
                                conf_06 = str(conf_06[0]) + str(conf_06[1]) + str(conf_06[2])
                                conf_06_list.append(conf_06)
                            if i == tmp + 17:
                                tmp_line = str(line).split(sep=" ")
                                spearman_07 = tmp_line[-1]
                                spearman_07 = spearman_07[0:-1]
                                spearman_07_list.append(spearman_07)
                            if i == tmp + 18:
                                tmp_line = str(line).split()
                                conf_07 = tmp_line[-3:]
                                conf_07 = str(conf_07[0]) + str(conf_07[1]) + str(conf_07[2])
                                conf_07_list.append(conf_07)
                            if i == tmp + 19:
                                tmp_line = str(line).split(sep=" ")
                                spearman_08 = tmp_line[-1]
                                spearman_08 = spearman_08[0:-1]
                                spearman_08_list.append(spearman_08)
                            if i == tmp + 20:
                                tmp_line = str(line).split()
                                conf_08 = tmp_line[-3:]
                                conf_08 = str(conf_08[0]) + str(conf_08[1]) + str(conf_08[2])
                                conf_08_list.append(conf_08)
                            if i == tmp + 21:
                                tmp_line = str(line).split(sep=" ")
                                spearman_09 = tmp_line[-1]
                                spearman_09 = spearman_09[0:-1]
                                spearman_09_list.append(spearman_09)
                            if i == tmp + 22:
                                tmp_line = str(line).split()
                                conf_09 = tmp_line[-3:]
                                conf_09 = str(conf_09[0]) + str(conf_09[1]) + str(conf_09[2])
                                conf_09_list.append(conf_09)
                            if i == tmp + 23:
                                tmp_line = str(line).split(sep=" ")
                                spearman_010 = tmp_line[-1]
                                spearman_010 = spearman_010[0:-1]
                                spearman_010_list.append(spearman_010)
                            if i == tmp + 24:
                                tmp_line = str(line).split()
                                conf_010 = tmp_line[-3:]
                                conf_010 = str(conf_010[0]) + str(conf_010[1]) + str(conf_010[2])
                                conf_010_list.append(conf_010)
                            
                    print(data,model,add_info,score,"Done")


    data = {
        "Datatype": data_list,
        "Modeltype": model_list,
        "Featuretype": add_info_list,
        "Scoretype": score_list,
        "Top1: # correct binding poses:": top1_number_list,
        "Top1: '%' correct binding poses:": top1_percent_list,
        "Top2: # correct binding poses:": top2_number_list,
        "Top2: '%' correct binding poses:": top2_percent_list,
        "Top3: # correct binding poses:": top3_number_list,
        "Top3: '%' correct binding poses:": top3_percent_list,
        "Spearman correlation coefficient in rmsd range [0-2]:": spearman_02_list,
        "90'%' confidence interval in rmsd range [0-2]:": conf_02_list,
        "Spearman correlation coefficient in rmsd range [0-3]:": spearman_03_list,
        "90'%' confidence interval in rmsd range [0-3]:": conf_03_list,
        "Spearman correlation coefficient in rmsd range [0-4]:": spearman_04_list,
        "90'%' confidence interval in rmsd range [0-4]:": conf_04_list,
        "Spearman correlation coefficient in rmsd range [0-5]:": spearman_05_list,
        "90'%' confidence interval in rmsd range [0-5]:": conf_05_list,
        "Spearman correlation coefficient in rmsd range [0-6]:": spearman_06_list,
        "90'%' confidence interval in rmsd range [0-6]:": conf_06_list,
        "Spearman correlation coefficient in rmsd range [0-7]:": spearman_07_list,
        "90'%' confidence interval in rmsd range [0-7]:": conf_07_list,
        "Spearman correlation coefficient in rmsd range [0-8]:": spearman_08_list,
        "90'%' confidence interval in rmsd range [0-8]:": conf_08_list,
        "Spearman correlation coefficient in rmsd range [0-9]:": spearman_09_list,
        "90'%' confidence interval in rmsd range [0-9]:": conf_09_list,
        "Spearman correlation coefficient in rmsd range [0-10]:": spearman_010_list,
        "90'%' confidence interval in rmsd range [0-10]:": conf_010_list,
    }
    df = pd.DataFrame(data)
    df.to_csv(f"results_power_docking_{sign}.csv")
