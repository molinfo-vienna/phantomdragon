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

data_list = []
model_list = []
add_info_list = []
score_list = []
start = False
SP_list = []
tau_list =[]
PI_list = []
SPscipy_list = []
conf_list = []

for data in datatypes:
    for model in modeltypes:
        for add_info in additional_information:
            for score in scoretypes:
                output = open(
                    f"results_power_ranking/{data}{model}{add_info}{score}.out", "r"
                )
                lines = output.readlines()
                for i, line in enumerate(lines):
                    if line[:29] == "Summary of the ranking power:":
                        data_list.append(data)
                        model_list.append(model)
                        add_info_list.append(add_info)
                        score_list.append(score)
                        tmp = i
                        start = True
                    elif start == True:
                        if i == tmp + 1:
                            tmp_line = str(line).split(sep="=")
                            blub = tmp_line[-1]
                            SP_list.append(blub)
                        if i == tmp + 2:
                            tmp_line = str(line).split(sep="=")
                            blub = tmp_line[-1]
                            SPscipy_list.append(blub)
                        if i == tmp + 3:
                            tmp_line = str(line).split(sep="=")
                            blub = tmp_line[-1]
                            conf_list.append(blub)
                        if i == tmp + 4:
                            tmp_line = str(line).split(sep="=")
                            blub = tmp_line[-1]
                            tau_list.append(blub)
                        if i == tmp + 5:
                            tmp_line = str(line).split(sep="=")
                            blub = tmp_line[-1]
                            PI_list.append(blub)
                        
                print(data,model,add_info,score,"Done")


data = {
    "Datatype": data_list,
    "Modeltype": model_list,
    "Featuretype": add_info_list,
    "Scoretype": score_list,
    "Spearman correlation coefficient (SP)": SP_list,
    "Scipy Spearman correlation coefficent (SP)": SPscipy_list,
    "Confidence interval of 90% Spearman SP": conf_list,
    "Kendall correlation coefficient (tau)": tau_list,
    "Predictive index (PI)": PI_list,
}
df = pd.DataFrame(data)
df.to_csv("results_power_ranking.csv")
