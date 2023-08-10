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
ligand1_list = []
succes1_list = []
ligand5_list = []
succes5_list = []
ligand10_list = []
succes10_list = []

'''
Summary of the reverse screening power: =========================================================
The best target is found among top 1% candidates for  5 ligand(s); success rate = 1.8%
The best target is found among top 5% candidates for 10 ligand(s); success rate = 3.5%
The best target is found among top 10% candidates for 22 ligand(s); success rate = 7.7%
=================================================================================================

Summary of the reverse screening power: =========================================================
The best target is found among top 1% candidates for  2 ligand(s); success rate = 0.7%
The best target is found among top 5% candidates for  5 ligand(s); success rate = 1.8%
The best target is found among top 10% candidates for 14 ligand(s); success rate = 4.9%
=================================================================================================
'''

for data in datatypes:
    for model in modeltypes:
        for add_info in additional_information:
            for score in scoretypes:
                output = open(
                    f"results_reverse_screening/{data}{model}{add_info}{score}.out", "r"
                )
                lines = output.readlines()
                for i, line in enumerate(lines):
                    if line[:39] == "Summary of the reverse screening power:":
                        data_list.append(data)
                        model_list.append(model)
                        add_info_list.append(add_info)
                        score_list.append(score)
                        tmp = i
                        start = True
                    elif start == True:
                        if i == tmp + 1:
                            tmp_line = str(line).split()
                            ligand1 = tmp_line[10]
                            ligand1_list.append(ligand1)
                            succes1 = tmp_line[-1]
                            succes1_list.append(succes1)
                        if i == tmp + 2:
                            tmp_line = str(line).split()
                            ligand5 = tmp_line[10]
                            ligand5_list.append(ligand5)
                            succes5 = tmp_line[-1]
                            succes5_list.append(succes5)
                        if i == tmp + 3:
                            tmp_line = str(line).split()
                            ligand10 = tmp_line[10]
                            ligand10_list.append(ligand10)
                            succes10 = tmp_line[-1]
                            succes10_list.append(succes10)
                print(data,model,add_info,score,"Done")

data = {
    "Datatype": data_list,
    "Modeltype": model_list,
    "Featuretype": add_info_list,
    "Scoretype": score_list,
    "1'%' ligand number:": ligand1_list,
    "1'%' success rate:": succes1_list,
    "5'%' ligand number:": ligand5_list,
    "5'%' success rate:": succes5_list,
    "10'%' ligand number:": ligand10_list,
    "10'%' success rate:": succes10_list,
}
df = pd.DataFrame(data)
df.to_csv("results_reverse_screening.csv")
