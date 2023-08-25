import pandas as pd

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest"]

# additional_information = ["basic", "-w", "basic-el", "-w-el", "basic-vdw", "-w-vdw", "basic-el-vdw", "-w-el-vdw","final"]
additional_information = ["final","GAP"]
scoretypes = ["delta_G","Affinity_Data_Value","pKd_pKi_pIC50"]

regression_eq_list = []
N_list = []
R_list = []
SD_list = []
data_list = []
model_list = []
add_info_list = []
score_list = []
confidence_list = []

for data in datatypes:
    for model in modeltypes:
        for add_info in additional_information:
            for score in scoretypes:
                output = open(f'results_power_scoring/{data}{model}{add_info}{score}.out', 'r')
                lines = output.readlines()
                for line in lines:
                    if line[:29] == "Summary of the scoring power:":
                        data_list.append(data)
                        model_list.append(model)
                        add_info_list.append(add_info)
                        score_list.append(score)
                    elif line[:24] == "The regression equation:":
                        equ,val1 = str(line).split(sep="=")
                        regression_eq_list.append(val1[1:-1])
                    elif line[:30] == "Number of favorable sample (N)":
                        equ,val2 = str(line).split(sep="=")
                        N_list.append(val2[1:-1])
                    elif line[:35] == "Pearson correlation coefficient (R)":
                        equ,val3 = str(line).split(sep="=")
                        R_list.append(val3[1:-1])
                    elif line[:10] == "Confidence":
                        equ, val5 = str(line).split(sep="=")
                        confidence_list.append(val5[1:-1])
                    elif line[:34] == "Standard deviation in fitting (SD)":
                        equ,val4 = str(line).split(sep="=")
                        SD_list.append(val4[1:-1])
                print(f"{data} {model} {add_info} {score} done!")

                    
data = {'Datatype':data_list,'Modeltype':model_list,'Featuretype':add_info_list,'Scoretype':score_list,'The regression equation: logKa =':regression_eq_list,'Number of favorable sample (N):':N_list,'Pearson correlation coefficient (R)':R_list,'Confidence interval of 90% Pearson R':confidence_list,'Standard deviation in fitting (SD)':SD_list}
df = pd.DataFrame(data)
df.to_csv("results_power_scoring.csv")