import os
import pandas as pd
import numpy as np
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest","XGBoost"]
scoretypes = ["delta G"]

modeltype_list = []
scoretype_list = []
datatype_list = []
mae_list = []
mse_list = []
sd_list = []
pear_list = []
spear_list = []
coef_list = []
confidence_interval_list =[]
add_info_list = []
system_list = []


for system in os.listdir("/data/shared/datasets/PL-REX/"):
    if os.path.isdir(f"/data/shared/datasets/PL-REX/{system}") and not system.startswith("."):
        print("-----------------------------------")
        print(f"Starting {system}")
        print("-----------------------------------")
        for k in datatypes:
            for modeltype in modeltypes:
                for score in scoretypes:

                    blub = ph.parameterCollector(add_information=f"X-GRADE",modeltype=modeltype,scoretype=score)
                    x_train,x_test,y_train,y_test = ph.prepare_data(score,"GAP/data/ref_set_grail_descr.csv","/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_GAP.csv",f"../data/PDBbind_refined_set_{k}.csv",f"../data/experiments/general_set.csv","X-GRADE")
                    blub.set_trainingdata(x_train,y_train)
                    
                    df_exp = pd.read_csv(f"/data/shared/datasets/PL-REX/{system}/experimental_dG.txt", delim_whitespace=True, comment='#', header=None, names=['ID', 'BindingFreeEnergy'])
                    df_exp = df_exp.sort_values(by='ID')
                    df_exp.reset_index(drop=True, inplace=True)
                    
                    df_desc = pd.read_csv(f"/data/shared/datasets/PL-REX/{system}/structures_pl-rex/X-GRADE_charged.csv")
                    df_desc = df_desc.sort_values(by='PDB code')
                    df_desc.reset_index(drop=True, inplace=True)
                    
                    df_exp = df_exp.drop("ID", axis=1)
                    df_desc = df_desc.drop("PDB code", axis=1)
                    
                    y_test = np.array(df_exp)
                    y_test = y_test.reshape(-1)
                    x_test = np.array(df_desc)
                    
                    blub.set_testingdata(x_test,y_test)
                    blub.set_datatype(f"{k}")
                    # blub.train_and_save_model(savepath="../models/")
                    blub.phantomtest(loadpath="../models/")
                    modeltype,scoret,datatype,mae,mse, sd,pearsonr,confidence_interval,r_2,spearman_r,add_info = blub.get_stats(spearman=True)
                    blub.plot_phantomtest("../plots/",name=f"{system}_{modeltype}_{scoret}_{datatype}_{add_info}_charged")
                    modeltype_list.append(modeltype)
                    scoretype_list.append(scoret)
                    datatype_list.append(datatype)
                    mae_list.append(mae)
                    mse_list.append(mse)
                    sd_list.append(sd)
                    pear_list.append(pearsonr)
                    confidence_interval_list.append(confidence_interval)
                    coef_list.append(r_2)
                    spear_list.append(spearman_r)
                    add_info_list.append(add_info)
                    system_list.append(system)
                    print(k,modeltype,score,"done (X-GRADE)")
                    print("Training Set size",len(x_train))
                    print("Testing Set size",len(x_test))

        for k in datatypes:
            for modeltype in modeltypes:
                for score in scoretypes:

                    blub = ph.parameterCollector(add_information="GRADE",modeltype=modeltype,scoretype=score)
                    x_train,x_test,y_train,y_test = ph.prepare_data(score,"../../GAP/data/ref_set_descrs.csv","/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_slim.csv",f"../data/PDBbind_refined_set_{k}.csv",f"../data/experiments/general_set.csv","GRADE")
                    blub.set_trainingdata(x_train,y_train)
                    
                    df_exp = pd.read_csv(f"/data/shared/datasets/PL-REX/{system}/experimental_dG.txt", delim_whitespace=True, comment='#', header=None, names=['ID', 'BindingFreeEnergy'])
                    df_exp = df_exp.sort_values(by='ID')
                    df_exp.reset_index(drop=True, inplace=True)
                    
                    df_desc = pd.read_csv(f"/data/shared/datasets/PL-REX/{system}/structures_pl-rex/GRADE_charged.csv")
                    df_desc = df_desc.sort_values(by='PDB code')
                    df_desc.reset_index(drop=True, inplace=True)
                    
                    df_exp = df_exp.drop("ID", axis=1)
                    df_desc = df_desc.drop("PDB code", axis=1)
                    
                    y_test = np.array(df_exp)
                    y_test = y_test.reshape(-1)
                    x_test = np.array(df_desc)
                    
                    blub.set_testingdata(x_test,y_test)
                    blub.set_datatype(f"{k}")
                    # blub.train_and_save_model(savepath="../models/")
                    blub.phantomtest(loadpath="../models/")
                    modeltype,scoret,datatype,mae,mse,sd,pearsonr,confidence_interval,r_2,spearman_r,add_info = blub.get_stats(spearman=True)
                    blub.plot_phantomtest("../plots/",name=f"{system}_{modeltype}_{scoret}_{datatype}_{add_info}_charged")
                    modeltype_list.append(modeltype)
                    scoretype_list.append(scoret)
                    datatype_list.append(datatype)
                    mae_list.append(mae)
                    mse_list.append(mse)
                    sd_list.append(sd)
                    pear_list.append(pearsonr)
                    confidence_interval_list.append(confidence_interval)
                    coef_list.append(r_2)
                    spear_list.append(spearman_r)
                    add_info_list.append(add_info)
                    system_list.append(system)
                    print(k,modeltype,score,"done (GRADE)")
                    print("Training Set size",len(x_train))
                    print("Testing Set size",len(x_test))
    else:
        continue

print(len(modeltype_list),len(scoretype_list),len(datatype_list),len(mae_list),len(mse_list),len(sd_list),len(pear_list),len(confidence_interval_list),len(coef_list),len(spear_list),len(add_info_list))

data = {'Modeltype':modeltype_list,'Scoretype':scoretype_list,'Datatype':datatype_list,'Mean absolute error (mae)':mae_list,'Mean squared error (mse)':mse_list,'Standard Diviation (SD)':sd_list,'Pearson correlation coefficient (r)':pear_list,'90% Confidence interval':confidence_interval_list,'Coefficient of determination (r²)':coef_list,'Spearman correlation coefficient':spear_list,'add. information':add_info_list,'System':system_list}
df = pd.DataFrame(data)
df.to_csv("../results/PL-REX_results_charged_tmp.csv")