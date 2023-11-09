import pandas as pd
import phantomdragon.functions as ph

modeltypes = [
    "linearRegression",
    "Ridge",
    "Lasso",
    "ElasticNet",
    "SVR",
    "DecisionTree",
    "RandomForest",
    "XGBoost",
]
featuretypes = ["final", "GAP", "slim"]
criterias = ["f_regression", "r_regression", "mutual_info_regression"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]
datatypes = ["all","ki","kd"]




modeltype_list = []
scoretype_list = []
datatype_list = []
mae_list = []
mse_list = []
sd_list = []
pear_list = []
coef_list = []
confidence_interval_list = []
add_info_list = []

for data in datatypes:
    for score in scoretypes:
        for i in range(5, 101, 5):
            for criteria in criterias:
                for featuretype in featuretypes:
                    if featuretype == "slim" and i > 36:
                        continue
                    else:
                        for modeltype in modeltypes:
                            try:
                                blub = ph.parameterCollector(
                                    add_information=f"{featuretype}_size{i}_{criteria}",
                                    modeltype=modeltype,
                                    scoretype=score,
                                )
                                x_train, x_test, y_train, y_test = ph.prepare_data(
                                    score,
                                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_{featuretype}_size{i}_{criteria}_{data}_{score}.csv",
                                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_{featuretype}_size{i}_{criteria}_{data}_{score}.csv",
                                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/PDBbind_ref_set_affinity_data.csv",
                                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/PDBbind_core_set_affinity_data.csv",
                                    f"{featuretype}_size{i}_{criteria}"
                                )
                                blub.set_trainingdata(x_train, y_train)
                                blub.set_testingdata(x_test, y_test)
                                blub.set_datatype(data)
                                blub.train_and_save_model(savepath="../models/")
                                blub.phantomtest(loadpath="../models/")
                                blub.plot_phantomtest("../plots/")
                                (
                                    modeltype,
                                    scoret,
                                    datatype,
                                    mae,
                                    mse,
                                    sd,
                                    pearsonr,
                                    confidence_interval,
                                    r_2,
                                    add_info,
                                ) = blub.get_stats()
                                modeltype_list.append(modeltype)
                                scoretype_list.append(scoret)
                                datatype_list.append(datatype)
                                mae_list.append(mae)
                                mse_list.append(mse)
                                sd_list.append(sd)
                                pear_list.append(pearsonr)
                                confidence_interval_list.append(confidence_interval)
                                coef_list.append(r_2)
                                add_info_list.append(add_info)
                                print(f"{data} {score} {modeltype} {featuretype} size{i} {criteria} Done")
                                print("Training Set size", len(x_train))
                                print("Testing Set size", len(x_test))
                            except:
                                print(f"Something went wrong with {data} {score} size{i} {criteria} {featuretype} {modeltype}")


        for featuretype in featuretypes:
            for i,modeltype in enumerate(modeltypes):
                if i == 0:
                    model = "linearRegression"
                elif i == 1:
                    model = "Ridge"
                elif i == 2:
                    model = "Lasso"
                elif i == 3:
                    model = "ElasticNet"
                elif i == 4:
                    model = "SVR"                    
                elif i == 5:
                    model = "DecisionTree"
                elif i == 6:
                    model = "RandomForest"
                elif i == 7:
                    model = "XGB"    
                try:
                    blub = ph.parameterCollector(
                        add_information=f"{featuretype}_RFECV_{modeltype}",
                        modeltype=modeltype,
                        scoretype=score,
                    )
                    x_train, x_test, y_train, y_test = ph.prepare_data(
                        score,
                        f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_{featuretype}_RFECV_{data}_{score}_{model}.csv",
                        f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_{featuretype}_RFECV_{data}_{score}_{model}.csv",
                        f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/PDBbind_ref_set_affinity_data.csv",
                        f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/PDBbind_core_set_affinity_data.csv",
                        f"{featuretype}_RFECV_{modeltype}"
                    )
                    blub.set_trainingdata(x_train, y_train)
                    blub.set_testingdata(x_test, y_test)
                    blub.set_datatype(data)
                    blub.train_and_save_model(savepath="../models/")
                    blub.phantomtest(loadpath="../models/")
                    blub.plot_phantomtest("../plots/")
                    (
                        modeltype,
                        scoret,
                        datatype,
                        mae,
                        mse,
                        sd,
                        pearsonr,
                        confidence_interval,
                        r_2,
                        add_info,
                    ) = blub.get_stats()
                    modeltype_list.append(modeltype)
                    scoretype_list.append(scoret)
                    datatype_list.append(datatype)
                    mae_list.append(mae)
                    mse_list.append(mse)
                    sd_list.append(sd)
                    pear_list.append(pearsonr)
                    confidence_interval_list.append(confidence_interval)
                    coef_list.append(r_2)
                    add_info_list.append(add_info)
                    print(f"{data} {score} {modeltype} {featuretype} RFECV Done")
                    print("Training Set size", len(x_train))
                    print("Testing Set size", len(x_test))
                except:
                    print(f"Something went wrong with RFECV {data} {score} {featuretype} {model}")

# print(len(modeltype_list),len(featuretype_list),len(datatype_list),len(mse_list),len(pear_list),len(coef_list),len(add_info_list))
data = {
    "Modeltype": modeltype_list,
    "Scoretype": scoretype_list,
    "Datatype": datatype_list,
    "Mean absolute error (mae)": mae_list,
    "Mean squared error (mse)": mse_list,
    "Standard Diviation (SD)": sd_list,
    "Pearson correlation coefficient (r)": pear_list,
    "90% Confidence interval": confidence_interval_list,
    "Coefficient of determination (r²)": coef_list,
    "add. information": add_info_list,
}
df = pd.DataFrame(data)
df.to_csv("../results/feature_selection_results2.csv")
