import pandas as pd
import phantomdragon.functions as ph

datatypes = ["all", "ki", "kd"]
# modeltypes = [
#     "linearRegression",
#     "Ridge",
#     "Lasso",
#     "ElasticNet",
#     "SVR",
#     "DecisionTree",
#     "RandomForest",
#     "XGBoost",
# ]
modeltypes = ["RandomForest"]
scoretypes = ["pKd pKi pIC50"]
classes = ["1", "2", "3", "4", "5", "6", "7"]

modeltype_list = []
scoretype_list = []
datatype_list = []
mae_list = []
mse_list = []
sd_list = []
pear_list = []
spear_list = []
coef_list = []
confidence_interval_list = []
add_info_list = []
classes_list = []
setlist = []
setsizes = []

for c in classes:
    for k in datatypes:
        for modeltype in modeltypes:
            for score in scoretypes:
                blub = ph.parameterCollector(
                    add_information=f"X-GRADE_class{c}_test", modeltype=modeltype, scoretype=score
                )
                x_train, x_test, y_train, y_test = ph.prepare_data(
                    score,
                    "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_refined_set_grail_scores_GAP.csv",
                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_GAP_class{c}.csv",
                    "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/refined_set.csv",
                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/general_set_class{c}.csv",
                    f"X-GRADE_class{c}_test",
                )
                if len(x_train) > 1 and len(y_train) > 1 and len(x_test) > 1 and len(y_test) > 1:
                    blub.set_trainingdata(x_train, y_train)
                    blub.set_testingdata(x_test, y_test)
                    blub.set_datatype(f"{k}")
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
                        spearman_r,
                        add_info,
                    ) = blub.get_stats(spearman=True)
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
                    classes_list.append(c)
                    setlist.append("test")
                    setsizes.append(len(x_test))
                    print(k, modeltype, score, "done (X-GRADE)","Class",c,"Test")
                    print("Training Set size", len(x_train))
                    print("Testing Set size", len(x_test))
                else:
                    print(f"No data for class {c} Test (X-GRADE)")
############################################################################################
                # blub = ph.parameterCollector(
                #     add_information=f"X-GRADE_class{c}_val", modeltype=modeltype, scoretype=score
                # )
                # x_train, x_test, y_train, y_test = ph.prepare_data(
                #     score,
                #     "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_refined_set_grail_scores_GAP.csv",
                #     f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_GAP_class{c}.csv",
                #     "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/refined_set.csv",
                #     f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/core_set_class{c}.csv",
                #     f"X-GRADE_class{c}_val",
                # )
                # if len(x_train) > 1 and len(y_train) > 1 and len(x_test) > 1 and len(y_test) > 1:
                #     blub.set_trainingdata(x_train, y_train)
                #     blub.set_testingdata(x_test, y_test)
                #     blub.set_datatype(f"{k}")
                #     blub.train_and_save_model(savepath="../models/")
                #     blub.phantomtest(loadpath="../models/")
                #     blub.plot_phantomtest("../plots/")
                #     (
                #         modeltype,
                #         scoret,
                #         datatype,
                #         mae,
                #         mse,
                #         sd,
                #         pearsonr,
                #         confidence_interval,
                #         r_2,
                #         spearman_r,
                #         add_info,
                #     ) = blub.get_stats(spearman=True)
                #     modeltype_list.append(modeltype)
                #     scoretype_list.append(scoret)
                #     datatype_list.append(datatype)
                #     mae_list.append(mae)
                #     mse_list.append(mse)
                #     sd_list.append(sd)
                #     pear_list.append(pearsonr)
                #     confidence_interval_list.append(confidence_interval)
                #     coef_list.append(r_2)
                #     spear_list.append(spearman_r)
                #     add_info_list.append(add_info)
                #     classes_list.append(c)
                #     setlist.append("val")
                #     setsizes.append(len(x_test))
                #     print(k, modeltype, score, "done (X-GRADE)","Class",c,"Validation")
                #     print("Training Set size", len(x_train))
                #     print("Testing Set size", len(x_test))
                # else:
                #     print(f"No data for class {c} Validation (X-GRADE)")

    for k in datatypes:
        for modeltype in modeltypes:
            for score in scoretypes:
                blub = ph.parameterCollector(
                    add_information=f"GRADE_class{c}_test", modeltype=modeltype, scoretype=score
                )
                x_train, x_test, y_train, y_test = ph.prepare_data(
                    score,
                    "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_refined_set_grail_scores_slim.csv",
                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_slim_class{c}.csv",
                    "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/refined_set.csv",
                    f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/general_set_class{c}.csv",
                    f"GRADE_class{c}_test",
                )
                if len(x_train) > 1 and len(y_train) > 1 and len(x_test) > 1 and len(y_test) > 1:
                    blub.set_trainingdata(x_train, y_train)
                    blub.set_testingdata(x_test, y_test)
                    blub.set_datatype(f"{k}")
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
                        spearman_r,
                        add_info,
                    ) = blub.get_stats(spearman=True)
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
                    classes_list.append(c)
                    setlist.append("test")
                    setsizes.append(len(x_test))
                    print(k, modeltype, score, "done (GRADE)","Class",c,"Test")
                    print("Training Set size", len(x_train))
                    print("Testing Set size", len(x_test))
                else:
                    print(f"No data for class {c} Test (GRADE)")
############################################################################################
                # blub = ph.parameterCollector(
                #     add_information=f"GRADE_class{c}_val", modeltype=modeltype, scoretype=score
                # )
                # x_train, x_test, y_train, y_test = ph.prepare_data(
                #     score,
                #     "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_refined_set_grail_scores_slim.csv",
                #     f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/GRAIL-scores/PDBbind_general_set_grail_scores_slim_class{c}.csv",
                #     "/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/refined_set.csv",
                #     f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/experiments/core_set_class{c}.csv",
                #     f"GRADE_class{c}_val",
                # )
                # if len(x_train) > 1 and len(y_train) > 1 and len(x_test) > 1 and len(y_test) > 1:
                #     blub.set_trainingdata(x_train, y_train)
                #     blub.set_testingdata(x_test, y_test)
                #     blub.set_datatype(f"{k}")
                #     blub.train_and_save_model(savepath="../models/")
                #     blub.phantomtest(loadpath="../models/")
                #     blub.plot_phantomtest("../plots/")
                #     (
                #         modeltype,
                #         scoret,
                #         datatype,
                #         mae,
                #         mse,
                #         sd,
                #         pearsonr,
                #         confidence_interval,
                #         r_2,
                #         spearman_r,
                #         add_info,
                #     ) = blub.get_stats(spearman=True)
                #     modeltype_list.append(modeltype)
                #     scoretype_list.append(scoret)
                #     datatype_list.append(datatype)
                #     mae_list.append(mae)
                #     mse_list.append(mse)
                #     sd_list.append(sd)
                #     pear_list.append(pearsonr)
                #     confidence_interval_list.append(confidence_interval)
                #     coef_list.append(r_2)
                #     spear_list.append(spearman_r)
                #     add_info_list.append(add_info)
                #     classes_list.append(c)
                #     setlist.append("val")
                #     setsizes.append(len(x_test))
                #     print(k, modeltype, score, "done (GRADE)","Class",c,"Validation")
                #     print("Training Set size", len(x_train))
                #     print("Testing Set size", len(x_test))
                # else:
                #     print(f"No data for class {c} Validation (GRADE)")


print(
    len(modeltype_list),
    len(scoretype_list),
    len(datatype_list),
    len(mae_list),
    len(mse_list),
    len(sd_list),
    len(pear_list),
    len(confidence_interval_list),
    len(coef_list),
    len(spear_list),
    len(add_info_list),
    len(classes_list),
    len(setlist),
)
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
    "Spearman correlation coefficient": spear_list,
    "add. information": add_info_list,
    "Classes": classes_list,
    "Set": setlist,
    "Set size": setsizes,
}
df = pd.DataFrame(data)
df.to_csv("../results/classes_results_tmp.csv")
