import pandas as pd
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest","XGBoost"]
scoretypes = ["delta G","pKd pKi pIC50"]

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

for k in datatypes:
    for modeltype in modeltypes:
        for score in scoretypes:

            blub = ph.parameterCollector(add_information="X-GRADE_2021",modeltype=modeltype,scoretype=score)
            x_train,x_test,y_train,y_test = ph.prepare_data(score,"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_2021/X-GRADE_1981_2010.csv","/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_2021/X-GRADE_CASF2016.csv",f"../data/PDBbind_2021_1981_2010_{k}.csv",f"../data/PDBbind_2021_CASF2016_all.csv","X-GRADE")
            blub.set_trainingdata(x_train,y_train)
            blub.set_testingdata(x_test,y_test)
            blub.set_datatype(f"{k}")
            blub.train_and_save_model(savepath="../models/")
            blub.phantomtest(loadpath="../models/")
            blub.plot_phantomtest("../plots/")
            modeltype,scoret,datatype,mae,mse, sd,pearsonr,confidence_interval,r_2,spearman_r,add_info = blub.get_stats(spearman=True)
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
            print(k,modeltype,score,"done (X-GRADE)")
            print("Training Set size",len(x_train))
            print("Testing Set size",len(x_test))

for k in datatypes:
    for modeltype in modeltypes:
        for score in scoretypes:

            blub = ph.parameterCollector(add_information="GRADE_2021",modeltype=modeltype,scoretype=score)
            x_train,x_test,y_train,y_test = ph.prepare_data(score,"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_2021/GRADE_1981_2010.csv","/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_2021/GRADE_CASF2016.csv",f"../data/PDBbind_2021_1981_2010_{k}.csv",f"../data/PDBbind_2021_CASF2016_all.csv","GRADE")
            blub.set_trainingdata(x_train,y_train)
            blub.set_testingdata(x_test,y_test)
            blub.set_datatype(f"{k}")
            blub.train_and_save_model(savepath="../models/")
            blub.phantomtest(loadpath="../models/")
            blub.plot_phantomtest("../plots/")
            modeltype,scoret,datatype,mae,mse,sd,pearsonr,confidence_interval,r_2,spearman_r,add_info = blub.get_stats(spearman=True)
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
            print(k,modeltype,score,"done (GRADE)")
            print("Training Set size",len(x_train))
            print("Testing Set size",len(x_test))

print(len(modeltype_list),len(scoretype_list),len(datatype_list),len(mae_list),len(mse_list),len(sd_list),len(pear_list),len(confidence_interval_list),len(coef_list),len(spear_list),len(add_info_list))
data = {'Modeltype':modeltype_list,'Scoretype':scoretype_list,'Datatype':datatype_list,'Mean absolute error (mae)':mae_list,'Mean squared error (mse)':mse_list,'Standard Diviation (SD)':sd_list,'Pearson correlation coefficient (r)':pear_list,'90% Confidence interval':confidence_interval_list,'Coefficient of determination (r²)':coef_list,'Spearman correlation coefficient':spear_list,'add. information':add_info_list}
df = pd.DataFrame(data)
df.to_csv("../results/validation_results_2021.csv")