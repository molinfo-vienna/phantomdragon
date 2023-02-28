import pandas as pd
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]

modeltype_list = []
scoretype_list = []
datatype_list = []
mse_list = []
pear_list = []
coef_list = []
add_info_list = []

for k in datatypes:
    for modeltype in modeltypes:
        for score in scoretypes:
            blub = ph.parameterCollector(add_information="final",modeltype=modeltype,scoretype=score)
            x_train,x_test,y_train,y_test = ph.prepare_data(score,"../data/grail_scores_final.csv","../data/CASF_grail_scores.csv",f"../data/PDBbind_refined_set_{k}.csv","final")
            blub.set_trainingdata(x_train,y_train)
            blub.set_testingdata(x_test,y_test)
            blub.set_datatype(f"{k}")
            #blub.train_and_save_model(savepath="../models/")
            blub.phantomtest(loadpath="../models/")
            #blub.plot_phantomtest("../plots/")
            modeltype,scoret,datatype,mse,pearsonr,r_2,add_info = blub.get_stats()
            modeltype_list.append(modeltype)
            scoretype_list.append(scoret)
            datatype_list.append(datatype)
            mse_list.append(mse)
            pear_list.append(pearsonr)
            coef_list.append(r_2)
            add_info_list.append(add_info)
            print(k,modeltype,score,"done")

#print(len(modeltype_list),len(featuretype_list),len(datatype_list),len(mse_list),len(pear_list),len(coef_list),len(add_info_list))
data = {'Modeltype':modeltype_list,'Scoretype':scoretype_list,'Datatype':datatype_list,'Mean squared error (mse)':mse_list,'Pearson correlation coefficient (r)':pear_list,'Coefficient of determination (r²)':coef_list,'add. information':add_info_list}
df = pd.DataFrame(data)
df.to_csv("../results/CASF_all_results.csv")