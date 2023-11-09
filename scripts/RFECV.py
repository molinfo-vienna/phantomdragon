import pandas as pd
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LinearRegression,RidgeCV,LassoCV,ElasticNetCV
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing
from xgboost import XGBRegressor

datatypes = ["all","ki","kd"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]
modeltypes = [LinearRegression(),RidgeCV(),LassoCV(),ElasticNetCV(),DecisionTreeRegressor(),RandomForestRegressor(),XGBRegressor()]

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
        model = "DecisionTree"
    elif i == 5:
        model = "RandomForest"
    elif i == 6:
        model = "XGB"    
    
    print(model)

    for k in datatypes:
        for score in scoretypes:

            featurepath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_{k}_grail_scores.csv" # final
            experimentpath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_{k}.csv"

            identifier = "PDB code"

            features = pd.read_csv(featurepath, dtype={identifier: str})
            experiment = pd.read_csv(experimentpath, dtype={identifier: str})
            experiment = experiment.sort_values(identifier)
            experiment = experiment.reset_index(drop=True)
            experiment = experiment[score]

            features_num = features.drop(labels=identifier,axis=1)
            scaler = preprocessing.StandardScaler().fit(features_num)

            X = scaler.transform(features_num)
            y = experiment

            estimator = modeltype
            selector = RFECV(estimator, step=1, cv=5,min_features_to_select=5)
            selector = selector.fit(X, y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]          

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_final_RFECV_{k}_{score}_{model}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_core_set_all_grail_scores.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_final_RFECV_{k}_{score}_{model}.csv")

            if "_" in score:
                score = score.replace("_"," ")

            print(f"{model} {k} {score} Done!")

    for k in datatypes:
        for score in scoretypes:
            featurepath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/ref_set_grail_descr_{k}.csv" # GAP
            experimentpath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_{k}.csv"

            identifier = "PDB code"

            features = pd.read_csv(featurepath, dtype={identifier: str})
            experiment = pd.read_csv(experimentpath, dtype={identifier: str})
            experiment = experiment.sort_values(identifier)
            experiment = experiment.reset_index(drop=True)
            experiment = experiment[score]

            features_num = features.drop(labels=identifier,axis=1)
            scaler = preprocessing.StandardScaler().fit(features_num)

            X = scaler.transform(features_num)
            y = experiment

            estimator = modeltype
            selector = RFECV(estimator, step=1, cv=5,min_features_to_select=5)
            selector = selector.fit(X, y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_GAP_RFECV_{k}_{score}_{model}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/core_set_grail_descr.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_GAP_RFECV_{k}_{score}_{model}.csv")

            if "_" in score:
                score = score.replace("_"," ")

            print(f"{model} {k} {score} Done!")

    for k in datatypes:
        for score in scoretypes:

            featurepath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/ref_set_descrs_{k}.csv" # slim
            experimentpath = f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_{k}.csv"

            identifier = "PDB code"

            features = pd.read_csv(featurepath, dtype={identifier: str})
            experiment = pd.read_csv(experimentpath, dtype={identifier: str})
            experiment = experiment.sort_values(identifier)
            experiment = experiment.reset_index(drop=True)
            experiment = experiment[score]

            features_num = features.drop(labels=identifier,axis=1)
            scaler = preprocessing.StandardScaler().fit(features_num)

            X = scaler.transform(features_num)
            y = experiment

            estimator = modeltype
            selector = RFECV(estimator, step=1, cv=5,min_features_to_select=5)
            selector = selector.fit(X, y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_slim_RFECV_{k}_{score}_{model}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/core_set_descrs.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_slim_RFECV_{k}_{score}_{model}.csv")

            if "_" in score:
                score = score.replace("_"," ")

            print(f"{model} {k} {score} Done!")
