import pandas as pd
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression, r_regression, mutual_info_regression
from sklearn import preprocessing

datatypes = ["all","ki","kd"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]

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

        for i in range(5,171,5):
            selector = SelectKBest(f_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_final_size{i}_f_regression_{k}_{score}.csv")
            
            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_core_set_all_grail_scores.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_final_size{i}_f_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,171,5):
            selector = SelectKBest(r_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_final_size{i}_r_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_core_set_all_grail_scores.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_final_size{i}_r_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,171,5):
            selector = SelectKBest(mutual_info_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_final_size{i}_mutual_info_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_core_set_all_grail_scores.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_final_size{i}_mutual_info_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

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

        for i in range(5,176,5):
            selector = SelectKBest(f_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_GAP_size{i}_f_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/core_set_grail_descr.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_GAP_size{i}_f_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,176,5):
            selector = SelectKBest(r_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_GAP_size{i}_r_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/core_set_grail_descr.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_GAP_size{i}_r_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,176,5):
            selector = SelectKBest(mutual_info_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]

            if " " in score:
                score = score.replace(" ","_")

            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]
            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_GAP_size{i}_mutual_info_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/core_set_grail_descr.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_GAP_size{i}_mutual_info_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

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

        for i in range(5,36,5):
            selector = SelectKBest(f_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_slim_size{i}_f_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/core_set_descrs.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_slim_size{i}_f_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,36,5):
            selector = SelectKBest(r_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_slim_size{i}_r_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/core_set_descrs.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_slim_size{i}_r_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")

        for i in range(5,36,5):
            selector = SelectKBest(mutual_info_regression, k=i)
            selector.fit(X,y)
            selection = [identifier]
            for j in selector.get_feature_names_out(input_features=features_num.keys()):
                selection.append(j)
            sub_features = features[selection]

            if " " in score:
                score = score.replace(" ","_")

            sub_features.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/ref_set_slim_size{i}_mutual_info_regression_{k}_{score}.csv")

            core = pd.read_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/core_set_descrs.csv", dtype={identifier: str})
            core_new = core[selection]
            core_new.to_csv(f"/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/slim/core_set_slim_size{i}_mutual_info_regression_{k}_{score}.csv")

            if "_" in score:
                score = score.replace("_"," ")
