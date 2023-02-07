import numpy as np
import pandas as pd
import phantomdragon.functions as ph
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

datatypes = ["ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet"]
additional_information = ["basic", "-w", "basic-el", "-w-el", "basic-vdw", "-w-vdw", "basic-el-vdw", "-w-el-vdw"]

experiment = pd.read_csv(f"data/PDBbind_refined_set_ki.csv")
experiment = ph.filter_and_sort_scores(experiment, experiment)

scores = pd.read_csv("data/grail_scores.csv")
scores_basic = ph.filter_and_sort_scores(experiment, scores)
scores_basic = scores_basic.sort_values("PDB code")
scores_basic = scores_basic.drop(["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"], axis=1)
scores_basic = scores_basic.to_numpy()

experiment = np.array(experiment["pKd pKi pIC50"])

scaler = preprocessing.StandardScaler().fit(scores_basic)
scores_basic = scaler.transform(scores_basic)
x_train, x_test, y_train, y_test = train_test_split(scores_basic, experiment, test_size = 0.2)
        

blub = ph.parameterCollection(add_information="basic",modeltype="linearRegression",featuretype="Affinity Data Value")
blub.load_all_data("ki",score_train=x_train,score_test=x_test,experiment_train=y_train,experiment_test=y_test,path=False)
blub.train_and_save_model()
blub.phantomtest()
q,w,e,r,t,z,u = blub.extract_stats()
print(q,w,e,r,t,z,u)