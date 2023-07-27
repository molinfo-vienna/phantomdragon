import pandas as pd
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]


for k in datatypes:
    for modeltype in modeltypes:
        for score in scoretypes:
            blub = ph.parameterCollector(modeltype=modeltype,scoretype=score,add_information="final")
            blub.set_datatype(f"{k}")
            PDBs, phantomscore = blub.phantomscore("../data/PDBbind_core_set_all_grail_scores.csv","../models/")
            data = {"#code":PDBs,"score":phantomscore}
            df = pd.DataFrame(data)
            
            if "/" in score:
                score = score.replace("/","div")
            if " " in score:
                score = score.replace(" ","_")

            df.to_csv(f"../results/CASF_input/{k}{modeltype}final{score}.csv",index=False)

            if "div" in score:
                score = score.replace("div","/")
            if "_" in score:
                score = score.replace("_"," ")
            
            print(k,modeltype,score,"done")