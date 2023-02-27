import pandas as pd
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest"]
additional_information = ["basic", "-w", "basic-el", "-w-el", "basic-vdw", "-w-vdw", "basic-el-vdw", "-w-el-vdw"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50","1/K"]


for k in datatypes:
    for modeltype in modeltypes:
        for info in additional_information:
            for score in scoretypes:

                blub = ph.parameterCollector(add_information=info,modeltype=modeltype,scoretype=score)
                blub.set_datatype(f"{k}")
                PDBs, phantomscore = blub.phantomscore("../data/CASF_grail_scores.csv","../models/")
                data = {"#code":PDBs,"score":phantomscore}
                df = pd.DataFrame(data)
                
                if "/" in score:
                    score = score.replace("/","div")
                if " " in score:
                    score = score.replace(" ","_")

                df.to_csv(f"../results/CASF_input/{k}{modeltype}{info}{score}.csv",index=False)

                if "div" in score:
                    score = score.replace("div","/")
                if "_" in score:
                    score = score.replace("_"," ")
                
                print(k,modeltype,info,score,"done")