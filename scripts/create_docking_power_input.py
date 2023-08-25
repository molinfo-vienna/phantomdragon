import pandas as pd
import os
import phantomdragon.functions as ph

df = pd.read_csv("../data/PDBbind_core_set_all.csv")

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest","XGBoost"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]

for PDBcode in df["PDB code"]:
    os.system(f"python GAP/calc_descr_pdb_ligands.py -p ../data/CASF-2016/coreset/{PDBcode}/{PDBcode}_protein.pdb -l ../data/CASF-2016/decoys_docking/{PDBcode}_decoys.mol2 -o ../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv")

    for k in datatypes:
        for modeltype in modeltypes:
            for score in scoretypes:
                blub = ph.parameterCollector(modeltype=modeltype,scoretype=score,add_information="GAP")
                blub.set_datatype(f"{k}")
                PDBs, phantomscore = blub.phantomscore(f"../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv","../models/",identifier="Ligand")
                data = {"#code":PDBs,"score":phantomscore}
                df = pd.DataFrame(data)
                
                if "/" in score:
                    score = score.replace("/","div")
                if " " in score:
                    score = score.replace(" ","_")

                df.to_csv(f"../data/CASF-2016/power_docking/GAP_Scores/{PDBcode}_{k}{modeltype}GAP{score}.csv",index=False)

                os.system(f"mkdir ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}GAP{score}")
                os.system(f"mv ../data/CASF-2016/power_docking/GAP_Scores/*{k}{modeltype}GAP{score}* ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}GAP{score}/.")

                os.system(f"mv ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}GAP{score}/{PDBcode}_{k}{modeltype}GAP{score}.csv ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}GAP{score}/{PDBcode}_score.dat")

                if "div" in score:
                    score = score.replace("div","/")
                if "_" in score:
                    score = score.replace("_"," ")
                
                print(k,modeltype,score,"done")