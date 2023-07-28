import pandas as pd
import os
import phantomdragon.functions as ph

df = pd.read_csv("../data/PDBbind_core_set_all.csv")

datatypes = ["all","ki","kd"]
modeltypes = ["linearRegression","Ridge","Lasso","ElasticNet","SVR","DecisionTree","RandomForest"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]

for PDBcode in df["PDB code"]:
    os.system(f"python ../../GAP/calc_descr_pdb_ligands.py -p ../data/CASF-2016/coreset/{PDBcode}/{PDBcode}_protein.pdb -l ../data/CASF-2016/decoys_docking/{PDBcode}_decoys.mol2 -o ../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv")
    
    # just to make it comparable for now (removing TPSA)
    ############################################################################################
    tmp = pd.read_csv(f"../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv")
    tmp = tmp.drop("TPSA",axis=1)
    tmp.to_csv(f"../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv",index=False)
    ############################################################################################

    for k in datatypes:
        for modeltype in modeltypes:
            for score in scoretypes:
                blub = ph.parameterCollector(modeltype=modeltype,scoretype=score,add_information="final")
                blub.set_datatype(f"{k}")
                PDBs, phantomscore = blub.phantomscore(f"../data/CASF-2016/power_docking/GAP_descriptors/{PDBcode}_grail_scores.csv","../models/",identifier="Ligand")
                data = {"#code":PDBs,"score":phantomscore}
                df = pd.DataFrame(data)
                
                if "/" in score:
                    score = score.replace("/","div")
                if " " in score:
                    score = score.replace(" ","_")

                df.to_csv(f"../data/CASF-2016/power_docking/GAP_Scores/{PDBcode}_{k}{modeltype}final{score}.csv",index=False)

                os.system(f"mkdir ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}final{score}")
                os.system(f"mv ../data/CASF-2016/power_docking/GAP_Scores/*{k}{modeltype}final{score}* ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}final{score}/.")

                os.system(f"mv ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}final{score}/{PDBcode}_{k}{modeltype}final{score}.csv ../data/CASF-2016/power_docking/GAP_Scores/{k}{modeltype}final{score}/{PDBcode}_score.dat")

                if "div" in score:
                    score = score.replace("div","/")
                if "_" in score:
                    score = score.replace("_"," ")
                
                print(k,modeltype,score,"done")