import pandas as pd
import os
import phantomdragon.functions as ph

df = pd.read_csv("../data/PDBbind_core_set_all_grail_scores.csv")

datatypes = ["all", "ki", "kd"]
modeltypes = [
    "linearRegression",
    "Ridge",
    "Lasso",
    "ElasticNet",
    "SVR",
    "DecisionTree",
    "RandomForest",
]
scoretypes = ["delta G", "Affinity Data Value", "pKd pKi pIC50"]

for PDBcode1 in next(os.walk("../data/CASF-2016/decoys_screening"))[1]:
    data = pd.DataFrame(columns=df.columns)
    for PDBcode2 in df["PDB code"]:
        if PDBcode1 == PDBcode2:
            continue
        else:
            os.system(
                f"python ../../GAP/calc_descr_pdb_ligands.py -p ../data/CASF-2016/coreset/{PDBcode1}/{PDBcode1}_protein.pdb -l ../data/CASF-2016/decoys_screening/{PDBcode1}/{PDBcode1}_{PDBcode2}.mol2 -o ../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_{PDBcode2}_grail_scores.csv"
            )

            # just to make it comparable for now (removing TPSA)
            ############################################################################################
            tmp = pd.read_csv(
                f"../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_{PDBcode2}_grail_scores.csv"
            )
            tmp = tmp.drop("TPSA", axis=1)
            tmp.to_csv(
                f"../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_{PDBcode2}_grail_scores.csv",
                index=False,
            )
            ############################################################################################

            data1 = pd.read_csv(
                f"../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_{PDBcode2}_grail_scores.csv"
            )
            data = pd.concat([data, data1])
    data.to_csv(
        f"../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_grail_scores.csv",
        index=False,
    )

print("------------------------------")
print("Creation of descriptor Done!!!")
print("------------------------------")

for PDBcode1 in next(os.walk("../data/CASF-2016/decoys_screening"))[1]:
    for k in datatypes:
        for modeltype in modeltypes:
            for score in scoretypes:
                blub = ph.parameterCollector(
                    modeltype=modeltype, scoretype=score, add_information="final"
                )
                blub.set_datatype(f"{k}")
                PDBs, phantomscore = blub.phantomscore(
                    f"../data/CASF-2016/power_screening/GAP_descriptors/{PDBcode1}_grail_scores.csv",
                    "../models/",
                    identifier="Ligand",
                )
                data = {"#code": PDBs, "score": phantomscore}
                df = pd.DataFrame(data)

                if "/" in score:
                    score = score.replace("/", "div")
                if " " in score:
                    score = score.replace(" ", "_")

                df.to_csv(
                    f"../data/CASF-2016/power_screening/GAP_Scores/{PDBcode1}_{k}{modeltype}final{score}.csv",
                    index=False,
                )

                os.system(
                    f"mkdir ../data/CASF-2016/power_screening/GAP_Scores/{k}{modeltype}final{score}"
                )
                os.system(
                    f"mv ../data/CASF-2016/power_screening/GAP_Scores/*{k}{modeltype}final{score}* ../data/CASF-2016/power_screening/GAP_Scores/{k}{modeltype}final{score}/."
                )

                os.system(
                    f"mv ../data/CASF-2016/power_screening/GAP_Scores/{k}{modeltype}final{score}/{PDBcode1}_{k}{modeltype}final{score}.csv ../data/CASF-2016/power_screening/GAP_Scores/{k}{modeltype}final{score}/{PDBcode1}_score.dat"
                )

                if "div" in score:
                    score = score.replace("div", "/")
                if "_" in score:
                    score = score.replace("_", " ")

                print(k, modeltype, score, "done")
