import pandas as pd
import phantomdragon.functions as ph

final = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_all_grail_scores.csv")
GAP = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/ref_set_grail_descr.csv")
slim = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/ref_set_descrs.csv")

all = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_all.csv")
ki = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_ki.csv")
kd = pd.read_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_kd.csv")

final.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_all_grail_scores.csv",index=False)
final_ki = ph.filter_and_sort_features(ki,final)
print(final_ki)
final_ki.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_ki_grail_scores.csv",index=False)
final_kd = ph.filter_and_sort_features(kd,final)
print(final_kd)
final_kd.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/data/PDBbind_refined_set_kd_grail_scores.csv",index=False)

GAP.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/ref_set_grail_descr_all.csv",index=False)
GAP_ki = ph.filter_and_sort_features(ki,GAP)
print(GAP_ki)
GAP_ki.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/ref_set_grail_descr_ki.csv",index=False)
GAP_kd = ph.filter_and_sort_features(kd,GAP)
print(GAP_kd)
GAP_kd.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/phantomdragon/scripts/GAP/data/ref_set_grail_descr_kd.csv",index=False)

slim.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/ref_set_descrs_all.csv",index=False)
slim_ki = ph.filter_and_sort_features(ki,slim)
print(slim_ki)
slim_ki.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/ref_set_descrs_ki.csv",index=False)
slim_kd = ph.filter_and_sort_features(kd,slim)
print(slim_kd)
slim_kd.to_csv("/data/shared/projects/pharmacophore_hot_spot_analysis/GAP/data/ref_set_descrs_kd.csv",index=False)

print("DONE")