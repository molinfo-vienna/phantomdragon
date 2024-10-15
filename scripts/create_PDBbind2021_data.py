import pandas as pd
import math

tmp = open('/data/shared/datasets/PDBbind_2021/index/INDEX_general_PL_data.2021', 'r')
lines = tmp.readlines()

PDB_code_list = []
Affinity_Data_Type_list = []
Affinity_Data_Value_list = []
pKi_pKd_IC50_list = []
Delta_G_list = []
K_reciprocal_list = []

for line in lines:
    test = str(line).split()
    print(test)
    if test[0] == "#":
        pass
    else:
        PDB_code_list.append(test[0])
        if "=" in test[4]:
            Affinity_Data_Type, Affinity_Data_Value = test[4].split(sep="=")
        elif "~" in test[4]:
            Affinity_Data_Type, Affinity_Data_Value = test[4].split(sep="~")
        elif ">" in test[4]:
            Affinity_Data_Type, Affinity_Data_Value = test[4].split(sep=">")
        elif "<" in test[4]:
            Affinity_Data_Type, Affinity_Data_Value = test[4].split(sep="<")
        else:
            raise ValueError("No split character found")
        
        Affinity_Data_Type_list.append(Affinity_Data_Type)
        
        if Affinity_Data_Value[-2:] == "nM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
        elif Affinity_Data_Value[-2:] == "uM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)*1000
            Affinity_Data_Value = Affinity_Data_Value
        elif Affinity_Data_Value[-2:] == "pM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)/1000
            Affinity_Data_Value = Affinity_Data_Value
        elif Affinity_Data_Value[-2:] == "mM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)*1000000
            Affinity_Data_Value = Affinity_Data_Value
        elif Affinity_Data_Value[-2:] == "fM":
            Affinity_Data_Value = Affinity_Data_Value[0:-2]
            Affinity_Data_Value = float(Affinity_Data_Value)/1000000
            Affinity_Data_Value = Affinity_Data_Value
        
        Affinity_Data_Value_list.append(Affinity_Data_Value)
        
        pKi_pKd_IC50_list.append(test[3])
        # print(Affinity_Data_Value)
        Delta_G = round(0.001987*298*math.log(float(Affinity_Data_Value)/1000000000),2)
        Delta_G_list.append(Delta_G)

        K_reciprocal = round(1/float(Affinity_Data_Value),5)
        K_reciprocal_list.append(K_reciprocal)

# print(len(PDB_code_list),len(Affinity_Data_Type_list),len(Affinity_Data_Value_list),len(pKi_pKd_IC50_list),len(Delta_G_list),len(K_reciprocal_list))
data = {'PDB code':PDB_code_list,'Affinity Data Type':Affinity_Data_Type_list,'Affinity Data Value':Affinity_Data_Value_list,'pKd pKi pIC50':pKi_pKd_IC50_list,'delta G':Delta_G_list,'1/K':K_reciprocal_list}
df = pd.DataFrame(data)
df.to_csv("../data/PDBbind_2021_all.csv")

df_kd = df[df['Affinity Data Type'] == "Kd"]
df_ki = df[df['Affinity Data Type'] == "Ki"]

df_kd.to_csv("../data/PDBbind_2021_kd.csv")
df_ki.to_csv("../data/PDBbind_2021_set_ki.csv")

df_refined = pd.read_csv("../data/PDBbind_2021/GRADE_1981_2010.csv")
df_core = pd.read_csv("../data/PDBbind_2021/GRADE_CASF2016.csv")
df_general = pd.read_csv("../data/PDBbind_2021/GRADE_2011_2020.csv")

df_refined_subset = df[df['PDB code'].isin(df_refined['PDB code'])]
df_core_subset = df[df['PDB code'].isin(df_core['PDB code'])]
df_general_subset = df[df['PDB code'].isin(df_general['PDB code'])]

df_refined_subset.to_csv("../data/PDBbind_2021_1981_2010_all.csv")
df_core_subset.to_csv("../data/PDBbind_2021_CASF2016_all.csv")
df_general_subset.to_csv("../data/PDBbind_2021_2011_2020_all.csv")

df_refined_subset_kd = df_kd[df_kd['PDB code'].isin(df_refined['PDB code'])]
df_core_subset_kd = df_kd[df_kd['PDB code'].isin(df_core['PDB code'])]
df_general_subset_kd = df_kd[df_kd['PDB code'].isin(df_general['PDB code'])]

df_refined_subset_kd.to_csv("../data/PDBbind_2021_1981_2010_kd.csv")
df_core_subset_kd.to_csv("../data/PDBbind_2021_CASF2016_kd.csv")
df_general_subset_kd.to_csv("../data/PDBbind_2021_2011_2020_kd.csv")

df_refined_subset_ki = df_ki[df_ki['PDB code'].isin(df_refined['PDB code'])]
df_core_subset_ki = df_ki[df_ki['PDB code'].isin(df_core['PDB code'])]
df_general_subset_ki = df_ki[df_ki['PDB code'].isin(df_general['PDB code'])]

df_refined_subset_ki.to_csv("../data/PDBbind_2021_1981_2010_ki.csv")
df_core_subset_ki.to_csv("../data/PDBbind_2021_CASF2016_ki.csv")
df_general_subset_ki.to_csv("../data/PDBbind_2021_2011_2020_ki.csv")
