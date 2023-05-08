# tmp = open('../data/CASF_grail_scores.csv', 'r')
# lines = tmp.readlines()
# numbers = ["0","1","2","3","4","5","6","7","8","9"]

# for line in lines:
#     test = str(line).split(sep=",")
#     if "e" == test[0][1] and test[0][2] in numbers and test[0][3] in numbers:
#         print(test[0])

# tmp = open('../data/PDBbind_refined_set_kd.csv', 'r')
# lines = tmp.readlines()
# numbers = ["0","1","2","3","4","5","6","7","8","9"]

# for line in lines:
#     test = str(line).split(sep=",")
#     if "e" == test[1][1] and test[1][2] in numbers and test[1][3] in numbers:
#         print(test[1])

tmp = open('../data/CASF_grail_scores.csv', 'r')
lines = tmp.readlines()

coreset = []

for line in lines:
    test = str(line).split(sep=",")
    if test[0] == "PDB code":
        pass
    else:
        #print(test[0])
        coreset.append(test[0])

tmp = open('../data/PDBbind_refined_set_all.csv', 'r')
lines = tmp.readlines()

refinedset = []

for line in lines:
    test = str(line).split(sep=",")
    if test[1] == "PDB code":
        pass
    else:
        #print(test[1])
        refinedset.append(test[1])

i = 0
for complex in coreset:
    if complex in refinedset:
        i += 1
    
print('Number of Complexex unique to Coreset',len(coreset)-i)

print(len(coreset))
print(len(refinedset))