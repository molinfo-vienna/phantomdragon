import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model
from sklearn import svm
from sklearn import tree
from sklearn import ensemble
from sklearn import preprocessing
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.utils import shuffle
from scipy.stats import pearsonr

def filter_and_sort_scores(exp_data, feature_scores):
    df = feature_scores[feature_scores["PDB code"].isin(exp_data["PDB code"])]
    df = df.sort_values("PDB code")
    df = df.reset_index()
    df = df.drop("index", axis=1)
    return df

def combine(list1,list2,list3):
    combinations = []
    for lst1 in list1:
        for lst2 in list2:
            for lst3 in list3:
                combinations.append([lst1,lst2,lst3])
    return combinations
