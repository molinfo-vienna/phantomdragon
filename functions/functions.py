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
import joblib

def filter_and_sort_scores(exp_data, feature_scores):
    '''
    Parameters
    ----------
    exp_data : TYPE PANDAS Dataframe
        A pandas Dataframe of the avavible experimental values.
    feature_scores : TYPE PANDAS Dataframe
        A pandas Dataframe of the avaviable scores. 
        It is assumed that exp_data is corresponding to a subset of feature_scores.

    Returns
    -------
    df : TYPE PANDAS Dataframe
        Gives back a sorted subset of feature_scores. 
        This subset corresponds to the subset of the experimental data.
        It is sorted and filtered according to the "PDB code" bracket.

    '''
    df = feature_scores[feature_scores["PDB code"].isin(exp_data["PDB code"])]
    df = df.sort_values("PDB code")
    df = df.reset_index()
    df = df.drop("index", axis=1)
    return df

def combine(list1,list2,list3):
    '''
    Parameters
    ----------
    list1 : TYPE List
    list2 : TYPE List
    list3 : TYPE List

    Returns
    -------
    combinations: TYPE List
        List that contains a list of all possible combinations between the three input lists
    '''
    combinations = []
    for lst1 in list1:
        for lst2 in list2:
            for lst3 in list3:
                combinations.append([lst1,lst2,lst3])
    return combinations

def load_scores(add_information = "basic"):
    '''
    Parameters
    ----------
    add_information : TYPE String
        String containing one of the following options: 
        basic, -w, basic -el, -w -el, basic -vdw, -w -vdw, basic -el -vdw, -w -el -vdw 

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    scores : TYPE
        DESCRIPTION.
    '''
    if add_information == "basic":
        drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
    elif add_information == "-w":
        drop = ["PDB code"," H-H_SUM"," H-H_MAX"," ES"," VDW_ATT"," VDW_REP"]
    elif add_information == "basic -el":
        drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," VDW_ATT"," VDW_REP"]
    elif add_information == "-w -el":
        drop = ["PDB code"," H-H_SUM"," H-H_MAX"," VDW_ATT"," VDW_REP"]
    elif add_information == "basic -el -vdw":
        drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"]
    elif add_information == "-w -el -vdw":
        drop = ["PDB code"," H-H_SUM"," H-H_MAX"]
    elif add_information == "basic -vdw":
        drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
    elif add_information == "-w -vdw":
        drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
    else:
        raise ValueError("Unexpected string in add_information")
        
    scores = filter_and_sort_scores()
    scores = scores.drop(drop, axis=1)
    scores = scores.to_numpy()
    return scores

def train_and_save_model(scores, features, testset_size = 0.2, modeltype = "linearRegression"):
    scaler = preprocessing.StandardScaler().fit(scores)
    scores = scaler.transform(scores)
    if modeltype == "linearRegression":
        reg = linear_model.LinearRegression()
    elif modeltype == "Ridge":
        reg = linear_model.RidgeCV()
    elif modeltype == "Lasso":
        reg = linear_model.LassoCV()
    elif modeltype == "ElasticNet":
        reg = linear_model.ElasticNetCV()
    elif modeltype == "SVR":
        reg = svm.SVR()
    elif modeltype == "DecisionTree":
        reg = tree.DecisionTreeRegressor()
    elif modeltype == "RandomForest":
        reg = ensemble.RandomForestRegressor()
    else:
        raise ValueError("Unexpected string in modeltype")
        
    x_train, x_test, y_train, y_test = train_test_split(scores, features, test_size = testset_size)
    
    if modeltype == "SVR":
        params = {'kernel': ['linear', 'poly', 'rbf', 'sigmoid'], 'C': [1,2.5,5], 'gamma': ['scale','auto'], 'degree' : [2,3,4], 'epsilon': [0.001, 0.01, 0.1, 1, 10, 100]}
        gs_reg = GridSearchCV(reg, params)
        gs_reg.fit(x_train, y_train)
        reg.set_params(**gs_reg.best_params_)
    elif modeltype == "DecisionTree":
        params = {'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'splitter':['best','random'],'max_features':['auto', 'sqrt','log2']}
        gs_reg = GridSearchCV(reg, params)
        gs_reg.fit(x_train, y_train)
        reg.set_params(**gs_reg.best_params_)
    elif modeltype == "RandomForest":
        params = {'n_estimators':[100,200,300],'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'max_features':['auto', 'sqrt','log2']}
        gs_reg = GridSearchCV(reg, params)
        gs_reg.fit(x_train, y_train)
        reg.set_params(**gs_reg.best_params_)
    
    reg.fit(x_train, y_train)
    joblib.dump(reg, f"../models/{modeltype}.sav")
    
def phantomprediction():
    

    
    