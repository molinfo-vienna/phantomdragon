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
import pickle

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

class parameterCollection:
    def __init__(self, modeltype, add_information, featuretype, datatype = 0, scores_train = 0, scores_test = 0, features_train = 0, features_test = 0):
        self.modeltype = modeltype
        self.add_information = add_information
        self.scores_train = scores_train
        self.features_train = features_train
        self.scores_test = scores_test
        self.features_test = features_test
        self.featuretype = featuretype
        self.datatype = datatype
        self.mse = None
        self.r = None
        self.r_2 = None
    
    def load_all_data(self,datatype,score_train, experiment_train,score_test, experiment_test, path = True):
        self.load_trainingdata(datatype,score_train,experiment_train,path)
        self.load_testingdata(score_test,experiment_test,path)

    def load_trainingdata(self, datatype, scorepath, experimentpath, path = True):

        if path == False:
            self.scores_train = scorepath
            self.features_train = experimentpath
            return

        scores = pd.read_csv(scorepath)
        experiment = pd.read_csv(experimentpath)
        self.features_train =np.array(experiment[self.featuretype])

        if self.add_information == "basic":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "basic-el":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w-el":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "basic-el-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"]
        elif self.add_information == "-w -el -vdw":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"]
        elif self.add_information == "basic-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        else:
            raise ValueError("Unexpected string in add_information")
        

        scores = filter_and_sort_scores(experiment,scores)
        scores = scores.drop(drop, axis=1)
        scores = scores.to_numpy()
        self.scores_train = scores
        self.trainingscores_type = "All+"
        self.datatype = datatype
    
    def load_testingdata(self, scorepath, experimentpath, path = True):
        if path == False:
            self.scores_test = scorepath
            self.features_test = experimentpath
            return
        
        scores = pd.read_csv(scorepath)
        experiment = pd.read_csv(experimentpath)
        self.features_test =np.array(experiment[self.featuretype])

        if self.add_information == "basic":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "basic-el":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w-el":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "basic-el-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"]
        elif self.add_information == "-w -el -vdw":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"]
        elif self.add_information == "basic-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        else:
            raise ValueError("Unexpected string in add_information")
        

        scores = filter_and_sort_scores(experiment,scores)
        scores = scores.drop(drop, axis=1)
        scores = scores.to_numpy()
        self.scores_test = scores
    
    def train_and_save_model(self, training_scores = 0, training_features = 0):
        if training_scores != 0:
            self.scores_train = training_scores
        if training_features != 0:
            self.features_train = training_features
        self.scores_train, self.features_train = shuffle(self.scores_train, self.features_train)
        scaler = preprocessing.StandardScaler().fit(self.scores_train)
        self.scores_train = scaler.transform(self.scores_train)
        if self.modeltype == "linearRegression":
            reg = linear_model.LinearRegression()
        elif self.modeltype == "Ridge":
            reg = linear_model.RidgeCV()
        elif self.modeltype == "Lasso":
            reg = linear_model.LassoCV()
        elif self.modeltype == "ElasticNet":
            reg = linear_model.ElasticNetCV()
        elif self.modeltype == "SVR":
            reg = svm.SVR()
        elif self.modeltype == "DecisionTree":
            reg = tree.DecisionTreeRegressor()
        elif self.modeltype == "RandomForest":
            reg = ensemble.RandomForestRegressor()
        else:
            raise ValueError("Unexpected string in modeltype")
    
        if self.modeltype == "SVR":
            params = {'kernel': ['linear', 'poly', 'rbf', 'sigmoid'], 'C': [1,2.5,5], 'gamma': ['scale','auto'], 'degree' : [2,3,4], 'epsilon': [0.001, 0.01, 0.1, 1, 10, 100]}
            gs_reg = GridSearchCV(reg, params)
            gs_reg.fit(self.scores_train, self.features_train)
            reg.set_params(**gs_reg.best_params_)
        elif self.modeltype == "DecisionTree":
            params = {'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'splitter':['best','random'],'max_features':['auto', 'sqrt','log2']}
            gs_reg = GridSearchCV(reg, params)
            gs_reg.fit(self.scores_train, self.features_train)
            reg.set_params(**gs_reg.best_params_)
        elif self.modeltype == "RandomForest":
            params = {'n_estimators':[100,200,300],'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'max_features':['auto', 'sqrt','log2']}
            gs_reg = GridSearchCV(reg, params)
            gs_reg.fit(self.scores_train, self.features_train)
            reg.set_params(**gs_reg.best_params_)
        
        reg.fit(self.scores_train, self.features_train)
        pickle.dump(reg, open(f"models/{self.modeltype}_{self.featuretype}_{self.datatype}_{self.add_information}.sav","wb"))

    def phantomtest(self, testing_scores = 0, testing_features = 0, stats = True, plot = False):
        if testing_scores != 0:
            self.scores_test = testing_scores
        if testing_features != 0:
            self.features_test = testing_features

        reg = pickle.load(open(f"models/{self.modeltype}_{self.featuretype}_{self.datatype}_{self.add_information}.sav","rb"))
        features_pre = reg.predict(self.features_test)
        
        if stats == True:
            self.mse = mean_squared_error(features_pre, self.features_test)
            self.r = round(pearsonr(features_pre,self.features_test).statistic, 6)
            self.r_2 = round(r2_score(features_pre,self.features_test), 6)
        
        if plot == True:
            k, d = np.polyfit(self.features_test,features_pre,deg=1)
            fig, ax = plt.subplots()
            ax.plot(self.featues_test, features_pre,'o')
            plt.axline(xy1=(0, d), slope=k, label=f'r\u00b2 = {self.r_2}', color="black",ls="--")
            plt.title(f"{self.trainingscores_type}, {self.featuretype}, {self.modeltype}")
            plt.legend()
            plt.xlabel("y_real")
            plt.ylabel("y_predicted")
            plt.savefig(f"../plots/{self.trainingscores_type}_{self.featuretype}_{self.modeltype}_{self.add_information}")
            plt.close(fig)
    
    def extract_stats(self):
        return self.modeltype, self.featuretype, self.datatype, self.mse, self.r, self.r_2, self.add_information
    