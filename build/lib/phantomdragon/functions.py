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

def filter_and_sort_features(exp_data, features):
    '''
    Parameters
    ----------
    exp_data : TYPE PANDAS Dataframe
        A pandas Dataframe of the avavible experimental scores.
    features : TYPE PANDAS Dataframe
        A pandas Dataframe of the avaviable features. 
        It is assumed that exp_data is corresponding to a subset of score_features.

    Returns
    -------
    df : TYPE PANDAS Dataframe
        Gives back a sorted subset of score_features. 
        This subset corresponds to the subset of the experimental data.
        It is sorted and filtered according to the "PDB code" bracket.

    '''
    df = features[features["PDB code"].isin(exp_data["PDB code"])]
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

class parameterCollector:
    def __init__(self, modeltype, add_information, scoretype):
        self.modeltype = modeltype
        self.add_information = add_information
        self.scoretype = scoretype
        self.mse = None
        self.r = None
        self.r_2 = None
    
    def get_data(self):
        x_train, y_train = self.get_trainingdata()
        x_test, y_test = self.get_testingdata()
        return x_train, x_test, y_train, y_test
    
    def get_trainingdata(self):
        return self.features_train, self.scores_train
    
    def get_testingdata(self):
        return self.features_test, self.scores_test

    def define_trainingdata(self, features_train, scores_train):
        self.features_train = features_train
        self.scores_train = scores_train

    def define_testingdata(self, features_test, scores_test):
        self.features_test = features_test
        self.scores_test = scores_test

    def define_datatype(self, datatype):
        self.datatype = datatype

    def define_trainingscoretype(self, trainingscores_type):
        self.trainingscores_type = trainingscores_type

    def load_data(self,featurepath,experimentpath,split=0.2,sh=True):
        features = pd.read_csv(featurepath)
        experiment = pd.read_csv(experimentpath)
        experiment = experiment.sort_values("PDB code")
        experiment = experiment.reset_index()
        experiment = experiment.drop("index", axis=1)
        scores = np.array(experiment[self.scoretype])

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
        elif self.add_information == "-w-el-vdw":
            drop = ["PDB code"," H-H_SUM"," H-H_MAX"]
        elif self.add_information == "basic-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        elif self.add_information == "-w-vdw":
            drop = ["PDB code"," HW-HW_SUM"," HW-HW_MAX"," ES"," VDW_ATT"," VDW_REP"]
        else:
            raise ValueError("Unexpected string in add_information")
        

        features = filter_and_sort_features(experiment,features)
        features = features.drop(drop, axis=1)
        features = features.to_numpy()
        
        if sh == True:
            features, scores = shuffle(features, scores)
        
        scaler = preprocessing.StandardScaler().fit(features)
        features = scaler.transform(features)

        x_train, x_test, y_train, y_test = train_test_split(features, scores, test_size = split)
        
        self.features_train = x_train
        self.features_test = x_test
        self.scores_train = y_train
        self.scores_test = y_test
    
    def train_and_save_model(self,savepath=""):
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
            gs_reg.fit(self.features_train, self.scores_train)
            reg.set_params(**gs_reg.best_params_)
        elif self.modeltype == "DecisionTree":
            params = {'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'splitter':['best','random'],'max_features':['auto', 'sqrt','log2']}
            gs_reg = GridSearchCV(reg, params)
            gs_reg.fit(self.features_train, self.scores_train)
            reg.set_params(**gs_reg.best_params_)
        elif self.modeltype == "RandomForest":
            params = {'n_estimators':[100,200,300],'criterion':['squared_error','friedman_mse','absolute_error','poisson'],'max_features':['auto', 'sqrt','log2']}
            gs_reg = GridSearchCV(reg, params)
            gs_reg.fit(self.features_train, self.scores_train)
            reg.set_params(**gs_reg.best_params_)
        
        reg.fit(self.features_train, self.scores_train)
        
        if "/" in self.scoretype:
            self.scoretype = self.scoretype.replace("/","div")

        pickle.dump(reg, open(f"{savepath}{self.modeltype}_{self.scoretype}_{self.datatype}_{self.add_information}.sav","wb"))
        
        if "_" in self.scoretype:
            self.scoretype = self.scoretype.replace("div","/")
        

    def phantomtest(self, testing_features = 0, testing_scores = 0,loadpath=""):
        if testing_features != 0:
            self.features_test = testing_features
        if testing_scores != 0:
            self.scores_test = testing_scores

        if "/" in self.scoretype:
            self.scoretype = self.scoretype.replace("/","div")

        reg = pickle.load(open(f"{loadpath}{self.modeltype}_{self.scoretype}_{self.datatype}_{self.add_information}.sav","rb"))

        if "_" in self.scoretype:
            self.scoretype = self.scoretype.replace("div","/")

        scores_pre = reg.predict(self.features_test)
        
        self.scores_pre = scores_pre

        self.mse = mean_squared_error(self.scores_test, self.scores_pre)
        self.r = round(pearsonr(self.scores_test,self.scores_pre).statistic, 6)
        self.r_2 = round(r2_score(self.scores_test,self.scores_pre), 6)
    
    def plot_phantomtest(self,savepath):
        k, d = np.polyfit(self.scores_test,self.scores_pre,deg=1)
        fig, ax = plt.subplots()
        ax.plot(self.scores_test, self.scores_pre,'o')
        plt.axline(xy1=(0, d), slope=k, label=f'r\u00b2 = {self.r_2}', color="black",ls="--")
        plt.title(f"{self.scoretype}, {self.modeltype}")
        plt.legend()
        plt.xlabel("y_real")
        plt.ylabel("y_predicted")
        
        if " " in self.scoretype or " " in self.modeltype or " " in self.add_information or "/" in self.scoretype:
            self.scoretype = self.scoretype.replace(" ","_")
            self.scoretype = self.scoretype.replace("/","div")
            self.modeltype = self.modeltype.replace(" ","_")
            self.add_information = self.add_information.replace(" ","_")

        plt.savefig(f"{savepath}{self.scoretype}_{self.modeltype}_{self.add_information}.png")
        
        if "_" in self.scoretype or "_" in self.modeltype or "_" in self.add_information or "div" in self.scoretype:
            self.scoretype = self.scoretype.replace("_"," ")
            self.scoretype = self.scoretype.replace("div","/")
            self.modeltype = self.modeltype.replace("_"," ")
            self.add_information = self.add_information.replace("_"," ")
        
        plt.close(fig)

    def get_stats(self):
        return self.modeltype, self.scoretype, self.datatype, self.mse, self.r, self.r_2, self.add_information
    