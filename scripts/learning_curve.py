from sklearn import linear_model
from sklearn import svm
from sklearn import tree
from sklearn import ensemble
from sklearn.model_selection import learning_curve
import matplotlib.pyplot as plt
import numpy as np
import phantomdragon.functions as ph

datatypes = ["all","ki","kd"]
scoretypes = ["delta G","Affinity Data Value","pKd pKi pIC50"]
modeltypes = [linear_model.LinearRegression(), 
             linear_model.RidgeCV(), 
             linear_model.LassoCV(), 
             linear_model.ElasticNetCV(), 
             svm.SVR(), 
             tree.DecisionTreeRegressor(), 
             ensemble.RandomForestRegressor()]

for k in datatypes:
    for modeltype in modeltypes:
        for score in scoretypes:
            #
            # Load Dataset
            # Create training and test split
            #
            X_train, X_test, y_train, y_test = ph.prepare_data(score,"../data/grail_scores_final.csv","../data/CASF_grail_scores.csv",f"../data/PDBbind_refined_set_{k}.csv",f"../data/PDBbind_core_set_{k}.csv","final")

            # Grid search?
            # Use learning curve to get training and test scores along with train sizes
            #
            train_sizes, train_scores, test_scores = learning_curve(estimator=modeltype, X=X_train, y=y_train,
                                                                cv=10, train_sizes=np.linspace(0.1, 1.0, 10),
                                                                n_jobs=1)
            #
            # Calculate training and test mean and std
            #
            train_mean = np.mean(train_scores, axis=1)
            train_std = np.std(train_scores, axis=1)
            test_mean = np.mean(test_scores, axis=1)
            test_std = np.std(test_scores, axis=1)
            #
            # Plot the learning curve
            #
            plt.plot(train_sizes, train_mean, color='blue', marker='o', markersize=5, label='Training Accuracy')
            plt.fill_between(train_sizes, train_mean + train_std, train_mean - train_std, alpha=0.15, color='blue')
            plt.plot(train_sizes, test_mean, color='green', marker='+', markersize=5, linestyle='--', label='Validation Accuracy')
            plt.fill_between(train_sizes, test_mean + test_std, test_mean - test_std, alpha=0.15, color='green')
            plt.title(f'Learning Curve for {k} {score} {modeltype}')
            plt.xlabel('Training Data Size')
            plt.ylabel('Model accuracy')
            plt.grid()
            plt.legend(loc='lower right')

            if " " in score or "/" in score:
                score = score.replace(" ","_")
                score = score.replace("/","div")
            
            plt.savefig(f'../results/learning_curves/{k}_{score}_{modeltype}.png')  

            if "_" in score or "div" in score:
                score = score.replace("_"," ")
                score = score.replace("div","/")
            
            plt.close()
            print(k, score, modeltype, "Done")