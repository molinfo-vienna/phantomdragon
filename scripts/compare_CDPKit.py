from sklearn.metrics import mean_squared_error, r2_score
from scipy import stats
import statistics
import phantomdragon.functions as ph
import pandas as pd

k_list = ["ki","kd","all"]

modeltype_list = []
scoretype_list = []
datatype_list = []
mse_list = []
sd_list = []
pear_list = []
coef_list = []
confidence_interval_list =[]
add_info_list = []




for k in k_list:
    scores_test = pd.read_csv(f"../data/PDBbind_core_set_{k}.csv")
    scores_test = ph.filter_and_sort_features(scores_test,scores_test)
    scores_pre = pd.read_csv(f"../data/core_set_preds.csv")
    scores_pre = ph.filter_and_sort_features(scores_test,scores_pre)
    scores_test = scores_test["pKd pKi pIC50"]
    if k == "ki":
        scores_pre = scores_pre[" pKi"]
    elif k == "kd":
        scores_pre = scores_pre[" pKd"]
    elif k == "all":
        scores_pre = scores_pre[" pKi/pKd"]
    else:
        print("something went wrong")

    mse = mean_squared_error(scores_test, scores_pre)
    sd = statistics.stdev(scores_pre)
    r = round(stats.pearsonr(scores_test,scores_pre).statistic, 6)
    conf_int = stats.pearsonr(scores_test,scores_pre).confidence_interval(confidence_level=0.9)
    conf_int_low = round(conf_int.low,3)
    conf_int_high = round(conf_int.high,3)
    conf_int = f"[{conf_int_low} ~ {conf_int_high}]"
    r_2 = round(r2_score(scores_test,scores_pre), 6)
    modeltype_list.append("polynomialRegression")
    scoretype_list.append("pKd pKi pIC50")
    datatype_list.append(f"{k}")
    mse_list.append(mse)
    sd_list.append(sd)
    pear_list.append(r)
    confidence_interval_list.append(conf_int)
    coef_list.append(r_2)
    add_info_list.append("CDPKit")

data = {'Modeltype':modeltype_list,'Scoretype':scoretype_list,'Datatype':datatype_list,'Mean squared error (mse)':mse_list,'Standard Diviation (SD)':sd_list,'Pearson correlation coefficient (r)':pear_list,'90% Confidence interval':confidence_interval_list,'Coefficient of determination (r²)':coef_list,'add. information':add_info_list}
df = pd.DataFrame(data)
df.to_csv("../results/CDPKit_results.csv")