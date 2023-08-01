#!/usr/bin/python
import numpy as np
import sys,os
import pandas as pd
import scipy
import getopt
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from decimal import *

if len(sys.argv) < 2:
	print ("Please input parameter files or use -h for help")
	sys.exit()

def usage():
    print("-c or --coreset: specify the location of 'CoreSet.dat' (or a subset data file) in the CASF-2016 package")
    print("-s or --score: input your scoring file name. Remember the 1st column name is #code and the 2nd column name is score. Supported file separators are comma(,), tabs(\\t) and space character( )")
    print("-p or --prefer: input 'negative' or 'positive' string, depend on your scoring funtion preference")
    print("-o or --output: input the prefix of the output processed scoring files. Default name is My_Scoring_Power")
    print("-h or --help: print help message")
    print("\nExample: python scoring_power.py -c CoreSet.dat -s ./examples/X-Score.dat -p 'positive' -o 'X-Score' > MyScoringPower.out")


try:
	options,args = getopt.getopt(sys.argv[1:],"h:c:s:p:o:",["help","coreset=","score=","preferi=","output="])
except getopt.GetoptError:
	sys.exit()

def dec(x,y):
	if y == 2:
		return Decimal(x).quantize(Decimal('0.01'),rounding=ROUND_HALF_UP)
	if y == 3:
		return Decimal(x).quantize(Decimal('0.001'),rounding=ROUND_HALF_UP)
	if y == 4:
		return Decimal(x).quantize(Decimal('0.0001'),rounding=ROUND_HALF_UP)
#Read the CoreSet.dat and scoring results file
out='My_Scoring_Power'
for name,value in options:
	if name in ("-h","--help"):
		usage()
		sys.exit()
	if name in ("-c","--coreset"):
		f=open(value,'r')
		f1=open('cstemp','w+')
		for i in f.readlines():
			if i.startswith('#'):
				if i.startswith('#code'):
					f1.writelines(i)
				else:
					continue
			else:
				f1.writelines(i)
		f.close()
		f1.close()
		aa=pd.read_csv('cstemp',sep='[,,\t, ]+',engine='python')
		aa=aa.drop_duplicates(subset=['#code'],keep='first')
	if name in ("-s","--score"):
		filename=value
		bb=pd.read_csv(value,sep='[,,\t, ]+',engine='python')
	if name in ("-p","--prefer"):
		fav=value
	if name in ("-o","--output"):
		out=value
#Process the data and remove the outliers
testdf1=pd.merge(aa,bb,on='#code')
if str(fav) == 'positive':
	testdf2=testdf1[testdf1.score > 0]
	testdf2.to_csv(out+'_processed_score',columns=['#code','logKa','score'],sep='\t',index=False)
	if len(testdf2['score']) < 0.5*len(testdf1['score']):
		raise ValueError("less then half of the scores are positive")
elif str(fav) == 'negative':
	testdf1['score']= testdf1['score'].apply(np.negative)
	testdf2=testdf1[testdf1.score > 0]
	testdf2.to_csv(out+'_processed_score',columns=['#code','logKa','score'],sep='\t',index=False)
	if len(testdf2['score']) < 0.5*len(testdf1['score']):
		raise ValueError("less then half of the scores are negative")
else:
	print('please input negative or positive')
	sys.exit()

#Calculate the Pearson correlation coefficient
regr=linear_model.LinearRegression()
regr.fit(testdf2[['score']],testdf2[['logKa']])
testpredy=regr.predict(testdf2[['score']])
testr=scipy.stats.pearsonr(testdf2['logKa'].values,testdf2['score'].values)[0]
testmse=mean_squared_error(testdf2[['logKa']],testpredy)
num=testdf2.shape[0]
testsd=np.sqrt((testmse*num)/(num-1))
conf_int = scipy.stats.pearsonr(
	testdf2['logKa'].values,testdf2['score'].values
).confidence_interval(confidence_level=0.9)
conf_int_low = round(conf_int.low, 3)
conf_int_high = round(conf_int.high, 3)
conf_int = f"[{conf_int_low} ~ {conf_int_high}]"
if os.path.exists('cstemp'):
	os.remove('cstemp')
#Print the output of scoring power evluation
def f(x):
	return x+1
testdf1.rename(columns={'#code':'code'},inplace=True)
testdf1.index=testdf1.index.map(f)
testdf1.style.set_properties(align="right")
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
print(testdf1[['code','logKa','score']])
print("\nSummary of the scoring power: ===================================")
print("The regression equation: logKa = %.2f + %.2f * Score"%(dec(float(regr.coef_),2), dec(float(regr.intercept_),2)))
print("Number of favorable sample (N) = %d"%(num))
print("Pearson correlation coefficient (R) = %0.3f"%(dec(testr,3)))
print(f"Confidence interval of 90% Pearson R = {conf_int} ")
print("Standard deviation in fitting (SD) = %0.2f"%(dec(testsd,2)))
print("=================================================================")
print("\nTemplate command for running the bootstrap in R program==========\n")
print("rm(list=ls());\nrequire(boot);\ndata_all<-read.table(\"%s_processed_score\",header=TRUE);\naa<-c(1:nrow(data_all));\n"%(out))
print("mycore<-function(x,indices)\n{\ndata_1<-matrix(NA,%d,2);\nfor(i in 1:%d);\n    {\n        data_1[i,1]=data_all[x[indices][i],2];\n        data_1[i,2]=data_all[x[indices][i],3];\n    }\n        data_2<-data.frame(data_1);\n        names(data_2)<-c(\"a\",\"b\");\n        cor(data_2$a,data_2$b);\n};\n"%(num,num))
print("data.boot<-boot(aa,mycore,R=10000,stype=\"i\",sim=\"ordinary\");\nsink(\"%s-ci.results\");\na<-boot.ci(data.boot,conf=0.9,type=c(\"bca\"));\nprint(a);\nsink();\n"%(out))
print("=================================================================")
