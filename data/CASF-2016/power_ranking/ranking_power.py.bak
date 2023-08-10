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
try:
        options,args = getopt.getopt(sys.argv[1:],"hc:s:p:o:",["help","coreset=","score=","prefer=","output"])
except getopt.GetoptError:
        sys.exit()

def usage():
	print "-c or --coreset: specify the location of 'CoreSet.dat' (or a subset data file) in the CASF-2016 package"
        print "-s or --score: input your scoring file name. Remember the 1st column name is #code and the 2nd column name is score. Supported file separators are comma(,), tabs(\\t) and space character( )"
        print "-p or --prefer: input 'negative' or 'positive' string, depend on your scoring funtion preference"
	print "-o or --output: input the prefix of output result files. Default is My_Ranking_Power"
        print "-h or --help: print help message"
        print "\nExample: python ranking_power.py -c CoreSet.dat -s ./examples/X-Score.dat -p 'positive' -o 'X-Score' > MyRankingPower.out"

#Define the Predictive Index function
def cal_PI(df):
	dfsorted=df.sort_values(['logKa'],ascending=True)
	W=[]
	WC=[]
	lst=list(dfsorted.index)
	for i in np.arange(0,5):
		xi=lst[i]
		score=float(dfsorted.ix[xi]['score'])
		bindaff=float(dfsorted.ix[xi]['logKa'])
		for j in np.arange(i+1,5):
			xj=lst[j]
			scoretemp=float(dfsorted.ix[xj]['score'])
			bindafftemp=float(dfsorted.ix[xj]['logKa'])
			w_ij=abs(bindaff-bindafftemp)
			W.append(w_ij)
			if score < scoretemp:
				WC.append(w_ij)
			elif score > scoretemp:
				WC.append(-w_ij)
			else:
				WC.append(0)
	pi=float(sum(WC))/float(sum(W))
	return pi
			
	
def dec(x,y):
        if y == 2:
                return Decimal(x).quantize(Decimal('0.01'),rounding=ROUND_HALF_UP)
        if y == 3:
                return Decimal(x).quantize(Decimal('0.001'),rounding=ROUND_HALF_UP)
        if y == 4:
                return Decimal(x).quantize(Decimal('0.0001'),rounding=ROUND_HALF_UP)
	
#Read the CoreSet.dat and scoring results file	
out='My_Ranking_Power'
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

#Process the data
testdf1=pd.merge(aa,bb,on='#code')
if str(fav) == 'negative':
	testdf1['score']=testdf1['score'].apply(np.negative)
	group=testdf1.groupby('target')	
elif str(fav) == 'positive':
	group=testdf1.groupby('target')
else:
	print 'please input negative or positive'
	sys.exit()

#Get the representative complex in each cluster
def top(df,n=1,column='logKa'):
	return df.sort_values(by=column)[-n:]
toptardf=testdf1.groupby('target').apply(top)
targetlst=toptardf['#code'].tolist()

#Calculate the Spearman correlation coefficient, Kendall correlation coefficient and Predictive index
spearman=pd.DataFrame(index=targetlst,columns=['spearman'])
kendall=pd.DataFrame(index=targetlst,columns=['kendall'])
PI=pd.DataFrame(index=targetlst,columns=['PI'])
rankresults=pd.DataFrame(index=range(1,len(targetlst)+1),columns=['Target','Rank1','Rank2','Rank3','Rank4','Rank5'])
tmp=1
for i,j in group.__iter__():
        testdf2=group.get_group(i)[['#code','logKa','score']]
	testdf2=testdf2.sort_values('score',ascending=False)
	tartemp=top(testdf2)['#code'].tolist()
	tar=''.join(tartemp)
	if len(testdf2) == 5:
		spearman.ix[tar]['spearman']=testdf2.corr('spearman')['logKa']['score']
	        kendall.ix[tar]['kendall']=testdf2.corr('kendall')['logKa']['score']
        	PI.ix[tar]['PI']=cal_PI(df=testdf2)
		rankresults.ix[tmp]['Rank1']=''.join(testdf2[0:1]['#code'].tolist())
		rankresults.ix[tmp]['Rank2']=''.join(testdf2[1:2]['#code'].tolist())
		rankresults.ix[tmp]['Rank3']=''.join(testdf2[2:3]['#code'].tolist())
                rankresults.ix[tmp]['Rank4']=''.join(testdf2[3:4]['#code'].tolist())
		rankresults.ix[tmp]['Rank5']=''.join(testdf2[4:5]['#code'].tolist())
		rankresults.ix[tmp]['Target']=tar
		tmp+=1
	else:
		spearman.drop(tar,inplace=True)
		kendall.drop(tar,inplace=True)
		PI.drop(tar,inplace=True)

#Print the output of ranking power evluation
spearmanmean=dec(float(spearman['spearman'].sum())/float(spearman.shape[0]),3)
kendallmean=dec(float(kendall['kendall'].sum())/float(kendall.shape[0]),3)
PImean=dec(float(PI['PI'].sum())/float(PI.shape[0]),3)
tmplen=len(PI)
spearman.to_csv(out+'_Spearman.results',sep='\t',index_label='#Target')
kendall.to_csv(out+'_Kendall.results',sep='\t',index_label='#Target')
PI.to_csv(out+'_PI.results',sep='\t',index_label='#Target')
if os.path.exists('cstemp'):
	os.remove('cstemp')
#Output results
rankresults.dropna(axis=0,inplace=True)
rankresults.style.set_properties(align="right")
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
print rankresults
print ("\nSummary of the ranking power: ===========================================")
print ("The Spearman correlation coefficient (SP) = %0.3f"%(dec(spearmanmean,3)))
print ("The Kendall correlation coefficient (tau) = %0.3f"%(dec(kendallmean,3)))
print ("The Predictive index (PI) = %0.3f"%(dec(PImean,3)))
print ("=========================================================================\n")
print ("\nTemplate command for running the bootstrap in R program==================\n")
print ("rm(list=ls());\nrequire(boot);\ndata_all<-read.table(\"%s_Spearman.results\",header=TRUE);\ndata<-as.matrix(data_all[,2]);"%(out))
print ("mymean<-function(x,indices) sum(x[indices])/%d;"%(tmplen))
print ("data.boot<-boot(aa,mymean,R=10000,stype=\"i\",sim=\"ordinary\");\nsink(\"%s_Spearman-ci.results\");\na<-boot.ci(data.boot,conf=0.9,type=c(\"bca\"));\nprint(a);\nsink();\n"%(out))
print ("=========================================================================\n")

