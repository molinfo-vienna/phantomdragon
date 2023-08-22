#!/bin/bash

datatypes="all ki kd"
modeltypes="linearRegression Ridge Lasso ElasticNet SVR DecisionTree RandomForest"
additional_information="final"
scoretypes="delta_G Affinity_Data_Value pKd_pKi_pIC50"

mkdir results_power_ranking_negative

for data in $datatypes
do
	for model in $modeltypes
	do
		for add_info in $additional_information
		do
			for score in $scoretypes
			do
				python ranking_power.py -c CoreSet.dat -s GAP_Scores/$data$model$add_info$score.csv -p 'negative' -o 'GAP' > results_power_ranking_negative/$data$model$add_info$score.out
				echo $data $model $add_info $score "done"
			done
		done
	done
done
