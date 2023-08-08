#!/bin/bash

datatypes="all ki kd"
modeltypes="linearRegression Ridge Lasso ElasticNet SVR DecisionTree RandomForest"
additional_information="final"
scoretypes="delta_G Affinity_Data_Value pKd_pKi_pIC50"

mkdir results_reverse_screening

for data in $datatypes
do
	for model in $modeltypes
	do
		for add_info in $additional_information
		do
			for score in $scoretypes
			do
				python reverse_screening_power.py -c CoreSet.dat -s GAP_Scores/$data$model$add_info$score -l LigandInfo.dat -p 'positive' -o 'GAP-Score' > results_reverse_screening/$data$model$add_info$score.out || python reverse_screening_power.py -c CoreSet.dat -s GAP_Scores/$data$model$add_info$score -l LigandInfo.dat -p 'negative' -o 'GAP-Score' > results_reverse_screening/$data$model$add_info$score.out
				echo $data $model $add_info $score "done"
			done
		done
	done
done
