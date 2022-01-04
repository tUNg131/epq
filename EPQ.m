load kalmanfull.txt ;
data = dlmread('kalmanfull.txt');
[result_data, var_ran_walk, i] = AllenVariance(data, 100);
Kalman(result_data, var_ran_walk, i)