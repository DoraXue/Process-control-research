<h1> README </h1>

<h2> File description </h2>

`StabilityRegion_K15.mat`  --- Linear controller with K = 15. For contour and 3D plots of stability regions. Generated from <strong> `StabilityRegion.m` </strong>. 

`StabilityRegion_K20.mat`

`StabilityRegion_K25.mat`

`StabilityRegion_K30.mat`

<strong> `StabilityCondition.m` </strong> --- MATLAB code used to generate h(za). The first chunk generates the Hza corresponding to a given K, and the second chunk generates the stability region plot `Hza_Ks.eps`. `StabilityRegion_K*.mat` is called.


<strong> `main_OPT_controller.m` </strong> 
--- MATLAB code for time-domain simulations. Tunable parameters include Wa, Wu,  Wh, and H. Workspace is saved as:

`Wa_1_Wu_1_Wh_1_H_2_K_20.mat`, --- compare Wa's

`Wa_10_Wu_1_Wh_1_H_2_K_20.mat`, 

`Wa_100_Wu_1_Wh_1_H_2_K_20.mat`, 

`Wa_1000_Wu_1_Wh_1_H_2_K_20.mat`,

`Wa_1_Wu_10_Wh_1_H_2_K_20.mat`, --- compare Wu's 

`Wa_1_Wu_100_Wh_1_H_2_K_20.mat`, 

`Wa_1_Wu_1000_Wh_1_H_2_K_20.mat`, 

`Wa_1_Wu_1_Wh_10_H_2_K_20.mat`, --- compare Wh's

`Wa_1_Wu_1_Wh_100_H_2_K_20.mat`, 

`Wa_1_Wu_1_Wh_1000_H_2_K_20.mat`, 

`Wa_1_Wu_1_Wh_1_H_07_K_20.mat`, --- comapre H's

`Wa_1_Wu_1_Wh_1_H_2_K_20_long.mat`,

`Wa_1_Wu_1_Wh_1_H_4_K_20.mat`, 

`Wa_1_Wu_1_Wh_1_H_6_K_20.mat`, 

`Wa_1_Wu_1_Wh_1_H_10_K_20.mat`

<strong> `ComparisonPlots.m` </strong> --- MATLAB code calling the above mat files and generate time-domain comparisons on X, a, u, update, period, Zopt, and costs.
