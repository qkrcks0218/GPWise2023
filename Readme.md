# Replication Files for Park and Kang (2019)

This Github repository contains replication files for [Park and Kang (2019)](https://arxiv.org/abs/1908.04427 "GPwise").


## Folders

### ECLSK

ECLSK folder contains 1.Analysis.R, and 3.Summary.R files, and replicates the data analysis in Section 6 of the paper. 

* 1.Analysis.R file cleans the raw data (available at https://nces.ed.gov/ecls/dataproducts.asp), implements our method, and generates Result\_BATCH####.csv files in NData folder. 

* 2.Summary.R file summarizes Result\_BATCH####.csv files and estimates the groupwise effects based on the causal forest using `grf` package.

* NData folder contains csv files titled as BATCH\_BATCH####.csv. Here B#### represents the index of sample splitting procedure ranging from B0001 to B0100.

### Simulation

Simulation folder contains plot, and Result, Summary\_Result folders, and 1.Estimation.R and 2.Summary.R files.

* 1.Estimation.R replicates the simulation analysis in Section 4 of the paper. The result of the simulation is saved in Result folder. 

* 2.Summary.R summarizes the results in Result folder.

* Result folder contains csv files titled as Result\_###\_COR\_###\_rho###\_B####\_SubB####.csv. The first ### takes either CPS or VPS according to the propensity score model. The second ### takes either Est or CL according to the i.i.d. model or the clustered model. The third #### varies from 0000 to 0003 according to the effect heterogeneity. The last #### varies from 0000 to 0100 according to the index of sample splitting procedure. 

* Summary\_Result folder contains csv files titled as Result\_###\_COR\_###\_rho###\_B####.csv which is the merged file of the csv files in Result folder with the format Result\_###\_COR\_###\_rho###\_B####\_SubB####.csv.

* plot folder contains graphical summary of the simulation results. 


## Code

* MySL.R contains functions used for implementing superlearner algorithm and estimating the nuisance functions.

* SSLS.R contains functions used for the estimation of the groupwise effect using the methods in Section 2 and 3 of the main paper. 


## References
Chan Park & Hyunseung Kang (2019) **A Groupwise Approach for Inferring Heterogeneous Treatment Effects in Causal Inference**, _arXiv:1908.04427_ [[link](https://arxiv.org/abs/1908.04427 "GPwise")]