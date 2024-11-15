# Clustering of Motile Cell in a Confluent Layer
Code repository for paper "Formation of motile cell clusters in heterogeneous model tumors: the role of cell-cell alignment" published in Physical Review E by Quirine J.S. Braat, Cornelis Storm and Liesbeth M.C. Janssen. 

## Simulation code
We simulated the primary tumor tissue with the Cellular Potts Model. The simulation code is written in CompuCell3D (version 4.2.5) [1]. 

The simulations are run for multiple situations: 
1) Active cells in a confluent layer (full dynamic data): "ConfluentLayer_Random_full"
2) Active cells in an empty layer (full dynamic data): "EmptyLayer_Random_full"
3) Active cells in a confluent layer (single export data): "ConfluentLayer_Random_empty"
4) Active cells in an empty layer (single export data): "EmptyLayer_Random_empty"
5) Active cells without alignment interactions (full dynamic data): "EmptySpace_gamma0" and "ConfluentLayer_gamma0"
6) Active cells starting as a single aligned cluster in the confluent layer (single export data): "ConfluentLayer_initialblock"

## Data 
Each data point is generated by running 200 simulations using the simulation code for a set of parameters. We ran the "XXX_full" to get the full dynamic evolution of all cells (positions, vector orientation, cluster id etc.) and we ran "XXX_empty" to get some of the relevant steady state values. Both data is stored at Zenodo [2]. 

For each combination of parameters (tau, gamma), we either got the full dynamic information (Full) or only the dynamic evolution of the mean cluster size and the steady state values for the mean cluster size and order parameter (Single).

#### Confluent layer ($\phi = 0.25$)

| tau / gamma   | 1.0  | 0.5 | 0.2 | 0.1  | 0.05 | 0.02 | 0.01 | 0.005 |0.002 | 0.001 |
| ------------- |------|-----|-----|------|------|------|------|-------|------|-------|
| 500 mcs       | Full | Single | Single | Full | Single | Single | Full | Full | Single | Full  |
| 2500 mcs      | Full | Single | Single | Full | Single | Single | Full | Full | Single | Full  |
| 4000 mcs      | Full | Single | Single | Full | Single | Single | Full | Full | Single | Full  |


#### Empty layer ($\phi = 0.25$)

| tau / gamma   | 1.0  | 0.5 | 0.2 | 0.01 | 0.05 | 0.02 | 0.01 | 0.005 |0.002 | 0.001 |
| ------------- |------|-----|-----|------|------|------|------|-------|------|-------|
| 500 mcs       | Full | Full | Full | Full | Full | Full | Full | Full | Full | Full |
| 2500 mcs      | Full | Full | Full | Full | Full | Full | Full | Full | Full | Full |
| 4000 mcs      | Single | Single | Single | Single | Single | Single | Single | Single | Single | Single |

#### Other simulations

Apart from the main results, we also ran simulations with different parameter settings to get more insights into the dynamic behavior. 

| Simulations | Settings | Storage type |
| ------------|----------| -------------|
| Different fraction active cells | Confluent layer, $\phi = 0.1$, $\tau, \gamma$ as before | Single |
| No alignment | Confluent Layer + Empty Layer, $\phi = 0.25$, $\gamma = 0$, $\tau = 2500$ mcs | Full |
| Initially aligned cluster | Fully aligned cluster, Confluent Layer + Empty Layer, $\gamma = 1.0, 0.01, 0.001$ and $\tau = 2500$ mcs | Full |
| Finite-size effects | Confluent Layer, number of cells = 100, 400, 900, 1600, 2500 cells | Single | 


## Sources
[1] Swat, M. H., Thomas, G. L., Belmonte, J. M., Shirinifard, A., Hmeljak, D., & Glazier, J. A. (2012). Multi-scale modeling of tissues using CompuCell3D. Methods in cell biology, 110, 325–366. https://doi.org/10.1016/B978-0-12-388403-9.00013-8

[2] Braat, Q., Storm, C., & Janssen, L. (2024). Data for paper 'Formation of motile cell clusters in heterogeneous model tumors: The role of cell-cell alignment' (PRE, 2024) [Data set]. In Physical Review E. Zenodo. https://doi.org/10.5281/zenodo.13949358
