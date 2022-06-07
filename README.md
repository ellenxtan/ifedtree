# ifedtree

R package `ifedtree` and replication codes for paper "**A Tree-based Model Averaging Approach for Personalized Treatment Effect Estimation from Heterogeneous Data Sources**" [[arxiv]](https://arxiv.org/abs/2103.06261), which has recently been accepted by **ICML 2022** for publication.   
An earlier version received the Student Research Award at the 35th New England Statistics Symposium and Honorable Mention Award in the ASA Student Paper Award competition (SLDS section) at JSM 2021.

## Package installation

To install this package in R, run the following commands:

```R
install.packages("devtools")
devtools::install_github("ellenxtan/ifedtree")
```

## Usage examples

```R
library(ifedtree)
data(SimDataLst)
K <- length(SimDataLst)
covars <- grep("^X", names(SimDataLst[[1]]), value=TRUE)

# coordinating site
coord_id <- 1  
coord_test <- GenSimData(coord_id)
coord_df <- SimDataLst[[coord_id]]

# local models (causal forest from `grf` package as an example)
fit_lst <- list()
for (k in 1:K) {
    df <- SimDataLst[[k]]
    fit_lst[[k]] <- grf::causal_forest(X=as.matrix(df[, covars, with=FALSE]), Y=df$Y, W=df$Z)
}

# augmented coordinating site data
aug_df <- GenAugData(coord_id, coord_df, fit_lst, covars)

# ensemble tree
et_fit <- EnsemTree(coord_id, aug_df, "site", covars)$myfit
PlotTree(et_fit)
BestLinearProj(et_fit, coord_df, coord_id, "site", "Z", "Y", covars)

# ensemble forest
ef_fit <- EnsemForest(coord_id, aug_df, "site", covars, importance="impurity")$myfit
PlotForestImp(ef_fit)
PlotForestPred(aug_df, coord_df, coord_id, ef_fit, "site", covars, "site", "X1")
BestLinearProj(ef_fit, coord_df, coord_id, "site", "Z", "Y", covars)
```

## Replication for paper

- Simulation codes are under folder code_for_paper
    - Run the following bash script.
        ```bash
        sbatch run_hetero.sh
        ```
- Real data access and preparation
    - The eICU Collaborative Research Database [[paper]](https://www.nature.com/articles/sdata2018178) [[website]](https://eicu-crd.mit.edu/).
    - Data are preprocessed following paper "The Search for Optimal Oxygen Saturation Targets in Critically Ill Patients: Observational Data from Large ICU Databases" [[paper]](https://doi.org/10.1016/j.chest.2019.09.015) [[github]](https://github.com/nus-mornin-lab/oxygenation_kc).


## Cite

If you use our code, please cite
```
@article{tan2021tree,
  title={A tree-based federated learning approach for personalized treatment effect estimation from heterogeneous data sources},
  author={Tan, Xiaoqing and Chang, Chung-Chou H and Tang, Lu},
  journal={arXiv preprint arXiv:2103.06261},
  year={2021}
}
```
