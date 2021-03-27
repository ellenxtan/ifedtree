rm(list=ls())
library(dplyr)
library(data.table)
library(grf)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(ggpubr)
library(ranger)
library(causalToolbox)

mypath <- "eICUdata/sp02_study/"
load(paste0(mypath, "eICUdata_spo2.RData"))

savepath <- paste0(mypath, "hosp_mortality_", Sys.Date(), "/")
dir.create(file.path(savepath))

global_seed <- 1234567

## 26723 pats
dataout <- pat_eICU_subset %>%
    mutate(hospital_id = as.character(hospital_id),
           opt_spo2 = ifelse(94<=median & median<=98, 1, 0),
           female = ifelse(gender=="F", 1, 0),
           mortality_in_Hospt = ifelse(mortality_in_Hospt, 1, 0),
           mortality_in_ICU = ifelse(mortality_in_ICU, 1, 0),
           Y = mortality_in_Hospt,
           Z = opt_spo2
    ) %>%
    filter(has_respiratory_disease == TRUE) %>%
    dplyr::select(hospital_id, patient_ID, Y, Z, female, age, bmi, sofatotal, vent_duration) %>%
    add_count(hospital_id, name="ni")

(K0 <- length(unique(dataout$hospital_id)))

tab <- data.table(table(dataout$hospital_id, dataout$Z))
names(tab)
trt_sites <- with(tab, V1[which(V2==1 & N>=50)])
ctr_sites <- with(tab, V1[which(V2==0 & N>=50)])
used_sites <- intersect(trt_sites, ctr_sites)
dataout <- dataout[dataout$hospital_id %in% used_sites, ]
table(dataout$hospital_id, dataout$Z)
# summary(dataout)
# median: age 65, bmi 28.34, sofatotal 6, duration 126.55, female

(K <- length(unique(dataout$hospital_id)))

stopifnot(anyNA(dataout) == FALSE)
stopifnot(length(unique(dataout$patient_ID)) == nrow(dataout))

dataout$duration <- dataout$vent_duration
dataout$sofa <- dataout$sofatotal
covars <- c("age", "bmi", "sofa", "duration", "female")
dataout <- data.table(dataout)

avg_n <- dataout %>%
    distinct(hospital_id, .keep_all = TRUE) %>%
    summarise(avg_n = mean(ni))
avg_n <- avg_n$avg_n

dataout <- dataout %>%
    mutate(n_wt = ni / avg_n)


################################################################################
################################# data analysis ################################
################################################################################

coord_site <- 16

### individual model - XBART
ind.lst <- list()
df1 <- c()
avgtau_df <- c()
n_wts <- c()
for (k in 1:K) {
    print(k)
    df <- dataout[hospital_id == unique(dataout$hospital_id)[k], ]
    print(table(df$Z))
    n_wts <- c(n_wts, unique(df$n_wt))

    x_train <- as.data.frame(df[, ..covars])
    y_train <- df$Y
    w_train <- df$Z

    set.seed(global_seed*k)
    myfit <- X_BART(feat = x_train, tr = w_train, yobs = y_train)
    ind.lst[[k]] <- myfit

    if (k == coord_site) {  # save site 1 data
        df1 <- copy(df)
        myfit1 <- myfit
    }

    tau_hat <- EstimateCate(myfit, x_train)

    avgtau_df <- data.table(rbind(avgtau_df, cbind(site=k, avg_tauhat=mean(tau_hat))))
}

### augmented data
aug_df <- c()
for(k in 1:K) {
    print(k)
    usefit <- ind.lst[[k]]
    n_wt <- n_wts[k]
    df1$leaves <- EstimateCate(usefit, as.matrix(df1[, ..covars]))
    aug_df <- data.table(rbind(aug_df, cbind(df1, site=k, n_wt=n_wt)))
}

### ensemble tree
aug_df$site <- factor(aug_df$site)
aug_df$female <- factor(aug_df$female)
fml <- as.formula(paste("leaves ~ site +", paste(covars, collapse="+")))
set.seed(global_seed)
ensenble <- rpart(fml, data = aug_df, weights=aug_df$n_wt)
ensem.prune <- prune(ensenble, cp=ensenble$cptable[which.min(ensenble$cptable[,"xerror"]),"CP"])
rpart.plot(ensem.prune)

pdf(paste0(savepath,"hosp_eICU_xbart_ct_spo2_wt.pdf"), width=8, height=10)
rpart.plot(ensem.prune)
dev.off()

ensenble <- rpart(fml, data = aug_df)
ensem.prune <- prune(ensenble, cp=ensenble$cptable[which.min(ensenble$cptable[,"xerror"]),"CP"])
rpart.plot(ensem.prune)

pdf(paste0(savepath,"hosp_eICU_xbart_ct_spo2_no_wt.pdf"), width=8, height=10)
rpart.plot(ensem.prune)
dev.off()


### ensemble forest
rf_wt <- ranger(fml, data=na.omit(aug_df), respect.unordered.factors=TRUE,
             importance='impurity', case.weights=aug_df$n_wt)

rf <- ranger(fml, data=na.omit(aug_df), respect.unordered.factors=TRUE,
             importance='impurity')

(  rf.imp <- rf$variable.importance %>%
        as.data.frame() %>%
        dplyr::arrange(desc(.)) %>%
        dplyr::mutate(prop = . / sum(.)) %>%
        tibble::rownames_to_column(var = "names")  )

(  imp.plt <- rf.imp %>% ggplot(aes(reorder(names, prop), prop)) +
    geom_col() +
    ylab("Relative\n Importance") + xlab("Covariates") +
    scale_x_discrete(expand = c(0, 0), limits=c("site", "bmi", "age", "duration", "female", "sofa")) +
    scale_y_continuous(expand = c(0, 0), breaks=c(0,0.5,0.25)) +
    theme_classic() +
    theme(axis.text=element_text(size=13,face="bold"),
          axis.title=element_text(size=14,face="bold"))  )


pdf(paste0(savepath, "resp_hosp_eICU_xbart_imp.pdf"), width=10, height=1.8)
imp.plt
dev.off()


## rank site by average tau_hat using local model (for figures)
smry_aug_df <- avgtau_df %>%
    arrange(avg_tauhat, .by_group=TRUE) %>%  #desc(avg_tau)
    mutate(site_rank = seq(K)) %>%
    dplyr::select(site) %>%
    mutate(site = as.numeric(site))
siteRank <- smry_aug_df$site

################################################################################
########################### plots for ensemble models ##########################
################################################################################
library(scales)
library(gridExtra)
library(shades)
library(RColorBrewer)
source("R/PlotPredSiteInteract.R")

mygrid <- 1000
preds <- list()
i <- 1
for (covar in rf.imp$names[-1]) {
    stopifnot(covar != "site")
    preds[[i]] <- PlotPredSiteInteract(aug_df, rf, rf.imp$names, var1="site", var2=covar,
                                       set="real", siteRank=siteRank, grids=mygrid, myseed=global_seed*i)
    i <- i + 1
}
(  multi.preds <- ggarrange(preds[[1]], preds[[2]], preds[[3]], preds[[4]], preds[[5]],
          legend="top", common.legend=TRUE, ncol=1, align = "v") )

pdf(paste0(savepath, "resp_hosp_eICU_xbart_pred_fixed_median.pdf"), width=10, height=12)
multi.preds
dev.off()


################################################################################
############################ best linear projection ############################
################################################################################

## focus on coord site (df1)
source("R/BestLinearProj.R")

## For EF (federated model)
X <- as.matrix(df1[, ..covars])
W <- df1[, Z]
Y <- df1[, Y]
df1$Z.hat <- regression_forest(X, W)$predictions
df1$Y.hat <- regression_forest(X, Y)$predictions

df1$site <- coord_site
df1$site <- factor(df1$site)
levels(df1$site) = levels(aug_df$site)

BestLinearProj(rf, df1, A=as.matrix(df1[, ..covars]))


## For XL (localized model)
X <- as.matrix(df1[, ..covars])
W <- df1[, Z]
Y <- df1[, Y]
df1$Z.hat <- regression_forest(X, W)$predictions
df1$Y.hat <- regression_forest(X, Y)$predictions

BestLinearProj(myfit1, df1, A=as.matrix(df1[, ..covars]))


################################################################################
########################### Hospital-level info ################################
################################################################################

siteid <- as.data.frame(cbind(unique(dataout$hospital_id), seq(1,K)))
colnames(siteid) <- c("hospital_id", "site")
siteid <- siteid %>%
    left_join(smry_aug_df, by="site") %>%
    dplyr::select(hospital_id, site_rank) %>%
    rename(site=site_rank)

hospdata <- fread("eICUdata/hospital.csv") %>%
    rename(hospital_id = hospitalid) %>%
    mutate(hospital_id = as.character(hospital_id)) %>%
    filter(hospital_id %in% dataout$hospital_id)

nidata <- dataout %>%
    group_by(hospital_id) %>%
    slice(1) %>%
    ungroup() %>%
    rename(n = ni) %>%
    dplyr::select(hospital_id, n)

hosp_smry <- siteid %>%
    left_join(nidata, by="hospital_id") %>%
    left_join(hospdata, by="hospital_id") %>%
    arrange(site)

fwrite(hosp_smry, file=paste0(savepath, "hosp_smry.csv"))


################################################################################
################################################################################
################################################################################
################################################################################


