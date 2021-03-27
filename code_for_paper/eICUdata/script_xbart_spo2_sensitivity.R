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


################################################################################
################################# data analysis ################################
################################################################################

### sensitivity analysis: fix siteRank for plot
avgtau_df <- c()
for (k in 1:K) {
    print(k)
    df <- dataout[hospital_id == unique(dataout$hospital_id)[k], ]

    x_train <- as.data.frame(df[, ..covars])
    y_train <- df$Y
    w_train <- df$Z

    set.seed(global_seed*k)
    myfit <- X_BART(feat = x_train, tr = w_train, yobs = y_train)

    tau_hat <- EstimateCate(myfit, x_train)
    avgtau_df <- data.table(rbind(avgtau_df, cbind(site=k, avg_tauhat=mean(tau_hat))))
}

## rank site by average tauhat using local model (for figures)
smry_aug_df <- avgtau_df %>%
    arrange(avg_tauhat, .by_group=TRUE) %>%
    mutate(site_rank = seq(K)) %>%
    dplyr::select(site) %>%
    mutate(site = as.numeric(site))
siteRank <- smry_aug_df$site  # fix it for all coord_site (in sensitivity analysis)
cat("siteRank:", siteRank, "\n")
# [1]  4  9 11  8 16 18 17 15  3 12 20
# [12]  7 13 10  6 19  2 14  1  5


### begin analysis
## local XBART
for (s in 1:K) {  # try all as coordinating site
    coord_site <- s

    ### individual model - XBART
    ind.lst <- list()
    df1 <- c()
    avgtau_df <- c()
    for (k in 1:K) {
        df <- dataout[hospital_id == unique(dataout$hospital_id)[k], ]
        if (k == coord_site) { df1 <- copy(df) }  # save site 1 data

        x_train <- as.data.frame(df[, ..covars])
        y_train <- df$Y
        w_train <- df$Z

        set.seed(global_seed*k)
        myfit <- X_BART(feat = x_train, tr = w_train, yobs = y_train)
        ind.lst[[k]] <- myfit
    }

    ### augmented data
    aug_df <- c()
    for(k in 1:K) {
        cat("coord:", coord_site, "; k:", k, "\n")
        usefit <- ind.lst[[k]]
        df1$leaves <- EstimateCate(usefit, as.matrix(df1[, ..covars]))
        aug_df <- data.table(rbind(aug_df, cbind(df1, site = k)))
    }

    ### ensemble tree
    aug_df$site <- factor(aug_df$site)
    aug_df$female <- factor(aug_df$female)
    fml <- as.formula(paste("leaves ~ site +", paste(covars, collapse="+")))
    set.seed(global_seed)
    ensenble <- rpart(fml, data = aug_df)
    ensem.prune <- prune(ensenble, cp=ensenble$cptable[which.min(ensenble$cptable[,"xerror"]),"CP"])
    # rpart.plot(ensem.prune)

    # pdf(paste0(savepath,"hosp_eICU_ct_spo2.pdf"), width=8, height=10)#, units="in", res=800)
    # rpart.plot(ensem.prune)
    # dev.off()

    ### ensemble forest
    rf <- ranger(fml, data=na.omit(aug_df), respect.unordered.factors=TRUE, importance='impurity')

    rf.imp <- rf$variable.importance %>%
        broom::tidy() %>%
        dplyr::arrange(desc(x)) %>%
        mutate(prop = x / sum(x))  # site: 64.2%

    (  imp.plt <- rf.imp %>% ggplot(aes(reorder(names, prop), prop)) +
        geom_col() +
        coord_flip() +
        ylab("Relative importance") + xlab("Covariates") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() +
        theme(legend.position="top",
              legend.title=element_text(size=13,face="bold"),
              legend.text=element_text(size=11,face="bold"),
              axis.text=element_text(size=12,face="bold"),
              axis.title=element_text(size=14,face="bold"))  )

    pdf(paste0(savepath, "sensitivity/resp_hosp_eICU_xbart_imp_", coord_site, ".pdf"), width=10, height=12)
    print(imp.plt)
    dev.off()


    ################################################################################
    ########################### plots for ensemble models ##########################
    ################################################################################
    source("R/PlotPredSiteInteract.R")
    library(scales)
    library(gridExtra)
    library(shades)
    library(RColorBrewer)

    preds <- list()
    i <- 1
    for (covar in rf.imp$names[-1]) {
        stopifnot(covar != "site")
        preds[[i]] <- PlotPredSiteInteract(aug_df, rf, rf.imp$names, var1="site", var2=covar,
                                           set="real", siteRank=siteRank, grids=500, myseed=global_seed*i)
        i <- i + 1
    }
    multi.preds <- ggarrange(preds[[1]], preds[[2]], preds[[3]], preds[[4]], preds[[5]],
                             legend="top", common.legend=TRUE, ncol=1, align = "v")

    pdf(paste0(savepath, "sensitivity/resp_hosp_eICU_xbart_pred_coord", coord_site, ".pdf"), width=10, height=12)
    print(multi.preds)
    dev.off()

    print("Figure saved!")
}


