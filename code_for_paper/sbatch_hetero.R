files.sources <- list.files(c("R"), pattern="*.R$", full.names=TRUE)
suppressMessages(invisible(sapply(files.sources, source)))

####
loc.mod       <- "CT"
ni_tau        <- 500
mydate        <- paste0("tau-",loc.mod,"-ni",ni_tau)
####
outcomes      <- c("tau")
grp_types     <- c("2grp", "cont")
program_name1 <- "run_hetero.sh"
program_name2 <- "sbatch_hetero.R"

R_tau <- 1000
scale_tau  <- c(0, 0.6, 1, 2)

tau_fn_tau <- parse(text = "X1*(X1>0) -3*U*scale_c + X1*U*scale_c")
m_fn_tau   <- parse(text = "1/2*(X1) + X2+X3+X4 -3*U*scale_c + X1*U*scale_c")


K        <- 20
xpdf     <- "norm"
xdim     <- 5
coord    <- 1
n_test   <- 5000
honest1  <- TRUE
honest2  <- FALSE


################################### Default ####################################

## save folders
suppressWarnings(dir.create(file.path("results")))
date_folder <- paste0("results/res_", mydate)
suppressWarnings(dir.create(file.path(date_folder)))
message(paste0("date_folder: ", date_folder))
code_folder <- paste0(date_folder, "/codes")
suppressWarnings(dir.create(file.path(code_folder)))
message(paste0("code_folder: ", code_folder))
system(paste0("cp ", program_name1, " ", code_folder, "/"))
system(paste0("cp ", program_name2, " ", code_folder, "/"))
system(paste0("cp -r R ", code_folder, "/"))
log_folder <- paste0(date_folder, "/logs")
suppressWarnings(dir.create(file.path(log_folder)))
message(paste0("log_folder: ", log_folder))
for (out in outcomes) {
    for (grp in grp_types) {
        set_folder <- paste0(date_folder, "/", out, "_", grp)
        suppressWarnings(dir.create(file.path(set_folder)))
        message(paste0("set_folder: ", set_folder))
    }
}

## job tasks
task_lst <- list()
ii <- 1
for (out in outcomes) {
    if (out == "tau") { 
        scales <- scale_tau
        R <- R_tau
    }
    else if (out == "y") { 
        scales <- scale_y 
        R <- R_y
    }
    
    for (grp in grp_types) {
        for (scale_c in scales) {
            task_lst[[ii]] <- list()
            task_lst[[ii]]$outcome <- out
            task_lst[[ii]]$grp_type <- grp
            task_lst[[ii]]$scale_c <- scale_c
            task_lst[[ii]]$R <- R
            ii <- ii + 1
        }
    }
}
num_tasks <- ii - 1
cat("total # of tasks:", num_tasks, "\n")


#################################### Input #####################################

## read inputs
(ii      <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))
outcome  <- task_lst[[ii]]$outcome
grp_type <- task_lst[[ii]]$grp_type
scale_c  <- task_lst[[ii]]$scale_c
R        <- task_lst[[ii]]$R
note     <- ""

design_name <- paste0(outcome,"_",grp_type,"_scale",scale_c,"_R",R,note)

ni     <- ni_tau
R      <- R_tau
scales <- scale_tau
tau_fn <- tau_fn_tau
m_fn   <- m_fn_tau

cat("=========\n=========\n")
cat("date_folder:", date_folder, "\n")
cat("code_folder:", code_folder, "\n")
cat("log_folder:", log_folder, "\n")
cat("total # of tasks:", num_tasks, "\n")  # 42 tasks / 10 tasks
print(c("scales",scales))
cat("task", ii, "\n")
cat("outcome", outcome, "\n")
cat("R", R, "\n")
cat("ni", ni, "\n")
cat("grp_type", grp_type, "\n")
cat("scale_c", scale_c, "\n")
print(c("tau_fn",tau_fn))
print(c("m_fn",m_fn))
cat("loc.mod", loc.mod, "\n")
cat("=========\n=========\n")


#################################### Main ######################################

sink_file <- paste0(log_folder, "/log_", design_name, ".txt")
sink(sink_file)

cat("=========\n=========\n")
cat("date_folder:", date_folder, "\n")
cat("code_folder:", code_folder, "\n")
cat("log_folder:", log_folder, "\n")
cat("total # of tasks:", num_tasks, "\n")
print(c("scales",scales))
cat("task", ii, "\n")
cat("outcome", outcome, "\n")
cat("R", R, "\n")
cat("ni", ni, "\n")
cat("grp_type", grp_type, "\n")
cat("scale_c", scale_c, "\n")
print(c("tau_fn",tau_fn))
print(c("m_fn",m_fn))
cat("loc.mod", loc.mod, "\n")
cat("=========\n=========\n")

#### main
res.out <- replicate(R, 
                     Iter1(K, ni, xdim, xpdf, honest1, honest2, loc.mod,
                           tau_fn, m_fn, scale_c, coord, n_test, grp_type, outcome), 
                     simplify=FALSE)


#### print results
method_names <- names(res.out[[1]])
res.df <- list()
for (name in method_names) {
    print(name)
    res.df[[name]] <- LstToDF(res.out, name)
    res.df[[name]] %>% 
        apply(2, mean,na.rm=TRUE) %>% 
        print()
}


#### save results
res.save <- list(res.out=res.out, res.df=res.df, method_names=method_names,
                 outcome=outcome, grp_type=grp_type, scale_c=scale_c, R=R)

saveRDS(res.save, file=paste0(date_folder,"/",outcome,"_",grp_type,"/res_",design_name,".rds"))
cat("Save results!\n=========\n")


#################################### End #######################################

# remove null_CT.txt
if (file.exists("null_CT.txt")) {
    system('rm null_CT.txt')
}

### end
warnings()
cat("Finished scale", scale_c, "\n=========\n")

sink()

warnings()
cat("Finished scale", scale_c, "\n=========\n")
