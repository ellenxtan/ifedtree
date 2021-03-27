rm(list=ls())

## load functions
files.sources <- list.files(c("R"), pattern="*.R$", full.names=TRUE)
suppressMessages(invisible(sapply(files.sources, source)))

## result folder path
dir.create(file.path("results"))
out_folder <- paste0("results/res_", Sys.Date())
dir.create(file.path(out_folder))
message(paste0("out_folder: ", out_folder))

## global settings
R        <- 500
K        <- 20
ni       <- 100
is_ens   <- TRUE
is_meta  <- TRUE
is_agg   <- TRUE
is_save  <- TRUE
xpdf     <- "norm"

## parallel computing
core_num <- round(detectCores()*0.65)
cores <- makeCluster(core_num, type='PSOCK')
options('mc.cores' = cores)
registerDoParallel(cores)

for (honest in c(FALSE, TRUE)) {
    ## log file
    sink_file <- paste0("hetero","_hon",substr(honest,1,1),"_Ens",substr(is_ens,1,1),
                        "_M",substr(is_meta,1,1),"_B",substr(is_agg,1,1),"_xpdf",xpdf)
    sink(paste0(out_folder, "/", sink_file, ".txt"))
    cat(paste0("Creating result folder: ", out_folder, "\n"))
    
    cat("core:", core_num, "\n")
    cat("K:", K, "\n")
    cat("n:", ni, "\n")
    cat("xpdf:", xpdf, "\n")
    cat("honest:", honest, "; is_ens:", is_ens, "\n")
    cat("is_meta:", is_meta, "; is_agg", is_agg, "\n")
    cat("is_save", is_save, "\n")
    
    for (ind_model in c("X-BART", "CF", "CT")) {
        
        ######################## heterogeneous settings ########################
        message(sink_file)
        design = "hetero1"
        xdim <- 3
        tau_fn <- parse(text = "-4*(k%%2==0)")
        m_fn <- parse(text = "0")
        Main(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
             honest, is_meta, is_agg, is_save, out_path, core_num)
        
        message(sink_file)
        design = "hetero2"
        xdim <- 4
        tau_fn <- parse(text = "(2/(1+exp(-12*(X1-1/2)))) * (2/(1+exp(-12*(X2-1/2)))) -4*(k%%2==0)+X1*(k%%2==0)")
        m_fn <- parse(text = "0")
        Main(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
             honest, is_meta, is_agg, is_save, out_path, core_num)
        
        message(sink_file)
        design = "hetero3"
        xdim <- 5
        tau_fn <- parse(text = "X1*(X1>0) -4*(k%%2==0)+X1*(k%%2==0)")
        m_fn <- parse(text = "1/2*(X1) + X2+X3+X4 -4*(k%%2==0)+X1*(k%%2==0)")
        Main(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
             honest, is_meta, is_agg, is_save, out_path, core_num)
        
        message(sink_file)
        design = "hetero4"
        xdim <- 8
        tau_fn <- parse(text = "X1*(X1>0) + X2*(X2>0) -4*(k%%2==0)+X1*(k%%2==0)")
        m_fn <- parse(text = "1/2*(X1+X2) + X3+X4+X5+X6 -4*(k%%2==0)+X1*(k%%2==0)")
        Main(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
             honest, is_meta, is_agg, is_save, out_path, core_num)
        
    }
    
    sink()
    
}

print("Finished all!!!")
stopCluster(cores)

