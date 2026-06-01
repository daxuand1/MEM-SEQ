rm(list=ls())
suppressPackageStartupMessages({
  library(stringr)
  library(sp)
  library(spdep)
  library(dirmult) 
  library(splines)
  library(spdep)
  library(adespatial)
  library(MASS)
  library(psych)
  library(mvtnorm)
  library(vegan)
  library(Matrix)
  library(tidyr)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(parallel)
})
source("functions.R")

args = commandArgs(trailingOnly = TRUE)
fn = args[1]
filename = paste0("rdata/", fn, ".RData")
print(fn)

# ---- simulation config ----

ad = grepl("ad", fn)
bd_type = ifelse(grepl("bray", fn), "bray", "jaccard") 
p_mix = str_extract(fn, "(?<=mix)\\d{1}")
p_mix = as.numeric(ifelse(is.na(p_mix), 10, p_mix)) / 10
if(grepl("fwdseq", fn)){
  select_fun = select_fwdseq
}else if(grepl("seq", fn)){
  select_fun = select_seq
}else if(grepl("fwd", fn)){
  select_fun = select_fwd
}else if(grepl("glb", fn)){
  select_fun = select_glb
}else select_fun = NULL

# ---- parameter setup ----
a1 = 30 # top/bottom rectangle width
b1 = 5 # top/bottom rectangle height
a21 = 4 # mid rectangle normal width
a22 = 6 # mid rectangle tumor width
a23 = a1 - 2 * a21 - 2 * a22 # mid rectangle hole width
b2 = 20 # mid rectangle height

dt = rbind(cbind(x = rep(1:a1, b1),
                 y = rep(1:b1, each = a1),
                 fac = 0), # bottom
           cbind(x = rep(1:a21, b2),
                 y = rep(1:b2 + b1, each = a21),
                 fac = 0), # mid left normal
           cbind(x = rep(1:a22 + a21, b2),
                 y = rep(1:b2 + b1, each = a22),
                 fac = 1), # mid left tumor
           cbind(x = rep(1:a22 + a21 + a22 + a23, b2),
                 y = rep(1:b2 + b1, each = a22),
                 fac = 1), # mid right tumor
           cbind(x = rep(1:a21 + a21 + a22 + a23 + a22, b2),
                 y = rep(1:b2 + b1, each = a21),
                 fac = 0), # mid right normal
           cbind(x = rep(1:a1, b1),
                 y = rep(1:b1 + b1 + b2, each = a1),
                 fac = 0) # upper
)
# ggplot(dt, aes(x, y, color = factor(fac))) + geom_point()
# mean(dt[ ,3])
n = nrow(dt) 
coords = dt[, 1:2]
region = dt[, 3]

# Squared-exponential kernel:
dmat = as.matrix(dist(coords))
sigma = exp(- (dmat^2) / (2 * (2^2)))   # lengthscale 1.5

nsim = 500
seed_start = 0
p = 50
gamma = 0
gamma_alt = 1
gamma_b = 0
gamma_b_alt = ifelse(bd_type == "bray", 0.25, 0.5)
alpha0 = 1
L = 2000
nperm = 999
set.seed(0)
w = rnorm(p, mean = 0, sd = 1)
w2 = rnorm(p, mean = 0, sd = 1)

nb = dnearneigh(coords, 0, 1)
lw = nb2listw(nb)
mem = mem(lw, MEM.autocor = "positive")
# mem2 = dbmem(coords)

res_df = res_df_alt = NULL

# ---- parallel computation ----
ncores_p = 25
clp = makeCluster(ncores_p)
registerDoParallel(clp)

plist = c("dirmult", 
          "mvtnorm",
          "splines",
          "spdep",
          "adespatial",
          "MASS", 
          "psych",
          "vegan",
          "Matrix",
          "tidyr",
          "dplyr",
          "stringr")

# ---- Null case ----

cat("Null case:\n")

timestart = Sys.time()
res_df = foreach(i = 1:nsim + seed_start, .packages = plist, 
                 .combine = rbind) %dopar% 
  {
    source("functions.R")
    # i = 1
    seed = i
    out = simulate_one(region, sigma, ad,
                       bd_type, p_mix, 
                       select_fun,
                       p, alpha0, 
                       gamma, gamma_b,
                       w, w2,
                       lw, mem,
                       nperm, seed)
    df_i = data.frame(p_emp = out$p_emp,
                      num_mem = out$num_mem,
                      select_time = out$select_time,
                      analysis_time = out$analysis_time,
                      total_time = out$total_time)
    return(df_i)
  }
cat("Monte Carlo simulation time:\n")
Sys.time() - timestart

cat("P-value:\n")
summary(c(res_df$p_emp))
mean(res_df$p_emp < 0.05)  # empirical type I error estimate

# Kolmogorov-Smirnov test vs uniform
suppressWarnings({
  ks.test(res_df$p_emp, "punif")
})
cat("Number of selected MEM:\n")
summary(c(res_df$num_mem))
cat("Selection time:\n")
summary(c(res_df$select_time))
cat("Analysis time:\n")
summary(c(res_df$analysis_time))
cat("Total (selection + ANOVA) time:\n")
summary(c(res_df$total_time))

save(res_df, res_df_alt, file = filename)

# ---- Alternative case ----
cat("Alternative case:\n")

if(ad){
  gamma = gamma_alt
}else{
  gamma_b = gamma_b_alt
}

timestart = Sys.time()
res_df_alt = foreach(i = 1:nsim + seed_start, .packages = plist, 
                     .combine = rbind) %dopar% 
  {
    source("functions.R")
    seed = i
    out = simulate_one(region, sigma, ad,
                       bd_type, p_mix, 
                       select_fun,
                       p, alpha0, 
                       gamma, gamma_b,
                       w, w2,
                       lw, mem,
                       nperm, seed)
    df_i = data.frame(p_emp = out$p_emp,
                      num_mem = out$num_mem,
                      select_time = out$select_time,
                      analysis_time = out$analysis_time,
                      total_time = out$total_time)
  }
cat("Monte Carlo simulation time:\n")
Sys.time() - timestart

cat("P-value:\n")
summary(c(res_df_alt$p_emp))
mean(res_df_alt$p_emp < 0.05)  # empirical type I error estimate

# Kolmogorov-Smirnov test vs uniform
suppressWarnings({
  ks.test(res_df_alt$p_emp, "punif")
})
cat("Number of selected MEM:\n")
summary(c(res_df_alt$num_mem))
cat("Selection time:\n")
summary(c(res_df_alt$select_time))
cat("Analysis time:\n")
summary(c(res_df_alt$analysis_time))
cat("Total (selection + ANOVA) time:\n")
summary(c(res_df_alt$total_time))

save(res_df, res_df_alt, file = filename)

stopCluster(clp)

