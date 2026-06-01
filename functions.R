
# ---- Seq selection ----
select_seq = function(y, x, lw, mem, nperm = 999){
  
  ms = mem.select(y, lw, "positive", "global",
                  nperm.global = nperm)
  if(ms$global.test$pvalue <= 0.05){
    mod_obs = rda(y, mem)
    R2adj_full = RsquareAdj(mod_obs)$adj.r.squared
    
    f_ind = sapply(1:ncol(mem), function(i){
      mod_i = rda(y, mem[, i])
      an_i = anova(mod_i, permutations = 0)
      an_i$'F'[1]
    })
    
    ord = order(f_ind, decreasing = T)
    selected_mem = rep()
    
    # step 1.1
    step_pool = round(seq(1, ncol(mem), length.out = 10))
    R2adj_pool = rep()
    
    for(i in 1:length(step_pool)){
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i > R2adj_full) break
    }
    
    # step 1.2
    if(i > 1){
      step_pool = seq(step_pool[i - 1], step_pool[i], by = 1)
      R2adj_pool = rep()
      
      for(i in 1:length(step_pool)){
        
        mem_temp = mem[, ord[1:step_pool[i]], drop = F]
        
        R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
        R2adj_pool = c(R2adj_pool, R2adj_i)
        
        if(R2adj_i > R2adj_full) break
      }
    }
    
    mem_names_1 = colnames(mem_temp)
    mem_num_y_r2adjfull = step_pool[i]
    
    # step 2.1 - search until max R2adj
    step_pool = round(seq(mem_num_y_r2adjfull, ncol(mem), length.out = 10))
    R2adj_pool = R2adj_i
    
    for(i in 2:length(step_pool)){
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i < R2adj_pool[i - 1]) break
    }
    
    # step 2.2
    step_pool = seq(step_pool[i], 1, by = -1)
    # selected_mem = mem[, ord[1:step_pool[1]], drop = F]
    R2adj_pool = R2adj_i
    
    for(i in 2:length(step_pool)){
      
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i < R2adj_pool[i - 1]) break
    }
    
    mem_num_y_r2adjmax = step_pool[i - 1]
    mem_names_2 = 
      colnames(mem)[ord[(mem_num_y_r2adjfull + 1):mem_num_y_r2adjmax]]
    
    # check MEM for x
    
    mem_remain = mem[, !(colnames(mem) %in% mem_names_1)]
    
    mod_obs = rda(x, mem_remain)
    R2adj_full = RsquareAdj(mod_obs)$adj.r.squared
    
    f_ind = sapply(1:ncol(mem_remain), function(i){
      mod_i = rda(x, mem_remain[, i])
      an_i = anova(mod_i, permutations = 0)
      an_i$'F'[1]
    })
    
    ord = order(f_ind, decreasing = T)
    
    # step 3.1
    step_pool = round(seq(1, ncol(mem_remain), length.out = 10))
    R2adj_pool = rep()
    
    for(i in 1:length(step_pool)){
      mem_temp = mem_remain[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(x, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i > R2adj_full) break
    }
    
    # step 3.2
    if(i > 1){
      step_pool = seq(step_pool[i - 1], step_pool[i], by = 1)
      R2adj_pool = R2adj_i
      
      for(i in 1:length(step_pool)){
        
        mem_temp = mem_remain[, ord[1:step_pool[i]], drop = F]
        
        R2adj_i =  RsquareAdj(rda(x, mem_temp))$adj.r.squared
        R2adj_pool = c(R2adj_pool, R2adj_i)
        
        if(R2adj_i > R2adj_full) break
      }
    }
    
    mem_names_3 = colnames(mem_temp)
    ind = mem_names_3 %in% mem_names_2
    
    mem_names = unique(c(mem_names_1, mem_names_3[ind]))
    num_mem = length(mem_names)
  }else{
    num_mem = 0
    mem_names = colnames(mem)
  }
  
  return(list(num_mem = num_mem,
              mem_names = mem_names))
}

# ---- Forward selection + SEQ ----
select_fwdseq = function(y, x, lw, mem, nperm = 999){
  
  ms = mem.select(y, lw, "positive", "global",
                  nperm.global = nperm,
                  nperm = nperm)
  selected_mem = ms$MEM.select
  
  if(!is.null(selected_mem)){
    
    R2adj_full = RsquareAdj(rda(y, mem))$adj.r.squared
    
    f_ind = sapply(1:ncol(mem), function(i){
      mod_i = rda(y, mem[, i])
      an_i = anova(mod_i, permutations = 0)
      an_i$'F'[1]
    })
    
    ord = order(f_ind, decreasing = T)
    selected_mem = rep()
    
    # step 1.1
    step_pool = round(seq(1, ncol(mem), length.out = 10))
    R2adj_pool = rep()
    
    for(i in 1:length(step_pool)){
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i > R2adj_full) break
    }
    
    # step 1.2
    if(i > 1){
      step_pool = seq(step_pool[i - 1], step_pool[i], by = 1)
      R2adj_pool = rep()
      
      for(i in 1:length(step_pool)){
        
        mem_temp = mem[, ord[1:step_pool[i]], drop = F]
        
        R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
        R2adj_pool = c(R2adj_pool, R2adj_i)
        
        if(R2adj_i > R2adj_full) break
      }
    }
    
    mem_names_1 = colnames(mem_temp)
    mem_num_y_r2adjfull = step_pool[i]
    
    ms = forward.sel(y, mem[, mem_names_1],
                     adjR2thresh = R2adj_full,
                     nperm = nperm,
                     R2more = 0)
    
    mem_names_fwd = ms$variables
    
    # step 2.1 - search until max R2adj
    step_pool = round(seq(mem_num_y_r2adjfull, ncol(mem), length.out = 10))
    R2adj_pool = R2adj_i
    
    for(i in 2:length(step_pool)){
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i < R2adj_pool[i - 1]) break
    }
    
    # step 2.2
    step_pool = seq(step_pool[i], 1, by = -1)
    # selected_mem = mem[, ord[1:step_pool[1]], drop = F]
    R2adj_pool = R2adj_i
    
    for(i in 2:length(step_pool)){
      
      mem_temp = mem[, ord[1:step_pool[i]], drop = F]
      
      R2adj_i =  RsquareAdj(rda(y, mem_temp))$adj.r.squared
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i < R2adj_pool[i - 1]) break
    }
    
    mem_num_y_r2adjmax = step_pool[i - 1]
    mem_names_2 = colnames(mem)[ord[1:mem_num_y_r2adjmax]]
    
    # check MEM for x
    mem_remain = mem[, !(colnames(mem) %in% mem_names_fwd)]
    R2adj_full = RsquareAdj(rda(y, mem_remain))$adj.r.squared
    
    f_ind = sapply(1:ncol(mem_remain), function(i){
      mod_i = rda(x, mem_remain[, i])
      an_i = anova(mod_i, permutations = 0)
      an_i$'F'[1]
    })
    
    ord = order(f_ind, decreasing = T)
    
    # step 3.1
    step_pool = round(seq(1, ncol(mem_remain), length.out = 10))
    R2adj_pool = rep()
    
    for(i in 1:length(step_pool)){
      mem_temp = mem_remain[, ord[1:step_pool[i]], drop = F]
      
      fm = paste0("x~", 
                  paste0(colnames(mem_temp), collapse = "+"))
      mod = rda(formula = as.formula(fm), data = mem)
      R2adj_i = RsquareAdj(mod)$adj.r.squared
      
      R2adj_pool = c(R2adj_pool, R2adj_i)
      
      if(R2adj_i > R2adj_full) break
    }
    
    # step 3.2
    if(i > 1){
      step_pool = seq(step_pool[i - 1], step_pool[i], by = 1)
      R2adj_pool = R2adj_i
      
      for(i in 1:length(step_pool)){
        
        mem_temp = mem_remain[, ord[1:step_pool[i]], drop = F]
        
        fm = paste0("x~", 
                    paste0(colnames(mem_temp), collapse = "+"))
        mod = rda(formula = as.formula(fm), data = mem)
        R2adj_i = RsquareAdj(mod)$adj.r.squared
        
        R2adj_pool = c(R2adj_pool, R2adj_i)
        
        if(R2adj_i > R2adj_full) break
      }
    }
    
    mem_names_3 = colnames(mem_temp)
    ind = mem_names_3 %in% mem_names_2
    
    mem_names = unique(c(mem_names_fwd, mem_names_3[ind]))
    num_mem = length(mem_names)
    
  }else{
    num_mem = 0
    mem_names = colnames(mem)
  }
  
  return(list(num_mem = num_mem,
              mem_names = mem_names))
}


# ---- Forward selection ----
select_fwd = function(y, x, lw, mem, nperm = 999){
  
  ms = mem.select(y, lw, "positive", "FWD",
                  nperm.global = nperm,
                  nperm = nperm)
  selected_mem = ms$MEM.select
  
  if(!is.null(selected_mem)){
    num_mem = ncol(selected_mem)
    mem_names = colnames(selected_mem)
  }else{
    num_mem = 0
    mem_names = colnames(mem)
  }
  
  return(list(num_mem = num_mem,
              mem_names = mem_names))
}

# ---- Global selection ----
select_glb = function(y, x, lw, mem, nperm = 999){
  
  ms = mem.select(y, lw, "positive", "global",
                  nperm.global = nperm,
                  nperm = nperm)
  selected_mem = ms$MEM.select
  
  if(!is.null(selected_mem)){
    num_mem = ncol(selected_mem)
    mem_names = colnames(selected_mem)
  }else{
    num_mem = 0
    mem_names = colnames(mem)
  }
  
  return(list(num_mem = num_mem,
              mem_names = mem_names))
}

# ---- Single simulation ----
simulate_one = function(region, sigma, ad = T,
                        bd_type = "bray", p_mix = 1, 
                        select_fun = NULL,
                        p = 50, alpha0 = 1, 
                        gamma = 0, gamma_b = 0,
                        w = 1, w2 = 1,
                        lw = NULL, mem = NULL,
                        nperm = 999, seed = 1, ...) {
  
  # seed = 15
  set.seed(seed)
  
  # number of sample
  n = length(region)
  
  # base concentration parameters for taxa
  if(length(alpha0) == 1) alpha_vec = rep(alpha0, p) else alpha_vec = alpha0
  
  # mixture of tumor and normal tissue
  region_final = region
  if(p_mix != 1){
    region_final[region == 0] = rbinom(sum(region == 0), 1, 1 - p_mix)
    region_final[region == 1] = rbinom(sum(region == 1), 1, p_mix)
  }
  
  # latent sample effects b ~ MVN(0, sigma)
  b = as.numeric(rmvnorm(1, sigma = sigma))

  # library size
  lib_size = round(L + region_final * 100 * gamma +
                     b * 100 + rnorm(n, 0, 100))
  
  if(ad){
    y = lib_size
    y = y - mean(y)
    outcome = y
  }else{
    counts = matrix(NA, nrow = n, ncol = p)
    for(i in seq_len(n)){
      alpha_i = alpha_vec * exp(w * b[i] + gamma_b * w2 * region_final[i])
      p_i = rdirichlet(1, alpha_i)
      counts[i, ] = as.numeric(rmultinom(1, size = lib_size[i], prob = p_i))
    }
    rownames(counts) = paste0("s", seq_len(n))
    colnames(counts) = paste0("t", seq_len(p))
    
    if(bd_type == "jaccard"){
      outcome = vegdist(counts, "canberra", T)
    }else{
      outcome = vegdist(counts, bd_type)
    }
    
    # y = decostand(counts, "hellinger")
    # y = t(t(y) - colMeans(y))
    
    D = as.matrix(outcome)
    centerM = diag(n) - 1/n
    K = -0.5 * centerM %*% (D * D) %*% centerM
    eK = eigen(K, symmetric = TRUE)
    ord = order(abs(eK$values), decreasing = T)
    ind = which(cumsum(abs(eK$values)[ord]) >=
                      sum(abs(eK$values)) * 0.9)[1]
    val = sqrt(abs(eK$values[ord[1:ind]]))
    vec = eK$vectors[, ord[1:ind]]
    y = t(t(vec) * val)
  }
  
  if(is.null(select_fun)){
    fm = as.formula("outcome ~ region_final")
    num_mem = 0
    selected_mem = mem
    select_time = 0
  }else{
    timestart = Sys.time()
    ms = select_fun(y, region_final, lw, mem, nperm)
    num_mem = ms$num_mem
    
    if(num_mem == 0){
      fm = as.formula("outcome ~ region_final")
      selected_mem = mem
    }else{
      selected_mem = mem[, ms$mem_names]
      cond_terms = paste(colnames(selected_mem), collapse = " + ")
      fm = as.formula(paste("outcome ~ region_final + Condition(", 
                            cond_terms, ")", sep = ""))
    }
    select_time = difftime(Sys.time(), timestart, units = "mins")
  }
  
  timestart = Sys.time()
  if(ad){
    mod_obs = rda(fm, data = selected_mem)
  }else{
    mod_obs = dbrda(fm, data = selected_mem)
  }
  mod_obs = dbrda(fm, data = selected_mem)
  an_obs = anova(mod_obs, permutations = nperm)
  # an_obs = anova(mod_obs, permutations = 0)
  # an_obs
  p_emp = an_obs$`Pr(>F)`[1]
  p_emp = ifelse(is.na(p_emp), 1, p_emp)
  
  analysis_time = difftime(Sys.time(), timestart, units = "mins")
  
  total_time = select_time + analysis_time
  
  list(p_emp = p_emp,
       num_mem = num_mem,
       select_time = select_time,
       analysis_time = analysis_time,
       total_time = total_time)
}
