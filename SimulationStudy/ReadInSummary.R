Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
# First read in HMM results
# from summary file
for(tract in Tract.list){
  hmm.tract.summary <- NULL
  for(sim in 1:100){
    hmm.summary <- paste("./summary/Tract_", toString(tract), '.0/sim_', 
                         toString(sim), '/HMM_YDR418W_YEL054C_MG94_nonclock_sim_',
                         toString(sim), '_1D_summary.txt', sep = "")
    if (file.exists(hmm.summary)){
      all <- readLines(hmm.summary, n = -1)
      col.names <- paste("sim_", toString(sim), sep = "")
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- as.matrix(read.table(hmm.summary, 
                                          row.names = row.names, 
                                          col.names = col.names))
      hmm.tract.summary <- cbind(hmm.tract.summary, summary_mat)     
    }
    
  }
  assign(paste("HMM_Tract_", toString(tract), "_summary", sep = ""), hmm.tract.summary)
}

# from plots
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/Tract_", toString(tract), '.0/sim_', 
                      toString(sim), '/HMM_YDR418W_YEL054C_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in actual mean tract length in each simulated dataset
for (tract in Tract.list){
  sim.tract <- NULL
  for(sim in 1:100){
    sim_log <- paste("./Tract_", toString(tract), ".0_HKY/sim_", toString(sim), 
                     "/YDR418W_YEL054C_sim_", toString(sim), "_IGC.log", sep = "")
    # now read in log file
    log_info <- read.table(sim_log, header = TRUE)
    realized.tract.length <- log_info[, "stop_pos"] - log_info[, "start_pos"] + 1
    potential.tract.length <- log_info[, "tract_length"]
    # Now get longest length of variant subtract of each tract
    diff.tracts <- log_info[, "num_diff"] > 0
    if (sum(diff.tracts)){
      subtract.length.list <- NULL
      for(row.num in (1:dim(log_info)[1])[diff.tracts]){
        donor.seq <- strsplit(toString(log_info[row.num, "template_seq"]), "")[[1]]
        recipient.seq <- strsplit(toString(log_info[row.num, "overide_seq"]), "")[[1]]
        first.pos <- FALSE
        for(seq.pos in 1:length(donor.seq)){
          if(donor.seq[seq.pos] != recipient.seq[seq.pos]){
            last.pos <- seq.pos
            if(!first.pos){
              first.pos <- seq.pos
            }
          }
        }
        subtract.length.list <- c(subtract.length.list, last.pos - first.pos + 1)
      }
    }
    new.info <- matrix(c(dim(log_info)[1], mean(potential.tract.length), sd(potential.tract.length), 
                         mean(realized.tract.length), sd(realized.tract.length),
                         mean(potential.tract.length[diff.tracts]), 
                         mean(realized.tract.length[diff.tracts]),
                         mean(subtract.length.list),
                         sum(log_info[, "num_diff"] > 1)),
                       9, 1)
    rownames(new.info) <- c("num IGC", "mean potential tract length", "sd potential tract length", 
                            "mean realized tract length", "sd realized tract length", 
                            "mean potential nonidentical tract length", "mean realized nonidentical tract length",
                            "mean subtract length", "num IGC with at least two variant sites")
    colnames(new.info) <- paste("sim_", toString(sim), sep = "")
    sim.tract <- cbind(sim.tract, new.info)
  }
  assign(paste("sim.tract.", toString(tract), sep = ""), sim.tract)
}

guess.list <- c(50.0, 100.0, 250.0, 500.0)
# Now read in PSJS summary results 
for(tract in Tract.list){
  for(guess in guess.list){
    PSJS.tract.summary <- NULL
    for(sim in 1:100){
      
      PSJS.summary <- paste("./summary/Tract_", toString(tract), '.0_HKY/sim_', 
                            toString(sim), '/PSJS_HKY_rv_sim_',
                            toString(sim), "_Tract_", toString(tract), '.0_guess_',
                            toString(guess),'.0_nt_summary.txt', sep = "")
      if (file.exists(PSJS.summary)){
        all <- readLines(PSJS.summary, n = -1)
        col.names <- paste("sim_", toString(sim), sep = "")
        row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
        summary_mat <- as.matrix(read.table(PSJS.summary, 
                                            row.names = row.names, 
                                            col.names = col.names))
        PSJS.tract.summary <- cbind(PSJS.tract.summary, summary_mat)
      }
    }
    assign(paste("PSJS_HKY_Tract_", toString(tract), "_guess_", 
                 toString(guess), "_summary", sep = ""), PSJS.tract.summary)
  }
}

# Now combine all initial guess results
for(tract in Tract.list){
  combined.PSJS.tract.summary <- NULL
  col.list <- NULL
  for ( sim_num in 1:100){
    sim_col <- paste("sim_", toString(sim_num), sep = "")
    best.lnL <- -Inf
    best.guess <- NULL
    for(guess in guess.list){
      target_summary <- get(paste("PSJS_HKY_Tract_", toString(tract), "_guess_", toString(guess), "_summary", sep = ""))
      if(sim_col %in% colnames(target_summary) ){
        if (target_summary["ll", sim_col] > best.lnL){
          best.lnL <- target_summary["ll", sim_col]
          best.guess <- guess        
        }
      }
    }
    if(! is.null(best.guess)){
      combined.PSJS.tract.summary <- cbind(combined.PSJS.tract.summary, 
                                           get(paste("PSJS_HKY_Tract_", toString(tract), "_guess_", toString(best.guess), "_summary", sep = ""))[, sim_col]) 
      col.list <- c(col.list, sim_col)
    }
    
  }
  colnames(combined.PSJS.tract.summary) <- col.list
  assign(paste("PSJS_HKY_Tract_", toString(tract), "_combined_summary", sep = ""), combined.PSJS.tract.summary)
}

# Now read in PSJS summary results of Old Simulated Datasets 
for(tract in Tract.list){
  for(guess in guess.list){
    PSJS.tract.summary <- NULL
    for(sim in 1:100){
      
      PSJS.summary <- paste("./summary/IGCgeo_", toString(tract), '.0/sim_', 
                            toString(sim), '/PSJS_HKY_rv_sim_',
                            toString(sim), "_IGCgeo_", toString(tract), '.0_guess_',
                            toString(guess),'.0_summary.txt', sep = "")
      if (file.exists(PSJS.summary)){
        all <- readLines(PSJS.summary, n = -1)
        col.names <- paste("sim_", toString(sim), sep = "")
        row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
        summary_mat <- as.matrix(read.table(PSJS.summary, 
                                            row.names = row.names, 
                                            col.names = col.names))
        PSJS.tract.summary <- cbind(PSJS.tract.summary, summary_mat)
      }
    }
    assign(paste("Old_PSJS_HKY_Tract_", toString(tract), "_guess_", 
                 toString(guess), "_summary", sep = ""), PSJS.tract.summary)
  }
}

# Now combine all initial guess results using old simulated data
for(tract in Tract.list){
  combined.PSJS.tract.summary <- NULL
  col.list <- NULL
  for ( sim_num in 1:100){
    sim_col <- paste("sim_", toString(sim_num), sep = "")
    best.lnL <- -Inf
    best.guess <- NULL
    for(guess in guess.list){
      target_summary <- get(paste("Old_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(guess), "_summary", sep = ""))
      if(sim_col %in% colnames(target_summary) ){
        if (target_summary["ll", sim_col] > best.lnL){
          best.lnL <- target_summary["ll", sim_col]
          best.guess <- guess        
        }
      }
    }
    if(! is.null(best.guess)){
      combined.PSJS.tract.summary <- cbind(combined.PSJS.tract.summary, 
                                           get(paste("Old_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(best.guess), "_summary", sep = ""))[, sim_col]) 
      col.list <- c(col.list, sim_col)
    }
    
  }
  colnames(combined.PSJS.tract.summary) <- col.list
  assign(paste("Old_PSJS_HKY_Tract_", toString(tract), "_combined_summary", sep = ""), combined.PSJS.tract.summary)
}

# Now read in PSJS summary results for HalfTau simulated datasets 
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 300.0)
for(tract in Tract.list){
  for(guess in guess.list){
    PSJS.tract.summary <- NULL
    for(sim in 1:100){
      
      PSJS.summary <- paste("./summary/HalfTau/Tract_", toString(tract), '.0_HKY/sim_', 
                            toString(sim), '/PSJS_HKY_rv_sim_',
                            toString(sim), "_Tract_", toString(tract), '.0_guess_',
                            toString(guess),'.0_nt_summary.txt', sep = "")
      if (file.exists(PSJS.summary)){
        all <- readLines(PSJS.summary, n = -1)
        col.names <- paste("sim_", toString(sim), sep = "")
        row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
        summary_mat <- as.matrix(read.table(PSJS.summary, 
                                            row.names = row.names, 
                                            col.names = col.names))
        PSJS.tract.summary <- cbind(PSJS.tract.summary, summary_mat)
      }
    }
    assign(paste("HalfTau_PSJS_HKY_Tract_", toString(tract), "_guess_", 
                 toString(guess), "_summary", sep = ""), PSJS.tract.summary)
  }
}

# Now combine all initial guess results for HalfTau simulated datasets 
for(tract in Tract.list){
  combined.PSJS.tract.summary <- NULL
  col.list <- NULL
  for ( sim_num in 1:100){
    sim_col <- paste("sim_", toString(sim_num), sep = "")
    best.lnL <- -Inf
    best.guess <- NULL
    for(guess in guess.list){
      target_summary <- get(paste("HalfTau_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(guess), "_summary", sep = ""))
      if(sim_col %in% colnames(target_summary) ){
        if (target_summary["ll", sim_col] > best.lnL){
          best.lnL <- target_summary["ll", sim_col]
          best.guess <- guess        
        }
      }
    }
    if(! is.null(best.guess)){
      combined.PSJS.tract.summary <- cbind(combined.PSJS.tract.summary, 
                                           get(paste("HalfTau_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(best.guess), "_summary", sep = ""))[, sim_col]) 
      col.list <- c(col.list, sim_col)
    }
    
  }
  colnames(combined.PSJS.tract.summary) <- col.list
  assign(paste("HalfTau_PSJS_HKY_Tract_", toString(tract), "_combined_summary", sep = ""), combined.PSJS.tract.summary)
}

# Now read in actual mean tract length in each simulated dataset
for (tract in Tract.list){
  sim.tract <- NULL
  for(sim in 1:100){
    sim_log <- paste("./HalfTau/Tract_", toString(tract), ".0_HKY/sim_", toString(sim), 
                     "/YDR418W_YEL054C_sim_", toString(sim), "_IGC.log", sep = "")
    # now read in log file
    log_info <- read.table(sim_log, header = TRUE)
    realized.tract.length <- log_info[, "stop_pos"] - log_info[, "start_pos"] + 1
    potential.tract.length <- log_info[, "tract_length"]
    # Now get longest length of variant subtract of each tract
    diff.tracts <- log_info[, "num_diff"] > 0
    if (sum(diff.tracts)){
      subtract.length.list <- NULL
      for(row.num in (1:dim(log_info)[1])[diff.tracts]){
        donor.seq <- strsplit(toString(log_info[row.num, "template_seq"]), "")[[1]]
        recipient.seq <- strsplit(toString(log_info[row.num, "overide_seq"]), "")[[1]]
        first.pos <- FALSE
        for(seq.pos in 1:length(donor.seq)){
          if(donor.seq[seq.pos] != recipient.seq[seq.pos]){
            last.pos <- seq.pos
            if(!first.pos){
              first.pos <- seq.pos
            }
          }
        }
        subtract.length.list <- c(subtract.length.list, last.pos - first.pos + 1)
      }
    }
    new.info <- matrix(c(dim(log_info)[1], mean(potential.tract.length), sd(potential.tract.length), 
                         mean(realized.tract.length), sd(realized.tract.length),
                         mean(potential.tract.length[diff.tracts]), 
                         mean(realized.tract.length[diff.tracts]),
                         mean(subtract.length.list),
                         sum(log_info[, "num_diff"] > 1)),
                       9, 1)
    rownames(new.info) <- c("num IGC", "mean potential tract length", "sd potential tract length", 
                            "mean realized tract length", "sd realized tract length", 
                            "mean potential nonidentical tract length", "mean realized nonidentical tract length",
                            "mean subtract length", "num IGC with at least two variant sites")
    colnames(new.info) <- paste("sim_", toString(sim), sep = "")
    sim.tract <- cbind(sim.tract, new.info)
  }
  assign(paste("HalfTau.sim.tract.", toString(tract), sep = ""), sim.tract)
}

# Now read in PSJS summary results for TenthTau simulated datasets 
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 300.0)
for(tract in Tract.list){
  for(guess in guess.list){
    PSJS.tract.summary <- NULL
    for(sim in 1:100){
      
      PSJS.summary <- paste("./summary/TenthTau/Tract_", toString(tract), '.0_HKY/sim_', 
                            toString(sim), '/PSJS_HKY_rv_sim_',
                            toString(sim), "_Tract_", toString(tract), '.0_guess_',
                            toString(guess),'.0_nt_summary.txt', sep = "")
      if (file.exists(PSJS.summary)){
        all <- readLines(PSJS.summary, n = -1)
        col.names <- paste("sim_", toString(sim), sep = "")
        row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
        summary_mat <- as.matrix(read.table(PSJS.summary, 
                                            row.names = row.names, 
                                            col.names = col.names))
        PSJS.tract.summary <- cbind(PSJS.tract.summary, summary_mat)
      }
    }
    assign(paste("TenthTau_PSJS_HKY_Tract_", toString(tract), "_guess_", 
                 toString(guess), "_summary", sep = ""), PSJS.tract.summary)
  }
}

# Now combine all initial guess results for TenthTau simulated datasets 
for(tract in Tract.list){
  combined.PSJS.tract.summary <- NULL
  col.list <- NULL
  for ( sim_num in 1:100){
    sim_col <- paste("sim_", toString(sim_num), sep = "")
    best.lnL <- -Inf
    best.guess <- NULL
    for(guess in guess.list){
      target_summary <- get(paste("TenthTau_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(guess), "_summary", sep = ""))
      if(sim_col %in% colnames(target_summary) ){
        if (target_summary["ll", sim_col] > best.lnL){
          best.lnL <- target_summary["ll", sim_col]
          best.guess <- guess        
        }
      }
    }
    if(! is.null(best.guess)){
      combined.PSJS.tract.summary <- cbind(combined.PSJS.tract.summary, 
                                           get(paste("TenthTau_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(best.guess), "_summary", sep = ""))[, sim_col]) 
      col.list <- c(col.list, sim_col)
    }
    
  }
  colnames(combined.PSJS.tract.summary) <- col.list
  assign(paste("TenthTau_PSJS_HKY_Tract_", toString(tract), "_combined_summary", sep = ""), combined.PSJS.tract.summary)
}

# Now read in actual mean tract length in each simulated dataset
for (tract in Tract.list){
  sim.tract <- NULL
  for(sim in 1:100){
    sim_log <- paste("./TenthTau/Tract_", toString(tract), ".0_HKY/sim_", toString(sim), 
                     "/YDR418W_YEL054C_sim_", toString(sim), "_IGC.log", sep = "")
    # now read in log file
    log_info <- read.table(sim_log, header = TRUE)
    realized.tract.length <- log_info[, "stop_pos"] - log_info[, "start_pos"] + 1
    potential.tract.length <- log_info[, "tract_length"]
    # Now get longest length of variant subtract of each tract
    diff.tracts <- log_info[, "num_diff"] > 0
    if (sum(diff.tracts)){
      subtract.length.list <- NULL
      for(row.num in (1:dim(log_info)[1])[diff.tracts]){
        donor.seq <- strsplit(toString(log_info[row.num, "template_seq"]), "")[[1]]
        recipient.seq <- strsplit(toString(log_info[row.num, "overide_seq"]), "")[[1]]
        first.pos <- FALSE
        for(seq.pos in 1:length(donor.seq)){
          if(donor.seq[seq.pos] != recipient.seq[seq.pos]){
            last.pos <- seq.pos
            if(!first.pos){
              first.pos <- seq.pos
            }
          }
        }
        subtract.length.list <- c(subtract.length.list, last.pos - first.pos + 1)
      }
    }
    new.info <- matrix(c(dim(log_info)[1], mean(potential.tract.length), sd(potential.tract.length), 
                         mean(realized.tract.length), sd(realized.tract.length),
                         mean(potential.tract.length[diff.tracts]), 
                         mean(realized.tract.length[diff.tracts]),
                         mean(subtract.length.list),
                         sum(log_info[, "num_diff"] > 1)),
                       9, 1)
    rownames(new.info) <- c("num IGC", "mean potential tract length", "sd potential tract length", 
                            "mean realized tract length", "sd realized tract length", 
                            "mean potential nonidentical tract length", "mean realized nonidentical tract length",
                            "mean subtract length", "num IGC with at least two variant sites")
    colnames(new.info) <- paste("sim_", toString(sim), sep = "")
    sim.tract <- cbind(sim.tract, new.info)
  }
  assign(paste("TenthTau.sim.tract.", toString(tract), sep = ""), sim.tract)
}

# Now read in HMM results from plots
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_YDR418W_YEL054C_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in HMM results from plots for HalfTau case
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/HalfTau/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_YDR418W_YEL054C_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("HalfTau_HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in HMM results from plots for TenthTau case
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/TenthTau/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_YDR418W_YEL054C_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("TenthTau_HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in PSJS summary results for simulated datasets with each half for bias correction
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 300.0)
half.list <- c("first_half", "second_half")
for(tract in Tract.list){
  for(guess in guess.list){
    for(half_case in half.list){
      PSJS.tract.summary <- NULL
      for(sim in 1:100){
        
        PSJS.summary <- paste("./summary/Tract_", toString(tract), '.0_HKY/sim_', 
                              toString(sim), '/PSJS_HKY_rv_sim_',
                              toString(sim), "_Tract_", toString(tract), '.0_guess_',
                              toString(guess),'.0_nt_summary_',half_case,'.txt', sep = "")
        if (file.exists(PSJS.summary)){
          all <- readLines(PSJS.summary, n = -1)
          col.names <- paste("sim_", toString(sim), sep = "")
          row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
          summary_mat <- as.matrix(read.table(PSJS.summary, 
                                              row.names = row.names, 
                                              col.names = col.names))
          PSJS.tract.summary <- cbind(PSJS.tract.summary, summary_mat)
        }
      }
      assign(paste(half_case, "_PSJS_HKY_Tract_", toString(tract), "_guess_", 
                   toString(guess), "_summary", sep = ""), PSJS.tract.summary)
    }      
  }
}

# Now combine all initial guess results for TenthTau simulated datasets 
for(tract in Tract.list){
  for(half_case in half.list){
    combined.PSJS.tract.summary <- NULL
    col.list <- NULL
    for ( sim_num in 1:100){
      sim_col <- paste("sim_", toString(sim_num), sep = "")
      best.lnL <- -Inf
      best.guess <- NULL
      for(guess in guess.list){
        target_summary <- get(paste(half_case, "_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(guess), "_summary", sep = ""))
        if(sim_col %in% colnames(target_summary) ){
          if (target_summary["ll", sim_col] > best.lnL){
            best.lnL <- target_summary["ll", sim_col]
            best.guess <- guess        
          }
        }
      }
      if(! is.null(best.guess)){
        combined.PSJS.tract.summary <- cbind(combined.PSJS.tract.summary, 
                                             get(paste(half_case, "_PSJS_HKY_Tract_", toString(tract), "_guess_", toString(best.guess), "_summary", sep = ""))[, sim_col]) 
        col.list <- c(col.list, sim_col)
      }
      
    }
    colnames(combined.PSJS.tract.summary) <- col.list
    assign(paste(half_case, "_PSJS_HKY_Tract_", toString(tract), "_combined_summary", sep = ""), combined.PSJS.tract.summary)
  }
}