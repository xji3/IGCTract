# This R file is for Jeff to play with the data
rm(list=ls())  # clean up workspace
# First, load in RData file
load("./TractSummary.RData")

# Show the estimated tract length and IGC initiation rates.
# Estimated Tract length (unit: nucleotide)
seq.length <- JS.HKY.nonclock.summary["length", ]
PSJS.HKY.dim.1.nonclock.eff.lnL <- PSJS.HKY.dim.1.nonclock.summary["ll", ] / (seq.length - 1)
PSJS.HKY.dim.2.nonclock.eff.lnL <- PSJS.HKY.dim.2.nonclock.summary["ll", ] / (seq.length - 1)
PSJS.HKY.rv.NOSC.dim.1.nonclock.eff.lnL <- PSJS.HKY.rv.NOSC.dim.1.nonclock.summary["ll", ]/ (seq.length - 3)
PSJS.HKY.rv.NOSC.dim.2.nonclock.eff.lnL <- PSJS.HKY.rv.NOSC.dim.2.nonclock.summary["ll", ]/ (seq.length - 3)
PSJS.HKY.rv.SCOK.dim.1.nonclock.eff.lnL <- PSJS.HKY.rv.SCOK.dim.1.nonclock.summary["ll", ]/ (seq.length - 1)
PSJS.HKY.rv.SCOK.dim.2.nonclock.eff.lnL <- PSJS.HKY.rv.SCOK.dim.2.nonclock.summary["ll", ]/ (seq.length - 1)

show.mat <- rbind(PSJS.HKY.dim.1.nonclock.summary["tract_length", ],
                  PSJS.HKY.dim.2.nonclock.summary["tract_length", ],
                  PSJS.HKY.dim.1.nonclock.eff.lnL - PSJS.HKY.dim.2.nonclock.eff.lnL,
                  PSJS.HKY.rv.NOSC.dim.1.nonclock.summary["tract_length",],
                  PSJS.HKY.rv.NOSC.dim.2.nonclock.summary["tract_length",],
                  PSJS.HKY.rv.NOSC.dim.1.nonclock.eff.lnL - PSJS.HKY.rv.NOSC.dim.2.nonclock.eff.lnL,
                  PSJS.HKY.rv.SCOK.dim.1.nonclock.summary["tract_length",],
                  PSJS.HKY.rv.SCOK.dim.2.nonclock.summary["tract_length",],
                  PSJS.HKY.rv.SCOK.dim.1.nonclock.eff.lnL - PSJS.HKY.rv.SCOK.dim.2.nonclock.eff.lnL
)
row.names(show.mat) <- c("Homo D1", "Homo D2", "lnL (D1 - D2)" ,
                         "Heter NOSC D1", "Heter NOSC D2", "lnL (D1 - D2)", 
                         "Heter SCOK D1", "Heter SCOK D2", "lnL (D1 - D2)"
)
show.mat

# Now show equivalent lnL matrix
show.mat <- rbind(JS.HKY.nonclock.summary["ll", ], 
                  PSJS.HKY.dim.1.nonclock.eff.lnL, PSJS.HKY.dim.2.nonclock.eff.lnL, 
                  JS.HKY.rv.nonclock.summary["ll", ], 
                  PSJS.HKY.rv.NOSC.dim.1.nonclock.eff.lnL, PSJS.HKY.rv.NOSC.dim.2.nonclock.eff.lnL, 
                  PSJS.HKY.rv.SCOK.dim.1.nonclock.eff.lnL, PSJS.HKY.rv.SCOK.dim.1.nonclock.eff.lnL)
row.names(show.mat) <- c("Homo JS", "Homo PSJS D1", "Homo PSJS D2", 
                         "Heter JS", 
                         "Heter NOSC D1", "Heter NOSC D2", "Heter SCOK D1", "Heter SCOK D2")
show.mat
# Compare estimated Tau value
# Estimated Tract length (unit: nucleotide) 
show.mat <- rbind(
  JS.HKY.nonclock.summary["Tau",],
  PSJS.HKY.dim.2.nonclock.summary["tract_length", ] * PSJS.HKY.dim.2.nonclock.summary["init_rate",],
  JS.HKY.rv.nonclock.summary["Tau",] * 3.0 / colSums(rbind(1, JS.HKY.rv.nonclock.summary[c("r2", "r3"),])),
  (PSJS.HKY.rv.NOSC.dim.2.nonclock.summary["tract_length",] * PSJS.HKY.rv.NOSC.dim.2.nonclock.summary["init_rate", ]
   * 3.0 / colSums(rbind(1,  PSJS.HKY.rv.NOSC.dim.2.nonclock.summary[c("r2", "r3"), ]))), 
  PSJS.HKY.rv.SCOK.dim.2.nonclock.summary["tract_length",] * PSJS.HKY.rv.SCOK.dim.2.nonclock.summary["init_rate", ]
  * 3.0 / colSums(rbind(1, PSJS.HKY.rv.SCOK.dim.2.nonclock.summary[c("r2", "r3"), ]))
)
row.names(show.mat) <- c("Homo JS Tau", "Homo PSJS D2 Tau", "Heter JS Tau", "Heter PSJS NOSC Tau", "Heter PSJS SCOK Tau")
show.mat

##########################################################
##########################################################
##########################################################

#Plot lnL increase for pairs
library(fields)
for (pair in finished.pairs){
  for (dim in 1:2){
    # This code is for PSJS RV SCOK model only
    
    JS.rv.lnL <- get(paste(pair, "JS_HKY_rv_nonclock_lnL", sep = "_"))
    row.lnL <- JS.rv.lnL[, 2] %*% matrix(1, 1, dim(JS.rv.lnL)[1])
    col.lnL <- t(row.lnL)
    JS.rv.mat <- row.lnL + col.lnL
    diag(JS.rv.mat) <- 0.0
    
    PSJS.rv.lnL <- get(paste(pair, "PSJS_dim", toString(dim), "HKY_rv_SCOK_nonclock_lnL", sep = "_"))
    PSJS.rv.mat <- matrix(0, dim(JS.rv.lnL)[1], dim(JS.rv.lnL)[1])
    for(i in 1:dim(PSJS.rv.lnL)[1]){
      PSJS.rv.mat[PSJS.rv.lnL[i, 1] + 1, PSJS.rv.lnL[i, 2] + 1] <- PSJS.rv.lnL[i, 3]
    }
    PSJS.rv.mat <- PSJS.rv.mat + t(PSJS.rv.mat)
    
    diff.lnL.mat <- PSJS.rv.mat - JS.rv.mat
    to.plot.mat <- diff.lnL.mat
    to.plot.mat[to.plot.mat < 0.01] <- 0
    brk = quantile( c(diff.lnL.mat), c(0., 0.98, 1.0))
    #lab.brk = paste(names(brk), round(brk, digits = 2), sep = ":")
    lab.brk = names(brk)
    image.plot(1:dim(PSJS.rv.mat)[1], 1:dim(PSJS.rv.mat)[2], diff.lnL.mat, breaks = brk, col = rainbow(2), 
               lab.breaks = lab.brk, main = paste(pair, "PSJS_dim", toString(dim), "HKY_rv_SCOK", sep = "_"))
    print(brk)
    #hist(diff.lnL.mat)
    #write.table(diff.lnL.mat, paste("./", pair, "_matlab_test.txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}