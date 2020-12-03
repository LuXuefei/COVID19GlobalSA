# * Codes for Figure3, ICE plot
#  'Are Interventions in the COVID-19 Outbreak Really Important? A Global Sensitivity Approach', 
#   by Xuefei Lu and Emanuele Borgonovo, 2020
# 
# * Author: Xuefei Lu, xuefei.lu@ed.ac.uk
# * Date: Dec, 2020
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rm(list=ls())
library('R.matlab')
library('ICEbox')
library('ALEPlot')
library('randomForest')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sname = 'Italy'

# import Monte Carlo sample
f <- readMat('SEIRdata20200420_Italy.mat')

X <- data.frame(f$x)
ipnames <- c('alpha','beta','1/gamma', 'delta', 'I0', 'Interv. on March') 
names(X) <-ipnames
y <- f$y
remove('f')

# use subsample to train a randomForest
# number of Monte Carlo run, set to 1000 for a test, in the paper n = 10000
n = 1000
set.seed(1234)
ind <- sample(1:length(y), n)
X <- X[ind,]
#X[,6] <- X[,6] - 5 #adjust date shown in the graph
y <- y[ind]

#train randomForest
ptm <- proc.time()
seir_rf_mod = randomForest(X, y)
proc.time() - ptm

####### ICE plot, may take few seconds to run
ptm <- proc.time()
par(mfrow = c(2, 3))
par(mar = c(4.5, 4.5, 1, 1)) # margin c(bottom, left, top, right) make the plots be closer together
for (i in 1:length(X)) {
  seir.ice = ice(object = seir_rf_mod, X = X, y = y, predictor = ipnames[i], frac_to_build = .1)
  plot(seir.ice, x_quantile = FALSE, plot_pdp = TRUE, frac_to_plot = 1,
       cex.lab = 1.7,cex.axis = 1.5, 
       ylim = c(round(0.9*min(y)),round(1.1*max(y))), ylab = 'Number of total infections') #cex.axis = 1
  }
proc.time() - ptm
