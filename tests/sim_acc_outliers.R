# set Working directory
setwd("your/path")

library(MASS)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(ggpubr)
library(mixtools)
library(RColorBrewer)
library(fossil)
library(dplyr)
library(mclust)
library(tidyr)
library(nwREM)

num = 12700
set.seed(num)

# Set some of the variables
d = 2                   # Number of dimensions
n = 500                 # Total number of observations
k = 3                   # Number of clusters
rand_outliers = NULL    # Initialize this list
num_sim = 100           # Number of simulations per percent outlier
start = Sys.time()

# We're running different percent outliers - start with 5% 
# Also running with lambda = 10
for (i in 1:5) { # for 1:5 we range from 5% to 25%

  # Set the percent outliers
  perc_out = (i-1)*5
  print(i-1)
  
  for (j in 1:num_sim) {
    # for each k, we want the repeated samples to be the same
    seed_num = num + j
    set.seed(seed_num)
    # print(seed_num)         # Print seed number if preferred
    
    ################ Create the simulated data ################
    covs = createCovs(k, d)
    mus = createMus(k, d)
    mix_prop = createMixProp(k)
    simul_mvGauss = simMultiGauss(n, mix_prop, mus, covs, perc_out)
    sampleMat = simul_mvGauss$x
    
    sampleMat1 = as.data.frame(cbind(sampleMat, simul_mvGauss$z))
    colnames(sampleMat1) = c(1:d, "Cluster")
    sampleMat1$Cluster = as.character(sampleMat1$Cluster)
    
    ################ Run both clustering algorithms on sample #################
    # Run Mclust - should be similar to the standard EM results
    sim_EMfit_Mclust = Mclust(sampleMat, verbose=F, G=k, modelNames = "VVV", 
                              initialization = list(hcPairs = hc(sampleMat)), 
                              control = emControl(tol=c(1.e-5, 1.e-6), itmax=30))
    # Run standard EM - no penalty included
    sim_EMfit_stand = em.mv.gmm(x=sampleMat, k=k)
    # Run robust EM - L0 penalty
    sim_EMfit_robust = em_alg_GMM(x=sampleMat, k=k, lambda = 10, max_it = 100, tol = 1e-8)

    # calculate the Rand index for each
    true = as.numeric(sampleMat1$Cluster)
    rand_ind_robust = rand.index(true, sim_EMfit_robust$classification)
    rand_ind_stand = rand.index(true, sim_EMfit_stand$classification)
    rand_ind_Mclust = rand.index(true, sim_EMfit_Mclust$classification)
    
    curr_stand = c(perc_out, seed_num, "Standard", rand_ind_stand)
    curr_rob = c(perc_out, seed_num, "Robust", rand_ind_robust)
    curr_Mclust = c(perc_out, seed_num, "Mclust", rand_ind_Mclust)
    
    rand_outliers = rbind(rand_outliers, curr_stand, curr_rob, curr_Mclust)
  }
}
# Make into a data frame
rand_outliers = as.data.frame(rand_outliers)
colnames(rand_outliers) = c("Perc_outliers", "Seed_number", "Method", "Rand_index")
rand_outliers = rand_outliers %>% mutate(Perc_outliers = as.numeric(as.character(Perc_outliers)), 
                                         Seed_number = as.numeric(as.character(Seed_number)),
                                         Rand_index = as.numeric(as.character(Rand_index))) #%>%

# Create data frame for easy plotting
rand_outliers_mean = rand_outliers %>% group_by(Perc_outliers, Method) %>% 
  summarise(mean=mean(Rand_index), 
            lower = mean - 1.96*sd(Rand_index)/num_sim,
            upper = min(mean + 1.96*sd(Rand_index)/num_sim,1))
end=Sys.time()
run_time = end-start

## Plot the accuracy over the percent outliers
out_plot = ggplot(rand_outliers_mean, aes(x = Perc_outliers, y = mean, colour = Method)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, colour = NULL, group = Method), 
                width=2, position = position_dodge(0.6)) +
  geom_line(position = position_dodge(0.6), stat = "identity") + 
  geom_point(position = position_dodge(0.6), size=2.5) +
  labs(y = "Rand Index", x = "Percent Outliers (%)", colour = "Method") +
  ggtitle(paste("Clustering accuracy over percent outliers")) +
  scale_color_brewer(palette="Dark2") +
  theme(plot.title = element_text(hjust = 0.5))

out_plot

print(run_time)

