# Set working directory
setwd("your/path")

library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
library(mclust)
library(ggplot2)
library(mvtnorm)
library(ggpubr)
library(mixtools)
library(RColorBrewer)
library(fossil)
library(nwREM)

num = 12329
set.seed(num)

# Set number of dimensions
d = 2
# Set total number of observed points
n = 600
# Set number of clusters
k = 3
# Set percent outliers in observed data
perc_out = 15

# Simulate data
covs = createCovs(k, d)
mus = createMus(k, d)
mix_prop = createMixProp(k)
simul_mvGauss = simMultiGauss(n, mix_prop, mus, covs, perc_out)
sampleMat = simul_mvGauss$x

sampleMat1 = as.data.frame(cbind(sampleMat, simul_mvGauss$z))
colnames(sampleMat1) = c(1:d, "Cluster")
sampleMat1$Cluster = as.character(sampleMat1$Cluster)

################ Run both clustering methods on sample #################
sim_kmeans = kmeans(x = sampleMat, centers = k, iter.max = 10, nstart = 10)

sim_EMfit_stand = em.mv.gmm(x=sampleMat, k=k)

sim_EMfit_robust = em_alg_GMM(x=sampleMat, k=k, lambda = 10, max_it = 100, tol = 1e-8)


# calculate the Rand index for each
true = as.numeric(sampleMat1$Cluster)
rand_ind_robust = rand.index(true, sim_EMfit_robust$classification)
rand_ind_stand = rand.index(true, sim_EMfit_stand$classification)
rand_ind_Mclust = rand.index(true, sim_EMfit_Mclust$classification)

##### Plot each of these!
plot.new()
ellipse0 = lapply(1:k, function(j) mixtools::ellipse(mu=mus[,j], 
                                                     sigma=covs[,,j], npoints = 300, newplot=F))
ellipse1 = as.data.frame(ellipse0[[1]]); colnames(ellipse1) = c("X1","X2")
ellipse2 = as.data.frame(ellipse0[[2]]); colnames(ellipse2) = c("X1","X2")
ellipse3 = as.data.frame(ellipse0[[3]]); colnames(ellipse3) = c("X1","X2")
sim_mu = as.data.frame(t(mus)); colnames(sim_mu) = c("V1", "V2")

sampleMat1 = as.data.frame(sampleMat1)
mypalette = brewer.pal(n=8,name="Dark2")
orig_plot = ggplot(sampleMat1, aes(x=`1`, y=`2`, colour=Cluster)) + geom_point() + 
  scale_colour_brewer(palette="Dark2") + 
  geom_point(data = sim_mu, aes(x=sim_mu$V1, y=sim_mu$V2), color="black") +
  geom_point(data=ellipse1, aes(x=ellipse1$X1, y=ellipse1$X2), color = mypalette[1], size=0.3) +
  geom_point(data=ellipse2, aes(x=ellipse2$X1, y=ellipse2$X2), color = mypalette[2], size=0.3) +
  geom_point(data=ellipse3, aes(x=ellipse3$X1, y=ellipse3$X2), color = mypalette[3], size=0.3) +
  ggtitle("Simulated Data")
#dev.off()

stand_plot = plot_standard_EM(sim_EMfit_stand, sampleMat1)
#dev.off()
robust_plot = plot_robust_EM(sim_EMfit_robust, sampleMat1)
#dev.off()
kmeans_plot = plot_kmeans(sim_kmeans, sampleMat1)
#dev.off()

pdf(file="your/path/sim_visual_comp.pdf", width = 10, height = 10)
ggarrange(orig_plot, stand_plot, robust_plot, kmeans_plot, ncol=2, nrow=2,
          legend = "bottom", common.legend = T)
dev.off()