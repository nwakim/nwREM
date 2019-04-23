#' Function to plot the cluster assignment from the robust EM algorithm 
#'
#' @param sim_EMfit Output from the EM algorithm
#' @param sampleMat An (n x d+1) matrix with the observed data and true cluster assignments
#' @return A plot containing the observed data, colored by their clusters according to the EM algorithm. Estimated covariance and mean of clusters also plotted.
#' @export
plot_robust_EM = function(sim_EMfit, sampleMat){
  
  swap <- function(vec, from, to) {
    tmp <- to[ match(vec, from) ]
    tmp[is.na(tmp)] <- vec[is.na(tmp)]
    return(tmp)
  }
  # Assign values from returned list to mu, sigma, T_mat, tau, e
  mu = as.data.frame(t(sim_EMfit$mus))
  sigma = sim_EMfit$covs
  #error = sim_EMfit$e
  c = nrow(mu)
  d = ncol(mu)
  n = nrow(sim_EMfit$probs)
  lambda = 10
  
  
  # Assign the cluster number to whichever one has the highest probability
  hard_assign = sim_EMfit$classification
  hard_assign = swap(hard_assign, c(1,2), c(2,1))
  
  sampleMat2 = as.data.frame(sampleMat)
  sampleMat2$Cluster = as.factor(hard_assign)
  
  
  plot.new()
  ellipse0 = lapply(1:3, function(j) mixtools::ellipse(mu=sim_EMfit$mus[,j], 
                                                       sigma=sim_EMfit$covs[,,j], npoints = 300, newplot=F))
  ellipse1 = as.data.frame(ellipse0[1]); colnames(ellipse1) = c("X1","X2")
  ellipse2 = as.data.frame(ellipse0[2]); colnames(ellipse2) = c("X1","X2")
  ellipse3 = as.data.frame(ellipse0[3]); colnames(ellipse3) = c("X1","X2")
  
  mypalette = brewer.pal(c,"Dark2")
  
  y = ggplot(sampleMat2, aes(x=`1`, y=`2`, colour=Cluster)) + geom_point() +
    scale_colour_brewer(palette="Dark2") +
    geom_point(data = mu, aes(x=mu$V1, y=mu$V2), color="black") +
    geom_point(data=ellipse1, aes(x=ellipse1$X1, y=ellipse1$X2), color = mypalette[2], size=0.3) +
    geom_point(data=ellipse2, aes(x=ellipse2$X1, y=ellipse2$X2), color = mypalette[1], size=0.3) +
    geom_point(data=ellipse3, aes(x=ellipse3$X1, y=ellipse3$X2), color = mypalette[3], size=0.3) +
    ggtitle(paste("Robust EM Algorithm with lamba =", lambda))
  # ggsave(paste("~/Box Sync/Research - Hui/Clustering Simulations/Plots of EM vs Actual/2017-11-03/EM_robust_L",lambda, "v3.pdf", sep =""))
  y
}
#' Function to plot the cluster assignment from the standard EM algorithm 
#'
#' @param sim_EMfit Output from the standard EM algorithm
#' @param sampleMat An (n x d+1) matrix with the observed data and true cluster assignments
#' @return A plot containing the observed data, colored by their clusters according to the EM algorithm. Estimated covariance and mean of clusters also plotted.
#' @export
plot_standard_EM = function(sim_EMfit, sampleMat){
  # Assign values from returned list to mu, sigma, T_mat, tau, e
  swap <- function(vec, from, to) {
    tmp <- to[ match(vec, from) ]
    tmp[is.na(tmp)] <- vec[is.na(tmp)]
    return(tmp)
  }
  # Assign values from returned list to mu, sigma, T_mat, tau, e
  mu = as.data.frame(t(sim_EMfit$mus))
  sigma = sim_EMfit$covs
  #error = sim_EMfit$e
  c = nrow(mu)
  d = ncol(mu)
  n = nrow(sim_EMfit$probs)
  lambda = sim_EMfit$lambda
  
  
  # Assign the cluster number to whichever one has the highest probability
  hard_assign = sim_EMfit$classification
  hard_assign = swap(hard_assign, c(1,2), c(2,1))
  
  sampleMat2 = as.data.frame(sampleMat)
  sampleMat2$Cluster = as.factor(hard_assign)
  
  
  plot.new()
  ellipse0 = lapply(1:3, function(j) mixtools::ellipse(mu=sim_EMfit$mus[,j], 
                                                       sigma=sim_EMfit$covs[,,j], npoints = 300, newplot=F))
  ellipse1 = as.data.frame(ellipse0[1]); colnames(ellipse1) = c("X1","X2")
  ellipse2 = as.data.frame(ellipse0[2]); colnames(ellipse2) = c("X1","X2")
  ellipse3 = as.data.frame(ellipse0[3]); colnames(ellipse3) = c("X1","X2")
  
  mypalette = brewer.pal(c,"Dark2")
  
  y = ggplot(sampleMat2, aes(x=`1`, y=`2`, colour=Cluster)) + geom_point() +
    scale_colour_brewer(palette="Dark2") +
    geom_point(data = mu, aes(x=mu$V1, y=mu$V2), color="black") +
    geom_point(data=ellipse1, aes(x=ellipse1$X1, y=ellipse1$X2), color = mypalette[2], size=0.3) +
    geom_point(data=ellipse2, aes(x=ellipse2$X1, y=ellipse2$X2), color = mypalette[1], size=0.3) +
    geom_point(data=ellipse3, aes(x=ellipse3$X1, y=ellipse3$X2), color = mypalette[3], size=0.3) +
    ggtitle(paste("Standard EM Algorithm"))
  # ggsave(paste("~/Box Sync/Research - Hui/Clustering Simulations/Plots of EM vs Actual/2017-11-03/EM_robust_L",lambda, "v3.pdf", sep =""))
  y
}
#' Function to plot the cluster assignment from the k-means algorithm.  
#'
#' @param sim_kmeans Output from the k-means function
#' @param sampleMat An (n x d+1) matrix with the observed data and true cluster assignments
#' @return A plot containing the observed data, colored by their clusters according to k-means. Estimated mean of clusters also plotted.
#' @export
plot_kmeans = function(sim_kmeans, sampleMat){
  
  sampleMat$Cluster = as.character(sim_kmeans$cluster)
  mu = as.data.frame(sim_kmeans$centers)
  
  # plot.new()
  # ellipse0 = lapply(1:3, function(j) mixtools::ellipse(mu=sim_EMfit$mu[j,], 
  #                                                      sigma=sim_EMfit$sigma[[j]], npoints = 300, newplot=F))
  # ellipse1 = as.data.frame(ellipse0[1]); colnames(ellipse1) = c("X1","X2")
  # ellipse2 = as.data.frame(ellipse0[2]); colnames(ellipse2) = c("X1","X2")
  # ellipse3 = as.data.frame(ellipse0[3]); colnames(ellipse3) = c("X1","X2")
  # 
  mypalette = brewer.pal(nrow(mu),"Dark2")
  
  y = ggplot(sampleMat, aes(x=`1`, y=`2`, colour=Cluster)) + geom_point() +
    scale_colour_brewer(palette="Dark2") +
    geom_point(data = mu, aes(x=mu$V1, y=mu$V2), color="black") +
    #geom_point(data=ellipse1, aes(x=ellipse1$X1, y=ellipse1$X2), color = mypalette[2], size=0.3) +
    #geom_point(data=ellipse2, aes(x=ellipse2$X1, y=ellipse2$X2), color = mypalette[1], size=0.3) +
    #geom_point(data=ellipse3, aes(x=ellipse3$X1, y=ellipse3$X2), color = mypalette[3], size=0.3) +
    ggtitle("K-means")
  # ggsave(paste("~/Box Sync/Research - Hui/Clustering Simulations/Plots of EM vs Actual/2017-11-03/EM_robust_L",lambda, "v3.pdf", sep =""))
  y
}