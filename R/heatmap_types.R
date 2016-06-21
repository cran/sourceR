
heatmap_types <- function(object, iters, hclust_method = "complete", hm_colors = c("blue", "white")) {

  if ("q" %in% names(object$posterior)) {
      q_data <- object$posterior$q
      q_names <- colnames(q_data)
      
      # the cluster list in the posterior contains the group that each
      # bacterial type has been assigned to at each iteration of the mcmc
      if (missing(iters)) {
        groups <- as.data.frame(t(object$posterior$cluster))
      } else {
        groups <- as.data.frame(t(object$posterior$cluster[iters,]))
      }
      groups <- as.data.frame(apply(groups, 2, function(x) as.factor(x)))
      rownames(groups) <- q_names
      
      # compute dissimilarity matrix for the type effect clusters
      disim_clust_g <- cluster::daisy(groups)
      clu <- stats::hclust(disim_clust_g, method = hclust_method)
      dend <- stats::as.dendrogram(clu)
      
      # OPTIONAL: change the colour of the heatmap. The lighter the colour
      # (when using the default white blue colour scheme),
      # the higher the dissimilarity between the 2 types (i.e. the less
      #often two type effects are assigned to the same group in the mcmc)
      hmcols <- colorRampPalette(hm_colors)(299)
      
      heatmap_data <- as.matrix(disim_clust_g)
      if (sapply(q_names, function(x) substring(as.character(x), 1, last = 4))[1] != "type") {
        stop("The type names should beging with the word type. Please contact package maintainer unless you have manually changed the type names in the posterior.")
      }
      q_names <- sapply(q_names, function(x) substring(as.character(x), 5, last = 1000000L))
      rownames(heatmap_data) <- colnames(heatmap_data) <- q_names

        gplots::heatmap.2(heatmap_data,
                          density.info = "none",# turns off density plot in the legend
                          trace = "none",       # turns off trace lines in the heat map
                          col = hmcols,         # use color palette defined earlier
                          dendrogram = "col",   # only draw a row dendrogram
                          Colv = dend, 
                          Rowv = dend, 
                          symm = TRUE,
                          key = F)
  } else {
    stop("There are no type effects in the posterior. \nPerhaps the type effects were fixed (e.g. a negative binomial likelihood was used.")
  }

}
