# 1. Turns vector into matrix
vector2matrix <- function(vec, p, diag = FALSE, bycolumn = FALSE) {
  m <- matrix(0, p, p)

  if(!bycolumn){
    m[lower.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
  } else {
    m[upper.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
  }
  return(m)
}

# 2. Transform precision into partial correlations for interpretation
pr2pc <- function(K) {
  D.Prec = diag(diag(K)^(-.5))
  R <- diag(2,dim(K)[1])-D.Prec%*%K%*%D.Prec
  colnames(R) <- colnames(K)
  rownames(R) <- rownames(K)
  return(R)
}

# 3. BDgraph stores graphs as byte strings for efficiency
string2graph <- function(Gchar, p) {
  Gvec <- rep(0, p*(p-1)/2)
  edges <- which(unlist(strsplit(as.character(Gchar), "")) == 1)
  Gvec[edges] = 1
  G <- matrix(0, p, p)
  G[upper.tri(G)] <- Gvec
  G <- G + t(G)
  diag(G) <- 0
  return(G)
}

# 4. BDgraph extract posterior distribution for estimates
extract_posterior <- function(fit, data, method = c("ggm", "gcgm"), not_cont, no_samples = 10000){
  m <- length(fit$all_graphs)
  k <- no_samples
  n <- nrow(data)
  p <- ncol(data)
  j <- 1
  densities <- rep(0, k)
  #Rs = array(0, dim=c(k, p, p))
  Rs = matrix(0, nrow = k, ncol = (p*(p-1))/2)
  if(method == "gcgm") {
    S <- BDgraph::get_S_n_p(data, method = method, n = n, not.cont = not_cont)$S
  } else {
    S <- t(data) %*% data
  }
  for (i in seq(1, m, length.out=k)) {
    graph_ix <- fit$all_graphs[i]
    G <- string2graph(fit$sample_graphs[graph_ix], p)
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)

    Rs[j,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    densities[j] <- sum(sum(G)) / (p*(p-1))
    j <- j + 1
  }
  return(list(Rs, densities))
}

# 5. Samples from the G-wishart distribution
gwish_samples <- function(G, S, nsamples=1000) {
  n <- nrow(S)
  p <- ncol(S)
  #Rs <- array(0, dim=c(nsamples, p, p))
  Rs = matrix(0, nrow = nsamples, ncol = (p*(p-1))/2)

  for (i in 1:nsamples) {
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)*(G + diag(p))
    Rs[i,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    #Rs[i,,] <- .pr2pc(K)
  }
  return(Rs)
}


# 6. Centrality of weighted graphs

# Strength centrality only ## FASTER CODE
centrality <- function(res){
  Nsamples <- nrow(res$samples_posterior)
  p <- nrow(res$parameters)
  strength_samples <- matrix(0, nrow = Nsamples, ncol = p)
  for(i in 1:Nsamples){
    strength_samples[i, ] <- rowSums(abs(vector2matrix(res$samples_posterior[i,], p, bycolumn = T)))
  }
  return(strength_samples)
}

# Strength, betweenness and closeness centrality ## SLOWER CODE
centrality_all <- function(res){
  Nsamples <- nrow(res$samples_posterior)
  p <- as.numeric(nrow(res$parameters))
  samples <- res$samples_posterior
  for(i in 1:Nsamples){

    #Strength
    strength_samples <- rowSums(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE)))
    #EI
    influence_samples <- rowSums(vector2matrix(samples[i, ], p, bycolumn = TRUE))

    DistMat <- 1/(ifelse(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE))==0,0,abs(vector2matrix(samples[i, ], p, bycolumn = T))))
    igraphObject <- igraph::graph.adjacency(DistMat, weighted = TRUE, mode = "undirected")
    # Closeness
    closeness_samples <- igraph::closeness(igraphObject)
    # Betweenness
    betweenness_samples <- igraph::estimate_betweenness(igraphObject,cutoff = 1/1e-10)

    if(i > 1){
      centrality_samples <- rbind(centrality_samples, cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                                                            rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples)))
    } else {
      centrality_samples <- cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                                  rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples))
    }
  }
  colnames(centrality_samples) <- c("Centrality", colnames(res$parameters))
  centrality_samples <- as.data.frame(centrality_samples)
  centrality_samples[, 2:(p+1)] <- sapply(centrality_samples[, 2:(p+1)], as.numeric)
  return(centrality_samples)
}


# 7. turn list into matrix
list2matrix <- function(obj, p) {
  nlist <- length(obj)/(p*p)
  m <- obj[, , 1]
  nest <- sum(lower.tri(m))
  res <- matrix(0, nrow = nlist, ncol = nest)
  for(i in 1:nlist){
    m <- obj[, , i]
    res[i, ] <- as.vector(m[lower.tri(m)])
  }
  return(res)
}

# 8. Set defaults of a function
set_defaults <- function(args, ...) {
  dots <- list(...)
  def_args <- setdiff(names(args), names(dots))
  dots[def_args] <- args[def_args]

  return(dots)
}

# 9. Check for empty ...

dots_check <- function(...){
  if(...length() > 0){
    warning("Arguments specified with ... are unused. ")
  }
}
