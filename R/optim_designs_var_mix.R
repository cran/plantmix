##' Plot design for varietal mixtures
##'
##' Plots design for varietal mixtures, as a graph and as a diallel.
##' @param graph graph
##' @param levGenos character vector with the names of the genotypes
##' @param main main title for the graph plot
##' @param subplots character vector indicating which object(s) to plot
##' @return nothing
##' @author Timothee Flutre
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' plotDesignVarMix(design$graph, levGenos)
##' @export
plotDesignVarMix <- function(graph, levGenos, main = NULL, subplots = c("graph", "diallel")) {
  combs <- convGraph2Combs(graph)
  if ("diallel" %in% subplots) {
    diallel <- getDiallelVarMix(levGenos, combs)
  }
  if (length(subplots) == 2) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 2))
  }
  if (is.null(main)) {
    idxS <- getShannonIndexFromGraph(graph)
    main <- paste0("entropy = ", round(idxS, 4))
  }
  if ("graph" %in% subplots) {
    plot(graph, main = main)
  }
  if ("diallel" %in% subplots) {
    plotDiallel(diallel)
  }
}

##' @noRd
getInitGraphVarMix <- function(levGenos, nbMixes) {
  nbGenos <- length(levGenos)
  out <- sample_gnm(n = nbGenos, m = nbMixes)
  stopifnot(all(!which_multiple(out)))
  out <- set_vertex_attr(out, "name", value = levGenos)
  return(out)
}

##' @noRd
getDiallelVarMix <- function(levGenos, mixtures, fillDiag = 0) {
  stopifnot(ncol(mixtures) >= 2)
  nbGenos <- length(levGenos)
  diallel <- matrix(0,
    nrow = nbGenos, ncol = nbGenos,
    dimnames = list(sort(levGenos), sort(levGenos))
  )
  diag(diallel) <- fillDiag
  for (i in 1:nrow(mixtures)) {
    genos <- sort(unlist(mixtures[i, 1:2]))
    diallel[genos[1], genos[2]] <- 1
  }
  return(diallel)
}

##' @noRd
rmvEdgeBtwMostConnectedNodesVarMix <- function(graph, verbose = FALSE) {
  degs <- degree(graph)
  degs <- sort(degs, decreasing = TRUE)
  combs <- convGraph2Combs(graph)
  combs$deg1 <- degs[combs$comp1]
  combs$deg2 <- degs[combs$comp2]
  combs$degs <- rowSums(combs[, c("deg1", "deg2")])
  combs$weight <- combs$degs / sum(combs$degs)
  idx <- which.max(combs$weight)
  ## idx <- sample(1:nrow(combs), 1, prob=combs$weight)
  edge_tormv <- paste0(combs$comp1[idx], "|", combs$comp2[idx])
  if (verbose) {
    print(paste0("edge to remove: ", edge_tormv))
  }
  out <- delete_edges(graph, edge_tormv)
  stopifnot(sum(degree(out)) == sum(degree(graph)) - 2)
  return(out)
}

##' @noRd
addEdgeBtwLessConnectedNodesVarMix <- function(graph, verbose) {
  degs <- degree(graph)
  degs <- sort(degs, decreasing = FALSE)
  combs <- convGraph2Combs(graph)
  existingMixes <- c(
    paste(combs$comp1, combs$comp2, sep = "|"),
    paste(combs$comp2, combs$comp1, sep = "|")
  )
  i <- 1
  j <- 2
  new_edge <- names(degs)[c(i, j)]
  while (TRUE) {
    if (!paste(new_edge, collapse = "|") %in% existingMixes) {
      break
    }
    if (j < length(degs)) {
      j <- j + 1
    } else {
      i <- i + 1
      j <- i + 1
    }
    new_edge <- names(degs)[c(i, j)]
    if (verbose) {
      print(paste0("i=", i, " j=", j, ": ", paste(new_edge, collapse = "|")))
    }
  }
  if (verbose) {
    print(paste0("new edge: ", paste(new_edge, collapse = "|")))
  }
  out <- add_edges(graph, new_edge)
  stopifnot(
    sum(degree(out)) == sum(degree(graph)) + 2,
    all(!which_multiple(out))
  )
  return(out)
}

##' Optimize design for binary varietal mixtures
##'
##' Optimizes a design for binary varietal mixtures
##' @param levGenos character vector with the names of the genotypes
##' @param nbMixes number of mixtures
##' @param seed seed for the pseudo-random number generator
##' @param showPlots logical indicating if plots should be showed
##' @param verbose verbosity level
##' @return list
##' @author Timothee Flutre
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' plotDesignVarMix(design$graph, levGenos)
##' ggplot(design$entropies) + aes(x=iteration, y=entropy) + geom_point() + geom_line() +
##'   geom_hline(yintercept = design$max_ent) + theme_bw()
##' @export
getDesignBinaryVarMix <- function(levGenos, nbMixes, seed = NULL, showPlots = FALSE, verbose = FALSE) {
  stopifnot(
    is.character(levGenos),
    nbMixes <= choose(length(levGenos), 2)
  )

  out <- list()

  tmp <- getOptimalProps(levGenos, nbMixes)
  max_idxS <- getShannonIndex(tmp)
  if (verbose) {
    print(paste0("max entropy = ", max_idxS))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  entropies <- c()
  g_curr <- getInitGraphVarMix(levGenos, nbMixes)
  idxS <- getShannonIndexFromGraph(g_curr)
  if (showPlots) {
    plotDesignVarMix(g_curr, levGenos)
  }
  entropies <- c(entropies, idxS)

  i <- 1
  while (idxS < max_idxS) {
    if (verbose) {
      print(paste0("iter ", i, "; entropy = ", idxS))
    }
    g_tmp <- rmvEdgeBtwMostConnectedNodesVarMix(g_curr, verbose)
    g_new <- addEdgeBtwLessConnectedNodesVarMix(g_tmp, verbose)
    idxS <- getShannonIndexFromGraph(g_new)
    entropies <- c(entropies, idxS)
    if (showPlots) {
      plotDesignVarMix(g_new, levGenos)
    }
    g_curr <- g_new
    i <- i + 1
  }
  stopifnot(all(!which_multiple(g_curr)))
  out$graph <- g_curr

  combs <- convGraph2Combs(g_curr)
  out$combs <- combs

  diallel <- getDiallelVarMix(levGenos, combs)
  out$diallel <- diallel

  out$entropies <- data.frame(
    iteration = seq_along(entropies),
    entropy = entropies
  )
  out$max_ent <- max_idxS

  return(out)
}
