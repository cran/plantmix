##' Plot design for crop mixtures
##'
##' Plots design for crop mixtures, as a graph and as a diallel.
##' @param graph graph
##' @param levGenos1 character vector with the names of the genotypes of the first species
##' @param levGenos2 character vector with the names of the genotypes of the second species
##' @param main main title for the graph plot
##' @return nothing
##' @author Timothee Flutre
##' @examples
##' nbGenos1 <- 30
##' levGenos1 <- sprintf(fmt=paste0("S1_%0", floor(log10(nbGenos1))+1, "i"),
##'                      1:nbGenos1)
##' nbGenos2 <- 8
##' levGenos2 <- sprintf(fmt=paste0("S2_%0", floor(log10(nbGenos2))+1, "i"),
##'                      1:nbGenos2)
##' nbMixes <- 60 # only binary and balanced
##' design <- getDesignBinaryCropMix(levGenos1, levGenos2, nbMixes, seed=12345)
##' plotDesignCropMix(design$graph, levGenos1, levGenos2)
##' @export
plotDesignCropMix <- function(graph, levGenos1, levGenos2, main = NULL) {
  combs <- convGraph2Combs(graph)
  diallel <- getDiallelCropMix(levGenos1, levGenos2, combs)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, 2))
  tmp <- layout_as_bipartite(graph) # , hgap=2)
  colnames(tmp) <- c("y", "x_init") # change two-row into two-column layout
  tmp <- cbind(tmp, "x" = 0)
  tmp[tmp[, "x_init"] == 0, "x"] <- 1
  if (is.null(main)) {
    idxS <- getShannonIndexFromGraph(graph)
    main <- paste0(
      "entropy 1 = ", round(idxS[1], 4),
      "\nentropy 2 = ", round(idxS[2], 4)
    )
  }
  plot(graph, main = main, layout = tmp[, c("x", "y")])
  plotDiallel(diallel)
}

##' @noRd
getInitGraphCropMix <- function(levGenos1, levGenos2, nbMixes) {
  nbGenos1 <- length(levGenos1)
  nbGenos2 <- length(levGenos2)
  out <- sample_bipartite(n1 = nbGenos1, n2 = nbGenos2, type = "gnm", m = nbMixes)
  stopifnot(all(!which_multiple(out)))
  out <- set_vertex_attr(out, "name", value = c(levGenos1, levGenos2))
  out <- set_vertex_attr(out, "color", value = c(rep("orange", nbGenos1), rep("green", nbGenos2)))
  return(out)
}

##' @noRd
getDiallelCropMix <- function(levGenos1, levGenos2, mixtures) {
  stopifnot(ncol(mixtures) >= 2)
  nbGenos1 <- length(levGenos1)
  nbGenos2 <- length(levGenos2)
  diallel <- matrix(0,
    nrow = nbGenos1, ncol = nbGenos2,
    dimnames = list(sort(levGenos1), sort(levGenos2))
  )
  for (i in 1:nrow(mixtures)) {
    genos <- sort(unlist(mixtures[i, 1:2]))
    diallel[genos[1], genos[2]] <- 1
  }
  return(diallel)
}

##' @noRd
rmvEdgeBtwMostConnectedNodesCropMix <- function(graph, verbose = FALSE) {
  degs <- degree(graph)
  degs <- sort(degs, decreasing = TRUE)
  combs <- convGraph2Combs(graph, retAllNodes = TRUE)
  combs$deg1 <- degs[combs$comp1]
  combs$deg2 <- degs[combs$comp2]
  combs$degs <- rowSums(combs[, c("deg1", "deg2")], na.rm = TRUE)
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
addEdgeBtwLessConnectedNodesCropMix <- function(graph, max_idxS, verbose) {
  degs <- degree(graph)
  degs <- sort(degs, decreasing = FALSE)
  combs <- convGraph2Combs(graph, retAllNodes = TRUE)
  existingMixes <- c(
    paste(combs$comp1, combs$comp2, sep = "|"),
    paste(combs$comp2, combs$comp1, sep = "|")
  )
  tmp <- data.frame(
    node = vertex_attr(graph)$name,
    deg = degs[vertex_attr(graph)$name],
    type = vertex_attr(graph)$type
  )
  tmp <- list(
    droplevels(tmp[!tmp$type, ]),
    droplevels(tmp[tmp$type, ])
  )
  tmp <- lapply(tmp, function(x) {
    x[order(x$deg, decreasing = FALSE), ]
  })
  i <- 1
  j <- 1
  new_edge <- c(
    rownames(tmp[[1]])[i],
    rownames(tmp[[2]])[j]
  )
  while (TRUE) {
    if (!paste(new_edge, collapse = "|") %in% existingMixes) {
      break
    }
    idxS <- sapply(tmp, function(x) {
      deg <- x$deg[x$deg > 0]
      getShannonIndex(deg / sum(deg))
    })
    prop_idxS <- idxS / max_idxS
    if (prop_idxS[1] <= prop_idxS[2]) {
      i <- i + 1
    } else {
      j <- j + 1
    }
    new_edge <- c(
      rownames(tmp[[1]])[i],
      rownames(tmp[[2]])[j]
    )
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

##' Optimize design for binary crop mixtures
##'
##' Optimizes a design for binary crop mixtures
##' @param levGenos1 character vector with the names of the genotypes of the first species
##' @param levGenos2 character vector with the names of the genotypes of the second species
##' @param nbMixes number of mixtures
##' @param seed seed for the pseudo-random number generator
##' @param showPlots logical indicating if plots should be showed
##' @param verbose verbosity level
##' @return list
##' @author Timothee Flutre
##' @examples
##' nbGenos1 <- 30
##' levGenos1 <- sprintf(fmt=paste0("S1_%0", floor(log10(nbGenos1))+1, "i"),
##'                      1:nbGenos1)
##' nbGenos2 <- 8
##' levGenos2 <- sprintf(fmt=paste0("S2_%0", floor(log10(nbGenos2))+1, "i"),
##'                      1:nbGenos2)
##' nbMixes <- 60 # only binary and balanced
##' design <- getDesignBinaryCropMix(levGenos1, levGenos2, nbMixes, seed=12345)
##' plotDesignCropMix(design$graph, levGenos1, levGenos2)
##' ggplot(design$entropies) + aes(x=iteration, y=entropy, group=type) + geom_point() + geom_line() +
##'   geom_hline(yintercept = design$max_ent) + theme_bw() + facet_wrap(~ type)
##' @export
getDesignBinaryCropMix <- function(levGenos1, levGenos2, nbMixes, seed = NULL,
                                   showPlots = FALSE, verbose = FALSE) {
  stopifnot(
    is.character(levGenos1),
    is.character(levGenos2),
    length(levGenos1) >= length(levGenos2),
    nbMixes <= length(levGenos1) * length(levGenos2)
  )

  out <- list()

  max_idxS <- c(
    getShannonIndex(getOptimalProps(levGenos1, nbMixes)),
    getShannonIndex(getOptimalProps(levGenos2, nbMixes))
  )
  if (verbose) {
    print(paste0("max entropy 1 = ", max_idxS[1]))
    print(paste0("max entropy 2 = ", max_idxS[2]))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  entropies <- matrix(nrow = 0, ncol = 2)
  g_curr <- getInitGraphCropMix(levGenos1, levGenos2, nbMixes)
  idxS <- getShannonIndexFromGraph(g_curr)
  if (showPlots) {
    plotDesignCropMix(g_curr, levGenos1, levGenos2)
  }
  entropies <- rbind(entropies, idxS)

  i <- 1
  while (all(idxS[1] < max_idxS[1] | idxS[2] < max_idxS[2])) {
    if (verbose) {
      print(paste0(
        "iter ", i, "; entropy1 = ", idxS[1],
        "; entropy2 = ", idxS[2]
      ))
    }
    g_tmp <- rmvEdgeBtwMostConnectedNodesCropMix(g_curr, verbose)
    g_new <- addEdgeBtwLessConnectedNodesCropMix(g_tmp, max_idxS, verbose)
    idxS <- getShannonIndexFromGraph(g_new)
    entropies <- rbind(entropies, idxS)
    if (showPlots) {
      plotDesignCropMix(g_new, levGenos1, levGenos2)
    }
    g_curr <- g_new
    i <- i + 1
  }
  stopifnot(all(!which_multiple(g_curr)))
  out$graph <- g_curr

  combs <- convGraph2Combs(g_curr)
  out$combs <- combs

  diallel <- getDiallelCropMix(levGenos1, levGenos2, combs)
  out$diallel <- diallel

  out$entropies <- data.frame(
    iteration = rep(1:nrow(entropies), times = 2),
    type = rep(c("1", "2"), each = nrow(entropies)),
    entropy = c(entropies[, 1], entropies[, 2])
  )
  out$entropies <- out$entropies[order(out$entropies$iteration, out$entropies$type), ]
  out$max_ent <- max_idxS

  return(out)
}
