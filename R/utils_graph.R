##' Combinations to graph
##'
##' Returns a weighted graph from a list of combinations (mixtures), the distance between two mixtures being their Jaccard index.
##' @param combs matrix, data.frame or list with one component per mixture, each one being a vector of mixture elements
##' @param method weighing method (only "jaccard" for the moment)
##' @return igraph object
##' @author Timothee Flutre
##' @examples
##' combs <- data.frame(comp1=c("comp1","comp1"),
##'                     comp2=c("comp2","comp3"),
##'                     row.names=c("mix1","mix2"))
##' (out <- convCombs2Graph(combs))
##' combs <- list(mix1=c("comp1","comp2","comp3"),
##'               mix2=c("comp1","comp3"),
##'               mix3=c("comp2","comp3"))
##' (out <- convCombs2Graph(combs))
##' @noRd
convCombs2Graph <- function(combs, method = "jaccard") {
  if (is.matrix(combs)) {
    combs <- as.data.frame(combs)
  }
  if (is.data.frame(combs)) {
    colnames(combs) <- NULL
    mix_names <- rownames(combs)
    combs <- lapply(1:nrow(combs), function(i) {
      unlist(combs[i, ])
    })
    names(combs) <- mix_names
  }
  stopifnot(method %in% c("jaccard"))

  dist_btw_mixes <- NA
  if (method == "jaccard") {
    dist_btw_mixes <- jaccardBtwCombs(combs)
  }
  out <- graph_from_adjacency_matrix(dist_btw_mixes, "undirected", weighted = TRUE)
  ## out <- set_vertex_attr(out, "name", value=names(combs))

  return(out)
}

##' Graph to mixtures
##'
##' Returns a list of mixtures from a graph.
##' @param graph igraph object
##' @param retAllNodes logical
##' @noRd
convGraph2Combs <- function(graph, retAllNodes = FALSE) {
  stopifnot(is(graph, "igraph"))

  out <- as_data_frame(graph)
  colnames(out) <- c("comp1", "comp2")

  if (retAllNodes) {
    degs <- degree(graph)

    if (any(degs == 0)) {
      if (is_bipartite(graph)) {
        tmp <- data.frame(
          node = vertex_attr(graph)$name,
          deg = degs[vertex_attr(graph)$name],
          type = vertex_attr(graph)$type
        )
        nodes0 <- names(degs)[which(degs == 0)]
        out2 <- rbind(out, data.frame(
          comp1 = rep(NA, length(nodes0)),
          comp2 = NA
        ))
        for (i in seq_along(nodes0)) {
          node <- nodes0[i]
          comp <- ifelse(!tmp[node, "type"], "comp1", "comp2")
          idx <- i + nrow(out)
          out2[idx, comp] <- node
        }
        out <- out2
      } else { # graph is not bipartite
        msg <- "retAllNodes=TRUE is not (yet) implemented when graph is not bipartite"
        warning(msg, immediate. = TRUE)
        ## TODO
      }
    } # end if any(degs == 0)
  } # end if retAllNodes

  return(out)
}

##' @noRd
getShannonIndexFromGraph <- function(graph) {
  out <- NA
  if (!is_bipartite(graph)) {
    combs <- convGraph2Combs(graph)
    props <- getGenoPropsFromCombs(combs)
    out <- getShannonIndex(props)
  } else {
    degs <- degree(graph)
    tmp <- data.frame(
      node = vertex_attr(graph)$name,
      deg = degs[vertex_attr(graph)$name],
      type = vertex_attr(graph)$type
    )
    tmp <- list(
      droplevels(tmp[!tmp$type, ]),
      droplevels(tmp[tmp$type, ])
    )
    out <- sapply(tmp, function(x) {
      deg <- x$deg[x$deg > 0]
      getShannonIndex(deg / sum(deg))
    })
  }
  return(out)
}
