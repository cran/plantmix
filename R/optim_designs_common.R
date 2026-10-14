## Functions used in both optim_designs_var_mix.R and optim_designs_crop_mix.R

##' @noRd
convDiallel2Combs <- function(diallel) {
  ind <- which(upper.tri(diallel, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(diallel)
  out <- data.frame(
    row = nn[[1]][ind[, 1]],
    col = nn[[2]][ind[, 2]],
    val = diallel[ind],
    stringsAsFactors = FALSE
  )
  out <- out[out$val == 1, ]
  out$val <- NULL
  colnames(out) <- c("comp1", "comp2")
  return(out)
}

##' Plot diallel
##'
##' Plots a diallel matrix.
##' @param diallel matrix
##' @param main title
##' @param colLabs optional vector of colors for each label (black otherwise)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{getDesignBinaryVarMix}}, \code{\link{plotDesignVarMix}}, \code{\link{plotDesignCropMix}}
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' plotDiallel(design$diallel)
##' @noRd
plotDiallel_OLD <- function(diallel, main = NULL, colLabs = NULL) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(1, 4, 6, 1) + 0.1)
  graphics::image(t(diallel)[, nrow(diallel):1],
    axes = FALSE,
    col = c("white", "black")
  )
  if (is.null(main)) {
    combs <- convDiallel2Combs(diallel)
    idxS <- getShannonIndex(getGenoPropsFromCombs(combs))
    main <- ""
    if (nrow(diallel) == ncol(diallel)) {
      main <- paste0(nrow(diallel), " genotypes")
    } else {
      main <- paste0(nrow(diallel), " x ", ncol(diallel), " genotypes")
    }
    main <- paste0(
      main, "; ", sum(diallel), " binary mixtures",
      "; entropy = ", round(idxS, 4)
    )
  }
  title(main = main, line = 4)

  y_labs <- rev(rownames(diallel))
  for (i in 1:length(y_labs)) {
    c_col <- if (!is.null(colLabs) && y_labs[i] %in% names(colLabs)) colLabs[y_labs[i]] else "black"
    axis(side = 2, at = seq(0, 1, length.out = nrow(diallel))[i], labels = y_labs[i], las = 1, col.axis = c_col)
  }

  x_labs <- colnames(diallel)
  for (i in 1:length(x_labs)) {
    c_col <- if (!is.null(colLabs) && x_labs[i] %in% names(colLabs)) colLabs[x_labs[i]] else "black"
    axis(side = 3, at = seq(0, 1, length.out = ncol(diallel))[i], labels = x_labs[i], col.axis = c_col)
  }
}

##' Plot diallel
##'
##' Plots a diallel matrix.
##' @param diallel matrix
##' @param main title
##' @param colLabs optional vector of colors for each label (black otherwise)
##' @param fill_values values with which to fill the tiles
##' @param margins named numeric vector with margins around the plot (unit with the sizes of the top, right, bottom, and left margins)
##' @param xlab_angle rotation angle for the x-axis labels
##' @return ggplot object
##' @author Timothee Flutre
##' @seealso \code{\link{getDesignBinaryVarMix}}, \code{\link{plotDesignVarMix}}, \code{\link{plotDesignCropMix}}
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' plotDiallel(design$diallel)
##' @export
plotDiallel <- function(diallel, main = NULL, colLabs = NULL,
                        fill_values = c("0" = "white", "1" = "black"),
                        margins = c("t" = 10, "r" = 25, "b" = 20, "l" = 40),
                        xlab_angle = 45) {
  stopifnot(
    all(as.character(c(diallel)) %in% names(fill_values)),
    !is.null(names(margins))
  )
  if (is.null(main)) {
    combs <- convDiallel2Combs(diallel)
    idxS <- getShannonIndex(getGenoPropsFromCombs(combs))
    if (nrow(diallel) == ncol(diallel)) {
      main_text <- paste0(nrow(diallel), " genotypes")
    } else {
      main_text <- paste0(nrow(diallel), " x ", ncol(diallel), " genotypes")
    }
    main <- paste0(
      main_text, "; ", sum(diallel), " binary mixtures",
      "; entropy = ", round(idxS, 4)
    )
  }
  if (!"t" %in% names(margins)) {
    margins["t"] <- 20
  }
  if (!"r" %in% names(margins)) {
    margins["r"] <- 20
  }
  if (!"b" %in% names(margins)) {
    margins["b"] <- 20
  }
  if (!"l" %in% names(margins)) {
    margins["l"] <- 40
  }

  ## reshape data matrix into long format, reversing row order
  row_names <- rownames(diallel)
  col_names <- colnames(diallel)
  plot_data <- expand.grid(
    x_lab = col_names,
    y_lab = rev(row_names),
    stringsAsFactors = FALSE
  )

  ## extract values corresponding to the grid pairs
  ## with the layout as when printing the input matrix
  matrix_flipped <- t(diallel)[, nrow(diallel):1]
  plot_data$value <- as.vector(matrix_flipped)

  ## convert labels to factors to preserve the specific matrix ordering
  plot_data$x_lab <- factor(plot_data$x_lab, levels = col_names)
  plot_data$y_lab <- factor(plot_data$y_lab, levels = rev(row_names))

  ## process custom colors for axis labels (default is black)
  x_colors <- sapply(col_names, function(lbl) if (!is.null(colLabs) && lbl %in% names(colLabs)) colLabs[lbl] else "black")
  x_label_df <- data.frame(
    x_lab = factor(col_names, levels = col_names),
    y_pos = nrow(diallel) + 0.55, # position just above the top cell boundary
    color = x_colors,
    stringsAsFactors = FALSE
  )

  y_colors <- sapply(rev(row_names), function(lbl) if (!is.null(colLabs) && lbl %in% names(colLabs)) colLabs[lbl] else "black")
  y_label_df <- data.frame(
    x_pos = 0.4, # position just to the left of the first cell boundary
    y_lab = factor(rev(row_names), levels = rev(row_names)),
    color = y_colors,
    stringsAsFactors = FALSE
  )

  ## build the ggplot object
  p <- ggplot() +
    geom_tile(
      data = plot_data,
      aes(x = .data[["x_lab"]], y = .data[["y_lab"]], fill = as.factor(.data[["value"]])),
      color = "gray90"
    ) + # matrix cells
    ## custom x axis (on top)
    geom_text(
      data = x_label_df,
      aes(x = .data[["x_lab"]], y = .data[["y_pos"]], label = .data[["x_lab"]], color = .data[["color"]]),
      angle = xlab_angle, hjust = 0, vjust = 0, fontface = "plain", inherit.aes = FALSE
    ) +
    ## custom y axis (on the left)
    geom_text(
      data = y_label_df,
      aes(x = .data[["x_pos"]], y = .data[["y_lab"]], label = .data[["y_lab"]], color = .data[["color"]]),
      hjust = 1, fontface = "plain", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fill_values) +
    scale_color_identity() +
    labs(title = main, x = NULL, y = NULL) +
    theme_minimal() +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none", # hide binary fill legend
      plot.title = element_text(
        hjust = 0.5, face = "bold", margin = margin(b = 45)
      ),
      ## hide native axis texts
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(), # remove background grid lines
      plot.margin = margin(
        t = margins["t"], r = margins["r"],
        b = margins["b"], l = margins["l"]
      )
    )

  return(p)
}
