library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

# Simple plots use base plot, others are ggplot based

## ------------------------------------------------------------------------
# Simple plotting theme for ggplot using arial family font
theme.plain <- function(p, base_size = 11, base_family = "ArialMT") {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.75),
          axis.ticks = element_line(size=0.75),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}

## ------------------------------------------------------------------------
## Make a PCA plot based on a prcomp object with the % of variance on the axis labels
pca.plot <- function(prcomp.obj, comp.1=1, comp.2=2, labels=NULL, fill=NULL, color=NULL, shape=NULL, size=1, main="PCA plot", ...){
  
  var.prop <- prcomp.obj$sdev^2/sum(prcomp.obj$sdev^2)
  df.plot  <- as.data.frame(prcomp.obj$x[,c(comp.1, comp.2)])
  
  if (is.null(fill)) {
    #df.plot$fill <- rep("orange", nrow(df.plot))
  } else {
    df.plot$fill <- fill
  }

  if (is.null(shape)) {
   # df.plot$shape <- rep("20", nrow(df.plot))
  } else {
    df.plot$shape <- shape
  }
  
  if (is.null(color)) {
   # df.plot$color <- rep("black", nrow(df.plot))
  } else {
    df.plot$color <- color
  }
  
  
  colnames(df.plot)[1:2] <- c("C1", "C2")
  
  p <- ggplot(df.plot, aes(y=C2, x=C1, shape=shape, color=color, fill=fill)) +
    geom_point(size=size, ...) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_line("#E3E3E3", 0.5, 2),
          panel.grid.minor=element_line("#E3E3E3", 0.25, 2)) +
    labs(x=paste0("PC", comp.1, " :", round(var.prop[comp.1], 3)*100, "% variance"),
         y=paste0("PC", comp.2, " :", round(var.prop[comp.2], 3)*100, "% variance")) 

  if (!is.null(labels)) {
    p <- p + geom_text_repel(aes(label=labels))
  }
  
  return(p)
  
}

## ------------------------------------------------------------------------
# PCA plot based on a DEseq object, thanks Aaron :)
pca.plot.v2 <- function (object, intgroup = "condition", pc1=1, pc2=2,  ntop = 10000, returnData = FALSE, ...) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument ‘intgroup’ should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse =" : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group,
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(pc1, pc2)]
    return(d)
  }
  p <- ggplot2::ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 3, ...) +
    xlab(paste0("PC", pc1 , ": ", round(percentVar[pc1] * 100), "% variance")) +
    ylab(paste0("PC", pc2 , ": ", round(percentVar[pc2] * 100), "% variance")) +
    coord_fixed() 
  return(p)
}

## ------------------------------------------------------------------------
# Simple heatmap with auto labels
simple.hm <- function(data, cellwidth=12, cellheight=12, limit=NULL, ...) {
  break.list <- seq(-max(abs(data)), max(abs(data)), by=max(abs(data))/100)
  
  pheatmap(data,
           breaks=break.list,
           col=colorRampPalette(rev(brewer.pal(n=7, name ="RdBu")))(length(break.list)),
           cellwidth=cellwidth,
           cellheight=cellheight,
           ...)
  
}

## ------------------------------------------------------------------------
# QQplot ONLY FOR PVALUES!
simple.qq.plot <- function (observedPValues) {
  observedPValues <- na.omit(observedPValues)
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}

## ------------------------------------------------------------------------
# The same function as in the TCseq package, but commented out the multiplot line at the end
# As I just want a list of the plots to edit.
timeclustplot.return <- function (object = NULL, categories = "timepoint", value = "expression", 
                                  cols = NULL, cl.color = "gray50", membership.color = rainbow(30, 
                                                                                               s = 3/4, v = 1, start = 1/6), title.size = 18, axis.line.size = 0.6, 
                                  axis.title.size = 18, axis.text.size = 16, legend.title.size = 14, 
                                  legend.text.size = 14) 
{
  if (class(object) != "clust" && class(object) != "TCA") {
    stop("object should be a 'timeclust' object or a 'TCA' object")
  }
  if (class(object) == "clust") {
    data <- object@data
    cluster <- object@cluster
    membership <- object@membership
  }
  if (class(object) == "TCA") {
    data <- object@clusterRes@data
    cluster <- object@clusterRes@cluster
    membership <- object@clusterRes@membership
  }
  ncl <- max(cluster)
  membercolor <- vector(length = length(cluster))
  membervalue <- list()
  counter <- 0
  if (!sum(dim(membership) == 0) == 2) {
    color <- membership.color
    colorseq <- seq(0, 1, length = length(color))
    for (i in seq_len(ncl)) {
      mtmp <- membership[cluster == i, i]
      membervalue[[i]] <- mtmp
      for (j in seq_len(length(mtmp))) {
        counter <- counter + 1
        ind <- which(abs(colorseq - mtmp[j]) == min(abs(colorseq - 
                                                          mtmp[j])))
        membercolor[counter] <- color[ind]
      }
    }
    membervalue <- unlist(membervalue)
    names(membercolor) <- membervalue
  }
  plotlist <- list()
  for (i in seq_len(ncl)) {
    title <- paste0("Cluster ", i)
    dtmp <- data[cluster == i, ]
    a <- which(cluster == i)
    if (length(a) == 1) {
      dtmp <- data.frame(time = 1:length(dtmp), value = dtmp)
      if (!sum(dim(membership) == 0) == 2) {
        m <- membership[cluster == i, i]
        colorname = toString(m)
        plotlist[[i]] <- ggplot(dtmp, aes(x = time, y = value)) + 
          geom_line(colour = membercolor[colorname]) + 
          theme_bw() + ggtitle(title) + scale_x_continuous(breaks = dtmp$time, 
                                                           labels = row.names(dtmp)) + labs(x = categories, 
                                                                                            y = value) + theme(plot.title = element_text(size = title.size), 
                                                                                                               axis.line.x = element_line(color = "black", 
                                                                                                                                          size = axis.line.size), axis.line.y = element_line(color = "black", 
                                                                                                                                                                                             size = axis.line.size), axis.title = element_text(size = axis.title.size), 
                                                                                                               axis.text = element_text(size = axis.text.size), 
                                                                                                               legend.position = "none", panel.border = element_blank(), 
                                                                                                               panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      else {
        plotlist[[i]] <- ggplot(dtmp, aes(x = time, y = value)) + 
          geom_line(colour = cl.color) + theme_bw() + 
          ggtitle(title) + scale_x_continuous(breaks = dtmp$time, 
                                              labels = row.names(dtmp)) + labs(x = categories, 
                                                                               y = value) + theme(plot.title = element_text(size = title.size), 
                                                                                                  axis.line.x = element_line(color = "black", 
                                                                                                                             size = axis.line.size), axis.line.y = element_line(color = "black", 
                                                                                                                                                                                size = axis.line.size), axis.title = element_text(size = axis.title.size), 
                                                                                                  axis.text = element_text(size = axis.text.size), 
                                                                                                  legend.position = "none", panel.border = element_blank(), 
                                                                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
    }
    else {
      dtmp_m <- melt(dtmp)
      colnames(dtmp_m) <- c("group", "time", "value")
      if (sum(dim(membership) == 0) == 2) {
        plotlist[[i]] <- ggplot(dtmp_m, aes(x = time, 
                                            y = value)) + geom_line(aes(group = group), 
                                                                    colour = cl.color) + theme_bw() + ggtitle(title) + 
          labs(x = categories, y = value) + theme(plot.title = element_text(size = title.size), 
                                                  axis.line.x = element_line(color = "black", 
                                                                             size = axis.line.size), axis.line.y = element_line(color = "black", 
                                                                                                                                size = axis.line.size), axis.title = element_text(size = axis.title.size), 
                                                  axis.text = element_text(size = axis.text.size), 
                                                  legend.position = "none", panel.border = element_blank(), 
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      if (!sum(dim(membership) == 0) == 2) {
        mem <- membership[cluster == i, i]
        mem1 <- data.frame(group = names(mem), member = mem)
        dtmp_m1 <- merge(dtmp_m, mem1, by = "group")
        colnames(dtmp_m1) <- c("group", "time", "value", 
                               "membership")
        dtmp_m1 <- dtmp_m1[order(dtmp_m1[, 4]), ]
        new.factor <- unique(as.vector(dtmp_m1$group))
        dtmp_m1$group <- factor(dtmp_m1$group, levels = new.factor)
        plotlist[[i]] <- ggplot(dtmp_m1, aes(x = time, 
                                             y = value, colour = membership)) + geom_line(aes(group = group)) + 
          scale_colour_gradientn(colours = membership.color) + 
          guides(colour = guide_colourbar()) + theme_bw() + 
          ggtitle(title) + labs(x = categories, y = value) + 
          theme(plot.title = element_text(size = title.size), 
                axis.line.x = element_line(color = "black", 
                                           size = axis.line.size), axis.line.y = element_line(color = "black", 
                                                                                              size = axis.line.size), axis.title = element_text(size = axis.title.size), 
                axis.text = element_text(size = axis.text.size), 
                legend.title = element_text(size = legend.title.size), 
                legend.text = element_text(size = legend.title.size), 
                panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
      }
    }
  }
  #suppressWarnings(multiplot(plotlist = plotlist, cols = cols))
  plotlist
}

## ------------------------------------------------------------------------
# Scatter plot, takes two DEseq resuts table as input
simple.replication.plot <- function(de.res.a, de.res.b, col.to.match="log2FoldChange", xlab="Fold change A", ylab="Fold change B", main="") {
  ol  <- intersect(rownames(de.res.a), rownames(de.res.b))
  min <- min(c(de.res.a[ol, col.to.match], de.res.b[ol, col.to.match]), na.rm=T)
  max <- max(c(de.res.a[ol, col.to.match], de.res.b[ol, col.to.match]), na.rm=T)
  plot(de.res.a[ol, col.to.match],
       de.res.b[ol, col.to.match],
       xlab=xlab,
       ylab=ylab,
       main=main,
       xlim=c(min, max),
       ylim=c(min, max))
  abline(a=0, b=1, col="red")
  abline(h=0, col="grey", lty=2)
  abline(v=0, col="grey", lty=2)
}