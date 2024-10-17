#' Generate color from a customed color palette
#'
#' @param n number of color
#' @importFrom grDevices colorRampPalette
#' @return
#' @export
#'
#' @examples
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}
#' Chord diagram for visualizing QS communication from a weighted adjacency matrix
#'
#' @param net weighted matrix with columns defining the QS network
#' @param color.use color use for groups
#' @param group group labels the sector name should be used as the names in the vector
#' @param cell.order defining the bacteria type orders
#' @param sources.use vector giving the index or the name of source bacterial groups
#' @param targets.use vector giving the index or the name of target bacterial groups
#' @param lab.cex font size for the text
#' @param small.gap small gap between sectors
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack height
#' @param remove.isolate whether remove sectors without any links
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied.
#' @param scale scale each sector to same width
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param transparency transparency of link colors
#' @param link.border border for links, single vector or a matrix with names or a data frame with source, target , prob, three column information
#' @param title.name title name of the plot
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom grDevices recordPlot colorRampPalette
#' @return
#' @export
#'
#' @examples
netVisual_chord_cell_internal <- function(net, color.use = NULL, group = NULL, cell.order = NULL,
                                          sources.use = NULL, targets.use = NULL,
                                          lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                          remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                          transparency = 0.4, link.border = NA,
                                          title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source","target")
  } else if (is.data.frame(net)) {
    if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source,net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob > 0)
  if(dim(net)[1]<=0){message("No interaction between those cells")}
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if(nrow(net) > nrow(net.fake)){
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      scale = TRUE
    }
  }

  df <- net
  cells.use <- union(df$source,df$target)

  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]

  # define grid color
  if (is.null(color.use)){
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector

  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }

  # define edge color
  edge.color <- color.use[as.character(df$source)]

  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type, # link.border = "white",
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               group = group,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  if(!is.null(title.name)){
    # title(title.name, cex = 1)
    text(-0, 1.02, title.name, cex=1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}

#' Chord diagram of QS communication for single or multiple communication pairs
#'
#' @param object a BACON object
#' @param interaction_use QS communication indexes used for aggregation ("all" means use all communication pairs); for a single communication pair,set interaction_use as the communication name
#' @param synthesis_subpopulation synthesis group used for plot
#' @param reception_subpopulation reception group used for plot
#' @param lab.cex
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom dplyr select %>% group_by summarize
#' @return
#' @export
#'
#' @examples
netVisual_chord_BACON <- function(object, interaction_use='all', synthesis_subpopulation=NULL,reception_subpopulation=NULL,lab.cex=1.3){
  if(is.null(synthesis_subpopulation)){synthesis_subpopulation=rownames(object@net[[1]])}
  if(is.null(reception_subpopulation)){reception_subpopulation=colnames(object@net[[1]])}
  if(interaction_use=='all') {
    net_aggregated <- net_aggregation(object@net)}
  else {net_aggregated <- net_aggregation(object@net[interaction_use])}
  netVisual_chord_cell_internal(net_aggregated,lab.cex=lab.cex)
}

#' Circle plot for QS communication network
#'
#' @param net_ori A weighted matrix representing the connections
#' @param title.name the name of the title
#' @param synthesis_sub vector giving the index or the name of the ligand bacterial groups
#' @param reception_sub vector giving the index or the name of the target bacterial groups
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of communications to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight the weight of vertex:either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param alpha.edge the transprency of edge
#' @param label.edge Whether or not shows the label of edges
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' @param vertex.size
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#'
#' @return
#' @export
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot colorRampPalette
#' @import grDevices
#' @examples
netVisual_circle_BACON <-function(net_ori, color.use = NULL,title.name = NULL, synthesis_sub = NULL, reception_sub = NULL, remove.isolate = FALSE, top = 1,
                                  weight.scale = TRUE, vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1.2,vertex.label.color= "black",
                                  edge.weight.max = NULL, edge.width.max=5, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                  edge.curved=0.2,shape=NULL,layout=in_circle(), margin=0.2, vertex.size = NULL,
                                  arrow.width=1,arrow.size = 0.8){
  net.names <- unique(c(rownames(net_ori),colnames(net_ori)))
  net <- matrix(0,nrow=length(net.names),ncol=length(net.names))
  dimnames(net) <- list(net.names,net.names)
  net[rownames(net_ori),colnames(net_ori)] <- net_ori
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  vertex.size.max <- if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) 5 else 15
  } else {
    vertex.size.max
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(synthesis_sub)) | (!is.null(reception_sub))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(synthesis_sub)){
      if (is.numeric(synthesis_sub)) {
        synthesis_sub <- cells.level[synthesis_sub]
      }
      df.net <- subset(df.net, source %in% synthesis_sub)
    }
    if (!is.null(reception_sub)){
      if (is.numeric(reception_sub)) {
        reception_sub <- cells.level[reception_sub]
      }
      df.net <- subset(df.net, target %in% reception_sub)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g1 <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
    g <-igraph::permute(g1,match(V(g1)$name,net.names))
    group <- structure(rep(1,length(V(g))),names=names(V(g)))
    group <- group[names(V(g))]

  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  # color.use.label  <- c('red','black')
  color.use.label <- c('black',scPalette(length(unique(group))))[1:length(unique(group))]
  color.use.label <- structure(color.use.label,names=unique(group))
  # color.use.label <- color.use.label[-length(color.use.label)]
  # if(length(color.use.label)==1){ color.use.label='black'}
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  if(length(vertex.weight)==1){ igraph::V(g)$size<-vertex.weight} else {igraph::V(g)$size<-vertex.weight[names(V(g))]}
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$shape <- c('circle')#, 'square', 'csquare', 'rectangle', 'crectangle', 'vrectangle', 'pie', 'raster','sphere')[group[igraph::V(g)]]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  #igraph::V(g)$frame.color <- border.color.use[igraph::V(g)]
  igraph::V(g)$label.color <-color.use.label[group[names(V(g))]]#unlist(lapply(V(g), FUN= function(x){color.use.label[group[x]]}))
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+4
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"


  if (!is.null(title.name)) {
    text(0,1.6,title.name, cex = 1.5)
  }

  gg <- recordPlot()
  return(gg)
}
