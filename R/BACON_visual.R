#' Chord diagram of QS communication for single or multiple communication pairs
#'
#' @param object a BACON object
#' @param interaction_use QS communication indexes used for aggregation ("all" means use all communication pairs); for a single communication pair,set interaction_use as the communication name
#' @param synthesis_subpopulation synthesis group used for plot
#' @param reception_subpopulation reception group used for plot
#' @param lab.cex
#' @importFrom CellChat netVisual_chord_cell_internal
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
  CellChat::netVisual_chord_cell_internal(net_aggregated,lab.cex=lab.cex)
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
#' @importFrom grDevices recordPlot
#' @import CellChat
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
    color.use = CellChat::scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  # color.use.label  <- c('red','black')
  color.use.label <- c('black',CellChat::scPalette(length(unique(group))))[1:length(unique(group))]
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
