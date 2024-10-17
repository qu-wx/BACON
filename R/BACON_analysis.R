#' The BACON Class
#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom data.table data.table
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))
setClassUnion(name = 'AnyDF', members = c("data.frame"))
#' The key slots used in the BACON object are described below
#'
#' @slot data.raw raw count data matrix
#' @slot data gene count matrix for BACON analysis (Genes should be in rows and bacteria in columns)
#' @slot data_qs a data frame only containing QS related genes and bacteria sub-populations or species labels
#' @slot net0 a list containing the original communication strength matrix without filtered by p-value(0.05)
#' @slot pvalue a list containing the p-value (calculated by permutation test) matrix corresponding to net0
#' @slot net a list of containing the communication strength matrix (number of synthesis bacteria groups X number of reception bacteria groups) filted by p-value. Each row of the communication strength matrix indicates the communication probability originating from the synthesis bacteria group to other groups.
#' @slot fc deprecated
#' @slot info information of each synthesis-reception pair
#' @slot ligand.abundance ligand abundance of synthesis bacteria groups(row) for synthesis-reception pair (column)
#' @slot target.abundance target abundance of reception bacteria groups(row) for synthesis-reception pair (column)
#' @slot meta data frame containing information associated with each bacteria
#' @slot idents a factor defining the cell identity used for all analysis. It becomes a list for a merged BACON object
#' @slot DB synthesis-reception communication database used in the analysis
#' @slot LR a list of information related with synthesis-reception communication pairs
#'
#' @return
#' @export
#'
#' @examples
BACON <- methods::setClass("BACON",
                           slots = c(data.raw = 'AnyMatrix',
                                     data = 'AnyMatrix',
                                     data_qs = "AnyDF",
                                     net0 = "list",
                                     pvalue = "list",
                                     net = "list",
                                     fc = 'numeric',
                                     info ='numeric',
                                     ligand.abundance = 'matrix',
                                     target.abundance = 'matrix',
                                     meta = "data.frame",
                                     idents = "AnyFactor",
                                     DB = "list",
                                     LR = "character")
)
#' Create a new BACON object from a data matrix
#'
#' @param object a gene count data matrix (rownames = gene names; colnames = bacteria)
#' @param meta a data frame (rownames = bacteria barcode) consisting of bacteria information, which will be used for defining baceria groups
#' @param group.by a vector to indicate group annotation of bacteria
#' @param database BACON database selected by users
#'
#' @return
#' @export
#'
#' @examples
createBACON <- function(object, meta = NULL, group.by = NULL,database = NULL) {
  # data matrix as input
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    message("Create a BACON object from a data matrix")
    data <- object
    group.by <- if(is.null(group.by)){'labels'}else{group.by}
  }
  #database
  QSsystem_db <- database
  system_gene <- unlist(lapply(QSsystem_db, function(db){c(db$ligand_gene, db$target_gene)}))
  system_gene <- unique(system_gene)
  toupper(system_gene) -> system_gene
  toupper(rownames(data)) -> rownames(data)
  system_gene <- system_gene[system_gene %in% rownames(data)]
  data_qs <- as.data.frame(Matrix::t(data[system_gene,]))
  data_qs$bacteria_type <- group.by

  if (!is.null(meta)) {
    if (inherits(x = meta, what = c("matrix", "Matrix"))) {
      meta <- as.data.frame(x = meta)
    }
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame")
    }
    if (!identical(rownames(meta), colnames(data))) {
      cat("The bacteria barcodes in 'meta' is ", head(rownames(meta)),'\n')
      warning("The bacteria barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }

  object <- methods::new(Class = "BACON",
                         data = data,
                         meta = meta,
                         data_qs = data_qs,
                         DB = QSsystem_db,
                         LR = names(QSsystem_db),
                         idents = factor(group.by))
  return(object)
}


#' Calculate the average gene expression by bacterial sub-populations or species
#'
#' @param qs_df data frame containing QS-related gene expression (rownames = bacteria; colnames = genes;the last column is bacteria group labels)
#' @param gene_used gene symbols used to calculate average expression for bacteria groups
#' @param calculate_method  set "mean" as default
#' @import data.table
#' @return
#' @export
#'
#' @examples
expression_calculate_by_select_group <- function(qs_df,gene_used,calculate_method=NULL) {
  qs_df_used <- qs_df[,c(gene_used,'bacteria_type'),drop = F]
  if(is.null(calculate_method)){
    getx <- setDT(qs_df_used)[, lapply(.SD, mean), keyby = bacteria_type]
    gety <- getx[,gene_used,with=FALSE]
    rownames(gety) <- getx$bacteria_type
    getx <- gety
  } else if(calculate_method=='quantile') {
    quantile_1 <- setDT(qs_df_used)[, lapply(.SD, quantile,0.25), keyby = bacteria_type]
    quantile_2 <- setDT(qs_df_used)[, lapply(.SD, quantile,0.5), keyby = bacteria_type]
    quantile_3 <- setDT(qs_df_used)[, lapply(.SD, quantile,0.75), keyby = bacteria_type]

    getx <- 0.25*quantile_1[,gene_used,with=FALSE] +
      0.5*quantile_2[,gene_used,with=FALSE] +
      0.25*quantile_3[,gene_used,with=FALSE]
    rownames(getx) <- quantile_1$bacteria_type
  }
  return(getx)
}

#' Calculate the communication strength matrix for a single synthesis-reception pair (without permutation test)
#'
#' @param qs_df data frame containing QS-related gene expression (rownames = bacteria; colnames = genes;the last column is bacteria group labels)
#' @param ligand_cell synthesis bacteria groups
#' @param target_cell reception bacteria groups
#' @param synthsis_complex_new vector with updated gene symbols (removing those excluded from the expression data) that are regarded as synthesis complex
#' @param ligand_complex vector indicating the groups of synthesis complex
#' @param ligand_complex_number vector indicating the stoichiometry of each group of synthesis complex for calculating ligand abundance
#' @param receptor_complex_new vector with updated gene symbols (removing those excluded from the expression data) that are regarded as receptor complex
#' @param target_complex vector indicating the groups of receptor complex
#' @param target_complex_number vector indicating the stoichiometry of each group of receptor complex for calculating target abundance
#'
#' @return
#' @export
#'
#' @examples
get_orig_strength_matrix <- function(qs_df,ligand_cell,target_cell,synthsis_complex_new,ligand_complex,ligand_complex_number,receptor_complex_new,target_complex,target_complex_number,calculate_method=NULL){
  gene_used <- c(synthsis_complex_new,receptor_complex_new)
  true_exp_genes <- expression_calculate_by_select_group(qs_df,gene_used,calculate_method = calculate_method) #use mean method to measure gene expression level
  bacteria_rownames <- rownames(true_exp_genes)
  true_exp_genes <- as.data.frame(true_exp_genes)
  cellgroup_number <- nrow(true_exp_genes)
  ligand_group_expr_mtx <- matrix(0,nrow = cellgroup_number,ncol = length(ligand_complex_number))
  for(i in 1:length(ligand_complex_number)){
    ind_i <- which(ligand_complex==i)
    ligand_group_expr_mtx[1:cellgroup_number,i] <- true_exp_genes[1:cellgroup_number,ind_i]
  }
  if(length(ligand_complex_number)>1){
    rep_number <- rep(1:length(ligand_complex_number),ligand_complex_number)
    ligand_group_expr <- exp(rowMeans(log(ligand_group_expr_mtx[, rep_number])))} else {
      ligand_group_expr <- ligand_group_expr_mtx
    }
  names(ligand_group_expr) <- bacteria_rownames
  # receptor
  target_group_expr_mtx <- matrix(0,cellgroup_number,length(target_complex_number))
  for(i in 1:length(target_complex_number)){
    ind_i <- which(target_complex==i)
    target_group_expr_mtx[1:cellgroup_number,i] <- true_exp_genes[1:cellgroup_number,length(synthsis_complex_new)+ind_i]

  }

  if(length(target_complex_number)>1){
    rep_number <- rep(1:length(target_complex_number),target_complex_number)
    target_group_expr <- exp(rowMeans(log(target_group_expr_mtx[, rep_number])))} else {
      target_group_expr <- target_group_expr_mtx
    }
  names(target_group_expr) <- bacteria_rownames
  lig_tar_mtx <- ((as.matrix(ligand_group_expr))%*%t(as.matrix(target_group_expr)))
  colnames(lig_tar_mtx) <- bacteria_rownames
  rownames(lig_tar_mtx) <- bacteria_rownames
  lig_tar_mtx_n <- lig_tar_mtx
  rownames(lig_tar_mtx_n) <- bacteria_rownames
  colnames(lig_tar_mtx_n) <- bacteria_rownames
  lig_tar_mtx_n <- lig_tar_mtx_n[ligand_cell,target_cell]
  lig_tar_mtx <- lig_tar_mtx[ligand_cell,target_cell]
  result_list <- list("lig_tar_mtx_n" = lig_tar_mtx_n,
                      "ligand_score" = ligand_group_expr,
                      "target_score" = target_group_expr,
                      "bacteria_group"= bacteria_rownames,
                      "lig_tar_mtx" = lig_tar_mtx)
  return(result_list)
}

#' Calculate the communication strength matrix for a single synthesis-reception pair (with permutation test)
#'
#' @param qs_df data frame containing QS-related gene expression (rownames = bacteria; colnames = genes;the last column is bacteria group labels)
#' @param ligand_cell synthesis bacteria groups
#' @param target_cell reception bacteria groups
#' @param interaction_name name of the synthesis-reception groups communication
#' @param ligand_gene vector with gene symbols that are regarded as synthesis groups
#' @param ligand_complex vector indicating the groups of synthesis complex
#' @param ligand_complex_number vector indicating the stoichiometry of each group of synthesis complex for calculating ligand abundance
#' @param target_gene vector indicating the groups of reception complex
#' @param target_complex vector indicating the groups of receptor complex
#' @param target_complex_number vector indicating the stoichiometry of each group of receptor complex for calculating target abundance
#' @param M number of permutation tests
#' @param calculate_method
#'
#' @return
#' @export
#'
#' @examples
get_permutation_sig_matrix <- function(qs_df,ligand_cell,target_cell,interaction_name,ligand_gene,ligand_complex,ligand_complex_number,
                                       target_gene,target_complex,target_complex_number,M,calculate_method=NULL){
  ind_lig <- which(ligand_gene %in% colnames(qs_df));synthsis_complex_new <- ligand_gene[ind_lig]
  ind_tar <- which(target_gene %in% colnames(qs_df));receptor_complex_new <- target_gene[ind_tar]
  ligand_complex <- ligand_complex[ind_lig]
  target_complex <- target_complex[ind_tar]
  result_list <- get_orig_strength_matrix(qs_df,ligand_cell,target_cell,synthsis_complex_new,ligand_complex,ligand_complex_number, receptor_complex_new,target_complex,target_complex_number,calculate_method)

  lig_tar_mtx_n <- result_list$lig_tar_mtx_n
  FC.lig <- max(result_list$ligand_score) - min(result_list$ligand_score)
  FC.tar <- max(result_list$target_score) - min(result_list$target_score)
  FC <- max(FC.lig,FC.tar)
  ## permutation
  if(M==0 | max(lig_tar_mtx_n)==0){net <- lig_tar_mtx_n;pvalue <- NA*net} else {
    pvalue <- matrix(0, dim(lig_tar_mtx_n)[1],dim(lig_tar_mtx_n)[2]);
    df_j <- qs_df;
    prob_mtx_permutation <- sapply(seq_along(1:M), function(j) {
      df_j$bacteria_type[qs_df$bacteria_type %in% ligand_cell] <- sample(qs_df$bacteria_type[qs_df$bacteria_type %in% ligand_cell],length(qs_df$bacteria_type[qs_df$bacteria_type %in% ligand_cell]), FALSE)
      df_j$bacteria_type[qs_df$bacteria_type %in% target_cell] <- sample(qs_df$bacteria_type[qs_df$bacteria_type %in% target_cell],length(qs_df$bacteria_type[qs_df$bacteria_type %in% target_cell]), FALSE)
      result_list_n <- get_orig_strength_matrix(df_j,ligand_cell,target_cell,synthsis_complex_new,ligand_complex,ligand_complex_number,receptor_complex_new,target_complex,target_complex_number,calculate_method = calculate_method)
      prob_mtx_permutation_j <- (result_list_n$lig_tar_mtx_n > lig_tar_mtx_n)*1
      prob_mtx_permutation_j
    },simplify = 'array')
    pvalue <- apply(prob_mtx_permutation,1:2,sum)
    pvalue <- pvalue/M
    pvalue_v <- c(pvalue)
    pvalue_v_sort <- sort.int(pvalue_v, decreasing = FALSE, index.return = TRUE)
    k <- which(pvalue_v_sort$x<0.05)
    prob_mtx_sig <- 0*pvalue
    prob_mtx_sig[pvalue_v_sort$ix[k]]<-lig_tar_mtx_n[pvalue_v_sort$ix[k]]
    net <- prob_mtx_sig;
  }
  list_return <- list(net=net,FC=FC,pvalue=pvalue,net0=result_list$lig_tar_mtx_n, info=sum(net),ligand.abundance = result_list$ligand_score, target.abundance = result_list$target_score)
  return(list_return)

}

#' Calculate the communication strength matrix for all synthesis-reception pair
#'
#' @param object a BACON object
#' @param ligand_cell synthesis bacteria groups
#' @param target_cell reception bacteria groups
#' @param M number of permutation tests
#'
#' @return
#' @export
#'
#' @examples
run_BACON <- function(object=object,ligand_cell=NULL,target_cell=NULL,M=100,calculate_method=NULL){
  start.time <- Sys.time()
  QSsystem_db <- object@DB
  net0_all <- list() # original probability mtx
  net_all <- list() # probability mtx filtered by pvalue with cutoff 0.05
  pvalue_all <- list()
  FC_all <- rep(0,length(QSsystem_db))
  info_all <- rep(0,length(QSsystem_db))
  net.rownames <- sort(unique(object@data_qs$bacteria_type),method='radix')
  if(is.null(ligand_cell)){ligand_cell <- net.rownames};if(is.null(target_cell)){target_cell <- net.rownames}
  ligand.abundance_all <- matrix(0,nrow=length(net.rownames),ncol=length(QSsystem_db))
  target.abundance_all <- ligand.abundance_all
  for(j in 1:length(QSsystem_db)){
    interaction_name <- names(QSsystem_db)[j]
    ligand_gene <- QSsystem_db[[j]]$ligand_gene
    ligand_gene <- toupper(ligand_gene)
    target_gene <- QSsystem_db[[j]]$target_gene
    target_gene <- toupper(target_gene)
    ligand_complex <- QSsystem_db[[j]]$ligand_complex
    ligand_complex_number <- QSsystem_db[[j]]$ligand_complex_number
    target_complex <- QSsystem_db[[j]]$target_complex
    target_complex_number <- QSsystem_db[[j]]$target_complex_number
    lig_boolean_group <- (ligand_gene %in% names(object@data_qs))*ligand_complex;
    rec_boolean_group <- (target_gene %in% names(object@data_qs))*target_complex;
    lig_boolean <- prod(unique(ligand_complex) %in% lig_boolean_group)
    rec_boolean <- prod(unique(target_complex) %in% rec_boolean_group)

    if(prod(lig_boolean,rec_boolean)==0){net_all[[j]]=matrix(0,nrow=length(ligand_cell),ncol=length(target_cell),
                                                             dimnames=list(ligand_cell,target_cell))

    net0_all[[j]] <- net_all[[j]]
    pvalue_all[[j]] <- net_all[[j]]
    } else
    {tmp <- get_permutation_sig_matrix(object@data_qs,ligand_cell,target_cell,interaction_name,ligand_gene,ligand_complex,ligand_complex_number,target_gene,target_complex,M=100,target_complex_number,calculate_method = calculate_method)
    net_all[[j]] <- tmp$net
    net0_all[[j]] <- tmp$net0
    pvalue_all[[j]] <- tmp$pvalue
    FC_all[j] <- tmp$FC
    info_all[j] <- tmp$info
    ligand.abundance_all[1:length(net.rownames),j] <- tmp$ligand.abundance
    target.abundance_all[1:length(net.rownames),j] <- tmp$target.abundance
    }
  }
  rownames(ligand.abundance_all) <- names(tmp$ligand.abundance);colnames(ligand.abundance_all) <- names(QSsystem_db)
  rownames(target.abundance_all) <- names(tmp$target.abundance);colnames(target.abundance_all) <- names(QSsystem_db)
  names(net_all) <- names(QSsystem_db)
  names(net0_all) <- names(QSsystem_db)

  end.time <- Sys.time();time.taken <- end.time - start.time;
  print(time.taken)

  object@net0 <- net0_all
  object@pvalue <- pvalue_all
  object@net <- net_all
  object@fc <- FC_all
  object@info <- info_all
  object@ligand.abundance <- ligand.abundance_all
  object@target.abundance <- target.abundance_all
  return(object)
}

#' Aggregation of communication networks over all synthesis-reception pairs
#'
#' @param net_list list containing communication strength matrix for all synthesis-reception pairs
#'
#' @return
#' @export
#'
#' @examples
net_aggregation <- function(net_list){
  net_aggregated <- 0*net_list[[1]]
  for(jj in 1:length(net_list)){
    net_aggregated <- net_aggregated+sum(net_list[[jj]])*(net_list[[jj]]>0);
  }
  return(net_aggregated)
}
