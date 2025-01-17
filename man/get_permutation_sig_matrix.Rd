% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BACON_analysis.R
\name{get_permutation_sig_matrix}
\alias{get_permutation_sig_matrix}
\title{Calculate the communication strength matrix for a single synthesis-reception pair (with permutation test)}
\usage{
get_permutation_sig_matrix(
  qs_df,
  ligand_cell,
  target_cell,
  interaction_name,
  ligand_gene,
  ligand_complex,
  ligand_complex_number,
  target_gene,
  target_complex,
  target_complex_number,
  M,
  calculate_method = NULL
)
}
\arguments{
\item{qs_df}{data frame containing QS-related gene expression (rownames = bacteria; colnames = genes;the last column is bacteria group labels)}

\item{ligand_cell}{synthesis bacteria groups}

\item{target_cell}{reception bacteria groups}

\item{interaction_name}{name of the synthesis-reception groups communication}

\item{ligand_gene}{vector with gene symbols that are regarded as synthesis groups}

\item{ligand_complex}{vector indicating the groups of synthesis complex}

\item{ligand_complex_number}{vector indicating the stoichiometry of each group of synthesis complex for calculating ligand abundance}

\item{target_gene}{vector indicating the groups of reception complex}

\item{target_complex}{vector indicating the groups of receptor complex}

\item{target_complex_number}{vector indicating the stoichiometry of each group of receptor complex for calculating target abundance}

\item{M}{number of permutation tests}

\item{calculate_method}{}
}
\description{
Calculate the communication strength matrix for a single synthesis-reception pair (with permutation test)
}
