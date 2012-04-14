% File src/library/QuasiSeq/man/NBDev.Rd
     \name{NBDev}
     \alias{NBDev}
     \title{Compute negative binomial deviances (up to a constant) for given design matrix}
     \description{A function called within \code{QL.fit} to compute quasi-negative binomial deviances of each gene for given design matrix}
     \usage{
     NBDev(counts,design,log.offset,nb.disp,print.progress=TRUE)
     }
     \arguments{
      \item{counts}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
      \item{design}{A single element from the \code{design.list} argument given to \code{QL.fit}.}
	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include,\code{log.offset=log(colSums(counts))} or \code{log.offset=log(apply(counts,2,quantile,.75))}.  The default setting  in \code{QLfit} makes no adjustment for library sizes (i.e. log.offset=0).}
	\item{nb.disp}{estimated negative binomial dispersion parameters obtained from either \code{estimateGLMTrendedDisp} or \code{estimateGLMCommonDisp} in package \code{edgeR}.  These estimates are treated as known and are used to compute deviances.}
	\item{print.progress}{logical. If TRUE, the function will provide an update on what gene (row number) is being analyzed.  Updates occur frequently to start then eventually occur every 5000 genes.}
}
	
\value{list containing:
\item{dev}{vector containing the deviance for each gene under a Poisson model fit to design matrix (or vector, for one-factor experiments) specified by \code{design}. This vector is passed along within the \code{QL.fit} function.}
  \item{means}{matrix of fitted mean values for each gene}
\item{parms}{matrix of estimated coefficients for each gene}	
}

\author{Steve Lund \email{lunds@iastate.edu}}

     \keyword{RNA-seq, quasi-likelihood, differential expression}




