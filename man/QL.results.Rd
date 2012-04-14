% File src/library/QuasiSeq/man/QL.results.Rd
     \name{QL.results}
     \alias{QL.results}
     \title{Obtain p-values and q-values using results from \code{QL.fit}}
     \description{
       Obtain significance results for quasi-likelihood model fits to RNA-seq expression count data using the methods detailed in Lund, Nettleton, McCarthy, and Smyth (2012).}
     \usage{	
	QL.results(fit,Dispersion="Deviance",spline.df=NULL,Plot=FALSE)
	     }
     \arguments{
      \item{fit}{The list returned by the function QL.fit}
	\item{Dispersion}{Must be one of "Deviance" or "Pearson", specifying which type of estimator should be used for estimating quasi-likelihood dispersion parameter.}
	\item{spline.df}{Optional.  User may specify the degrees of freedom to use when fitting a cubic spline to log-scale(estimated dispersion) versus the log(average count).  Default uses cross-validation in \code{sreg} function to pick appropriate degrees of freedom.} 	
	\item{Plot}{logical.  If TRUE, the estimated dispersion versus the average count are plotted on a log-scale with the corresponding cubic spline fit overlaid. }
	}

\value{list containing:
	\item{"P.values"}{list of matrices providing p-values.  The i^th element of \code{pvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
If there is only one test for each gene, \code{P.values} is returned as a matrix, rather than a list.  Each matrix contains three columns, one for each of the following approaches: QL, QLShrink and QLSpline.  Note that column names for QLShrink and QLSpline methods contain the number of additional denominator d.f. gained from shrinking dispersion estimates.}
	\item{"Q.values"}{list of matrices providing q-values.  The i^th element of \code{qvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
If there is only one test for each gene, \code{Q.values} is returned as a matrix, rather than a list.  Each matrix contains three columns, one for each of the following approaches: QL, QLShrink and QLSpline.  Q-values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
	\item{"m0"}{matrix providing estimated number of true null hypotheses for each test(arranged by row) under each of the three methods(arranged by column). \code{m0} values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
}

\author{Steve Lund \email{lunds@iastate.edu}}

\examples{  
### Create example data set from a gammma-Poisson model.  
trt<-rep(1:2,each=5);n<-length(trt); n.genes<-2000

#### Simulate means for 2000 genes
sim.mn<-matrix(exp(rnorm(n.genes,2.5,2)),n.genes,2)

#### make first 1000 genes DE
sim.mn[1:(.5*n.genes),1]<-exp(rnorm(.5*n.genes,2.5,2))

### Simulate dispersion parameters
phi<-rgamma(n.genes,2,.7)+1

b<-1/(phi-1); a<-sim.mn[,trt]*b

### Simulate library size factors
offset<-2^(rnorm(n,0,.5))
a<-t(t(a)*offset)
lambda<-matrix(rgamma(n.genes*n,a,rep(b,n)),n.genes,n)
simdat<-matrix(rpois(n.genes*n,lambda),n.genes,n)

### Keep genes with at least 10 total counts
simdat<-simdat[rowSums(simdat)>9,]

### Create list of designs describing model under null and alternative hypotheses
design.list<-vector("list",2)
design.list[[1]]<-model.matrix(~as.factor(trt))  #This also could have just been ``trt''.
design.list[[2]]<-rep(1,length(trt))

log.offset<-log(apply(simdat,2,quantile,.75))

### Analyze using QL, QLShrink, and QLSpline methods applied to quasi-Poisson model
fit<-QL.fit(simdat, design.list,log.offset=log.offset, Model="Poisson")
results<-QL.results(fit)

### How many significant genes at FDR=.05?
apply(results$Q.values<.05,2,sum)

### Indexes for Top 100 most significant genes from QLSpline method
order(results$P.values[,3])[1:100]

### Analyze using QL, QLShrink, and QLSpline methods 
### applied to quasi-negative binomial model
fit2<-QL.fit(simdat, design.list,log.offset=log.offset, Model="NegBin")
results2<-QL.results(fit2)

### How many significant genes at FDR=.05?
apply(results2$Q.values<.05,2,sum)

### Indexes for Top 100 most significant genes from QLShrink method
order(results2$P.values[,2])[1:100]

}
     \keyword{RNA-seq, quasi-likelihood, differential expression}



