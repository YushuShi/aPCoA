\name{aPCoA}
\alias{aPCoA}
\docType{package}
\title{
 Covariate Adjusted PCoA Plot
}
\description{
  Adjusted confounding covariates to show the effect of the primary covariate in a PCoA plot. This method is designed for non-Euclidean distance. This function will plot the original PCoA plot along with the covariate adjusted PCoA plot.
  }
\usage{
aPCoA(formula,data,maincov,drawEllipse=TRUE,drawCenter=TRUE)
}

\arguments{
  \item{formula}{A typical formula such as Y~ A, but here Y is a dissimilarity distance. The formula has the same requirements as in adonis function of the vegan package.}
  \item{data}{A dataset with the rownames the same as the rownames in distance. This dataset should include both the confounding covariate and the primary covariate.}
  \item{maincov}{the covariate of interest in the dataset, must be a factor}
  \item{drawEllipse}{Do you want to draw the 95\% confidence elipse for each cluster?}
  \item{drawCenter}{Do you want to show the connection between cluster center (medoid) and cluster members?}
}

\value{
 Two PCoA plots. One is the original one, while the other is the PCoA plot after adjusting for the confounding covariate.
 \item{plotMatrix}{The matrix for plotting the adjusted PCoA plot.}
}

\examples{
library(mvabund)
library(vegan)
library(aPCoA)
options(stringsAsFactors = FALSE)
data("Tasmania")
data<-data.frame(treatment=Tasmania$treatment,block=Tasmania$block)
bray<-vegdist(Tasmania$abund, method="bray")
rownames(data)<-rownames(as.matrix(bray))
opar<-par(mfrow=c(1,2),
    mar=c(3.1, 3.1, 3.1, 5.1),
    mgp=c(2, 0.5, 0),
    oma=c(0, 0, 0, 4))
result<-aPCoA(bray~block,data,treatment)
par(opar)
}

