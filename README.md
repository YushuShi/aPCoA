# aPCoA
# aPCoA: Covariate Adjusted Principal Coordinates Analysis

In fields such as ecology, microbiology, and genomics, non-Euclidean distances are widely applied to describe pairwise dissimilarity between samples. Given these pairwise distances, principal coordinates analysis (PCoA) is commonly used to construct a visualization of the data. However, confounding covariates can make patterns related to the scientific question of interest difficult to observe. We provide 'aPCoA' as an easy-to-use tool to improve data visualization in this context, enabling enhanced presentation of the effects of interest.
## Main Function
```
aPCoA(formula,data,maincov,drawEllipse=TRUE,drawCenter=TRUE)
```
+ **formula** A typical formula such as Y~ A, but here Y is a dissimilarity distance. The formula has the same requirements as in adonis function of the vegan package.
+ **data** A dataset with the rownames the same as the rownames in distance. This dataset should include both the confounding covariate and the primary covariate.
+ **maincov** the covariate of interest in the dataset, must be a factor
+ **drawEllipse** Do you want to draw the 95% confidence elipse for each cluster?
+ **drawCenter** Do you want to show the connection between cluster center (medoid) and cluster members?

## Example
```
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
```
## Result
![](aPCoA.png)

## Authors

**Yushu Shi**, Liangliang Zhang, Kim-Anh Do, Christine Peterson and Robert Jenq

Department of Biostatistics and Department of Genomic Medicine, the University of Texas, MD Anderson Cancer Center
