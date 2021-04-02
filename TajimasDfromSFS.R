### sfs is the site frequency spectrum for 4-fold degenerate sites
### sfs excludes polymorphic sites.  sfs is vector.
### It doesn't matter whether you are using the folded or unfolded sfs
### First entry is number of sites that are singletons, 2nd is number of doubletons, etc.
### sfs should look something like this c(25, 15, 9, 4, 4, 3, 2, 2, 0, 2, 1, 0, 0, 1, 0)
### n is the sample size.  the length of the sfs vector should never be greater than (n-1)
### L is the is the number of 4-fold sites (both polymorphic and non-polymorphic)

GetTajimasD<-function(sfs, L, n){  ## following p. 303 of Walsh and Lynch book.
  S<- sum(sfs)/L
  an <- sum(sapply(1:(n-1), function(x) 1/x))  ## p. 301
  bn <-  sum(sapply(1:(n-1), function(x) 1/(x^2))) ## p. 301
  thetaW = S/an
  betaD <- (1/(an^2 + bn))*( 2*(n^2 + n + 3)/(9*n*(n-1)) - (n+2)/(an * n) + bn/(an^2) ) ## p. 304
  alphaD<- (1/an)*( (n+1)/(3*(n-1)) - (1/an)) - betaD ## p. 304

  PI = sum(sapply(1:length(sfs), function(x) sfs[x]*2*(x/n)*((n - x)/(n - 1)))) / L

  TajD = (PI - thetaW)/sqrt(alphaD*S + betaD*(S^2))
  return(TajD)
}

#exampleSFS <- c(25, 15, 9, 4, 4, 3, 2, 2, 0, 2, 1, 0, 0, 1, 0)
#GetTajimasD(exampleSFS, 15)


GetTajimasD.NotSizeAdjusted<-function(sfs, n){  ## following p. 303 of Walsh and Lynch book.
  S<- sum(sfs)
  an <- sum(sapply(1:(n-1), function(x) 1/x))  ## p. 301
  bn <-  sum(sapply(1:(n-1), function(x) 1/(x^2))) ## p. 301
  thetaW = S/an
  betaD <- (1/(an^2 + bn))*( 2*(n^2 + n + 3)/(9*n*(n-1)) - (n+2)/(an * n) + bn/(an^2) ) ## p. 304
  alphaD<- (1/an)*( (n+1)/(3*(n-1)) - (1/an)) - betaD ## p. 304

  PI = sum(sapply(1:length(sfs), function(x) sfs[x]*2*(x/n)*((n - x)/(n - 1))))

  TajD = (PI - thetaW)/sqrt(alphaD*S + betaD*(S^2))
  return(TajD)
}


#vcf1<-read.csv2("/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexDifferences_IsoFormUsage/gene1.txt", sep = "\t")
#vcf2<-read.csv2("/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexDifferences_IsoFormUsage/gene2.txt", sep = "\t")
#vcf3<-read.csv2("/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexDifferences_IsoFormUsage/gene3.txt", sep = "\t")
#L.FakeValueForTesting<-100

### The function below does the following
### 1) Excludes sites with fewer some threshold level of samples (used 100)
### and rescales L accordingly, i.e., if 5% of sites are dropped then L_adjusted = 0.95*L.
### 2) Of retained sites, get minimum n across sites (nmin)
### 3) Subsample:For each retained site, sample without replacement nmin times
### 4) Make SFS from this this subsampled data and calculate Tajima's D
### 5) Repeat steps 3-4 100 times and returns the median D and the corresponding S


GetTajimasDFromVCF<-function(vcf, L){
  vcf = as.data.frame(vcf)
  sampleSizePerSite<- sapply(1:(dim(vcf)[1]), function(x) sum(!is.na(vcf[x,])))
  minAllowableSampleSize<-75
  numExcludedSites<-sum(sampleSizePerSite < minAllowableSampleSize)

  vcf.Filtered <-vcf[sampleSizePerSite>=minAllowableSampleSize, ]
  if (nrow(vcf.Filtered) > 0){
    sampleSizePerSite.UsableSites<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(!is.na(vcf.Filtered[x,])))
    numNonRefSamplesPerSite.UsableSites<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(vcf.Filtered[x,], na.rm = TRUE))
    minSampleSize<-min(sampleSizePerSite.UsableSites)

    nSubSamples = 100; listOfTajDValues<-rep(NA, nSubSamples); listOfSValues<-rep(NA, nSubSamples)
    for(i in 1:nSubSamples){
      nNonRefSNPsPerSite<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(sample(c(rep(0, sampleSizePerSite.UsableSites[x]- numNonRefSamplesPerSite.UsableSites[x]),
                                                              rep(1, numNonRefSamplesPerSite.UsableSites[x])), size = minSampleSize, replace = FALSE)))
      sfs<-rep(NA, dim(vcf.Filtered)[2])
      for(j in 1:length(sfs)) sfs[j] = sum(nNonRefSNPsPerSite == j)
      listOfTajDValues[i]<-GetTajimasD(sfs, L*(1- (numExcludedSites/(dim(vcf)[1]))), minSampleSize)
      listOfSValues[i]<-sum(sfs)
    }
    returnThisTajD<-median(listOfTajDValues)
    indexOfMedian<-which(abs(listOfTajDValues - median(listOfTajDValues)) == min(abs(listOfTajDValues - median(listOfTajDValues))))
    returnThisS<- mean(listOfSValues[indexOfMedian])
    return(c(returnThisTajD, returnThisS))
  } else {
    return(c(NA,nrow(vcf)))
  }
}

# test it
#GetTajimasDFromVCF(vcf1, L.FakeValueForTesting)
#GetTajimasDFromVCF(vcf2, L.FakeValueForTesting)
#GetTajimasDFromVCF(vcf3, L.FakeValueForTesting)








#####################################################################################################
########################################## DEPRECATED CODE ##########################################
#####################################################################################################


### sfs is the site frequency spectrum for 4-fold degenerate sites
### sfs excludes polymorphic sites.  sfs is vector.
### It doesn't matter whether you are using the folded or unfolded sfs
### First entry is number of sites that are singletons, 2nd is number of doubletons, etc.
### sfs should look something like this c(25, 15, 9, 4, 4, 3, 2, 2, 0, 2, 1, 0, 0, 1, 0)
### n is the sample size.  the length of the sfs vector should never be greater than (n-1)
### L is the is the number of 4-fold sites (both polymorphic and non-polymorphic)

#GetTajimasD<-function(sfs, L, n){  ## following p. 303 of Walsh and Lynch book.
#  S<- sum(sfs)/L
#  an <- sum(sapply(1:(n-1), function(x) 1/x))  ## p. 301
#  bn <-  sum(sapply(1:(n-1), function(x) 1/(x^2))) ## p. 301
#  thetaW = S/an

#  betaD <- ( (1/(an^2 + bn)) * (2*(n^2 + n + 3))) / ( (9*n*(n-1)) - ((n+2)/(an * n)) + (bn/(an^2)) ) ## p. 304
#  alphaD<- ( (1/an)*( ((n+1)/(3*(n-1))) - (1/an) ) ) - betaD ## p. 304
#
#  PI = sum(sapply(1:length(sfs), function(x) sfs[x]*2*(x/n)*((n - x)/(n - 1)))) / L
#
#  TajD = (PI - thetaW)/sqrt(alphaD*S + betaD*(S^2))
#
#  return(TajD)
#}

#exampleSFS <- c(2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#GetTajimasD(exampleSFS,15,20)


#GetTajimasD.NotSizeAdjusted<-function(sfs, n){  ## following p. 303 of Walsh and Lynch book.
#  S<- sum(sfs)
#  an <- sum(sapply(1:(n-1), function(x) 1/x))  ## p. 301
#  bn <-  sum(sapply(1:(n-1), function(x) 1/(x^2))) ## p. 301
#  thetaW = S/an
#  betaD =
#  betaD <- (1/(an^2 + bn))*( 2*(n^2 + n + 3)/(9*n*(n-1)) - (n+2)/(an * n) + bn/(an^2) ) ## p. 304
#  alphaD<- (1/an)*( (n+1)/(3*(n-1)) - (1/an)) - betaD ## p. 304

#  PI = sum(sapply(1:length(sfs), function(x) sfs[x]*2*(x/n)*((n - x)/(n - 1))))
#
#  TajD = (PI - thetaW)/sqrt(alphaD*S + betaD*(S^2))
#  return(TajD)
#}
