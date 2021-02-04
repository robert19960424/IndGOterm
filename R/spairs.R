#' get stable pairs from normal samples
#'
#' the function is used to get stable gene pairs from a host of normal samples.
#'
#' @param data1 a (normal sample) dataset which is used to screen the stable_pair.Notable that rownames is gene identifier(geneID/genesymbol/ensembleID).
#' @param gid1 a vector of gene indentifier which can completely match data1
#' @param threshold a threshold value is used to stipulate that a gene pair maintain the circumstance of > the threshold so that the gene pair should be retained.The threshold should between 0-1, and the default is 0.99.
#'
#' @return a list which contains a dataframe of stable pair and a runtime
#' @examples
#' geneid<-1:10;
#' samples<-runif(100,min = 0,max = 50);
#' data1<-matrix(c(geneid,samples),nrow = 10,byrow=F);
#' stable_pairs<-spairs(data1 = data1,gid1=geneid,threshold = 0.99)
#'
#' @export
#'
spairs=function(data1,gid1,threshold=0.99){
  timestart<-Sys.time()
  data1=as.matrix(data1)
  if (missing(threshold)){
    threshold=0.99
  }
  freqs=list()
  for(i in 1:(length(gid1)-1)){
    gid11=gid1[-c(1:i)]
    pair1=cbind(gid1[i],gid11)
    coms=data1[match(pair1[,1],gid1),,drop=F]-data1[match(pair1[,2],gid1),,drop=F]
    freq1=rowMeans(coms>0)
    freqs[[i]]=freq1
    rm(gid11,pair1,coms,freq1)
    gc()
  }
  pairs=t(combn(gid1,2))
  freqs=unlist(freqs)
  stablepair=rbind(pairs[freqs>threshold,],pairs[freqs<(1-threshold),c(2,1)])
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  get_spairs=list(stablepair=stablepair,runningtime=runningtime)
  return(get_spairs)
}
