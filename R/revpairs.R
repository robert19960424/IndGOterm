#' get reverse pairs from all the tumor samples
#'
#' the function is used to get reverse gene pairs from a host of disease samples.
#' this function support parallel computation . You need to set thread numbers to make sure how many threads do you need to use .
#'
#' @param stable.pair a matrix or data.frame of stable pairs from normal samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param threshold a numeric value which is used to control false discovery rate under the p_value of the chip-square test or fisher's exact probability test , default is 0.05.
#' @param threads an integer value to make sure how many threads you will use to complete the computation
#' @param patients a matrix of data.frame of tumor samples , the first column must be the geneID ,and tumor samples start with the second column.
#' @param spairs_threshold a threshold same with the "threshold" in function "spairs".
#' @param capacity an integer value to depict the computation capacity, ruling how many lines of stable pairs would be computed within one time.the default is 300000
#'
#' @return a matrix containing four columns respectively represent higher expression gene , lower expression gene , p_value under binomial distribution , false discovery rate under p.adjust
#'
#' @import parallel
#' @import data.table
#' @examples
#' stable.pair<-t(combn(sample(1:10,10),2));
#' geneid<-1:10;
#' samples<-runif(100,min = 0,max = 50);
#' patients<-matrix(c(geneid,samples),nrow = 10,byrow=F);
#' reverse_pairs<-revpairs(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=1L,capacity=300000L)
#' #compute with parallel
#' reverse_pairs<-revpairs(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=10L,capacity=300000L)
#' @export
#'
revpairs<-function(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=1L,capacity=300000L){
  if(missing(stable.pair)){
    stop("argument 'stable.pair' of two columns is missing.")
  }
  if(missing(patients)){
    stop("argument 'patients' is missing.")
  }
  if(dim(stable.pair)[2]!=2){
    stop("stable.pair should only with two columns")
  }
  if(missing(threshold)){
    threshold<-0.05
  }
  if(missing(capacity)){
    warning("capacity is set to 300000")
  }
  stable.pair<-data.table(stable.pair)
  setnames(stable.pair,c('higher',"lower"))
  patients<-data.table(patients)
  colnames(patients)[1]<-'GeneID'
  stab_gene1<-stable.pair$higher
  stab_gene2<-stable.pair$lower
  index<-patients$GeneID
  patients<-patients[,!'GeneID']
  reverse_num<-c()
  i=1
  while (i<=ceiling(length(stab_gene1)/capacity)) {
    if((capacity*i)>length(stab_gene1)){
      gene1<-patients[match(stab_gene1[((i-1)*capacity+1):length(stab_gene1)],index),]
      gene2<-patients[match(stab_gene2[((i-1)*capacity+1):length(stab_gene2)],index),]
      diif<-gene1-gene2
      reverse_num<-c(reverse_num,rowSums(diif<0))
      rm(gene1,gene2,diif)
      gc()
    }else{
      gene1<-patients[match(stab_gene1[((i-1)*capacity+1):(capacity*i)],index),]
      gene2<-patients[match(stab_gene2[((i-1)*capacity+1):(capacity*i)],index),]
      diif<-gene1-gene2
      reverse_num<-c(reverse_num,rowSums(diif<0))
      rm(gene1,gene2,diif)
      gc()
    }
    print(paste(i,"/",ceiling(length(stab_gene1)/capacity),sep = ""))
    i=i+1
  }
  rm(stab_gene1,stab_gene2,i)
  gc()
  binom_value<-function(x,sample_num,spairs_threshold){
    y<-pbinom(x-1,sample_num,1-spairs_threshold,lower.tail = F)
  }
  if(threads==1){
    binom_value1<-c()
    for (i in 1:length(reverse_num)) {
      binom_value1<- c(binom_value1,pbinom(reverse_num[i]-1,dim(patients)[2],1-spairs_threshold,lower.tail = F))
    }
  }else{
    c1<-makeCluster(threads)
    binom_value1<-parSapply(c1,reverse_num,binom_value,sample_num=dim(patients)[2],spairs_threshold=spairs_threshold)
    stopCluster(c1)
  }
  stable.pair<-as.matrix(stable.pair)
  stable_pair_pvalue_unadjusted<-cbind(stable.pair,binom_value1,reverse_num)
  stable_pair_pvalue_unadjusted<-stable_pair_pvalue_unadjusted[order(stable_pair_pvalue_unadjusted[,3]),]
  p_value_adjusted<-p.adjust(stable_pair_pvalue_unadjusted[,3],method = "fdr")
  stable_pair_pvalue_adjusted<-cbind(stable_pair_pvalue_unadjusted,p_value_adjusted)
  index<-which(stable_pair_pvalue_adjusted[,5]<threshold)
  reverse_pair<-stable_pair_pvalue_adjusted[index,]
  reverse_pair<-reverse_pair[,c(2,1,3,5,4)]
  colnames(reverse_pair)<-c("higher","lower","p_value","fdr","reverse_num")
  return(reverse_pair)
}
