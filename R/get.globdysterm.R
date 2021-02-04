#' get holistic dyregulated pathway
#'
#' the function is used to get global dysregulation terms from a host of disease samples
#' this function support parallel computation . You need to set thread numbers to make sure how many threads do you need to use.
#'
#' @param stable.pair a matrix or data.frame of stable pairs from normal samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param reverse.pair a matrix or data.frame of reverse pairs from disease samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param patients a matrix of data.frame of tumor samples , the first column must be the geneID ,and tumor samples start with the second column.
#' @param threshold a numeric value which is used to control false discovery rate under the p_value of the chip-square test or fisher exact probability test , default is 0.05.
#' @param threads an integer value to make sure how many threads you will use to complete the computation.
#' @param gene.minimum an integer value to control the minimum amount of the genes of the terms. The deault is 1.
#' @param term_gene_list a list of term genes. list type. each list is a vector type, containing all the genes of a term. By default the variables are "BPtermgene_list" and taken from environment(formula).
#'
#' @return a data.frame containing three columns represent "termID" , "p value" , "fdr value" respectively
#' @import parallel
#' @import data.table
#' @examples
#' stable.pair<-t(combn(sample(1:10,10),2));
#' geneid<-1:10;
#' samples<-runif(100,min = 0,max = 50);
#' patients<-matrix(c(geneid,samples),nrow = 10,byrow=F);
#' reverse_pairs<-revpairs(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=1L,capacity=300000L)
#' global_dysterm<-get.globdysterm(stable.pair,reverse_pairs,threads = 1,gene.minimum = 2,term_gene_list = BPtermgene_list,patients = patients,threshold = 0.05)
#' @export
#'
get.globdysterm<-function(stable.pair,reverse.pair,threads=1L, gene.minimum=1, term_gene_list=BPtermgene_list,patients,threshold=0.05){
  if(missing(stable.pair)){
    stop("argument 'stable.pair' of two columns is missing.")
  }
  if(missing(reverse.pair)){
    stop("argument 'reverse.pair' of two coloums is missing.")
  }
  if(missing(threads)){
    threads=1L
  }
  if(missing(threshold)){
    threshold<-0.05
  }
  if(missing(patients)){
    stop("argument 'patients' is missing. This argument should be completely same with the argument 'patients' used in function 'revpairs'.")
  }
  if(missing(term_gene_list)){
    warning("argument 'term_gene_list' is missing, default is 'BPtermgene_list'")
  }
  if(class(term_gene_list)!="list"){
    warning("argument 'term_gene_list' must be list type")
  }
  for(i in 1:length(term_gene_list)){
    if(length(intersect(patients[,1],term_gene_list[[i]]))!=length(term_gene_list[[i]])){
      term_gene_list[[i]]<-intersect(patients[,1],term_gene_list[[i]])
    }
  }
  index<-unlist(lapply(term_gene_list,function(x) length(x)))
  index1<-which(index>=gene.minimum)
  term_gene_list<-term_gene_list[index1]
  reverse.pair<-as.matrix(reverse.pair[,1:2])
  rm(list = setdiff(ls(),c("term_gene_list","stable.pair","reverse.pair","threshold","threads")))
  gc()
  path.dyregulated.in.entirety<-function(path.gene,stable.pair,reverse.pair){
    library(data.table)
    stable.pair<-data.table(stable.pair)
    reverse.pair<-data.table(reverse.pair)
    setnames(stable.pair,c("higher","lower"))
    setnames(reverse.pair,c("higher","lower"))
    pathgene_stable_pair<-stable.pair[higher%in%path.gene | lower%in%path.gene][!(higher%in%path.gene & lower%in%path.gene)]
    if(dim(pathgene_stable_pair)[1]!=0){
      a<-dim(pathgene_stable_pair[higher%in%path.gene])[1]
      b<-dim(pathgene_stable_pair[lower%in%path.gene])[1]
      pathgene_reverse_pair<-reverse.pair[higher%in%path.gene | lower%in%path.gene][!(higher%in%path.gene & lower%in%path.gene)]
      if(dim(pathgene_reverse_pair)[1]!=0){
        c1<-dim(pathgene_reverse_pair[higher%in%path.gene])[1]
        d1<-dim(pathgene_reverse_pair[lower%in%path.gene])[1]
        c<-a-d1+c1
        d<-b-c1+d1
        if(a>=0 & b>=0 & c>=0 & d>=0 & a+b+c+d!=0){
          p_value<-tryCatch(chisq.test(matrix(c(a,b,c,d),nrow = 2),correct=T)$p.value,warning=function(w){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value})
          ifelse((a/b)<(c/d),return(c(p_value,"upregulation")),return(c(p_value,"downregulation")))
        }else{return(c("negative value in a b c d",a,b,c,d))}
      }else{return("no reverse pair in this path")}
    }else{
      return("no stable pair in this path")
    }
  }
  if(threads==1){
    aa<-lapply(term_gene_list, path.dyregulated.in.entirety,stable.pair=stable.pair,reverse.pair=reverse.pair)
  }else{
    c1<-makeCluster(threads)
    clusterExport(c1,"data.table",envir = environment())
    aa<-parLapply(c1,term_gene_list, path.dyregulated.in.entirety,stable.pair=stable.pair,reverse.pair=reverse.pair)
    stopCluster(c1)
  }
  nonrandom_path<-as.data.frame(do.call(rbind,aa),stringsAsFactors=F)
  colnames(nonrandom_path)<-c("p_value","direction")
  nonrandom_path<-cbind(as.integer(row.names(nonrandom_path)),nonrandom_path)
  index<-which(nonrandom_path[,3]=="no reverse pair in this path" | nonrandom_path[,3]=="no stable pair in this path")
  if(length(index)>0){
    nonrandom_path<-nonrandom_path[-index,]
  }
  nonrandom_path[,2]<-as.numeric(nonrandom_path[,2])
  nonrandom_path<-nonrandom_path[,1:2]
  nonrandom_path<-nonrandom_path[order(nonrandom_path[,2]),]
  nonrandom_path_fdr<-p.adjust(nonrandom_path[,2],method = "fdr")
  nonrandom_path<-cbind(nonrandom_path,nonrandom_path_fdr)
  colnames(nonrandom_path)<-c("termID","p_value","fdr")
  result<-nonrandom_path[nonrandom_path[,3]<threshold,]
  return(result)
}

