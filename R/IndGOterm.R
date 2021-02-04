#' The main function of individPath
#'
#' The function is used to screen personalized dysregulation terms.
#' this function support parallel computation . You need to set thread numbers to make sure how many threads do you need to use.
#'
#' @param stable.pair a matrix or data.frame of stable pairs from normal samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param patients a matrix of data.frame of tumor samples , the first column must be the geneID ,and tumor samples start with the second column.
#' @param threads an integer value to make sure how many threads you will use to complete the computation
#' @param threshold a numeric value which is used to control false discovery rate under the p_value of the chip-square test or fisher exact probability test , default is 0.05.
#' @param threads an integer value to make sure how many threads you will use to complete the computation
#' @param term a matrix or data.frame contains the informations of the terms in the Gene Ontology database. Download from http://archive.geneontology.org/lite/2019-07-06/. By default the variables are taken from environment(formula).
#' @param term_gene_list a list of term genes. list type. each list is a vector type, containing all the genes of a term.By default the variables are "BPtermgene_list" and taken from environment(formula).
#' @param non_redun a vector of termID used to the screening of the personalized dysregulation terms.
#' @param gene.minimum an integer value to control the minimum amount of the genes of the terms. The deault is 1.
#'
#' @return a list contained two data.frame , individual dyregulated pathway with term ID and individual dyregulated pathway with term name
#' @import parallel
#' @import data.table
#' @examples
#' stable.pair<-t(combn(sample(1:10,10),2));
#' geneid<-1:10;
#' samples<-runif(100,min = 0,max = 50);
#' patients<-matrix(c(geneid,samples),nrow = 10,byrow=F);
#' reverse_pairs<-revpairs(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=1L,capacity=300000L)
#' global_dysterm<-get.globdysterm(stable.pair=stable.pair,reverse.pair=reverse_pairs,threads = 1,gene.minimum = 2,term_gene_list = BPtermgene_list,patients = patients,threshold = 0.05)
#' nonredundance<-remove.redundance(stable.pair = stable.pair,reverse.pair = reverse_pairs,gene.minimum = 2,patients = patients,re_graph = bp_graph,holidyregpath = global_dysterm,term_gene_list = BPtermgene_list,threads = 1)
#' samples2<-runif(50,min = 0,max = 50);
#' patients2<-matrix(c(geneid,samples2),nrow = 10,byrow=F);
#' result<-IndGOterm(stable.pair = stable.pair,term = term,patients = patients2,threads = 1,gene.minimum = 2,threshold = 0.05,term_gene_list = BPtermgene_list,non_redun = nonredundance$`non-redundant terms`)
#' @export
#'
IndGOterm<-function(stable.pair,term,patients,threads=1L,gene.minimum=1L,threshold=0.05,term_gene_list,non_redun){
  if(missing(stable.pair)){
    stop("argument 'stable.pair' of two columns is missing.")
  }
  if(missing(patients)){
    stop(" A dataset with disease samples to proceed the screening of individual dyregulated path is necessary")
  }
  if(missing(term_gene_list)){
    warning("argument 'term_gene_list' is missing, default is 'BPtermgene_list'")
    term_gene_list<-BPtermgene_list
  }
  if(class(term_gene_list)!='list'){
    warning("argument 'term_gene_list' isn't list class. The gene list of 'BPtermgene_list' is invoked.")
    term_gene_list<-BPtermgene_list
  }
  if(missing(non_redun)){
    stop("argument 'non_redun' is missing")
  }
  if(missing(term)){
    warning("argument 'term' is missing, default is 'term'")
    term=term
  }
  patients<-as.matrix(patients)
  stable.pair<-as.matrix(stable.pair[,1:2])
  colnames(stable.pair)<-c("higher","lower")
  for(i in 1:length(term_gene_list)){
    if(length(intersect(patients[,1],term_gene_list[[i]]))!=length(term_gene_list[[i]])){
      term_gene_list[[i]]<-intersect(patients[,1],term_gene_list[[i]])
    }
  }
  index<-unlist(lapply(term_gene_list,function(x) length(x)))
  index1<-which(index>=gene.minimum)
  term_gene_list<-term_gene_list[index1]
  dy_term_gene_list<-list()
  index<-match(non_redun,as.integer(names(term_gene_list)))
  dy_term_gene_list<-term_gene_list[index]
  dy_term<-data.frame()
  for (i in 1:length(dy_term_gene_list)) {
    dy_term[i,1]<-names(dy_term_gene_list)[i]
    dy_term[i,2]<-paste(dy_term_gene_list[[i]],collapse = "-")
  }
  colnames(dy_term)<-c("GOID","path_gene")
  if(sum(duplicated(dy_term[,2]))>0){
    dy_term<-data.table(dy_term)
    union_dy_term<-as.data.frame(dy_term[,.(newgoid=paste(GOID,collapse="-"),path_gene),by=path_gene][,2:3,with=T],stringsAsFactors=F)
    union_dy_term_gene_list<-list()
    for (i in 1:dim(union_dy_term)[1]) {
      union_dy_term_gene_list[[i]]<-as.integer(strsplit(union_dy_term[i,2],split="-")[[1]])
      names(union_dy_term_gene_list)[i]<-union_dy_term[i,1]
    }
  }else{
    union_dy_term_gene_list<-list()
    for (i in 1:dim(dy_term)[1]) {
      union_dy_term_gene_list[[i]]<-as.integer(strsplit(dy_term[i,2],split="-")[[1]])
      names(union_dy_term_gene_list)[i]<-dy_term[i,1]
    }
  }
  indGOterm<-function(stable.pair,patient,path.gene){
    library(data.table)
    stable.pair<-data.table(stable.pair)
    pathgene_all_stable_pair<-stable.pair[higher%in%path.gene | lower%in%path.gene][!(higher%in%path.gene & lower%in%path.gene)]
    a<-length(na.omit(match(pathgene_all_stable_pair$higher,path.gene)))
    b<-length(na.omit(match(pathgene_all_stable_pair$lower,path.gene)))
    patient<-data.table(patient)
    setnames(patient,c("higher","sample"))
    index<-pathgene_all_stable_pair[higher%in%path.gene]
    setkey(index,higher)
    setkey(patient,higher)
    index<-na.omit(index[patient])
    setnames(patient,c("lower","sample"))
    setkey(index,lower)
    setkey(patient,lower)
    index<-na.omit(index[patient])
    setnames(index,c("higher","lower","exp_h","exp_l"))
    index<-index[,comparison:=exp_h>exp_l]
    c1<-sum(index$comparison)
    d1<-sum(!index$comparison)


    setnames(patient,c("higher","sample"))
    index1<-pathgene_all_stable_pair[lower%in%path.gene]
    setkey(index1,higher)
    setkey(patient,higher)
    index1<-na.omit(index1[patient])
    setnames(patient,c("lower","sample"))
    setkey(index1,lower)
    setkey(patient,lower)
    index1<-na.omit(index1[patient])
    setnames(index1,c("higher","lower","exp_h","exp_l"))
    index1<-index1[,comparison:=exp_h>exp_l]
    c2<-sum(!index1$comparison)
    d2<-sum(index1$comparison)
    c<-c1+c2
    d<-d1+d2
    if(a+b+c+d==0){
      stop(" a+b+c+d=0")
    }
    if((a/b)>(c/d)){
      p_value<-tryCatch(chisq.test(matrix(c(a,b,c,d),nrow = 2),correct=T)$p.value,warning=function(w){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value})
      direction<-"downregulation"
    }else if((a/b)<(c/d)){
      p_value<-tryCatch(chisq.test(matrix(c(a,b,c,d),nrow = 2),correct=T)$p.value,warning=function(w){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value})
      direction<-"upregulation"
    }else if((a/b)==(c/d)){
      p_value<-tryCatch(chisq.test(matrix(c(a,b,c,d),nrow = 2),correct=T)$p.value,warning=function(w){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value})
      direction<-"nonregulation"
    }
    return(matrix(c(p_value,direction),ncol = 2))
  }
  patients_individpvalue<-list()
  rm(list=setdiff(ls(),c("stable.pair","patients","threads","threshold",'union_dy_term_gene_list',"patients_individpvalue","indGOterm")))
  if(threads==1){
    for (i in 1:(dim(patients)[2]-1)) {
      result<-lapply(union_dy_term_gene_list,indGOterm,stable.pair=stable.pair,patient=patients[,c(1,i+1)])
      patients_individpvalue[[i]]<-do.call(rbind,result)
    }
  }else{
    c1<-makeCluster(threads)
    for (i in 1:(dim(patients)[2]-1)) {
      result<-parLapply(c1,union_dy_term_gene_list,indGOterm,stable.pair=stable.pair,patient=patients[,c(1,i+1)])
      patients_individpvalue[[i]]<-do.call(rbind,result)
    }
    stopCluster(c1)
  }
  result<-do.call(cbind,patients_individpvalue)
  colnames(result)<-rep(colnames(patients)[2:dim(patients)[2]],each=2)
  row.names(result)<-names(union_dy_term_gene_list)
  index1<-grep("-",names(union_dy_term_gene_list))
  if(length(index1)>0){
    aa<-names(union_dy_term_gene_list)[index1]
    aa<-strsplit(aa,split = "-")
    aaa<-list()
    m=1
    n=1
    for (i in 1:length(aa)) {
      for (j in 1:length(aa[[i]])) {
        aaa[[n]]<-union_dy_term_gene_list[[index1[m]]]
        names(aaa)[n]<-aa[[i]][j]
        n=n+1
      }
      m=m+1
    }
    union_dy_term_gene_list1<-union_dy_term_gene_list[-index1]
    union_dy_term_gene_list1<-c(union_dy_term_gene_list1,aaa)
    index1<-grep("-",row.names(result))
    index1<-row.names(result)[index1]
    if(length(index1)>0){
      for (i in 1:length(aa)) {
        m<-matrix(rep(result[match(index1[i],row.names(result)),],times=length(aa[[i]])),nrow = length(aa[[i]]),byrow = T)
        row.names(m)<-aa[[i]]
        result<-result[-match(index1[i],row.names(result)),,drop=F]
        result<-rbind(result,m)
      }
      a<-matrix(0,nrow = length(union_dy_term_gene_list1))
      for (i in 1:(dim(result)[2]/2)) {
        a<-cbind(a,matrix(as.numeric(result[,2*i-1]),ncol = 1))
      }
      a<-as.data.frame(a[,-1,drop=F],stringsAsFactors=F)
      row.names(a)<-row.names(result)
      allsample_fdr<-matrix(0,nrow = length(union_dy_term_gene_list1))
      for (i in 1:(dim(result)[2]/2)) {
        index<-sort(a[,i],decreasing = F)
        index1<-rank(a[,i])
        fdr<-p.adjust(index,method = "fdr")
        fdr1<-fdr[index1]
        allsample_fdr<-cbind(allsample_fdr,fdr1)
      }
      allsample_fdr<-as.data.frame(allsample_fdr[,-1,drop=F],stringsAsFactors=F)
      colnames(allsample_fdr)=colnames(patients)[2:dim(patients)[2]]
      colnames(a)<-colnames(patients)[2:dim(patients)[2]]
      row.names(allsample_fdr)<-row.names(a)
      dyregulation<-allsample_fdr<threshold
      if( any(dim(dyregulation)==0)){
        dyregulation_result<-"no dysregulation terms in any patients"
      }else{
        for (i in 1:(dim(result)[2]/2)) {
          downregulation<-which(result[,2*i]=="downregulation")
          dyregulation[downregulation,i]<-dyregulation[downregulation,i]*(-1)
        }
        for (i in 1:(dim(result)[2]/2)) {
          downregulation<-which(result[,2*i]=="downregulation")
          allsample_fdr[downregulation,i]<-allsample_fdr[downregulation,i]*(-1)
        }
        dyregulation_result<-list()
        dyregulation_result[[1]]<-dyregulation
        index<-match(as.integer(row.names(dyregulation)),term[,1])
        row.names(dyregulation)<-term[index,2]
        dyregulation_result[[2]]<-dyregulation
        dyregulation_result[[3]]<-allsample_fdr
        names(dyregulation_result)<-c('result with term ID',"result with term name","result with FDR")
      }
    }
  }else{
    result<-as.data.frame(result,stringsAsFactors=F)
    result[,1]<-as.numeric(result[,1])
    a<-matrix(0,nrow = nrow(result))
    for (i in 1:(dim(result)[2]/2)) {
      a<-cbind(a,matrix(as.numeric(result[,2*i-1]),ncol = 1))
    }
    a<-as.data.frame(a[,-1],stringsAsFactors=F)
    row.names(a)<-row.names(result)
    colnames(a)<-colnames(patients)[2:dim(patients)[2]]
    allsample_fdr<-matrix(0,nrow = nrow(result))
    for (i in 1:(dim(result)[2]/2)) {
      index<-sort(a[,i],decreasing = F)
      index1<-rank(a[,i])
      fdr<-p.adjust(index,method = "fdr")
      fdr1<-fdr[index1]
      allsample_fdr<-cbind(allsample_fdr,fdr1)
    }
    allsample_fdr<-as.data.frame(allsample_fdr[,-1],stringsAsFactors=F)
    colnames(allsample_fdr)=colnames(patients)[2:dim(patients)[2]]
    colnames(a)<-colnames(patients)[2:dim(patients)[2]]
    row.names(allsample_fdr)<-row.names(a)
    dyregulation<-allsample_fdr<threshold
    if( any(dim(dyregulation)==0)){
      dyregulation_result<-"no dysregulation terms in any patients"
    }else{
      for (i in 1:(dim(result)[2]/2)) {
        downregulation<-which(result[,2*i]=="downregulation")
        dyregulation[downregulation,i]<-dyregulation[downregulation,i]*(-1)
      }
      for (i in 1:(dim(result)[2]/2)) {
        downregulation<-which(result[,2*i]=="downregulation")
        allsample_fdr[downregulation,i]<-allsample_fdr[downregulation,i]*(-1)
      }
      dyregulation_result<-list()
      dyregulation_result[[1]]<-dyregulation
      index<-match(as.integer(row.names(dyregulation)),term[,1])
      row.names(dyregulation)<-term[index,2]
      dyregulation_result[[2]]<-dyregulation
      dyregulation_result[[3]]<-allsample_fdr
      names(dyregulation_result)<-c('result with term ID',"result with term name","result with FDR")
    }
  }
  return(dyregulation_result)
}
