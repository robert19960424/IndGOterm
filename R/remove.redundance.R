#' remove redundant pathway
#'
#' this function is used to remove redundant parent or son terms from all of global dysregulation terms
#' this function support parallel computation . You need to set thread numbers to make sure how many threads do you need to use.
#'
#'
#' @param stable.pair a matrix or data.frame of stable pairs from normal samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param patients a matrix of data.frame of tumor samples , the first column must be the geneID ,and tumor samples start with the second column.
#' @param threads an integer value to make sure how many threads would be used to complete the computation.
#' @param reverse.pair a matrix or data.frame of reverse pairs from disease samples with two columns, the expression of the genes in the first column is higher than that in the second column.
#' @param re_graph a matrix or data.frame with three columns about the relationship between parent term and whose son term. The first to third column is parent termID, son termID and distance between parent term and son term separately.
#' @param holidyregpath a matrix or data.frame with three columns depicts the coorelation levels of terms and disease. Which computed by method get.holidyregpath. The first to third column is termID, p-value and FDR value.
#' @param term_gene_list a list of term genes. list type. each list is a vector type, containing all the genes of a term.By default the variables are "BPtermgene_list" and taken from environment(formula).
#' @param gene.minimum an integer value to control the minimum amount of the genes of the terms. The deault is 1.
#'
#' @return a list.The first list is all terms and their redundant parent terms. of all the nonredundant holistic-dyregulated pathway
#' @examples
#' stable.pair<-t(combn(sample(1:10,10),2));
#' geneid<-1:10;
#' samples<-runif(100,min = 0,max = 50);
#' patients<-matrix(c(geneid,samples),nrow = 10,byrow=F);
#' reverse_pairs<-revpairs(stable.pair,patients,threshold=0.05,spairs_threshold=0.99,threads=1L,capacity=300000L)
#' global_dysterm<-get.globdysterm(stable.pair=stable.pair,reverse.pair=reverse_pairs,threads = 1,gene.minimum = 2,term_gene_list = BPtermgene_list,patients = patients,threshold = 0.05)
#' nonredundance<-remove.redundance(stable.pair = stable.pair,reverse.pair = reverse_pairs,gene.minimum = 2,patients = patients,re_graph = bp_graph,holidyregpath = global_dysterm,term_gene_list = BPtermgene_list,threads = 1)
#' @export
#'
remove.redundance<-function(stable.pair,reverse.pair,gene.minimum=1,patients,re_graph,holidyregpath,term_gene_list,threads=1L){
  if(missing(stable.pair)){
    stop("argument 'stable.pair' of two columns is missing.")
  }
  if(missing(reverse.pair)){
    stop("argument 'reverse.pair' of two columns is missing.")
  }
  if(missing(patients)){
    stop("argument 'patients' is missing. This argument should be completely same with the argument 'patients' used in function 'revpairs'.")
  }
  if(missing(re_graph)){
    stop("The argument 're_graph' is missing, which should be one of 'bp_graph','cc_graph','mf_graph'.")
  }
  if(missing(holidyregpath)){
    stop("argument 'holidyregpath' is missing, with no default")
  }
  if(missing(term_gene_list)){
    warning("argument 'term_gene_list' is missing, default is 'BPtermgene_list'")
    term_gene_list<-BPtermgene_list
  }
  if(class(term_gene_list)!='list'){
    warning("argument 'term_gene_list' isn't list class. The gene list of 'BPtermgene_list' is invoked.")
    term_gene_list<-BPtermgene_list
  }
  if(missing(threads)){
    threads=1L
  }
  patients<-as.matrix(patients)
  stable.pair<-as.matrix(stable.pair[,1:2])
  reverse.pair<-as.matrix(reverse.pair[,1:2])
  colnames(stable.pair)<-c('higher','lower')
  for(i in 1:length(term_gene_list)){
    if(length(intersect(patients[,1],term_gene_list[[i]]))!=length(term_gene_list[[i]])){
      term_gene_list[[i]]<-intersect(patients[,1],term_gene_list[[i]])
    }
  }
  index<-unlist(lapply(term_gene_list,function(x) length(x)))
  index1<-which(index>=gene.minimum)
  term_gene_list<-term_gene_list[index1]
  term_gene<-list()
  for (i in 1:length(term_gene_list)) {
    term_gene[[i]]<-matrix(c(term_gene_list[[i]],rep(as.integer(names(term_gene_list)[i]),times=length(term_gene_list[[i]]))),ncol = 2)
  }
  term_gene<-do.call(rbind,term_gene)
  colnames(term_gene)<-c("geneID","termID")
  local.redundant<-function(son.termID,re_graph,term_gene,stable.pair,reverse.pair,holidyregpath){
    library(data.table)
    re_graph<-unique(data.table(re_graph[,1:2]))
    setnames(re_graph,c("termID","son_ID"))
    parent_termID<-re_graph[son_ID==son.termID]$termID
    parent_termID<-intersect(parent_termID,holidyregpath[,1])
    term_gene<-data.table(term_gene)
    setnames(term_gene,c("geneID","termID"))
    son_gene<-term_gene[termID==son.termID]$geneID
    if(length(son_gene)==0){
      stop(paste("there are no genes in son term : ",son.termID,sep = ""))
    }
    stable.pair<-data.table(stable.pair)
    setnames(stable.pair,c("higher","lower"))
    reverse.pair<-data.table(reverse.pair)
    setnames(reverse.pair,c("higher","lower"))
    redundant<-c()
    for (i in parent_termID) {
      parent_gene<-term_gene[termID==i]$geneID
      parent_self_gene<-setdiff(parent_gene,son_gene)
      if((length(parent_gene)-length(parent_self_gene))!=length(son_gene)){
        stop("parent term:",i," doesn't contain all his son term's gene","\n And the son term is ",son.termID,sep="")
      }
      parent_self_gene_stable_pair<-stable.pair[higher%in%parent_self_gene | lower%in%parent_self_gene][!(higher%in%parent_gene & lower%in%parent_gene)]
      stable_pair_parentself_num<-nrow(parent_self_gene_stable_pair)
      parent_self_gene_reverse_pair<-reverse.pair[higher%in%parent_self_gene | lower%in%parent_self_gene][!(higher%in%parent_gene & lower%in%parent_gene)]
      reverse_pair_parentself_num<-nrow(parent_self_gene_reverse_pair)
      a<-length(na.omit(match(parent_self_gene_stable_pair$higher,parent_self_gene)))
      b<-length(na.omit(match(parent_self_gene_stable_pair$lower,parent_self_gene)))
      c1<-dim(parent_self_gene_reverse_pair[higher%in%parent_self_gene])[1]
      d1<-dim(parent_self_gene_reverse_pair[lower%in%parent_self_gene])[1]
      c<-a-d1+c1
      d<-b-c1+d1
      p_value<-tryCatch(chisq.test(matrix(c(a,b,c,d),nrow = 2),correct=T)$p.value,warning=function(w){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value},error=function(e){fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value})
      if(p_value>0.05){
        stable_pair_parentall_num<-dim(stable.pair[higher%in%parent_gene | lower%in%parent_gene][!(higher%in%parent_gene & lower%in%parent_gene)])[1]
        reverse_pair_parentall_num<-dim(reverse.pair[higher%in%parent_gene | lower%in%parent_gene][!(higher%in%parent_gene & lower%in%parent_gene)])[1]
        if(stable_pair_parentself_num==0){
          redundant<-c(redundant,i)
        }else if((reverse_pair_parentself_num/stable_pair_parentself_num)<(reverse_pair_parentall_num/stable_pair_parentall_num)){
          redundant<-c(redundant,i)
        }
      }
      if(!(son.termID%in%redundant)){
        if(p_value<0.05){
          fdr_sonterm<-holidyregpath[holidyregpath[,1]==son.termID,3]
          if(length(fdr_sonterm)==0){
            stop(paste("there isn't a false dicovery rate which can match the son term : ",son.termID,"\n",sep=""))
          }
          fdr_parentterm<-holidyregpath[holidyregpath[,1]==i,3]
          if(length(fdr_parentterm)==0){
            stop(paste("there isn't a false dicovery rate which can match the parent term : ",parent_termID,"\n","And the son term ID is ",son.termID,sep=""))
          }
          if(fdr_parentterm<fdr_sonterm){
            redundant<-c(redundant,son.termID)
          }else if(fdr_parentterm>=fdr_sonterm){
            parent_gene_stable_pair<-stable.pair[higher%in%parent_gene | lower%in%parent_gene][!(higher%in%parent_gene & lower%in%parent_gene)]
            parent_gene_reverse_pair<-reverse.pair[higher%in%parent_gene | lower%in%parent_gene][!(higher%in%parent_gene & lower%in%parent_gene)]
            N<-nrow(parent_gene_stable_pair)
            M<-nrow(parent_gene_reverse_pair)
            n<-nrow(parent_gene_stable_pair[higher%in%son_gene | lower%in%son_gene][!(higher%in%son_gene & lower%in%son_gene)])
            m<-nrow(parent_gene_reverse_pair[higher%in%son_gene | lower%in%son_gene][!(higher%in%son_gene & lower%in%son_gene)])
            p_value_son_phyper<-phyper(m-1,M,N-M,n,lower.tail = F)
            if(p_value_son_phyper>0.05){
              redundant<-c(redundant,son.termID)
            }
          }
        }
      }
    }
    return(redundant)
  }
  termbp<-holidyregpath[,1]
  rm(list = setdiff(ls(),c("termbp","local.redundant","re_graph","term_gene","stable.pair","reverse.pair","holidyregpath","threads")))
  if(threads==1){
    redundant_path<-sapply(termbp, local.redundant,re_graph=re_graph,term_gene=term_gene,stable.pair=stable.pair,reverse.pair=reverse.pair,holidyregpath=holidyregpath)
  }else{
    c1<-makeCluster(threads)
    redundant_path<-parSapply(c1,termbp, local.redundant,re_graph=re_graph,term_gene=term_gene,stable.pair=stable.pair,reverse.pair=reverse.pair,holidyregpath=holidyregpath)
    stopCluster(c1)
  }
  all_redun<-unique(unlist(redundant_path))
  non_redun<-setdiff(termbp,all_redun)
  names(redundant_path)<-holidyregpath[,1]
  result<-list(redundant_path,non_redun)
  names(result)<-c("redundance terms of each term","non-redundant terms")
  return(result)
}



