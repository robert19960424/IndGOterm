#' get all the term genes under the term "biological process", "molecular function" or "cellular composition"
#'
#' the function can get all the offspring pathway and their entire genes under "biological_process" , "molecular_function" or "cellular_component"
#'
#' @param relation a character string specifying the research direction , must be one of "biological_process","molecular_function" or "cellular_component"
#' @param term a matrix or data.frame contains the informations of the terms in the Gene Ontology database. Download from http://archive.geneontology.org/lite/2019-07-06/. By default the variables are taken from environment(formula).
#' @param gene2go_human a matrix or data.frame contains the informations of the geneID matching to the GOID. Download from https://ftp.ncbi.nih.gov/gene/DATA/ and should be preprocessed to get the matrix or data.frame of geneID matching to GOID of human beings. By default the variables are taken from environment(formula).
#' @param graph_path a matrix or data.frame contains the informations of relations of parent term and son term. Doownload from http://archive.geneontology.org/lite/2019-07-06/. By default the variables are taken from environment(formula).
#'
#' @return a matrix whose row is geneID and col is termID
#'
#' @import data.table
#' @examples
#' term_gene<-getterm.gene<-function(term,gene2go_human,graph_path,relation="biological_process")
#' @export
getterm.gene<-function(term,gene2go_human,graph_path,relation=c("biological_process","molecular_function","cellular_component")){
  if(missing(relation) | length(relation)>1 | sum(!(relation%in%c("biological_process","molecular_function","cellular_component")))>0){
    stop("relation must be one of the 'biological_process','molecular_function','cellular_component'")
  }
  if(missing(term)){
    warning("argument 'term' is missing, default is 'term'")
  }
  if(missing(gene2go_human)){
    warning("argument 'gene2go_human' is missing, default is 'gene2go_human'")
  }
  if(missing(graph_path)){
    warning("argument 'graph_path' is missing, default is 'graph_path'")
  }
  relationid=term[term[,2]==relation,1]
  term_relation=term[which(term[,3]==relation),]
  term_hum=intersect(term_relation[,4],gene2go_human[,2])
  gene2term=c()
  terminfo=c()
  for (i in term_hum){
    term1=gene2go_human[gene2go_human[,2]==i,]
    term2=term_relation[term_relation[,4]==i,]
    term1[,2]=term2[1,1]
    gene2term=rbind(gene2term,term1)
    terminfo=rbind(terminfo,term2)
  }
  gene2term=unique(gene2term)
  colnames(gene2term)<-c("geneID","termID")
  terminfo=terminfo[,-3]
  graph_path<-graph_path[graph_path[,3] %in% c(1,20,22),]
  reterm=terminfo[,1]
  index1=match(graph_path[,1],reterm)
  index2=match(graph_path[,2],reterm)
  re_graph=na.omit(cbind(graph_path,index1,index2))
  re_graph=re_graph[,1:4]
  re_graph=re_graph[re_graph[,4]!=0,]
  re_graph=unique(re_graph)
  rm(list=setdiff(ls(),c("gene2term","terminfo","re_graph","relationid")))
  gc()
  colnames(gene2term)<-c("GeneID","GO_ID")
  re_graph1<-data.table(re_graph)
  setnames(re_graph1,c("parent_ID","termID","relation","distance"))
  setkey(re_graph1,termID)
  regene2term<-data.table(gene2term)
  setnames(regene2term,c("geneID","termID"))
  setkey(regene2term,termID)
  for (i in sort(unique(re_graph1$distance))) {
    setkey(regene2term,termID)
    distance_re<-re_graph1[distance==i]
    term_gene<-regene2term[distance_re,allow.cartesian=T,on="termID"]
    setnames(term_gene,c("termID","parent_ID"),c("son_ID","termID"))
    regene2term<-unique(rbind(regene2term,term_gene[,c(1,3)]))
  }
  for (i in sort(unique(re_graph1$distance),decreasing = T)) {
    setkey(regene2term,termID)
    distance_re<-re_graph1[distance==i]
    term_gene<-regene2term[distance_re,allow.cartesian=T,on="termID"]
    setnames(term_gene,c("termID","parent_ID"),c("son_ID","termID"))
    regene2term<-unique(rbind(regene2term,term_gene[,c(1,3)]))
  }
  GOID_geneID<-list()
  goid<-unique(regene2term$termID)
  for (i in 1:length(goid)) {
    GOID_geneID[[i]]<-regene2term[termID==goid[i]]$geneID
    names(GOID_geneID)[i]<-goid[i]
  }
  goid_genenum<-c()
  for (i in 1:length(GOID_geneID)) {
    goid_genenum<-c(goid_genenum,length(GOID_geneID[[i]]))
  }
  z<-matrix(ncol =length(GOID_geneID),nrow = max(goid_genenum) )
  for (i in 1:ncol(z)) {
    z[1:length(GOID_geneID[[i]]),i]<-GOID_geneID[[i]]
  }
  colnames(z)<-names(GOID_geneID)
  return(z)
}
