options(stringsAsFactors = FALSE)
library(dplyr)
library(stringr)
library(igraph)
library(visNetwork)



generate_cie10<-function(path="~/Documents/diseasomeCMBD2016/CIE10.txt"){
  cie10<-read.csv(path, sep=";", header=FALSE)
  cie10$V1<-NULL
  cie10$V4<-NULL
  names(cie10)[1]<-"id"
  names(cie10)[2]<-"str"
  return(cie10)
}


generate_cmbd<-function(sex,age_min,age_max,path="~/Downloads/cmbd_madrid/CMBD_HOS_ANONIMO_20160101_20161231.csv"){
  cmbd <- read.csv(path, sep=";")
  cmbd<-cmbd[!duplicated(cmbd$HISTORIA_Anonimo),]
  cmbd$age<-as.numeric(format(as.Date(cmbd$FECING,"%d/%m/%Y"),"%Y"))-as.numeric(format(as.Date(cmbd$FECNAC,"%d/%m/%Y"),"%Y"))
  cmbd<-cmbd[cmbd$SEXO==sex&cmbd$age>=age_min&cmbd$age<=age_max,]
  cmbd$diag<-paste(cmbd$C1,
                   cmbd$C2,
                   cmbd$C3,
                   cmbd$C4,
                   cmbd$C5,
                   cmbd$C6,
                   cmbd$C7,
                   cmbd$C8,
                   cmbd$C9,
                   cmbd$C10,
                   cmbd$C11,
                   cmbd$C12,
                   cmbd$C13,
                   cmbd$C14,
                   cmbd$C15,
                   cmbd$C16,
                   cmbd$C17,
                   cmbd$C18,
                   cmbd$C19,
                   cmbd$C20,
                   sep="|")
  
  reduced.cmbd<-data.frame(
    sex=cmbd$SEXO,
    age=cmbd$age,
    diag=cmbd$diag
  )
  return(reduced.cmbd)
}

generate_l<-function(cmbd){
  l<-cmbd$diag %>%
    str_split("\\|") %>%
    lapply(function(x){x[!x==""]})
  return(l)
}

generate_v<-function(l){
  v<-l %>%
    unlist %>%
    table %>%
    data.frame %>%
    arrange(-Freq)
  v$id<-v$.
  return(v)
}

generate_e<-function(l) {
  e <- l %>%
    lapply(function(x) {
      expand.grid(x, x, weight = 1, stringsAsFactors = FALSE)
    }) %>%
    bind_rows
  
  e <- apply(e[, -3], 1, str_sort) %>%
    t %>%
    data.frame(stringsAsFactors = FALSE) %>%
    mutate(weight = e$weight)
  
  e <- group_by(e, X1, X2) %>%
    summarise(weight = sum(weight)/2) %>%
    filter(X1 != X2)
  e<-e %>% filter(weight>1)
  e<-mutate(e,wobs=(weight/length(l)))
  e<-mutate(e,wexpA=(as.numeric(v[v['.']==X1][2])/length(l)))
  e$.<-e$X2
  e<-merge(e,v,by=".")
  e<-mutate(e,wexpB=Freq/length(l))
  e<-mutate(e,wexp=wexpA*wexpB)
  e<-mutate(e,weight=wobs/wexp)
  e$.<-NULL
  e$wexpA<-NULL
  e$wexpB<-NULL
  e$Freq<-NULL
  e<-mutate(e,pvalue=apply(e[, c(4,6)], 1, function(row) prop.test(x=c(row[1]*length(l), row[2]*length(l)), n=c(length(l), length(l)))$p.value))
  return(e)
}

build_edges<-function(level,v,e){
  e <- e %>% filter(weight>1)
  e <- e %>% filter(pvalue<0.01)
  ef <- e %>% filter(weight>quantile(weight,level))
  ef<-mutate(ef,width=weight/mean(weight))
  names(ef)[1]<-"from"
  names(ef)[2]<-"to"
  ef<-ef[ef$from %in% v$id | ef$to %in% v$id,]
  ef$id<-NULL
  ef$label<-round(ef$weight,1)
  return(ef)
}

build_nodes<-function(edges,original_nodes,cie10){
  v<-original_nodes
  ef<-edges
  v<-merge(v,cie10,by="id")
  v2<-v[v$id %in% ef$from,]
  v3<-v[v$id %in% ef$to,]
  v4<- rbind(v3,v2[!v2$id %in% v3$id,])
  v4$label<-substr(v4$str,1,30)
  v4$title<-v4$str
  v4$group<-substr(v4$id,1,1)
  v4$.<-NULL
  return(v4)
}



cie10<-generate_cie10("~/Documents/diseasomeCMBD2016/CIE10.txt")
cmbd<-generate_cmbd(1,30,40,"~/Documents/diseasomeCMBD2016/CMBD_HOS_ANONIMO_20160101_20161231.csv")
l<-generate_l(cmbd)
v<-generate_v(l)
e<-generate_e(l)
edges<-build_edges(0,v,e)
vertex<-build_nodes(edges,v,cie10)
g<-graph.data.frame(edges,vertices=vertex)
g.sym <- as.undirected(g, mode= "collapse")


visNetwork(vertex,edges,main="Red de comorbilidad (Autor: Julio Bonis)") %>% 
  visLayout(improvedLayout=TRUE) %>%
  visEdges(shadow=T,smooth=F,dashes=F) %>%
  visPhysics(solver="barnesHut",
             barnesHut=list(
               gravitationalConstant=-10000,
               centralGravity=0.1,
               springLength=95,
               springConstant=0.01,
               damping=0.5,
               avoidOverlap=0.8
             ),
             stabilization=T) %>%
  visNodes(shadow=T) %>%
  visConfigure(enabled=T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = F),
            selectedBy="group",
            collapse=FALSE)
  

sort(degree(g.sym))
edge_density(g.sym,loops=F)
transitivity(g,type="global")
transitivity(g,type="local")
triad_census(g)
diameter(g,directed=F,weights=NA)
get_diameter(g,directed=F)
deg<-degree(g,mode="all")
hist(deg)
centr_degree(g,mode="all",normalized="T")
largest_cliques(g.sym)

