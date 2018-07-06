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

generate_e<-function(l,v) {
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

build_edges<-function(v,e){
  e <- e %>% filter(weight>1)
  e <- e %>% filter(pvalue<0.01)
  e<-mutate(e,width=weight/mean(weight))
  names(e)[1]<-"from"
  names(e)[2]<-"to"
  e<-e[e$from %in% v$id | e$to %in% v$id,]
  e$id<-NULL
  e$label<-round(e$weight,1)
  return(e)
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


build_igraph<-function(sex,age_min,age_max,cie10,cmbd_file="~/Documents/diseasomeCMBD2016/CMBD_HOS_ANONIMO_20160101_20161231.csv"){
  cmbd<-generate_cmbd(sex,age_min,age_max,cmbd_file)
  l<-generate_l(cmbd)
  v<-generate_v(l)
  e<-generate_e(l,v)
  edges<-build_edges(v,e)
  vertex<-build_nodes(edges,v,cie10)
  g<-graph.data.frame(edges,vertices=vertex)
  g.sym <- as.undirected(g, mode= "collapse")
  return(g.sym)
}


build_summary_edges<-function(sex,age_min,age_max,g.sym){
  r.degree <- degree(g.sym)
  r.degree <- data.frame(id=names(r.degree),degree=r.degree)
  r.strength <- strength(g.sym)
  r.strength <- data.frame(id=names(r.strength),strength=r.strength)
  r.closeness <- closeness(g.sym,normalized=TRUE)
  r.closeness <- data.frame(id=names(r.closeness),closeness=r.closeness)
  r.betweenness <- betweenness(g.sym)
  r.betweenness <- data.frame(id=names(r.betweenness),betweenness=r.betweenness)
  r.summary <- merge(cie10,r.strength,by="id")
  r.summary <- merge(r.summary,r.degree,by="id")
  r.summary <- merge(r.summary,r.closeness,by="id")
  r.summary <- merge(r.summary,r.betweenness,by="id")
  return(r.summary)
}


build_summary_global<-function(sex,age_min,age_max,g.sym){
  return(data.frame(sex=sex,
                    age_min=age_min,
                    age_max=age_max,
                    num_vertex=length(V(g.sym)),
                    num_edges=length(E(g.sym)),
                    w_diameter=diameter(g.sym,directed=F),
                    diameter=diameter(g.sym,directed=F,weights=NA),
                    mean_distance=mean_distance(g.sym,directed=F),
                    edge_density=edge_density(g.sym),
                    transitivity=transitivity(g.sym)))
}

cie10_file<-"~/Documents/diseasomeCMBD2016/CIE10.txt"
cie10<-generate_cie10(cie10_file)

m0009.igraph<-build_igraph(1,0,9,cie10)
m1019.igraph<-build_igraph(1,10,19,cie10)
m2029.igraph<-build_igraph(1,20,29,cie10)
m3039.igraph<-build_igraph(1,30,39,cie10)
m4049.igraph<-build_igraph(1,40,49,cie10)
m5059.igraph<-build_igraph(1,50,59,cie10)
m6069.igraph<-build_igraph(1,60,69,cie10)
m7079.igraph<-build_igraph(1,70,79,cie10)
m8099.igraph<-build_igraph(1,80,120,cie10)

m0009.igraph.main<-decompose.graph(m0009.igraph)[[1]]
m1019.igraph.main<-decompose.graph(m1019.igraph)[[1]]
m2029.igraph.main<-decompose.graph(m2029.igraph)[[1]]
m3039.igraph.main<-decompose.graph(m3039.igraph)[[1]]
m4049.igraph.main<-decompose.graph(m4049.igraph)[[1]]
m5059.igraph.main<-decompose.graph(m5059.igraph)[[1]]
m6069.igraph.main<-decompose.graph(m6069.igraph)[[1]]
m7079.igraph.main<-decompose.graph(m7079.igraph)[[1]]
m8099.igraph.main<-decompose.graph(m8099.igraph)[[1]]


m0009.summary.edges<-build_summary_edges(1,0,9,m0009.igraph)
m1019.summary.edges<-build_summary_edges(1,10,19,m1019.igraph)
m2029.summary.edges<-build_summary_edges(1,20,29,m2029.igraph)
m3039.summary.edges<-build_summary_edges(1,30,39,m3039.igraph)
m4049.summary.edges<-build_summary_edges(1,40,49,m4049.igraph)
m5059.summary.edges<-build_summary_edges(1,50,59,m5059.igraph)
m6069.summary.edges<-build_summary_edges(1,60,69,m6069.igraph)
m7079.summary.edges<-build_summary_edges(1,70,79,m7079.igraph)
m8099.summary.edges<-build_summary_edges(1,80,89,m8099.igraph)

m7099.summary.edges<-merge(m7079.summary.edges,m8099.summary.edges,by=c("id","str"),suffixes = c(".m7079", ".m8099"))
m6099.summary.edges<-merge(m6069.summary.edges,m7099.summary.edges,by=c("id","str"))
m5099.summary.edges<-merge(m5059.summary.edges,m6099.summary.edges,by=c("id","str"),suffixes = c(".m5059", ".m6069"))


m0009.summary.global<-build_summary_global(1,0,9,m0009.igraph)
m1019.summary.global<-build_summary_global(1,10,19,m1019.igraph)
m2029.summary.global<-build_summary_global(1,20,29,m2029.igraph)
m3039.summary.global<-build_summary_global(1,30,39,m3039.igraph)
m4049.summary.global<-build_summary_global(1,40,49,m4049.igraph)
m5059.summary.global<-build_summary_global(1,50,59,m5059.igraph)
m6069.summary.global<-build_summary_global(1,60,69,m6069.igraph)
m7079.summary.global<-build_summary_global(1,70,79,m7079.igraph)
m8099.summary.global<-build_summary_global(1,80,120,m8099.igraph)

m0099.summary.global<-rbind(m0009.summary.global,
      m1019.summary.global,
      m2029.summary.global,
      m3039.summary.global,
      m4049.summary.global,
      m5059.summary.global,
      m6069.summary.global,
      m7079.summary.global,
      m8099.summary.global)

m0009.summary.main.global<-build_summary_global(1,0,9,m0009.igraph.main)
m1019.summary.main.global<-build_summary_global(1,10,19,m1019.igraph.main)
m2029.summary.main.global<-build_summary_global(1,20,29,m2029.igraph.main)
m3039.summary.main.global<-build_summary_global(1,30,39,m3039.igraph.main)
m4049.summary.main.global<-build_summary_global(1,40,49,m4049.igraph.main)
m5059.summary.main.global<-build_summary_global(1,50,59,m5059.igraph.main)
m6069.summary.main.global<-build_summary_global(1,60,69,m6069.igraph.main)
m7079.summary.main.global<-build_summary_global(1,70,79,m7079.igraph.main)
m8099.summary.main.global<-build_summary_global(1,80,120,m8099.igraph.main)

m0099.summary.main.global<-rbind(m0009.summary.main.global,
                            m1019.summary.main.global,
                            m2029.summary.main.global,
                            m3039.summary.main.global,
                            m4049.summary.main.global,
                            m5059.summary.main.global,
                            m6069.summary.main.global,
                            m7079.summary.main.global,
                            m8099.summary.main.global)

plot_commorbidity<-function(my_igraph,layout,size_by,min_size,max_size,physics,smooth) {
  if(size_by=="betweenness"){
    
    #V(my_igraph)$size<-betweenness(my_igraph)
    vector<-betweenness(my_igraph)
    V(my_igraph)$size<- min_size + ((max_size-min_size)*((vector-min(vector))/(max(vector)-min(vector))))
    # y = N + ((M-N)*((x-a)/(b-a)))
    # N = min_size
    # M = max_size
    # x = value
    # a = min(value)
    # b = max(value)
    print("betweenness")
  } else {
    print("no")
    #V(my_igraph)$size<-degree(my_igraph)
    vector<-degree(my_igraph)
    V(my_igraph)$size<- min_size + ((max_size-min_size)*((vector-min(vector))/(max(vector)-min(vector))))
  }
  V(my_igraph)$size[V(my_igraph)$size<min_size]<-min_size
  V(my_igraph)$size[V(my_igraph)$size>max_size]<-max_size
  visIgraph(my_igraph, layout=layout, physics=physics, smooth=smooth) %>%
    visPhysics(solver="barnesHut",
               barnesHut=list(
                 gravitationalConstant=-10000,
                 centralGravity=0.5,
                 springLength=95,
                 springConstant=0.01,
                 damping=0.05,
                 avoidOverlap=0.8
               ),
               stabilization=F)
}

plot_commorbidity(m5059.igraph.main,"layout_nicely","betweenness",10,100,T,T)



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
             stabilization=F) %>%
  visNodes(shadow=T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = F),
             selectedBy="group",
             collapse=FALSE)
  



