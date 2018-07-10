options(stringsAsFactors = FALSE)
if(!require("dplyr")) {install.packages("dplyr")}
library(dplyr)
if(!require("stringr")) {install.packages("stringr")}
library(stringr)
if(!require("igraph")) {install.packages("igraph")}
library(igraph)
if(!require("visNetwork")) {install.packages("visNetwork")}
library(visNetwork)



generate_cie10<-function(path="~/Documents/diseasomeCMBD2016/CIE10.txt"){
  cie10<-read.csv(path, sep=";", header=FALSE, encoding = "UTF-8")
  cie10$V1<-NULL
  cie10$V4<-NULL
  names(cie10)[1]<-"id"
  names(cie10)[2]<-"str"
  cie10$str<-enc2utf8(cie10$str)
  return(cie10)
}


generate_cmbd<-function(sex,age_min,age_max,path="~/Downloads/cmbd_madrid/CMBD_HOS_ANONIMO_20160101_20161231.csv"){
  cmbd <- read.csv(path, sep=";", encoding = "UTF-8")
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


build_igraph<-function(sex,age_min,age_max,cie10,cmbd_file="CMBD_HOS_ANONIMO_20160101_20161231.csv"){
  age_min<-max(age_min,0)
  age_max<-min(age_max,120)
  cmbd<-generate_cmbd(sex,age_min,age_max,cmbd_file)
  l<-generate_l(cmbd)
  v<-generate_v(l)
  e<-generate_e(l,v)
  edges<-build_edges(v,e)
  vertex<-build_nodes(edges,v,cie10)
  g<-graph.data.frame(edges,vertices=vertex)
  g.sym <- as.undirected(g, mode= "collapse")
  write_graph(g.sym,paste0(sex,formatC(age_min,width=3,format="d",flag="0"),formatC(age_max,width=3,format="d",flag="0"),".graphml"),"graphml")
  return(g.sym)
}

load_igraph<-function(sex,age_min,age_max) {
  return(read_graph(paste0(sex,formatC(age_min,width=3,format="d",flag="0"),formatC(age_max,width=3,format="d",flag="0"),".graphml"),"graphml"))
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

plot_commorbidity<-function(my_igraph,layout,size_by,min_size,max_size,physics,smooth) {
  if(size_by=="betweenness"){
    vector<-betweenness(my_igraph)
    V(my_igraph)$size<- min_size + ((max_size-min_size)*((vector-min(vector))/(max(vector)-min(vector))))
    # y = N + ((M-N)*((x-a)/(b-a)))
    # N = min_size
    # M = max_size
    # x = value
    # a = min(value)
    # b = max(value)
  } else {
    vector<-degree(my_igraph)
    V(my_igraph)$size<- min_size + ((max_size-min_size)*((vector-min(vector))/(max(vector)-min(vector))))
  }
  V(my_igraph)$size[V(my_igraph)$size<min_size]<-min_size
  V(my_igraph)$size[V(my_igraph)$size>max_size]<-max_size
  main<-"Red de comorbilidad (Autor: Julio Bonis)"
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

precalculate<-function(){
  cie10_file<-"CIE10.txt"
  cie10<-generate_cie10(cie10_file)
  
  for(sex in c(1,2)) {
    for(age_min in c(0,10,20,30,40,50,60,70,80)) {
      if(age_min<80){
        print(c(sex,age_min,age_min+9))
        build_igraph(sex,age_min,age_min+9,cie10)  
      } else {
        print(c(sex,age_min,age_min+9))
        build_igraph(sex,age_min,120,cie10)  
      }
    }
  }
}






setwd("C:/Users/jbonis_fcsai/Developer/trusty32/code/diseasomeCMBD2016")
cie10_file<-"CIE10.txt"
cie10<-generate_cie10(cie10_file)

for(sex in c(1,2)) {
  for(age_min in c(0,10,20,30,40,50,60,70,80)) {
    if(age_min<80){
      age_max<-age_min+9
    } else {
      age_max<-120
    }
    print(c(sex,age_min,120))
    var_name<-paste0(sex,formatC(age_min,width=3,format="d",flag="0"),formatC(age_max,width=3,format="d",flag="0"))
    mygraph<-read_graph(paste0(var_name,".graphml"),"graphml") 
    assign(paste0("igraph.",var_name),mygraph)
    assign(paste0("igraph.",var_name,".main"),decompose.graph(mygraph)[[which.max(components(mygraph)$csize)]])
  }
}



if(exists("summary.global.main")){remove("summary.global.main")}
for(sex in c(1,2)) {
  for(age_min in c(0,10,20,30,40,50,60,70,80)) {
    if(age_min<80){
      age_max<-age_min+9
    } else {
      age_max<-120
    }
    print(c(sex,age_min,age_max))
    var_name<-paste0(sex,formatC(age_min,width=3,format="d",flag="0"),formatC(age_max,width=3,format="d",flag="0"))
    mygraph<-get(paste0("igraph.",var_name))
    mygraph.main<-get(paste0("igraph.",var_name,".main"))
    assign(paste0("igraph.",var_name,".summary.edges"),build_summary_edges(sex,age_min,age_max,mygraph))
    assign(paste0("igraph.",var_name,".summary.global"),build_summary_global(sex,age_min,age_max,mygraph))
    
    
    if(exists("summary.global.main")){
      print(c("añado summary.global.main",sex,age_min,age_max))
      new_row<-build_summary_global(sex,age_min,age_max,mygraph.main)
      summary.global.main<-rbind(summary.global.main,new_row)  
    } else {
      print(c("Creo summary.global.main",sex,age_min,age_max))
      
      summary.global.main<-build_summary_global(sex,age_min,age_max,mygraph.main)
    }
    
    assign(paste0("igraph.",var_name,".main.summary.edges"),build_summary_edges(sex,age_min,age_max,mygraph.main))
    assign(paste0("igraph.",var_name,".main.summary.global"),build_summary_global(sex,age_min,age_max,mygraph.main))
  }
}










plot_commorbidity(igraph.1000009,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1010019,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1020029,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1030039,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1040049,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1050059,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1060069,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1070079,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.1080120,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2000009,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2010019,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2020029,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2030039,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2040049,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2050059,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2060069,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2070079,"layout_nicely","betweenness",10,100,T,T)
plot_commorbidity(igraph.2080120,"layout_nicely","betweenness",10,100,T,T)









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
  



