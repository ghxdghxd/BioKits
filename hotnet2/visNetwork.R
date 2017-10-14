#!/usr/bin/env Rscript
# @Date    : 2017-03-13 12:49:13
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

Args <- commandArgs(T)

if(length(Args)==0){
  print("visNetwork.R mat.txt heat.txt")
  q()
}

mat_file <- Args[1]
heat_file <- Args[2]

library(visNetwork)
mat <- read.table(mat_file ,header = T, stringsAsFactors = F, sep="\t")
heat <- read.table(heat_file ,header = T, stringsAsFactors = F, sep="\t")

mat$edges_color <- sapply(mat$networks,function(x){
  if(x=="hint+hi2012"){
    return("firebrick")
  }else if(x=="hprd"){
    return("green")
  }else if(x=="iref"){
    return("orange")
  }else if(x=="multinet"){
    return("blue")
  }
})

group <- as.data.frame.matrix(table(as.data.frame(unique(
  rbind(as.matrix(mat[,c(1,2)]), as.matrix(mat[,c(1,3)]))))))
group <- apply(group,2,function(x){sum(x)})
heat$size <- group[match(heat$gene,names(group))]
heat$color <- sapply(heat$type,function(x){
  if(x=="SNV"){
    return("steelblue")
  }else if(x=="SCNA"){
    return("mediumslateblue")
  }else if(x=="BOTH"){
    return("black")
  }
})
nodes <- data.frame(id = 1:nrow(heat),
                    label = heat$gene,
                    group = heat$type,
                    value = heat$size,
                    color = heat$color,
                    font.size = 100)

edges <- data.frame(
  from = as.numeric(sapply(mat$source,function(x){which(heat$gene==x)})),
  to = as.numeric(sapply(mat$target,function(x){which(heat$gene==x)})),
  color = mat$edges_color)

lnodes <- data.frame(label = c("SNV", "SCNA", "BOTH"),
  color = c("steelblue", "mediumslateblue","black"), shape="dot",
  font = list(size =  c("20", "20", "20"),
  align = c("horizontal","horizontal","horizontal")))
ledges <- data.frame(color = c("firebrick", "green", "orange", "blue"),
  label = c("hint+hi2012", "hprd", "iref", "multinet"),
  font = list(size =  c("20", "20", "20", "20"),
  align = c("bottom","bottom","bottom","bottom")),
  arrows =list(to = list(enabled = FALSE), to = list(enabled = FALSE),
  to = list(enabled = FALSE), to = list(enabled = FALSE)))

visNetwork(nodes, edges, height = "800px", width = "100%", main="gene interaction network") %>%
  visNodes(size=20) %>%
  addFontAwesome() %>%
  visEdges(length=100, width=2,smooth = T,
           font = list(align="bottom")) %>%
  visLegend(addNodes = lnodes, addEdges = ledges, width = 0.1,
            position = "right", useGroups=FALSE) %>%
  visIgraphLayout(smooth = TRUE, physics=FALSE, type="full") %>%
  visOptions(selectedBy = "group",
             highlightNearest = TRUE,
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE, solver = "barnesHut") %>%
  visSave(file = "network.html",selfcontained = FALSE)

