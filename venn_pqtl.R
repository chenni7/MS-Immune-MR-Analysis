library(VennDiagram)

setwd("C:/Users/Lenovo/Desktop/棕榈化代码/palmitoylation/1_get_palmitoylation_gene")

data1=read.table("pQTL.txt",header=F,sep = "\t")
dfset1=as.vector(data1[,1])
data2=read.table("palmitoylation.txt",header=F,sep = "\t")
dfset2=as.vector(data2[,1])


venn.plot <- venn.diagram(
  x = list(Set1 = dfset1, Set2 = dfset2),
  category.names = c("pQTL", "Palmitoylation"),
  fill = c("skyblue", "pink"),cat.pos=c(0,0),
  alpha = 0.5,cex = 1.5,
  cat.cex = 1.5,ext.dist = -0.1,
  ext.length = 0.85,
  filename = NULL
)
pdf("venn_pqtl.pdf",10,8);
grid.draw(venn.plot);
dev.off()
mylist= list(Set1 = dfset1, Set2 = dfset2)
vedata=Reduce(intersect,mylist)
write.table(vedata,"Venn_pqtl.txt",sep="\t",quote=F,col.names=F,row.names=F)
