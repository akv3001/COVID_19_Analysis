Create_Design <- function(input,condition,comparison){
  library("limma")
  # print(condition)
  input<- as.data.frame(input)
  colnames(input)<-c("Samples","target")
  input$target <- as.matrix(input$target)
  design <<- model.matrix(~0+ unlist(input$target))
  colnames(design) <<- condition
  contmatrix <<- makeContrasts(as.character(comparison),levels=design)
}


CallLimma <- function(Data,design,contmatrix){ 
  suppressMessages(library(limma))
  suppressMessages(library(edgeR))
  Data = Data[which(rowSums(Data) >0 ),]
  isexpr<-apply(Data,1,function(x){length(which(cpm(x)>1)) >= 2 })
  Data.norm <-  Data[isexpr,]
  Data.voom<- voom(Data.norm,plot=FALSE)
  fit<- lmFit(Data.voom,design)
  fit2<<- contrasts.fit(fit,contmatrix)
  fit2<<-eBayes(fit2)
  return(toptable(fit2,n=Inf))
}
