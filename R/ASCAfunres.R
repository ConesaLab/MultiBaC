#' @include PCA-GENES.R
NULL

ASCAfun.res<-function (X,Fac) {


 PCA<-PCA.GENES(X)
      sc<-PCA$scores[,1:Fac]
      ld<-PCA$loadings[,1:Fac]
      ssq<-PCA$var.exp

      if(Fac==1) {
      sc<-as.matrix(sc)
      ld<-as.matrix(ld)
      }
      TPres<-sc%*%t(ld)

      if(Fac==0){
      sc=0
      ld=0
      TPres<-matrix(0,nrow(X),ncol(X))
      }

  Eres<-X-TPres

output<-list(sc,ld,ssq,X,TPres,Eres)
names(output)<-c("scores","loadings","var.exp","X","TP","E")
output


}
