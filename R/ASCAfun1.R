NULL

ASCAfun1<-function (X,Design,Fac,Variability) {

n <- ncol(X) # number of genes
I <- ncol(Design) # number of levels in the factor

NK<-NULL
XK<-matrix(NA,nrow=I,ncol=n)

for (i in seq_len(I)) {
     sub<-X[Design[,i]==1,]
     NK[i]<-nrow(sub)
     XK[i,]<-apply(sub,2,mean)
}
  NK<-sqrt(NK)

# Weigh the data of the Submodel with the corresponding number of measurement occasions

  XKw<- NK*XK

  PCA<-PCA.GENES(XKw)
      scw<-PCA$scores[,1:Fac]
      ld<-PCA$loadings[,1:Fac]
      ssq<-PCA$var.exp
      if(Fac==1) {
      scw<-as.matrix(scw)
      ld<-as.matrix(ld)
      }

# Re-weigth the scores
    sc<-scw/NK

XKrec<-sc%*%t(ld)

Xa<-NULL
TPa<-NULL
for (i in seq_len(nrow(X))){
    position<-which(Design[i,]==1)
    Xa<-rbind(Xa,XK[position,])
    TPa<-rbind(TPa,XKrec[position,])
}

Ea<-Xa-TPa

output<-list(XK,sc,ld,ssq,Xa,TPa,Ea,Variability)
names(output)<-c("data","scores","loadings","var.exp","X","TP","E","Variability")
output

}

