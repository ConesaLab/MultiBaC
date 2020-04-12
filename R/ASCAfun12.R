ASCAfun12<-function (X,Desa,Desb,Fac,Variability) {

n <- ncol(X) # number of genes
I <- ncol(Desa) # number of levels in the factor TIME
J <- ncol(Desb) # number of levels in the other factor


XK1<-matrix(NA,nrow=I,ncol=n)

for (i in seq_len(I)) {
     sub<-X[Desa[,i]==1,]
     XK1[i,]<-apply(sub,2,mean)
}



NK<-matrix(NA,nrow=I,ncol=J)
XK<-matrix(NA,nrow=I*J,ncol=n)

k=1
for (j in seq_len(J)){
  for (i in seq_len(I)){
    sub<-X[(Desa[,i]+Desb[,j])==2,]
    NK[i,j]<-sqrt(nrow(sub))
    XK[k,]<-apply(sub,2,mean)-XK1[i,]
    k=k+1
  }
}
NK<-sqrt(NK) #######
XKw<-XK*(as.numeric(NK))

PCA<-PCA.GENES(XKw)
      scw<-PCA$scores[,1:Fac]
      ld<-PCA$loadings[,1:Fac]
      ssq<-PCA$var.exp
      if(Fac==1) {
      scw<-as.matrix(scw)
      ld<-as.matrix(ld)
      }
# Re-weigth the scores
sc<-scw/(as.numeric(NK))

XKrec<-sc%*%t(ld)

Xab<-NULL
TPab<-NULL
for (i in seq_len(nrow(X))){
     position1<-which(Desa[i,]==1)
     position2<-which(Desb[i,]==1)
     Xab<-rbind(Xab,XK[I*(position2-1)+position1,])
     TPab<-rbind(TPab,XKrec[I*(position2-1)+position1,])
}
Eab<-Xab-TPab

output<-list(XK,sc,ld,ssq,Xab,TPab,Eab,Variability)
names(output)<-c("data","scores","loadings","var.exp","X","TP","E", "Variability")
output


}
