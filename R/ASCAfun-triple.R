#' @include ASCAfun12.R
NULL

ASCAfun.triple<-function (X,Desa,Desb,Desc,Fac) {

n <- ncol(X) # number of genes
I <- ncol(Desa) # number of levels in the factor TIME
J <- ncol(Desb) # number of levels in the other factor
H <- ncol(Desc) # number of levels in the other factor

#Matrices con medias efectos individuales
XK1<-matrix(NA,nrow=I,ncol=n)
for (i in seq_len(I)) {
     sub<-X[Desa[,i]==1,]
     XK1[i,]<-apply(sub,2,mean)
}


XK2<-matrix(NA,nrow=J,ncol=n)
for (j in seq_len(J)) {
     sub<-X[Desb[,j]==1,]
     XK2[j,]<-apply(sub,2,mean)
}

XK3<-matrix(NA,nrow=H,ncol=n)
for (h in seq_len(H)) {
     sub<-X[Desc[,h]==1,]
     XK3[h,]<-apply(sub,2,mean)
}

#Matrices con medias de efectos simples

XK12<-matrix(NA,nrow=I*J,ncol=n)
k=1
for (j in seq_len(J)){
  for (i in seq_len(I)){
    sub<-X[(Desa[,i]+Desb[,j])==2,]
    XK12[k,]<-apply(sub,2,mean)
    k=k+1
  }
}

XK13<-matrix(NA,nrow=I*H,ncol=n)
k=1
for (h in seq_len(H)){
  for (i in seq_len(I)){
    sub<-X[(Desa[,i]+Desc[,h])==2,]
    XK13[k,]<-apply(sub,2,mean)
    k=k+1
  }
}


XK23<-matrix(NA,nrow=J*H,ncol=n)
k=1
for (h in seq_len(H)){
  for (j in seq_len(J)){
    sub<-X[(Desb[,j]+Desc[,h])==2,]
    XK23[k,]<-apply(sub,2,mean)
    k=k+1
  }
}


NK<-matrix(NA,nrow=I,ncol=J*H)
XK<-matrix(NA,nrow=I*J*H,ncol=n)

k=1
for (h in seq_len(H)){
 for (j in seq_len(J)){
  for (i in seq_len(I)){
    sub<-as.matrix(rbind(X[(Desa[,i]+Desb[,j]+Desc[,h])==3,]))
    NK[i,(h-1)*J+j]<-sqrt(nrow(sub))
    XK[k,]<-apply(sub,2,mean)+XK1[i,]+XK2[j,]+XK3[h,]-XK12[(j-1)*I+i,]-XK13[(h-1)*I+i,]-XK23[(h-1)*J+j,]
    k=k+1
  }
 }
}

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

Xabc<-NULL
TPabc<-NULL
for (i in seq_len(nrow(X))){
     position1<-which(Desa[i,]==1)
     position2<-which(Desb[i,]==1)
     position3<-which(Desc[i,]==1)
     Xabc<-rbind(Xabc,XK[I*(position2-1)+I*J*(position3-1)+position1,])
     TPabc<-rbind(TPabc,XKrec[I*(position2-1)+I*J*(position3-1)+position1,])
}
Eabc<-Xabc-TPabc

 #Leverage & SPE
    leverage<-apply(ld^2,1,sum)
    SPE<-apply(Eabc^2,2,sum)

output<-list(XK,sc,ld,ssq,Xabc,TPabc,Eabc,leverage,SPE)
names(output)<-c("data","scores","loadings","var.exp","X","TP","E","leverage","SPE")
output


}
