ARSyNcomponents<-function(asca = asca, Variability = 0.75, beta=beta)
  {
    # This program selects the number of components that explain more than the Variability%
    # For residuals model the number of components selected are beta*average-variability.

    MODEL = asca[-length(asca)]
    M = length(MODEL)-1
    output<-NULL

    for (i in seq_len(M))
    {
      t<-which(MODEL[[i]]$var.exp[,2]>=Variability)[1]
      names(t)<-names(MODEL)[i]
      output<-c(output,t)
    }

    ### Residuals model
    i=M+1
    lim <- beta*1/Matrix::rankMatrix(MODEL[[i]]$X)[1]
    t<-table(MODEL[[i]]$var.exp[,1]>lim)
    if(length(t)==1) {t[2]=0}
    t<-t[2]
    names(t)<-names(MODEL)[i]
    output<-c(output,t)

    output
  }
