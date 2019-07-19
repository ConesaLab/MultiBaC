## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5)

## ---- include=TRUE-------------------------------------------------------
library(MultiBaC)
data(multiyeast)

## ---- include=TRUE-------------------------------------------------------
head(A.rna)

## ---- include=TRUE-------------------------------------------------------
input <- inputData (A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, 
                    batches = c(1,1,2,2,3,3),
                    omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"), 
                    batchesNames = c("A", "B", "C"))

## ---- include=TRUE-------------------------------------------------------
batchEstPlot(input, commonOmic = 1)

## ---- include = TRUE-----------------------------------------------------
# Create a global matrix
inData <- t(cbind(A.rna, B.rna, C.rna, A.gro, B.ribo, C.par))

# Generate PCA analysis
pc <- ropls::opls(inData, predI = 2, scaleC = "center", 
                  fig.pdfC = NULL, info.txtC = NULL)

# Plot PCA
{mypar <- par(no.readonly = TRUE)
  mai <- c(par()$mai[1:3],par()$pin[2]-0.2)
  par(mfrow= c(1,1), xpd = TRUE, mai = mai)}
  
plot(pc@scoreMN,
     xlab = paste0("PC 1: ",
                   pc@modelDF$`R2X(cum)`[1]*100, " %"),
     ylab = paste0("PC 2: ",
                   pc@modelDF$`R2X(cum)`[2]*100, " %"),
     # pch = omic; fill = condition
     pch = c(19,19,19,1,1,1,19,19,1,1,19,19,1,1, # RNA
             15,15,15,0,0,0, # GRO
             17,17,2,2, # RIBO
             18,18,5,5),# PAR
     # col = batch
     col = rep(c(rep("brown", 6),
                 rep("dodgerblue2", 4),
                 rep("forestgreen", 4)),2),
     # other arguments
     asp = 1, cex.axis = 1.3, cex.lab = 1.2, cex = 1.7, bty = "L",
     main = "Scoreplot of corrected data sets", cex.main = 1.7,
     xlim = c(-30, 40), ylim = c(-20,20), font.lab = 2)
abline(v = 0, lty = 5, col = "gray", xpd = FALSE)
abline(h = 0, lty = 5, col = "gray", xpd = FALSE)
legend(60,40, legend = c("Shape: Omic",
                             "RNA", "GRO", "RIBO", "PAR",
                             "Fill: Condition",
                             "Glu +", "Glu -",
                             "Color: Batch",
                             "Batch A", "Batch B", "Batch C"),
       col = c("white", rep(1, 4),
               "white", 1, 1,
               "white", "brown", "dodgerblue2", "forestgreen"),
       pch = c(0, 19, 15, 17, 18,
               0, 19, 1,
               0, 15, 15, 15), cex = 1.4, bty = "n")

## ---- include=TRUE-------------------------------------------------------
cond.factor = list("A" = c("Glu+", "Glu+", "Glu+", "Glu-", "Glu-", "Glu+"),
                   "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
                   "C" = c("Glu+", "Glu+", "Glu-", "Glu-"))

## ---- include=TRUE-------------------------------------------------------
res <- MultiBaC(input, test.comp = 5, cond.factor = cond.factor,
                showplot = TRUE)

## ---- include=TRUE-------------------------------------------------------
# Extract corrected matrices
A.rnacor <- res$A@ExperimentList$RNA
B.rnacor <- res$B@ExperimentList$RNA
C.rnacor <- res$C@ExperimentList$RNA
A.grocor <- res$A@ExperimentList$GRO
B.ribocor <- res$B@ExperimentList$RIBO
C.parcor <- res$C@ExperimentList$PAR

# Create a global matrix
outData <- t(cbind(A.rnacor, B.rnacor, C.rnacor, A.grocor, B.ribocor, C.parcor))

# Generate PCA analysis
pc <- ropls::opls(outData, predI = 2, scaleC = "center", plotL = FALSE, printL = FALSE)

# Plot PCA
{mypar <- par(no.readonly = TRUE)
  mai <- c(par()$mai[1:3],par()$pin[2]-0.2)
  par(mfrow= c(1,1), xpd = TRUE, mai = mai)}

plot(pc@scoreMN,
     xlab = paste0("PC 1: ",
                   pc@modelDF$`R2X(cum)`[1]*100, " %"),
     ylab = paste0("PC 2: ",
                   pc@modelDF$`R2X(cum)`[2]*100, " %"),
     # pch = omic; fill = condition
     pch = c(19,19,19,1,1,1,19,19,1,1,19,19,1,1, # RNA
             15,15,15,0,0,0, # GRO
             17,17,2,2, # RIBO
             18,18,5,5),# PAR
     # col = batch
     col = rep(c(rep("brown", 6),
                 rep("dodgerblue2", 4),
                 rep("forestgreen", 4)),2),
     # other arguments
     asp = 1, cex.axis = 1.3, cex.lab = 1.2, cex = 1.7, bty = "L",
     main = "Scoreplot of corrected data sets", cex.main = 1.7,
     xlim = c(-30, 40), ylim = c(-30,40))
abline(v = 0, lty = 5, col = "gray", xpd = FALSE)
abline(h = 0, lty = 5, col = "gray", xpd = FALSE)
legend(60,40, legend = c("Shape: Omic",
                             "RNA", "GRO", "RIBO", "PAR",
                             "Fill: Condition",
                             "Glu +", "Glu -",
                             "Color: Batch",
                             "Batch A", "Batch B", "Batch C"),
       col = c("white", rep(1, 4),
               "white", 1, 1,
               "white", "brown", "dodgerblue2", "forestgreen"),
       pch = c(0, 19, 15, 17, 18,
               0, 19, 1,
               0, 15, 15, 15), cex = 1.4, bty = "n")

## ---- include=TRUE-------------------------------------------------------
modelList <- genModelList(input, test.comp = 5)$modelList

## ---- include=TRUE-------------------------------------------------------
missingOmics <- genMissingOmics(input, modelList, commonOmic = 1)

## ---- include=TRUE-------------------------------------------------------
res <- batchCorrection(missingOmics, cond.factor)$correctedOmics
res$A

