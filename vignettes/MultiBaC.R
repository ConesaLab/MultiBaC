## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 12, fig.height = 7)
devtools::load_all(".")
library("MultiAssayExperiment")
library("ggplot2")
library("ropls")
library("MultiBaC")

## ----init, eval = TRUE, echo = FALSE, fig.cap = "PCA plot of original gene expression data (before correction). Batches are completely separated from each other. Plot generated with MultiBaC package (see Visualization of results Section)."----
data_RNA <- createMbac (inputOmics = list(A.rna, B.rna, C.rna), 
                    batchFactor = c("A", "B", "C"),
                    experimentalDesign = list("A" =  c("Glu+", "Glu+", 
                                                "Glu+", "Glu-", 
                                               "Glu-", "Glu-"),
                                       "B" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-"),
                                       "C" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-")),
                    omicNames = "RNA")
custom_col <- c("brown2", "dodgerblue", "forestgreen")
custom_pch <- c(19,19,19,1,1,1,
                19,19,1,1,
                19,19,1,1)
plot(data_RNA, type="pca.org", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", levels(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 3),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 2))

## --------------------------------------------------------------------------
data("multiyeast")
head(A.rna)

## --------------------------------------------------------------------------
data_RNA<- createMbac (inputOmics = list(A.rna, B.rna, C.rna), 
                    batchFactor = c("A", "B", "C"),
                    experimentalDesign = list("A" =  c("Glu+", "Glu+", 
                                                "Glu+", "Glu-", 
                                               "Glu-", "Glu-"),
                                       "B" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-"),
                                       "C" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-")),
                    omicNames = "RNA")

## ---- eval = FALSE---------------------------------------------------------
#  ARSyNbac (mbac, batchEstimation = TRUE,
#            Interaction=FALSE, Variability = 0.90, beta = 2,
#            modelName = "Model 1",
#            showplot = TRUE)

## ----arsyn1, eval = TRUE, message = FALSE, fig.cap = "Batch correction when Interaction=FALSE. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the number of components (x axis) needed to explain a certain variability (y axis) of the effect of interest (batch effect). The number of components selected in the model is indicated with a triangle symbol. Gray dashed line indicates the threshold given by the Variability argument (in percentage). Right: PCA plot of corrected gene expression data when not considering the interaction batch-condition."----
par(mfrow = c(1,2))
arsyn_1 <- ARSyNbac(data_RNA, modelName = "RNA", Variability = 0.95, 
                 batchEstimation = TRUE, Interaction = FALSE)
plot(arsyn_1, type="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", levels(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 3),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))

## ----arsyn2, message = FALSE, fig.cap = "Batch correction when Interaction = TRUE. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the number of components (x axis) needed to explain a certain variability (y axis) of the effect of interest (batch effect). The number of components selected in the model is indicated with a triangle symbol. Gray dashed line indicates the threshold given by the Variability argument (in percentage). Right: PCA plot of corrected gene expression data considering the interaction batch-condition."----
par(mfrow = c(1,2))
arsyn_2 <- ARSyNbac(data_RNA, modelName = "RNA", Variability = 0.95, 
                 batchEstimation = TRUE, Interaction = TRUE)
plot(arsyn_2, type="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", levels(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 3),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))

## ----arsyn3, message = FALSE, fig.cap = "Noise reduction mode. Left: Explained variance plot. Default plot when showplot = TRUE. It represents the percentage of variability in the residuals (y axis) explained by a model with a given number of principal components (x axis). The number of selected components in the final model is indicated with a triangle symbol, and computed to explain beta times the average variability of the residuals. Right: PCA plot of corrected gene expression data."----
par(mfrow = c(1,2))
arsyn_3 <- ARSyNbac(data_RNA, modelName = "RNA", beta = 0.5, 
                 batchEstimation = FALSE)
plot(arsyn_3, type="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_RNA$ListOfBatches),
                     "Fill: Cond.", levels(colData(data_RNA$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 3),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))

## ----design, echo = FALSE, fig.cap = "Scheme of the yeast example data structure."----
knitr::include_graphics("designScheme.png", dpi = 30)

## --------------------------------------------------------------------------
my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, 
                                               B.rna, B.ribo, 
                                               C.rna, C.par), 
                    batchFactor = c("A", "A", 
                                    "B", "B", 
                                    "C", "C"),
                    experimentalDesign = list("A" =  c("Glu+", "Glu+", 
                                                "Glu+", "Glu-", 
                                               "Glu-", "Glu-"),
                                       "B" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-"),
                                       "C" = c("Glu+", "Glu+", 
                                               "Glu-", "Glu-")),
                    omicNames = c("RNA", "GRO", 
                                  "RNA", "RIBO", 
                                  "RNA", "PAR"))

## ---- eval = FALSE, echo = TRUE--------------------------------------------
#  MultiBaC (mbac,
#            test.comp = NULL, scale = FALSE,
#            center = TRUE, crossval = NULL,
#            Interaction = FALSE,
#            Variability = 0.90,
#            showplot = TRUE,
#            showinfo = TRUE)

## ----main-multibac2, include = TRUE, echo = TRUE, fig.cap = "Q2 and explained variance plots. Q2 plot (left) shows the number ob components (x) needed to reach a certain Q2 value (y). The number of components selected for each model is indicated with a triangle symbol. Gray dashed line indicates the 0.7 Q2 threshold. Explained variance plot (right) represents the number of components (x) needed to explain a certain varibility (y) of the effect of interest (batch effect). The number of components selected for each model is indicated with a triangle symbol. Gray dashed line indicates the Variability argument in percentage."----
my_final_mbac <- MultiBaC (my_mbac,
          				test.comp = NULL, scale = FALSE, 
          				center = TRUE, crossval = NULL, 
          				Interaction = TRUE,
          				Variability = 0.9,
          				showplot = TRUE,
          				showinfo = TRUE)

## --------------------------------------------------------------------------
my_mbac_2 <- genModelList (my_mbac, test.comp = NULL, 
              			  scale = FALSE, center = TRUE,
              			  crossval = NULL,
              			  showinfo = TRUE)

## --------------------------------------------------------------------------
multiBatchDesign <- genMissingOmics(my_mbac_2)

## --------------------------------------------------------------------------
my_finalwise_mbac <- batchCorrection(my_mbac_2, 
                                     multiBatchDesign = multiBatchDesign,
                                     Interaction = TRUE,
                                     Variability = 0.90)

## ---- eval = FALSE, echo = TRUE--------------------------------------------
#  plot (mbac, type = "def",
#        col.by.batch = TRUE,
#        col.per.group = NULL,
#        comp2plot = c(1,2),
#        legend.text = NULL,
#        args.legend = NULL, ...)

## ---- eval = TRUE, echo = FALSE--------------------------------------------
# just one plot in next chunk
my_aux <- my_final_mbac
my_final_mbac$PLSmodels <- my_final_mbac$PLSmodels[1]

## ----inner, fig.cap="Plot of inner relations of PLS components. Only results for batch 'A' are shown as example. Each panel represents the inner correlation of one component of the PCA model. Red line indicates the diagonal when the correlation is maximal (1:1)."----
plot (my_final_mbac, type = "inner", comp2plot = c(1,2))

## ---- eval = TRUE, echo = FALSE--------------------------------------------
# restore state
my_final_mbac <- my_aux

## ----batchest, fig.cap="Batch effect estimation plot. Dashed lines represent theoretical batch magnitudes. Violin plots represent the distribution of batch effect coefficents observed in data."----
plot (my_final_mbac, type = "batch")

## ----pca-org, fig.cap="Default PCA plot on the original data."-------------
plot (my_final_mbac, type = "pca.org",
      asp = 1, cex.axis = 1, cex.lab = 1, cex = 3, bty = "L", 
      cex.main = 1.2, pch = 19)

## ----pca-cor, fig.cap="Default PCA plot on the corrected data."------------
plot (my_final_mbac, type = "pca.cor", 
      asp = 1, cex.axis = 1, cex.lab = 1, cex = 3, bty = "L", 
      cex.main = 1.2, pch = 19)

## ---- include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.width = 15, fig.height = 12)

## ----pcaplot2, fig.cap="Customized PCA plots. Original (left panels) and Corrected (right panels) data. Upper panels show the second principal component (PC) against the first one while panels at the bottom show the third PC against the first one."----
custom_col <- c("brown2", "dodgerblue", "forestgreen")
custom_pch <- c(19,19,19,1,1,1,15,15,15,0,0,0, # batch A
                  19,19,1,1,17,17,2,2,  # batch B
                  19,19,1,1,18,18,5,5)  # batch C

plot(my_final_mbac, type = "pca.both", col.by.batch = TRUE, 
     col.per.group = custom_col, comp2plot = 1:3,
     asp = 1, cex.axis = 1.3, cex.lab = 1.2, cex = 3, bty = "L", 
     cex.main = 1.7, pch = custom_pch,
     legend.text = c("Color", names(my_final_mbac$ListOfBatches),
                     "Shape", c("RNA", "GRO", "RIBO", "PAR"),
                     "Fill", levels(colData(my_final_mbac$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 19, 15, 17, 18, 
                                  NA, 19, 1),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, rep(1, 4),
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 2))

