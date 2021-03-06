% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{plot_pca}
\alias{plot_pca}
\title{plot_pca}
\usage{
plot_pca(mbac, col.by.batch = TRUE, col.per.group = NULL,
  comp2plot = c(1, 2), typeP = "pca.both", legend.text = NULL,
  args.legend = NULL, ...)
}
\arguments{
\item{mbac}{Object of class mbac generated by *createMbac*, *ARSyNbac*, *MultiBaC*, *genModelList*, or *batchCorrection*.}

\item{col.by.batch}{Argument for PCA plots. TRUE or FALSE. If TRUE (default) samples are gruped according to the batch factor. If FALSE samples are gruped according to the experimental desing.}

\item{col.per.group}{Argument for PCA plot. Indicates the color for each group defined in "groups" argument. If NULL (default) the colors are taken from a predefined pallete.}

\item{comp2plot}{Indicates which components are plotted. Default is "c(1,2)", which means that component 1 is plotted in "x" axis and component 2 in "y" axis. If more components are indicated, the function will return as many plots as needed to show all the components.}

\item{typeP}{"pca.cor", "pca.org" or "pca.both". If inputOmics contains original matrices, set "pca.org". However, if inputOmics contains the corrected matrices, set "pca.cor".}

\item{legend.text}{a vector of text used to construct a legend for the plot.}

\item{args.legend}{list of additional arguments to pass to legend(); names of the list are used as argument names. Only used if legend.text is supplied.}

\item{...}{Other graphical arguments.}
}
\value{
A PCA plot is displayed.
}
\description{
plot_pca
}
\examples{
data('multiyeast')

my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par),
                       batchFactor = c("A", "A", "B", "B", "C", "C"),
                       experimentalDesign = list("A" =  c("Glu+", "Glu+",
                       "Glu+", "Glu-", "Glu-", "Glu-"),
                       "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
                       "C" = c("Glu+", "Glu+", "Glu-", "Glu-")),
                       omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"))

plot_pca(my_mbac, typeP = "pca.org")

my_final_mbac <- MultiBaC (my_mbac,
                           test.comp = NULL, scale = FALSE,
                           center = TRUE, crossval = NULL,
                           Variability = 0.90,
                           Interaction = TRUE ,
                           showplot = FALSE,
                           showinfo = FALSE)

plot_pca(my_final_mbac, typeP = "pca.cor")

}
