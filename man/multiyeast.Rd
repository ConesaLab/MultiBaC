\name{multiyeast}
\alias{multiyeast}
\docType{data}
\title{ A distributed yeast multiomic dataset }
\description{
The yeast expression data sets were collected from the Gene Expression Omnibus (GEO) database and from three different studies. All of them analyzed the effects of glucose starvation in yeast. Lab A is the Department of Biochemistry and Molecular Biology from Universitat de Valencia (accession number GSE11521) [@JE1, @JE2, @JE3]; Lab B is the Department of Molecular and Cellular Biology from Harvard University (accession number GSE56622) [@RIBO]; and Lab C is the Department of Biology from Johns Hopkins University (accession number GSE43747) [@PARCLIP]. These studies used equivalent yeast strains and experimental conditions but, as shown in Figure 1, the main effect on expression is due to data belonging to different labs, which are the batches in this case.
}
\usage{data("multiyeast")}
\format{
  6 matrices with 6 or 4 observations on 200 variables (genes).
  \describe{
    \item{\code{A.rna}}{Gene expression data matrix from lab A.}
    \item{\code{A.gro}}{Transcription rates data matrix from lab A.}
    \item{\code{B.rna}}{Gene expression data matrix from lab B.}
    \item{\code{B.ribo}}{Translation rates data matrix from lab B.}
    \item{\code{C.rna}}{Gene expression data matrix from lab C.}
    \item{\code{C.par}}{Protein binded RNA data matrix from lab C.}
  }
}
\details{
The yeast expression data sets were collected from the Gene Expression Omnibus (GEO) database and from three different studies. All of them analyzed the effects of glucose starvation in yeast. Lab A is the Department of Biochemistry and Molecular Biology from Universitat de Valencia (accession number GSE11521) [@JE1, @JE2, @JE3]; Lab B is the Department of Molecular and Cellular Biology from Harvard University (accession number GSE56622) [@RIBO]; and Lab C is the Department of Biology from Johns Hopkins University (accession number GSE43747) [@PARCLIP]. These studies used equivalent yeast strains and experimental conditions but, as shown in Figure 1, the main effect on expression is due to data belonging to different labs, which are the batches in this case.
After a proper data pre-processing for each case, a voom transformation (limma R package) was applied when necessary. Finally TMM normalization was performed on the whole set of samples from all labs. A reduced dataset was obtained by selecting 200 omic variables from each data matrix and just 3 samples from lab A.  This yeast multiomic reduced dataset is included in *MutiBaC* package to illustrate the usage of the package. The gene expression matrices can be loaded by using the *data("multiyeast")* instruction.
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
}

\examples{
data(multiyeast)
head(A.rna)
}
\keyword{datasets, yeast, multiomic}

\references{
[1] Pelechano, V. and Pérez‐Ortín, J.E. (2010), There is a steady‐state transcriptome in exponentially growing yeast cells. Yeast, 27: 413-422. https://doi.org/10.1002/yea.1768

[2] García-Martínez, J. and Aranda, A. and Pérez-Ortín, J. E. (2004). Genomic Run-On Evaluates Transcription Rates for All Yeast Genes and Identifies Gene Regulatory Mechanisms, Molecular Cell, 15. https://doi.org/10.1016/j.molcel.2004.06.004

[3] Pelechano, V., Chávez, S., & Pérez-Ortín, J. E. (2010). A complete set of nascent transcription rates for yeast genes. PloS one, 5(11), e15442. https://doi.org/10.1371/journal.pone.0015442

[4] Zid, B. M., & O'Shea, E. K. (2014). Promoter sequences direct cytoplasmic localization and translation of mRNAs during starvation in yeast. Nature, 514(7520), 117–121. https://doi.org/10.1038/nature13578

[5] Freeberg, M. A., Han, T., Moresco, J. J., Kong, A., Yang, Y. C., Lu, Z. J., Yates, J. R., & Kim, J. K. (2013). Pervasive and dynamic protein binding sites of the mRNA transcriptome in Saccharomyces cerevisiae. Genome biology, 14(2), R13. https://doi.org/10.1186/gb-2013-14-2-r13
}
