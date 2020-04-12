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
@article{JE1,
author = {Pelechano, Vicent and P{\'e}rez-Ort{\'i}n, Jos{\'e} E.},
title = {There is a steady-state transcriptome in exponentially growing yeast cells},
journal = {Yeast},
volume = {27},
number = {7},
pages = {413-422},
keywords = {transcription rate, steady state, transcriptome, Saccharomyces cerevisiae, mRNA stability},
doi = {10.1002/yea.1768},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/yea.1768},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/yea.1768},
year = {2010}
}

@article{JE2,
	Annote = {doi: 10.1016/j.molcel.2004.06.004},
	Author = {Garc{\'i}a-Mart{\'i}nez, Jos{\'e} and Aranda, Agust{\'i}n and P{\'e}rez-Ort{\'i}n, Jos{\'e}E},
	Booktitle = {Molecular Cell},
	Date = {2004/07/23},
	Date-Added = {2019-01-29 12:15:03 +0000},
	Date-Modified = {2019-01-29 12:15:03 +0000},
	Doi = {10.1016/j.molcel.2004.06.004},
	Isbn = {1097-2765},
	Journal = {Molecular Cell},
	M3 = {doi: 10.1016/j.molcel.2004.06.004},
	Month = {2019/01/29},
	Number = {2},
	Pages = {303--313},
	Publisher = {Elsevier},
	Title = {Genomic Run-On Evaluates Transcription Rates for All Yeast Genes and Identifies Gene Regulatory Mechanisms},
	Ty = {JOUR},
	Url = {https://doi.org/10.1016/j.molcel.2004.06.004},
	Volume = {15},
	Year = {2004},
	Year1 = {2004},
	Bdsk-Url-1 = {https://doi.org/10.1016/j.molcel.2004.06.004},
	Bdsk-Url-2 = {http://dx.doi.org/10.1016/j.molcel.2004.06.004}}

@article{JE3,
	An = {21103382},
	Author = {Pelechano, Vicent and Ch{\'a}vez, Sebasti{\'a}n and P{\'e}rez-Ort{\'i}n, Jos{\'e} E.},
	Date = {2010/11/16},
	Date-Added = {2019-01-29 12:15:29 +0000},
	Date-Modified = {2019-01-29 12:15:29 +0000},
	Db = {PubMed},
	Doi = {10.1371/journal.pone.0015442},
	Isbn = {1932-6203},
	J1 = {PLoS One},
	Journal = {PloS one},
	L2 = {https://www.ncbi.nlm.nih.gov/pmc/PMC2982843/},
	Month = {11},
	Number = {11},
	Pages = {e15442; e15442--e15442},
	Publisher = {Public Library of Science},
	Title = {A complete set of nascent transcription rates for yeast genes},
	Ty = {JOUR},
	U1 = {21103382{$[$}pmid{$]$}},
	U2 = {PMC2982843{$[$}pmcid{$]$}},
	Url = {https://www.ncbi.nlm.nih.gov/pubmed/21103382},
	Volume = {5},
	Year = {2010},
	Bdsk-Url-1 = {https://www.ncbi.nlm.nih.gov/pubmed/21103382},
	Bdsk-Url-2 = {http://dx.doi.org/10.1371/journal.pone.0015442}}


@article{RIBO,
	An = {25119046},
	Author = {Zid, Brian M and O'Shea, Erin K},
	Date = {2014/10/02},
	Date-Added = {2019-05-15 08:08:55 +0000},
	Date-Modified = {2019-05-15 08:08:55 +0000},
	Db = {PubMed},
	Doi = {10.1038/nature13578},
	Et = {2014/08/03},
	Isbn = {1476-4687; 0028-0836},
	J1 = {Nature},
	Journal = {Nature},
	L2 = {https://www.ncbi.nlm.nih.gov/pmc/PMC4184922/},
	La = {eng},
	Month = {10},
	Number = {7520},
	Pages = {117--121},
	Title = {Promoter sequences direct cytoplasmic localization and translation of mRNAs during starvation in yeast},
	Ty = {JOUR},
	U1 = {25119046{$[$}pmid{$]$}},
	U2 = {PMC4184922{$[$}pmcid{$]$}},
	U4 = {nature13578{$[$}PII{$]$}},
	Url = {https://www.ncbi.nlm.nih.gov/pubmed/25119046},
	Volume = {514},
	Year = {2014},
	Bdsk-Url-1 = {https://www.ncbi.nlm.nih.gov/pubmed/25119046},
	Bdsk-Url-2 = {http://dx.doi.org/10.1038/nature13578}}

@article{PARCLIP,
	An = {23409723},
	Author = {Freeberg, Mallory A and Han, Ting and Moresco, James J and Kong, Andy and Yang, Yu-Cheng and Lu, Zhi John and Yates, John R and Kim, John K},
	Date = {2013/02/14},
	Date-Added = {2019-01-29 12:26:46 +0000},
	Date-Modified = {2019-01-29 12:26:46 +0000},
	Db = {PubMed},
	Doi = {10.1186/gb-2013-14-2-r13},
	Isbn = {1474-760X; 1465-6906},
	J1 = {Genome Biol},
	Journal = {Genome biology},
	L2 = {https://www.ncbi.nlm.nih.gov/pmc/PMC4053964/},
	Month = {02},
	Number = {2},
	Pages = {R13--R13},
	Publisher = {BioMed Central},
	Title = {Pervasive and dynamic protein binding sites of the mRNA transcriptome in Saccharomyces cerevisiae},
	Ty = {JOUR},
	U1 = {23409723{$[$}pmid{$]$}},
	U2 = {PMC4053964{$[$}pmcid{$]$}},
	U4 = {gb-2013-14-2-r13{$[$}PII{$]$}},
	Url = {https://www.ncbi.nlm.nih.gov/pubmed/23409723},
	Volume = {14},
	Year = {2013},
	Bdsk-Url-1 = {https://www.ncbi.nlm.nih.gov/pubmed/23409723},
	Bdsk-Url-2 = {http://dx.doi.org/10.1186/gb-2013-14-2-r13}}


}
