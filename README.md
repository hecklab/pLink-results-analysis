# pLink-results-analysis

This jupyter notebook and functions were used in the data processing of the described paper below. The notebook reads in pLink output files which are deposited on Pride.

The notebook and functions reformat the data into flat tables that resemble an output one can get from the XlinkX node of proteome discoverer.

**Characterization of protein complexes in extracellular vesicles by intact extracellular vesicle crosslinking mass spectrometry (iEVXL)**

Julia Bauzá-Martinez<sup>1,2</sup>, Gad Armony<sup>1,2</sup>, Matti F. Pronker<sup>1,2</sup>, Wei Wu<sup>1,2,3*</sup>

<sup>1</sup> **Biomolecular Mass Spectrometry and Proteomics**, Bijvoet Center for Biomolecular Research and Utrecht Institute for Pharmaceutical Sciences, Utrecht University, Padualaan 8, 3584 CH Utrecht, The Netherlands\
<sup>2</sup> **Netherlands Proteomics Centre**, Padualaan 8, 3584 CH Utrecht, The Netherlands\
<sup>3</sup> **Singapore Immunology Network (SIgN)**, Agency for Science, Technology and Research (A*STAR), Singapore, Singapore.

\*correspondence: Wei Wu, wu_wei@immunol.a-star.edu.sg or w.wu1@uu.nl  



## Instructions for use

### Environment
The code is embedded in RMarkdown documents, one per analysis. RStudio (https://rstudio.com, no affiliation) is a convenient environment for "knitting" these documents, to create HTML or PDF output. The following packages need to be installed in R:

* tidyverse
* ggpubr
* reshape2
* stringr
* colorspace
* ggforce
* RColorbrewer
* VennDiagram
* psych

### Obtain the markdown documents
To obtain these documents, use git (available in RStudio as well) to download ('clone') the documents, or simply download the files as a zip file. The URL for cloning, or the link to download the zip, are availabe under the green "Code" button above the file listing.
We have added a project file for convenience, so you can double-click it to open the project in RStudio.

### Obtain the data
No data (MaxQuant output .txt files) is available in this repository, it needs to be downloaded from the Pride archive.
From there, obtain the MQ_output_txt.zip file and extract it in the same directory as the markdown documents. The scripts will locate the required .txt file in the subfolder and load all necessary libraries.
Second, open the .Rmd file you are interested in and 'knit' the document. The scripts will generate the plots used in the above mentioned paper. 
