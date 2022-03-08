# pLink-results-analysis

These R scripts and functions were used in the data processing of the described paper below. The notebook reads in pLink output files which are deposited on Pride.

The sctipts and functions reformat the data into flat tables that resemble an output one can get from the XlinkX node of proteome discoverer.\
The scripts clean the data and perform some processing which were useful in this project.

Some of the processing steps might not be suitable for some data sets. The functions are meant to be reused while the notebook should serve as a guide and an example.

**Characterization of protein complexes in extracellular vesicles by intact extracellular vesicle crosslinking mass spectrometry (iEVXL)**

Julia Bauzá-Martinez<sup>1,2</sup>, Gad Armony<sup>1,2</sup>, Matti F. Pronker<sup>1,2</sup>, Wei Wu<sup>1,2,3*</sup>

<sup>1</sup> **Biomolecular Mass Spectrometry and Proteomics**, Bijvoet Center for Biomolecular Research and Utrecht Institute for Pharmaceutical Sciences, Utrecht University, Padualaan 8, 3584 CH Utrecht, The Netherlands\
<sup>2</sup> **Netherlands Proteomics Centre**, Padualaan 8, 3584 CH Utrecht, The Netherlands\
<sup>3</sup> **Singapore Immunology Network (SIgN)**, Agency for Science, Technology and Research (A*STAR), Singapore, Singapore.

\*correspondence: Wei Wu, wu_wei@immunol.a-star.edu.sg or w.wu1@uu.nl  


## Instructions for use

### Environment
The code is embedded in a jupyter notebook which can be executed in several enviourments, in jupyter lab for example.\
Make sure that an R kernel is available to run the notebook.
The following packages need to be installed in R:
* tidyverse
* docstring (for rendring function documentation)

### Obtain the files
To obtain the functions, notebook, and data, use git to download ('clone') the files, or simply download the files as a zip file from the web interface. The URL for cloning, or the link to download the zip, are availabe under the green "Code" button above the file listing.


### Notes
* The frunctions assume that the protein names in the fasta file used to search in pLink, are in the uniprot format (sp|ACCESSION|NAME), while contaminants are not (do not start with sp| ).
* The renumber functions assumes that the fasta file used for the search in pLink was modified to remove signal peptides. These functions readjust the numers to keep the uniprot numbering in place.