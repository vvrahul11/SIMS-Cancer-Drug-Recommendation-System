# SIMS- Cancer-Drug-Recommendation-System

# Abstract

Non-small cell lung cancer (NSCLC) is a leading cause of death worldwide. Targeted monotherapies produce high regression rates, albeit for limited patient subgroups, who inevitably succumb. We present a novel strategy for identifying customized combinations of triplets of targeted agents, utilizing a simplified interventional mapping system (SIMS) that merges knowledge about existent drugs and their impact on the hallmarks of cancer. Based on interrogation of matched lung tumor and normal tissue using targeted genomic sequencing, copy number variation, transcriptomics, and miRNA expression, the activation status of 24 interventional nodes was elucidated. An algorithm was developed to create a scoring system that enables ranking of the activated interventional nodes for each patient. Based on the trends of co-activation at interventional points, combinations of drug triplets were defined in order to overcome resistance. This methodology will inform a prospective trial to be conducted by the WIN consortium, aiming to significantly impact survival in metastatic NSCLC and other malignancies.

Keywords: Tri-therapy, NSCLC, targeted therapies, algorithm, pathway

![alt text](https://github.com/vvrahul11/SIMS---Cancer-Drug-Recommendation-System/blob/master/sims_web_bgu/SIMS/WIN.png)


# Requirements
 * R
 * Shiny
 * R-studio

# Installing R and Shiny

Shiny works with R version >=3.0.0 (as of 15th-January- 2015). If you have an old version of R <3.0.0, please remove the version completely by executing the commands given below.

> sudo apt-get remove r-base-core
> sudo apt-get remove r-base sudo apt-get remove autoremove

Download the latest version of R from CRAN>= R-3.0.0 and install from source. Click on the below link to download the latest version of R.

http://cran.r-project.org/src/base/R-3/

- Download the latest version and unzip the file.
- Go to R directory and find the INSTALL file and follow the instructions. (Given below)

> ./configure --with-x=yes --enable-R-shlib

> make

> make check

> make pdf

> make info

> make install

> make install-info

> make install-pdf

Installing shiny
Execute the command inside the R console

> install.packages("shiny")

Installing shiny themes
> Install libcurl4-gnutls-dev first using terminal

> sudo apt-get -y build-dep libcurl4-gnutls-dev apt-get -y install
> sudo libcurl4-gnutls-dev

Install RStudio from [here](https://www.rstudio.com/)
# Usage
* Standalone version 

  > cd Standalone
  
  > Rscripts SIMS.standlone.v1.0.R
  
  or 
  
  Open R-studio project file 'Standalone.Rproj'
  
  * Web version
  
  > cd sims_web_bgu
  
  > R
  
  > library(shiny)
  
  > shiny::runApp("SIMS")
  
  