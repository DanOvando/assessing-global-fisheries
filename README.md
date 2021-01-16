#  Status of Global Unassessed Fisheries will Remain Highly Uncertain without Better Data
  
  - Daniel Ovando
  - Ray Hilborn
  - Cole Monnahan
  - Merrill Rudd
  - Rishi Sharma
  - James Thorson
  - Yannick Rousseau
  - Yimin Ye
  
  Assessments of the global state of fish populations play an important role in tracking the implementation of the United Nations Sustainable Development Goals. While we have reliable estimates of stock status for fish populations accounting for approximately half of recent global catch, our knowledge of the state of the majority of the word's 'unassessed' fish stocks remains highly uncertain. Numerous publications have produced estimates of the global status of these unassessed fisheries, but limited quantity and quality of data along with methodological differences have produced counterintuitive and conflicting results. Here, we show that despite numerous efforts, our understanding of the status of global fish stocks remains incomplete, even when new sources of broadly available data are added. Estimates of fish populations based primarily on catch histories alone frequently performed roughly as well as a random guess. Obtaining accurate estimates of stock status for the world's exploited fish populations depends on prioritizing the collection of high-priority
  
# Reproducing Results

All materials needed to reproduce our results and manuscript are contained in this repository. In order to reproduce

1. Fork the repository and clone to your machine

2. Open R and set your working directory of the cloned repository (or just use RStudio projects)

3. This project is set up with [`renv`](https://rstudio.github.io/renv/articles/renv.html) to manage package dependencies. Inside R (and with your working directory set correctly) run `renv::restore()`. This will install the correct versions of all the packages needed to replicate our results. Packages are installed in a stand-alone project library for this paper, and will not affect your installed R packages anywhere else. 


4. Run make-assessing-global-fisheries.R, setting all run_ options to TRUE (this will likey take over 48 hours to run). This will knit the manuscript for this paper automatically, generating ovando-etal-assessing-global-fisheries.docx in the documents folder