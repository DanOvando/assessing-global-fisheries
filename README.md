
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Improving Estimates of the State of Global Fisheries Depends on Better Data Data

## Fish and Fisheries (Link TBD)

-   Daniel Ovando
-   Ray Hilborn
-   Cole Monnahan
-   Merrill Rudd
-   Rishi Sharma
-   James Thorson
-   Yannick Rousseau
-   Yimin Ye

## Abstract

Implementation of the United Nations Sustainable Development Goals
requires assessments of the global state of fish populations. While we
have reliable estimates of stock status for fish populations accounting
for approximately half of recent global catch, our knowledge of the
state of the majority of the world’s ‘unassessed’ fish stocks remains
highly uncertain. Numerous publications have produced estimates of the
global status of these unassessed fisheries, but limited quantity and
quality of data along with methodological differences have produced
counterintuitive and conflicting results. Here, we show that despite
numerous efforts, our understanding of the status of global fish stocks
remains incomplete, even when new sources of broadly available data are
added. Estimates of fish populations based primarily on catch histories
on average performed 25% better than a random guess. But, on average
these methods assigned fisheries to the wrong FAO status category 57% of
the time. Within these broad summaries the performance of models trained
on our tested data sources varied widely across regions. Effective
improvement in estimates of the state of the world’s exploited fish
populations depends more on expanded collection of new information and
efficient use of existing data than development of new modeling
methods.s exploited fish populations depends on prioritizing the
collection of high-priority

*Mean classification accuracy (assignment to FAO stock status category)
by FAO statistical area arising from different data sources. Data source
panels are ordered in descending (starting from top left) mean accuracy
at the FAO region level. RLSADB Index refers to catch and abundance
index drawn from RLSADB. Effective CPUE refers to an index of abundance
based on reconstructed effort data. Effective CPUE+ uses CPUE along with
Fisheries Management Index (FMI) and/or swept area ratio (SAR) data. For
both CPUE series ‘nominal’ assumes a 0% technology creep, for
‘effective’ a 2.6% technology creep is assumed. FMI uses FMI scores to
develop a prior on recent fishing mortality rates, SAR does the same but
based on swept area ratio. CMSY uses the methods from Froese et al. 2017
\[@froese2017\]. Guess assigns a random recent B/B<sub>MSY</sub> of
0.4,1, or 1.6.* <img src="documents/figs/acc-map.png" width="2400" />

# Reproducing Results

All materials needed to reproduce our results and manuscript are
contained in this repository. In order to reproduce

1.  Fork the repository and clone to your machine

2.  Open R and set your working directory of the cloned repository (or
    just use RStudio projects)

3.  This project is set up with
    [`renv`](https://rstudio.github.io/renv/articles/renv.html) to
    manage package dependencies. Inside R (and with your working
    directory set correctly) run `renv::restore()`. This will install
    the correct versions of all the packages needed to replicate our
    results. Packages are installed in a stand-alone project library for
    this paper, and will not affect your installed R packages anywhere
    else.

4.  Run make-assessing-global-fisheries.R, setting all run\_ options to
    TRUE (this will likely take over 48 hours to run). This will knit
    the manuscript for this paper automatically, generating
    ovando-etal-assessing-global-fisheries.docx in the documents folder
