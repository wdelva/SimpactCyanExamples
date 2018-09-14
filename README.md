# SimpactCyanExamples

<!-- Created by Wim Delva and David Niyukuri, 14 September 2018 -->


This file accompanies the article “SimpactCyan: An Open-source Simulator for Individual-Based Models in HIV Epidemiology with R and Python Interfaces”, available at: <https://www.biorxiv.org/content/early/2018/09/12/411512>
<!-- This URL is a placeholder and must be replaced by the actual URL, once the paper has been uploaded. 
The original paper is available at: <http://www.journals.uchicago.edu/doi/full/10.1086/596510>
-->

## CONTENTS

This file contains the following information:

* [Code and data files](#code-and-data-files)
   * Code file
   * Data files
* [System and software requirements](#system-and-software-requirements)
* [Copyright and licensing information](#copyright-and-licensing-information)
* [Contact information](#contact-information)

## CODE AND DATA FILES 

Three types of files (not including this readme file) are stored on this repository: code files, data files and figure files.


### Code file

All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

  SimpactPaperPhyloExample.R -- This file contains the script to generate the data, run the phylogenetic analysis and produce the figures for the HIV transmission network, time-resolved phylogenetic tree, and the distribution of internal nodes in the reconstructed phylogenetic tree and simulated HIV transmission events.

### Data files

  his.seq.B.pol.j.fasta -- This file is a fasta file of the HIV subtype B pol gene sequence of the seed infection (root of the phylogenetic tree).

  SimpactPaperPhyloExample.RData -- This file is an RData object, containing the data needed for the reproduction of the figures in the phylo example. These data are: transNet.yrs.Old (HIV transmission network object), dater.tree (time-stamped phylogenetic tree), i.vec (vector of calendar years), int.node.vec (vector of the number of internal nodes in one-year calendar time intervals), and numb.tra (vector of the number of transmission events in one-year calendar time intervals).


### Figure files

  network_vsc.pdf -- Panel A of the figure shown in the phylo example.
  
  tree_vsc.pdf -- Panel B of the figure shown in the phylo example.
  
  events_vsc.pdf -- Panel C of the figure shown in the phylo example.  

 

## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system

  We have only tested this code on personal computers (OS X Version 10.11.6 and Linux Ubuntu Version 16.04), and on the XXXXXX cluster at the Cape Town Centre for High Performance Computing (CHPC) and the golett cluster of the Flemish Supercomputer Centre (VSC).

### Required software

  R version 3.4.4

  Seq-Gen version 1.3.4. <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network.

  FastTree version 2.1.10. <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset.

  SimpactCyan version 0.21 and RSimpactCyan. SimpactCyan is the core program that allows fast simulation of HIV transmission across a sexual network. RSimpactCyan is the R package that enables initiation and running of models built by SimpactCyan. Installation instructions for both are at: <https://github.com/j0r1/RSimpactCyan/blob/master/INSTALLATION.md>

  A long list of auxiliary R packages is required to run the post-simulation analysis for the MaxART and phylo examples in the paper”

install.packages("devtools")
install.packages("pacman")
library(devtools)
install_github("j0r1/readcsvcolumns/pkg")
install_github("wdelva/RSimpactHelp”, dependencies = TRUE)

p_load(Rcpp, ape, expoTree, data.table, readr, phangorn, dplyr, adephylo, treedater, geiger, picante, igraph, network, intergraph, ggtree, lubridate, ggplot2, ggnetwork, geomnet, RSimpactCyan, RSimpactHelper)
 

## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0 License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this work is suitable for commercial use, that licensees can modify the work, that they must release the source alongside with Derivative Work, and that Derivative Work must be released under the same terms.


## CONTACT INFORMATION

Please contact Prof Wim Delva with questions regarding the files provided and their authorized usage:

Wim Delva
Email: <DELVAW@sun.ac.za>


