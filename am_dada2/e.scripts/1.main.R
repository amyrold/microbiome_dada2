# General ----
# Working Title:
# Authors:
# Questions and Motivation ----

# R version
R.version.string # "R version 4.2.1 (2022-06-23)"
rm(list = ls())

#==============================================================================#

# Script Index ----
# 1.main.R        
# 2.data.manip.R
# 3.analysis.R
# 4.functions.R

#==============================================================================#

# To do ----
# remove primers (gilbert)
# skip merging step
# make metadata table (for end of dada2)
# pick questions from the paper (different if possible)
# 

#==============================================================================#

# Global Variables ----
# store current user's working directory
setwd('/cloud/project')
wk.dir <- getwd() #"/cloud/project"

#==============================================================================#

# Folder Management ----
# names of folders for output data (figures + data output)
folder.names <- c("a.data.raw",'b.data.trim',"c.data.filtered", "d.results", 'e.scripts', "f.figures",  'g.trash')
# if folder with name "i" does not exist, create it.
for(i in 1:length(folder.names)){
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i]) 
  }
}

# paths to the folders. The 'p.' indicates the variable is a path.
# make sure the variable names describe the folder.names
p.data.raw <- paste(wk.dir, "/", folder.names[1], "/", sep = "")
p.data.trim <- paste(wk.dir, "/", folder.names[2], "/", sep = "")
p.data.filtered <- paste(wk.dir, "/", folder.names[3], "/", sep = "")
p.results <- paste(wk.dir, "/", folder.names[4], "/", sep = "")
p.fig <- paste(wk.dir, "/", folder.names[5], "/", sep = "")
p.cutadapt <- '~/.local/bin/'

#==============================================================================#

# Libraries & Packages ----
# load libraries needed for analyses

# Install bioconductor controller, BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
BiocManager::install("phyloseq");
BiocManager::install("Biostrings");

library(dada2); packageVersion("dada2") # Load dada2 and check the version
library('phyloseq') #download and activate phyloseq
library('Biostrings') #download and activate biostrings
library(ggplot2); packageVersion("ggplot2") #activate gg_plot
theme_set(theme_bw()) #set the theme to black and white
library(ShortRead)


#==============================================================================#
# Download Data ----
# This creates a file that contains a series of bash scripts
# The scripts are used to download the raw fastq files
if(file.exists("sra_explorer_fastq_download.sh") == FALSE){
  source("e.scripts/5.download.create.R")
}

# Using the scripts from the above file, store file names and then create files
dl.cmds <- readLines('sra_explorer_fastq_download.sh')
f.names <- strsplit(dl.cmds, ' -o ');f.names <- sort(unlist(f.names)[c(FALSE,TRUE)]) #f.names <- gsub('.gz','',unlist(f.names)[c(FALSE,TRUE)])
setwd(p.data.raw) # Store the fastq files in the raw data directory
for(i in 1:length(dl.cmds)){
  if(file.exists(f.names[i]) == FALSE){ # If the file already exists, there is no need to recreate
    system(paste(dl.cmds[i])) #download the zipped files
    #system(paste('gzip -d', f.names[i], sep=' ')) #unzip and replace the zipped files
  }
};setwd(wk.dir) #restore working directory
metadata <- read.csv('FeFiFo_demographics_by_participant.csv') #Info on sample hosts

#==============================================================================#

# Run Scripts ----
#source("e.scripts/4.functions.R")
#source("e.scripts/3.analysis.R")
#source("e.scripts/5.figures.R")
# END OF MAIN ----
#==============================================================================#

