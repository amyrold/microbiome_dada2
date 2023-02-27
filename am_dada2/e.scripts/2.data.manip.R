# Getting Ready ----
# Create lists of the fastq files
sample.names <- gsub('.fastq.gz','', substring(sort(list.files(p.data.raw)),47,))
fq.names <- sort(list.files(p.data.raw, full.names = TRUE))
#================================================================================================#

# Trim primers ----
# Earth Microbiome Project Primers
# FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT (rc)
# ATTAGAWACCCBNGTAGTCC (5'-3')

#CutAdapt
for (i in 1:length(fq.names)){ #Loop through fastq files and trim primers
  #Trim EMBP primers from file i in fastq_names and output to b.data.trim
  system(paste0(p.cutadapt,'cutadapt -a ATTAGAWACCCBNGTAGTCC -g GTGYCAGCMGCCGCGGTAA -o /cloud/project/b.data.trim/t.', sample.names[i],'.fastq.gz ',fq.names[i])) #substring(fq.names[i],38,nchar(fq.names[i]))
} 
fq.names <- sort(list.files(p.data.trim, full.names = TRUE))
sample.names <- gsub('.fastq.gz','', substring(sort(list.files(p.data.trim)),3,))
#=============================================================================#

# Inspect read quality profiles ----
# Use heatmap plots to visualize the quality scores of the forward and reverse reads
plotQualityProfile(fq.names[1:2])
#=============================================================================#

# Filter and trim ----
#prep the file paths and names
filtFq <- paste0(p.data.filtered, 'f.', sample.names, '.fastq.gz')
setwd(p.data.raw) # Store the fastq files in the raw data directory
names(filtFq) <- sample.names
#use the filterAndTrim f(x) to cut the data at points determined above (240 and 160)
out <- filterAndTrim(fq.names, filtFq, truncLen=253, trimLeft=5,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out) #Displays a dataframe consisting of .fastq files and their respective reads in and out
#=============================================================================#

# Learn the Error Rates ----
#use a machine learning function to estimate the error rates
errF <- learnErrors(filtFq, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE) #plot to visualize and see that everything looks in order
# in this case, the black line seems congruent with the points and the red line seems to make sense
# (according to the tutorial. i've never seen plots like this before so im no expert)
#=============================================================================#

# Sample Inference ----
#Use dada to perform the sample inference algorithm
dadaFs <- dada(filtFq, err=errF, multithread=FALSE)
dadaFs[[1]] # 76 true sequence variants
#=============================================================================#

# Construct Sequence table ----
#this is the creation of the actual ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab) #311x7876
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) #248,7876
#=============================================================================#

# Remove chimeras ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim) # 311, 1985
sum(seqtab.nochim)/sum(seqtab) #.8113
#=============================================================================#

# Track reads through the pipeline ----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
#=============================================================================#