#==================================================#
# Assign taxonomy ----

# For alternative method, install decipher
BiocManager::install("DECIPHER") #used package viewer to turn on instead of library()
#Download respective training set
#system('wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData')
# Use DECIPHER to perform the assign taxonomy task
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("/cloud/project/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# use the IdTaxa function from DECIPHER to begin assigning taxonomy
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# specify taxonomic ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) { #create an in-line function
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
# manupulate taxid to fit correct format (set col names and row names)
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid

#=============================================================================#

# Handoff to phyloseq ----
# Structuring metadata
name_ids <- data.frame(Sample.name=sample.names,Participant=substring(sample.names,1,4))
metadata_trimmed <- metadata[,c(1:5,7,8,10)] #Pick out 'important' data from metadata that I want to use for analysis *could change? does it hurt to include all of it even if some of it might not seem relevant?*
metadata_by_sample <- merge(name_ids,metadata_trimmed,by='Participant') #Merge the trimmed metadata with the sample names by the participant ID 
rownames(metadata_by_sample) <- rownames(seqtab.nochim)

#Create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata_by_sample), 
               tax_table(taxid))

#Turn the names into ASV21 instead of full DNA strings. 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
# Store actualy DNA strings into the data frame incase needed later
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#==============================================================================#

# Plotting ----
#Plot alpha diversity by education level
plot_richness(ps, x="What.is.your.highest.level.of.education.", measures=c("Shannon", "Simpson"),color='What.is.your.highest.level.of.education.')

#Plot alpha diversity by diet
plot_richness(ps,x='Study.treatment.group',measures=c('Shannon','Simpson'),color='Study.treatment.group')




# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

# plot the Bray NMDS
plot_ordination(ps.prop, ord.nmds.bray,color='What.is.your.highest.level.of.education.',title="Bray NMDS")
plot_ordination(ps.prop, ord.nmds.bray,color='Study.treatment.group',title="Bray NMDS")



#Create a bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Participant", fill='family') + facet_wrap(~Study.treatment.group, scales="free_x")


#test
plot_bar(ps.top20, x="Study.treatment.group") + facet_wrap(~family, scales="free_x")
#==============================================================================#

