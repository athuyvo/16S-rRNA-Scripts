# Author: Anh Vo (email: athuyvo_at_uw.edu)
# Script incorporates DADA2 pipeline. See source below. 

# Source code info:
#   Title: DADA2 Pipeline Tutorial
#   Author: Benjamin Callahan
#   Code version: 1.12
#   Availability: https://benjjneb.github.io/dada2/tutorial.html

# This script runs the DADA2 pipeline for paired-end analysis on
# Illumina-sequenced demultiplexed fastq files.
# Requires DADA2 and ggplot2 packages.

# Import DADA2 and ggplot2 package
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")

# Sort files to ensure forward/reverse reads are in same order
# fastq files should be in .gz format
# Assign sample names using fastq file names
sortFastqs <- function() {

	fnFs <<- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
	fnRs <<- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
	sample.names <<- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

	cat("Forward files:\n")
	print(head(fnFs))
	print(head(fnRs))

	cat("Sample names:\n")
	print(head(sample.names))

	# Filter and trim reads, place files in filtered subdirectory.
	filtFs <<- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <<- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <<- sample.names
	names(filtRs) <<- sample.names
}



# Plot quality scores for forward and reverse reads using first two samples
# Saves plots in .pdf files to working directory
plotQScore <- function() {
	plotQualF <- plotQualityProfile(fnFs[1:2])
	ggsave("ForPlot.pdf", plotQualF, device="pdf")
	plotQualR <- plotQualityProfile(fnRs[1:2])
	ggsave("RevPlot.pdf", plotQualR, device="pdf")
}

# Plot estimated forward and reverse sequencing errors using Illumina quality scores
# Outputs error plots in .pdf files to working directory
plotErrorModel<- function() {
	plotErrF <- plotErrors(errF, nominalQ=TRUE)
	ggsave("ForErrPlot.pdf", plotErrF, device="pdf")
	plotErrR <- plotErrors(errR, nominalQ=TRUE)
	ggsave("RevErrPlot.pdf", plotErrR, device="pdf")
}

# Filters and trim reads according to user settings
# truncLen truncates reads at bp R1 and bp for R2
# maxN is the max number of ambiguous bases allowed
# maxEE is the number of estimated errors allowed for a read (default settings read1=15, read2=20)
# truncQ truncates reads with at specified Q score
# rm.phiX removes phiX control reads

# Outputs reads left after trimming
# Prompts user whether or not to rerun this filtering step
# Saves current RSession

filter <- function(read1, read2) {
	cat("Trimming reads... \n")
	out <<- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(read1, read2), minLen=50, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
	
	save.image(file='RSession.RData')
	savehistory(file='history.txt')

	cat("Trimmed reads: \n")
	print(out)
	cat("Percent Trimmed reads: \n")
	print(out[,2]/out[,1]*100)

	answer <- rerun("filtering?")
	return (answer)
}

# Learn forward and reverse errors using DADA2 algorithm
# Algorithm assumes max possible error rates.
# Apply DADA2 inference algorithm to filtered and trimmed reads using estimated error model  
# to denoise data by diciphering real vs spurious sequences.
# Merge denoised reads and output paired reads that overlap by at least 50 bases.
# Save error rates and current RSession 

denoiseAndMerge <- function() { 
	cat("Initializing error rates to the maximum possible estimate...\n")
	cat("R1: ")
	errF <<- learnErrors(filtFs, nbases=2.5e+08, multithread=TRUE)
	filepath <- file.path(wd, "errF.rds")
	saveRDS(errF, filepath)
	save.image(file='RSession.RData')

	cat("R2: ")
	errR <<- learnErrors(filtRs, nbases=2.5e+08, multithread=TRUE)
	filepath <- file.path(wd, "errR.rds")
	saveRDS(errR, filepath)
	save.image(file='RSession.RData')

	plotErrorModel()

	cat("Applying error model...\n")
	cat("R1: \n")
	#dadaFs <- dada(filtFs, err=errF, pool="pseudo", multithread=TRUE)
	dadaFs <<- dada(filtFs, err=errF, multithread=TRUE)
	filepath <- file.path(wd, "dadaFs.rds")
	saveRDS(dadaFs, filepath)
	save.image(file='RSession.RData')
	
	cat("R2: \n")
	#dadaRs <- dada(filtRs, err=errR, pool="pseudo", multithread=TRUE)
	dadaRs <<- dada(filtRs, err=errR, multithread=TRUE)
	#dadaFs[[1]]
	filepath <- file.path(wd, "dadaRs.rds")
	saveRDS(dadaRs, filepath)
	save.image(file='RSession.RData')
 
	mergers <<- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=50, verbose=TRUE)
	# Inspect the merger data.frame from the first sample
	head(mergers)
	save.image(file='RSession.RData')
}


# Remove chimeras from merged reads using input method analysis and 
# min input number of parent bimeras abundance to classify as a chimera.  
# Prompts user whether or not to rerun chimera removal step
# Default method= "per=sample"
# Default minParent = 700

# Saves sequence table, error rates, and current RSession 

removeChimeras <- function(method, minParent) {
	seqtab <<- makeSequenceTable(mergers)
	cat("Number of amplicon sequence variants: \n")
	print(dim(seqtab))
	filePath <- file.path(wd, "seqtab.rds")
	saveRDS(seqtab, filePath)

	# Inspect distribution of sequence lengths
	cat("Number of samples with distribution of sequence lengths: ")
	print(table(nchar(getSequences(seqtab))))

	# Remove chimerasFs <- dada(filtFs, err=errF, multithread=TRUE)
	seqtab.nochim <<- removeBimeraDenovo(seqtab, method=method, minParentAbundance=minParent, multithread=TRUE, verbose=TRUE)
	save.image(file='RSession.RData')
	savehistory(file='history.txt')

	cat("Number of samples and reads after chimera removal: ")
	print(dim(seqtab.nochim))

	cat("Proportion of reads remaining after chimera removal: ")
	print(sum(seqtab.nochim)/sum(seqtab))
	print(table(nchar(getSequences(seqtab.nochim))))
	filePath <- file.path(wd, "seqtabnochim.rds")
	saveRDS(seqtab.nochim, filePath)

	answer <- rerun("chimera removal?")
	return (answer)
}

# Filter nonchimeric reads by given read length	
# Prompts user whether or not to rerun length filtering
# Saves final sequences table

filterLen <- function() {
	cat("Set min length: ")
	MINLEN <- readline()
	cat("Set max length: ")
	MAXLEN <- readline()
	seqlens <- nchar(getSequences(seqtab.nochim))
	seqtab.filt <<- seqtab.nochim[,seqlens >= MINLEN & seqlens <= MAXLEN]
	print(sum(seqtab.filt)/sum(seqtab))
	filePath <- file.path(wd, "seqtabfilt.rds")
	saveRDS(seqtab.filt, filePath)

	answer <- rerun("filter read length step?")
	return (answer)
}

# Outputs number of reads that passed through each step in pipeline

trackReads <- function () {
	getN <- function(x) sum(getUniques(x))
	track <<- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100,1), rowSums(seqtab.filt), round(rowSums(seqtab.filt)/out[,1]*100,1))
	# ** If processing a single sample, remove the sapply calls:
	# e.g. replace sapply(dadaFs, getN) with getN(dadaFs) **
	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "perc_rds_retained", "len_filt", "final_per_rds_remain")
	rownames(track) <- sample.names
	head(track)
	filePath <- file.path(wd, "final_reads_filt.tsv")
	write.table(track, filePath, sep="\t", quote=F, col.names=NA)
}

# Assign taxonomy to sequence variants using Silva database
# Saves RSession and taxa as R objects
getTaxa <- function() {
	taxa <<- assignTaxonomy(seqtab.filt, "silva_nr_v132_train_set", taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), multithread=TRUE)
	save.image(file='RSession.RData')
	filePath <- file.path(wd, "taxa.rds")
	saveRDS(taxa,filePath)
}

# Write all ASV and taxonomy tables

writeFiles <- function() {
	taxa.print <- taxa #Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)
	asv_seqs <- colnames(seqtab.filt)
	asv_headers <- vector(dim(seqtab.filt)[2], mode="character")

	# Make read-friendly ASV and taxonomy tables 

	for (i in 1:dim(seqtab.filt)[2]) {
	    asv_headers[i] <- paste(">ASV", i, sep="_")
	}

	asv_fasta <- c(rbind(asv_headers, asv_seqs))
	asv<- t(seqtab.filt)
	row.names(asv) <- sub(">", "", asv_headers)
	asv_tax <- taxa
	row.names(asv_tax) <- sub(">", "", asv_headers)
	asv_tax_table <- merge(asv, asv_tax, by=0)
	fmla <- . ~ Kingdom + Phylum + Class + Order + Family + Genus
	asv_tax_table[] <- lapply(asv_tax_table, function(x) type.convert(as.character(x)))
	asv_tax_table_agg = aggregate(formula=fmla, data=subset(asv_tax_table, select=-Row.names), FUN=sum)

	# Writing out output files
	filePath <- file.path(wd, "ASV_fasta_filt.fa")
	write(asv_fasta, filePath)
	writeOut(asv, "ASV_filt.tsv")
	writeOut(asv_tax, "taxonomy_filt.tsv")
	writeOut(asv_tax_table, "ASV_tax_table.tsv")
	writeOut(asv_tax_table_agg, "ASV_tax_table_agg_filt.tsv")

	save.image(file='RSession.RData')
	savehistory(file='history.txt')
}

# Set path to fastq files and working directory
setPath <- function() {
	cat("Set fastq file path w/o quotes: ")
	path <<- readline()
	cat("Path is: ")
	print(path)
	cat("Confirm path: (y/n)")
	answer <- readline()
	return (answer)
}

# Set path to write files 
setWritePath <- function() {
	cat("Set directory to write files w/o quotes: ")
	wd <<- readline()
	cat("Write directory is: ")
	print(wd)
	cat("Confirm write directory: (y/n)")
	answer <- readline()
	return (answer)
}

# Write out tables 
writeOut <- function(table, filename) {
	filePath <- file.path(wd, filename)
	write.table(table, filePath, sep="\t", quote=F, col.names=NA)
}

# Prompt user whether or not to rerun functions
rerun <- function(string) {
	cat("Rerun", string, "(y/n): \n")
	answer <- readline()
	return (answer)
}


# Run all DADA2 filtering functions 
# Rerun any filtering steps when prompted by user
main <- function() {
	
	answer <- NULL
		
	while(answer != "y") {
		answer <- setPath()
	}

	setwd(path)
	answer <- setWritePath()

	while (answer != "y") {
		answer <- setWritePath()
	}

	cat("Fastq files:\n")
	print(list.files(path))
	cat("Proceed?: (y/n")
	answer <- readline()


	sortFastqs()
	plotQScore()
	answer <- filter(15,20)
	
	while(answer != "n") {
		cat("Set read1 maxEE: ")
		read1 <- readline() 
		cat("Set read2 maxEE: ")
		read2 <- readline()
		answer <- filter(read1, read2)
	}

	denoiseAndMerge()
	answer <- removeChimeras("per-sample", 700)

	while (answer != "n") {
		cat("Set chimera analysis method: ")
		method <- readline() 
		cat("Set min parent abundance: ")
		minParent <- readline()
		answer <- removeChimeras(method, minParent)
	}
	
	answer <- filterLen()

	while (answer != "n") {
		answer <- filterLen()
	}

	trackReads()
	cat("Assign taxonomy?(y/n): \n")
	answer <- readline()
	if (answer == "y") {
		getTaxa()
		writeFiles()
	}
	save.image(file='RSession.RData')
	savehistory(file='history.txt')
	
	cat("Done. Session saved.")
}

