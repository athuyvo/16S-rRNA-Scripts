# Author: Anh Vo (email: athuyvo_at_uw.edu)

# This source code performs alpha and beta diversity on
# 16S sequences using a given OTU, taxonomy, and sample table. 


install.packages(c("vegan", "metacoder", "taxa", "ggplot2", "dplyr", 
								"readr", "stringr", "agricolae", "ape"),
                 repos = "http://cran.rstudio.com", dependencies = TRUE)



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("DECIPHER")

install.packages("phangorn",repos="https://cloud.r-project.org",quiet=TRUE)

opts_chunk$set(cache = FALSE,fig.path="dadafigure/")

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

library(gridExtra)
library(knitr)
library(BiocStyle)
library(readr)
library(dplyr)
library(taxa)
library(phyloseq)
library(metacoder)
library(ape)
library(vegan)
library(ggplot2)
library(stringr)
library(agricolae)



# Import OTU, taxonomy, and sample tables
# Files must be in .txt format
# OTU data table should have "OTU_ID" as name
# Taxa data table should hve "OTU ID" as name 


read_tables(otu, taxa, sample) {
	org_otu_data <<- read_tsv(otu)
	tax_data <<- read_tsv(taxa)
	sample_data <<- read_tsv(sample)
}


print_tables <- function(otu_data, tax_data, sample_data) {
	cat("OTU table\n")
	print(otu_data)
	cat("Tax table\n")
	print(tax_data)
	cat(" table\n")
	print(sample_data)
}

# Combine taxonomy and OTU tables
# Create and return new taxmap object and filter out reads below 10
new_obj <- function(otu_data, tax_data, sample_data) {
	tax_data$`OTU ID` <- sub(tax_data$`OTU ID`, pattern = "OTU_", replacement = "")
	cat("tax table\n")
	print(tax_data)  

	# Combine OTU and tax table
	tax_data$`OTU ID` <- as.character(tax_data$`OTU ID`)
	otu_data$OTU_ID <- as.character(otu_data$OTU_ID) 
	otu_data <- left_join(otu_data, tax_data, by = c("OTU_ID" = "OTU ID"))
	cat("otu table\n")
	print(otu_data)
	
	cat("col names\n")
	tail(colnames(otu_data), n = 10)

	# Make new taxmap obj
	obj <- parse_tax_data(otu_data,
          class_cols = "taxonomy",
          class_sep = ";",
          class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
          class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
	cat("object\n")
	print(obj)
	cat("tax data\n")
	print(obj$data$tax_data)

	obj$data$class_data <- NULL
	names(obj$data) <- "otu_counts" #rename "tax_data" table
	cat("object\n")
	print(obj)

	# filter out associated taxonomy with less than 10
	obj$data$otu_counts <- zero_low_counts(obj, "otu_counts", min_count = 10,
                                       other_cols = TRUE)
	no_reads <- rowSums(obj$data$otu_counts[, sample_data$`Sample ID`]) == 0
	cat("num of OTUs without reads\n")
	sum(no_reads)

	# create new object 
	obj <- filter_obs(obj, "otu_counts", ! no_reads, drop_taxa = TRUE)
	print(obj)

	return (obj)
}


# Make histogram of read distribution 
make_hist <- function(obj, sample_data) {
	hist = hist(colSums(obj$data$otu_counts[, sample_data$`Sample ID`]))
	makePDF("ReadDistribution", hist, NULL, NULL)
}

# Plots given index type graph by given analysis type
# compare_alpha(sample_data, "Visit", "shannon")
compare_alpha <- function(sample_data, grouping_var, index_type) {
  # Calcuate alpha diversity
  sample_data$alpha <- diversity(obj$data$otu_counts[, sample_data$`Sample ID`],
                                 MARGIN = 2,
                                 index = index_type)

  # Write alpha values to table
  write.table(sample_data, paste(index_type, "Values.tsv", sep=""), sep="\t", quote=F, col.names=NA)

  # Do ANOVA
  sample_data$grouping <- sample_data[[grouping_var]] # needed for how `aov` works
  anova_result <- aov(alpha ~ grouping, sample_data)
  cat("ANOVA:\n")
  print(anova_result)
  sink(paste(index_type, "AnovaValues", sep=""))
  print(anova_result)
  sink()

  # Plot alpha result by grouping variable 
  pdf(paste(index_type, ".pdf", sep=""))
  print(ggplot(sample_data, aes(x = grouping, y = alpha)) + 
  geom_boxplot() + ggtitle("Alpha diversity") +
	xlab(grouping_var) +
	ylab("Alpha diversity index"))
  dev.off()
  
  # Do Tukey's HSD test
  tukey_result <- HSD.test(anova_result, "grouping", group = TRUE)
  cat("Turkey:\n")
  print(tukey_result)
  sink(paste(index_type, "TukeyValues", sep=""))
  print(tukey_result)
  sink()

  # Plot alpha result with Tukey test
  group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]

	pdf(paste(index_type, "Tukey.pdf", sep=""))
	print(ggplot(sample_data, aes(x = grouping, y = alpha)) +
	geom_text(data = data.frame(),
	          aes(x = rownames(group_data),
	              y = max(sample_data$alpha) + 1,
	              label = group_data$groups),
	          col = 'black',
	          size = 10) +
	geom_boxplot() +
	ggtitle("Alpha diversity") +
	xlab(grouping_var) +
	ylab("Alpha diversity index"))
	dev.off()
}


# Create and returns new phyloseq object with given taxmap object and sample data table
new_ps <- function (obj, sample_data, tree) {
	ps <- as_phyloseq(obj, otu_table = "otu_counts", otu_id_col = "OTU_ID",
	                sample_data = sample_data, sample_id_col = "Sample ID", phy_tree = tree)
	return (ps)
}


# Plot a variety of alpha diversity indexes per sample with given phyloseq object 
plot_alpha <- function(ps_obj, analysis_type, group_color, visit) {
	plot = plot_richness(ps_obj, color = group_color, x = analysis_type) 
	measures=c("Shannon", "Simpson")
	makePDF("AlphaPlot", plot, analysis_type, NULL)
}


# rarefy objects and remove OTUs that have no reads after rarefied 
# obj$data$otu_rarefied <- rarefy_obs(obj, "otu_counts", other_cols = TRUE)
# cat("count of rarefied samples")
# print(obj)
# no_reads <- rowSums(obj$data$otu_rarefied[, sample_data$SampleID]) == 0
# obj <- filter_obs(obj, "otu_rarefied", ! no_reads)
# print(obj)



# Plot phylogenetic tree
# Calculate distances and plot Bray and Unifract plots
plot_beta <- function (obj, sample_table, analysis_type, shape, visit, connect) {
	ps_obj <- new_ps(obj, sample_table, NULL)

	#Plot phylogentic trees 
	# Plot a naked tree
	naked_tree = rtree(ntaxa(ps_obj), rooted=TRUE, tip.label=taxa_names(ps_obj))
	makePDF("NakedTree", naked_tree, analysis_type, visit)

	# Create new phyloseq object with naked tree
	  ps_obj <- as_phyloseq(obj, otu_table = "otu_counts", otu_id_col = "OTU_ID",
	                      sample_data = sample_table, sample_id_col = "Sample ID", phy_tree =naked_tree)

	  phylotree = plot_tree(ps_obj, color=analysis_type, label.tips="taxa_names", ladderize="left", plot.margin=0.3)
	  makePDF("PhyloTree", phylotree, analysis_type, visit)

	  beta_dist(ps_obj, "bray", analysis_type, shape, visit, connect)
	  beta_dist(ps_obj, "uunifrac", analysis_type, shape, visit, connect)
	  beta_dist(ps_obj, "wunifrac", analysis_type, shape, visit, connect)
 	# beta_dist(ps_obj, "sorensen", analysis_type, shape, visit, connect)	
	# beta_dist(ps_obj, "jaccard", analysis_type, shape, visit, connect)
}

# Plots PCoA distances using input beta-index type
# If true, segment lines will connect between two paired vectors 
beta_dist <- function (ps_obj, distance, analysis_type, shape, visit, connect) {
	ps_ord <- ordinate(ps_obj, method = "PCoA", distance= distance)
	# print("ps_ord")
	# print(ps_ord$vectors) 
	shape <- as.factor(shape)
	plot = plot_ordination(ps_obj, ps_ord, type= "samples", color=analysis_type) + geom_point(show.legend=FALSE, aes(shape = shape, size= 3, alpha = 0.5)) # + geom_text(aes(label=sample$`Visit`))

	# Add segment lines for each pair of points in plot
	# Specific for 2 visits and sample treatment in sample data

	#length = dim(ps_ord$vectors)[1]
	if (connect == "TRUE") {
		plot <-	new_df(ps_ord, plot)
	}
	
	makePDF(distance, plot, analysis_type, visit)
}


# Create and return a new data frame for connecting points 
new_df <- function(ps_ord, plot) {
	length = dim(ps_ord$vectors)[1]
	pointA = 1
	pointB = length/2
	
	# If comparing different treatment limb, set second connecting point to matching visit limb
	if(sample$Treatment[1] == "Azithromycin" & sample$Treatment[length] == "Placebo") {
		pointB = length/4 - 1
	}

	axis1 = 1 
	axis2 = 2
	
		# Point A
	x1 <- ps_ord$vector[pointA, axis1] 
	y1 <- ps_ord$vector[pointA, axis2]

	# Point B
	x2 <- ps_ord$vector[pointB, axis1]
	y2 <- ps_ord$vector[pointB, axis2]

	#  Create data frame for paired points
	df <- data.frame(x1,y1, x2,y2)
	segment1 = NULL
	# print("df")
	# print(df)

	count = 0

	for (i in ps_ord$vector[,1]) {

		if (pointB != length) {
			pointA = pointA + 1
			pointB = pointB + 1

			# Reset A and B to next Treatment Limb
			if (sample$Treatment[pointB] == "Placebo" & count != 1) {
				count = 1
				pointA= pointB
				
				#Find matching 2nd connecting point 
				midpoint = (length - (pointA - 1))/2  
				pointB = pointA + midpoint

				# Set connecting lines for first treatment limb and reset data frame  
				segment1= geom_segment(data=df, aes(x=x1, y=y1, xend= x2, yend= y2), color= "red", size = 0.75, alpha=0.5) 
				df = NULL
			}	

			# print(sample$`Sample ID`[pointA])
			# print(sample$`Sample ID`[pointB])	

			# Point A
			x1 <- ps_ord$vector[pointA, axis1] 
			y1 <- ps_ord$vector[pointA, axis2]

			# Point B
			x2 <- ps_ord$vector[pointB, axis1]
			y2 <- ps_ord$vector[pointB, axis2]

			newdf <- data.frame(x1,y1,x2,y2)
			df <- rbind(df, newdf)

				# print("newdf")
				# print(df)
		}
	}
	segment2 = geom_segment(data=df, aes(x=x1, y=y1, xend= x2, yend= y2), color= "cyan3", size = 0.75, alpha=0.5) 
	plot = plot + segment1 + segment2
	return(plot)
}

# Create and return new sample and otu data.frames without input visit
# using given sample data sheet 
new_sample_visit <- function(sample, visit) {
	samp <- sample[sample$Visit!=visit, ] 
	return (samp)
}

# Create and return otu data.frames without input visit 
# using given out table
new_otu_visit <- function(otu, visit) {	
	elim_visit <- paste("*", visit, sep="")
	new_otu <- otu[,!grepl(elim_visit,names(otu))]
	return (new_otu)
}

new_sample_limb <- function(sample) {
	samp <- sample[sample$Limb!="UNK",]
	return (samp)
}

new_sample_trt <- function(sample, treatment) {
	samp <- sample[sample$Treatment==treatment,]
	return (samp)
}

#sample <- sample[!grepl("*\\*",sample$`Sample ID`),]

new_otu_limb <- function(otu, sample) {

	# Add samples from OTU table that are present in sample table
	new_otu <- otu[,grepl("OTU_ID",names(otu))]
	# new_otu <- otu[,1:6]
	for (j in sample$`Sample ID`) {
 		otu1 <- otu[,grepl(j,names(otu))]
 		new_otu <- bind_cols(new_otu,otu1)
	}

	return (new_otu) 
}


# Make .PDF for non-GUI interface 
makePDF <- function(name, item, analysis_type, visit) {
	pdf(paste(name,analysis_type, visit,".pdf", sep=""))
	plot(item)
	dev.off()
}


