#1#Plot a bar graph showing the absolute count of the 50 most abundant bacteria, grouped by time (day after weaning):

library(ggplot2)
library(dplyr)
library(vegan)
library(ShortRead)
library(edgeR)
library(Biostrings)
library(DESeq2)

#for loop to process and combine the data from multiple fastq files:

# Set the directory where the fastq files are located
dir <- "Path/to/files"

# Get a list of all the fastq files in the directory
files <- list.files(path = dir, pattern = "*.fastq")

# Create an empty dataframe to store the combined data
all_data <- data.frame()

# Loop through the files
for (file in files) {
  # Import the fastq file
  data <- import_function(paste0(dir, file))
  
  # Process the data to obtain the counts of the bacterial species
  counts <- process_function(data)
  
  # Append the counts data to the all_data dataframe
  all_data <- rbind(all_data, counts)
}

# Set the directory where the fastq files are located
dir <- "/home/amrgalal/Desktop/Teste_Amr/output"

# Get a list of all the fastq files in the directory
files <- list.files(path = dir, pattern = "*.fastq")

# Extract the information from the fastq files
sequences_data <- lapply(files, function(file){
    # Extract the time information from the file name
    time <- strsplit(basename(file), "_")[[1]][1]
    # Read the fastq file
    reads <- readFastq(file)
    data.frame(bacteria_name = names(reads), count = unlist(reads), time = time)
})
#Convert the list to a dataframe
sequencing_data <- do.call(rbind, sequences_data)

#Filter data to only include the top 50 bacteria by count
top50_data <- sequencing_data %>%
group_by(bacteria_name, time) %>%
summarize(total_count = sum(count)) %>%
ungroup() %>%
group_by(time) %>%
top_n(50, total_count)

#Create the bar graph
ggplot(top50_data, aes(x = bacteria_name, y = total_count, fill = time)) +
geom_bar(stat = "identity") +
xlab("Bacteria Name") +
ylab("Count") +
ggtitle("Absolute count of the 50 most abundant bacteria by day after weaning") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))


#2# Plot a PCoA plot showing the clustering profile among the samples by day after weaning:

# Assume your data is in a presence/absence matrix called "sequencing_data_matrix"
# and a vector called "time" indicating the day after weaning for each sample

# Perform the ordination
ordination_result <- metaMDS(sequencing_data_matrix, distance = "bray")

# Create the PCoA plot
plot(ordination_result, type = "n")
points(ordination_result, display = "sites", col = time)
legend("topright", legend = unique(time), col = unique(time), pch = 21)
title("PCoA plot of samples by day after weaning")


#3# Use edgeR to identify differentially abundant bacteria between weaning days:

# Assume your data is in a dataframe called "sequencing_data" with columns "bacteria_name", "time", and "count"

# Create the DGEList object
y <- DGEList(counts = sequencing_data[,c("bacteria_name","count")], group = sequencing_data$time)

# Estimate common dispersion
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

# Identify differentially expressed genes
de_genes <- topTags(y, n = Inf, adjust.method = "BH", sort.by = "logFC")

# Extract the results
de_genes_data <- data.frame(de_genes$table[,c("logFC","logCPM","PValue","FDR")])


#4# Use DESeq2 to identify differentially abundant bacteria between weaning days:

# Assume your data is in a dataframe called "sequencing_data" with columns "bacteria_name", "time", and "count"

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = sequencing_data[,c("bacteria_name","count")], colData = data
frame(condition = sequencing_data$time), design = ~ condition)

#Run the differential expression analysis
dds <- DESeq(dds)

#Extract the results
res <- results(dds)

#Extract the results of interest, such as log2 fold change, p-value and adjusted p-value
res_data <- data.frame(res[,c("log2FoldChange","pvalue","padj")])

#Filter the results by p-value and adjusted p-value
res_data_filtered <- res_data[which(res_data$pvalue < 0.05 & res_data$padj < 0.05),]



