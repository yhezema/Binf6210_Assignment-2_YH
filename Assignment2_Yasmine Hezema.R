#*****************************************
#              Assignment_2
#
#             Yasmine Hezema
#
#*****************************************

## -1-Packages used ---------
# Load necessary packages
library(rentrez) # For accessing NCBI databases.
library(stringr) # For string manipulation.
library(Biostrings) # For handling biological sequences.
library(ggplot2) # For plotting.
library(dplyr) # For data manipulation.
library(ape) # For sequence alignment.
library(phytools) # For phylogenetic analysis.
library(phangorn) # For additional phylogenetic tools.
library(muscle) # # For multiple sequence alignment.
library(styler) # For formats code for a consistent style


## 2-Code part1 - Data exploration and manipulation ------
# Search NCBI nucleotide database for "Cecidomyiidae"
#search_result_Cecidomyiidae <- entrez_search(db = "nucleotide", term = "Cecidomyiidae", retmax = 75)

# Check the class
#class(search_result_Cecidomyiidae) # It is an "esearch" object and a list.

# Check what it contains
#search_result_Cecidomyiidae # Entrez search result with 43361092 hits (object contains 75 IDs and no web_history object)

# Search for COI gene in 	Cecidomyiidae, sequences between 600-1000 bp
#COI_search <- entrez_search(db = "nucleotide", term = "Cecidomyiidae[ORGN] AND COI[Gene] AND 300:1000[SLEN]", retmax = 75)
#COI_search # Entrez search result with 156891 hits (object contains 75 IDs)

# Search for CAD gene in 	Cecidomyiidae , sequences between 600-1000 bp
#CAD_search <- entrez_search(db = "nucleotide", term = "Cecidomyiidae[ORGN] AND CAD[Gene] AND 300:1000[SLEN]", retmax = 75)
#CAD_search # Entrez search result with 237 hits (object contains 75 IDs)

# Check the class of the ID
#class(CAD_search$ids) # "character"
#class(COI_search$ids) # "character"
# Display IDs of COI and CAD
#COI_search$ids
#CAD_search$ids
# Search result count
#COI_search$count # count 156891 hits
#CAD_search$count # count = 237 hits
# check the search returned only 75 hits at it was set in the entrez_search
#length(COI_search$ids) # length = 75 as we set retmax to be 75
#length(CAD_search$ids) # length = 75 as we set retmax to be 75

# Get the data summaries about the nucleotide sequences using entrez_summary() to return metadata about the sequence records
## For COI
#COI_summ <- entrez_summary(db = "nucleotide", id = COI_search$ids)
#CAD_summ <- entrez_summary(db = "nucleotide", id = CAD_search$ids)

#COI_summ # List of  75 esummary records with 32 items
#CAD_summ ## List of  75 esummary records with 32 items
#class(COI_summ) # "esummary_list" "list"
#class(CAD_summ) # "esummary_list" "list"

# Extract and inspect organism information
#extract_from_esummary(COI_summ, "organism")
#extract_from_esummary(CAD_summ, "organism")
# will see at the console name of species at of each gene and their sequence id

# Web history for COI and CAD searches
#Cecidomyiidae_search_COI <- entrez_search(db = "nucleotide", term = "	Cecidomyiidae[ORGN] AND COI[Gene] AND 600:1000[SLEN]", retmax = 75, use_history = TRUE)
#Cecidomyiidae_search_CAD <- entrez_search(db = "nucleotide", term = "	Cecidomyiidae[ORGN] AND CAD[Gene] AND 600:1000[SLEN]", retmax = 75, use_history = TRUE)

# Fetch FASTA sequences for COI and CAD and save it in fasta files
#COI_fetch_Cecidomyiidae <- entrez_fetch(db = "nucleotide", id = COI_search$ids, rettype = "fasta")
#CAD_fetch_Cecidomyiidae <- entrez_fetch(db = "nucleotide", id = CAD_search$ids, rettype = "fasta")

# Write FASTA sequences to files at my work directory
#write(COI_fetch_Cecidomyiidae, "COI_fetch_Cecidomyiidae.fasta", sep = "\n")
#write(CAD_fetch_Cecidomyiidae, "CAD_fetch_Cecidomyiidae.fasta", sep = "\n")

# Read the DNA sequence back in as DNAStringSet from the FASTA file
COI_seqs <- readDNAStringSet("COI_fetch_Cecidomyiidae.fasta", format = "fasta")
CAD_seqs <- readDNAStringSet("CAD_fetch_Cecidomyiidae.fasta", format = "fasta")

class(COI_seqs)
class(CAD_seqs)
# it returned "DNAStringSet" attr(,"package"), "Biostrings", that means the COI_seqs and CAD_Seqs are a collection of DNA sequences managed by the Biostrings package
head(names(COI_seqs))
head(names(CAD_seqs))
# it shows the the title of head sequences of each gene

# Convert sequences in the fasta files into dataframes
dfCOI <- data.frame(COI_Title = names(COI_seqs), COI_Sequence = paste(COI_seqs))
dfCAD <- data.frame(CAD_Title = names(CAD_seqs), CAD_Sequence = paste(CAD_seqs))
# for each gene, from the objectives COI_seqs and CAD_seqs make a dataframe has a column with the original sequence titles called COI_Title or CAD_Title and column for the sequences called COI_Sequence or CAD_Sequence. So will get one dataframe for each genes with two columns original sequence title and the other is the sequence

# Create a new column contains the species name
dfCOI$Species_Name <- word(dfCOI$COI_Title, 2L, 3L)
dfCAD$Species_Name <- word(dfCAD$CAD_Title, 2L, 3L)
# Now each df has 75 rows and 3 columns in this order ("COI_Title", "COI_Sequence" and "Species_Name")

# Rearrange the columns
dfCOI <- dfCOI[, c("COI_Title", "Species_Name", "COI_Sequence")]
dfCAD <- dfCAD[, c("CAD_Title", "Species_Name", "CAD_Sequence")]

# view column names
colnames(dfCOI)
colnames(dfCAD)
# The columns order changed to be "COI_Title", "Species_Name", "COI_Sequence"

# Look at the dataframe after the modification
View(dfCOI)
View(dfCAD)

# Explore COI and CAD sequences:
# Explore the unique species
unique(dfCOI$Species_Name) # COI has 68 unique species shown at the console
unique(dfCAD$Species_Name) # CAD has 68 unique species shown at the console

# Check if the sequences have Na or "-"
sum(is.na(dfCOI$COI_Sequence)) # 0 NA was found
sum(is.na(dfCAD$CAD_Sequence)) # 0 NA was found
sum(str_count(dfCOI$COI_Sequence, "-")) # 0 NA was found
sum(str_count(dfCAD$CAD_Sequence, "-")) # 0 NA was found

# Explore the sequences lengths
summary(str_count(dfCOI$COI_Sequence))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#443.0   629.0   629.0   619.7   629.0   629.0
summary(str_count(dfCAD$CAD_Sequence))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#546.0   728.0   749.0   747.5   784.0   875.0 

# Calculate sequence lengths and add length column to the df
dfCOI$seq_length <- str_count(dfCOI$COI_Sequence)
dfCAD$seq_length <- str_count(dfCAD$CAD_Sequence)

# view column names
colnames(dfCOI)
colnames(dfCAD)
# Now each df contains "COI_Title", "Species_Name", "COI_Sequence" and "seq_length"

# Create a vector has sequence lengths
COI_lengths <- sapply(COI_seqs, length)
CAD_lengths <- sapply(CAD_seqs, length)
# This code calculates the lengths of each DNA sequence in the COI_seqs and CAD_seqs and store these lengths in the COI_lengths and CAD_lengths vectors, respectively for further data visualization.

# Plot histograms of sequence lengths
# Combine data into a single data frame
gene_data <- data.frame(
  Length = c(COI_lengths, CAD_lengths),
  Gene = rep(c("COI", "CAD"), each = length(COI_lengths))
)

# Plot histograms with each gene as a separate panel
ggplot(gene_data, aes(x = Length, fill = Gene)) +
  geom_histogram(bins = 20, color = "black", alpha = 0.7) +
  labs(title = "Sequence Lengths Frequency of CAD and COI", x = "Length (bp)", y = "Frequency") +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"), # Adjust font size and bold for facet labels
    plot.title = element_text(size = 20, face = "bold"), # Increase font size and bold for the title
    axis.title.x = element_text(size = 16), # Increase font size for x-axis label
    axis.title.y = element_text(size = 16) # Increase font size for y-axis label
  )

# Box plot for COI_lengths and CAD_lengths
boxplot(list(COI = COI_lengths, CAD = CAD_lengths),
        main = "COI and CAD Sequence Lengths",
        ylab = "Length (bp)",
        col = c("lightblue", "lightgreen")
)

# Filter rows where sequence length is between the 1st and 3rd quartile to remove the outliers
dfCOI_filtered <- dfCOI %>%
  filter(seq_length >= 629) # as the 3rd Qu. and the maximum = 629

dfCAD_filtered <- dfCAD %>%
  filter(seq_length >= 728 & seq_length <= 784)

# Box plot for filtered sequence lengths
# Extract filtered sequence lengths in a new vector
COI_filtered_lengths <- dfCOI_filtered$seq_length
CAD_filtered_lengths <- dfCAD_filtered$seq_length

# Box plot for filtered sequence lengths
boxplot(list(COI = COI_filtered_lengths, CAD = CAD_filtered_lengths),
        main = "Box Plot of Filtered Sequence Lengths",
        ylab = "Length (bp)",
        col = c("lightblue", "lightgreen")
)

# Combine non-filtered and filtered data for plotting in one graph
all_lengths <- list(
  COI = COI_lengths,
  CAD = CAD_lengths,
  COI_Filtered = COI_filtered_lengths,
  CAD_Filtered = CAD_filtered_lengths
)

# Plot boxplots
boxplot(all_lengths,
        main = "Box Plot of COI and CAD Sequence Lengths (Non-Filtered and Filtered)",
        ylab = "Length (bp)",
        col = c("lightblue", "lightgreen", "lightblue", "lightgreen"),
        names = c("Non-Filtered COI", "Non-Filtered CAD", "Filtered COI", "Filtered CAD"),
        cex.main = 1.5,  # Increase title size
        cex.lab = 1.6,   # Increase axis label size
        cex.axis = 1.4   # Increase axis tick label size
)

## 3-Code part2- Main Analysis-----
# Set a seed for reproducibility to select the same sequence every time run the code 
set.seed(123)

# Randomly sample one sequence per specie for COI and CAD for creating the tree
dfCOI_Subset <- dfCOI_filtered %>%
  group_by(Species_Name) %>%
  slice_sample(n = 1)

dfCAD_Subset <- dfCAD_filtered %>%
  group_by(Species_Name) %>%
  slice_sample(n = 1)

# explore the new subsets
dfCOI_Subset # 51 rows of unique species out of 75 before filtering
dfCAD_Subset # 41 rows of unique species out of 75 before filtering

# Merge COI and CAD dataframes to show only the species with sequences of the two genes
dfAllSeqs_subset <- merge(dfCOI_Subset, dfCAD_Subset, by = "Species_Name", all = FALSE)
View(dfAllSeqs_subset) # 31 species
colnames(dfAllSeqs_subset) # this df has 7 columns with names "Species_Name" "COI_Title" "COI_Sequence" "seq_length.x" "CAD_Title"  "CAD_Sequence" and  "seq_length.y" and 28 rows of the species

# Align the sequences using muscle
COI_seqs <- as.character(dfAllSeqs_subset$COI_Sequence)
CAD_seqs <- as.character(dfAllSeqs_subset$CAD_Sequence)

# Convert COI&CAD_seqs (character vector) to DNAStringSet
COI_seqs <- DNAStringSet(COI_seqs)
CAD_seqs <- DNAStringSet(CAD_seqs)

# Multiple sequence alignment
COI_aligned <- muscle::muscle(COI_seqs)
CAD_aligned <- muscle::muscle(CAD_seqs)
# The purpose of performing sequence alignment is to identify regions of similarity that may indicate evolutionary relationships among the sequences to construct the phylogenetic tree.

# Convert aligned sequences to DNAbin format for phylogenetic analysis
COI_aligned_bin <- as.DNAbin(COI_aligned)
CAD_aligned_bin <- as.DNAbin(CAD_aligned)
# Explore the DNAbin object
COI_aligned_bin
CAD_aligned_bin

# Save aligned sequences to files in FASTA format
write.dna(COI_aligned_bin, file = "COI_aligned.fasta", format = "fasta")
write.dna(CAD_aligned_bin, file = "CAD_aligned.fasta", format = "fasta")

# Construct phylogenetic trees using neighbor-joining on distance matrices
COI_tree <- nj(dist.dna(COI_aligned_bin))
CAD_tree <- nj(dist.dna(CAD_aligned_bin))

# Add species names to the tree
COI_tree$tip.label <- dfAllSeqs_subset$Species_Name
CAD_tree$tip.label <- dfAllSeqs_subset$Species_Name

# Check if both trees (COI_tree and CAD_tree) have matching tip labels
all(COI_tree$tip.label == CAD_tree$tip.label) # Returned TRUE which means both trees have the same species neames.

# Check if the trees are ultrametric
is.ultrametric(COI_tree) # Return FALSE
is.ultrametric(CAD_tree) # Return FALSE

# Modify the trees for consistency
# Set negative branch lengths to zero
COI_tree$edge.length[COI_tree$edge.length < 0] <- 0
CAD_tree$edge.length[CAD_tree$edge.length < 0] <- 0
# This code is to handle the negative branch lengths

# Force a trees to be ultrametric
if (!is.ultrametric(COI_tree)) COI_tree <- chronos(COI_tree)
if (!is.ultrametric(CAD_tree)) CAD_tree <- chronos(CAD_tree)
# checks if the trees are ultrametric. If they are not, it uses the chronos function from the ape package to convert them into ultrametric trees. This is important for this analyses as I assume a constant rate of evolution of the two genes.

# Recheck Ultrametric Status
is.ultrametric(COI_tree) # return TRUE
is.ultrametric(CAD_tree) # return TRUE

# Compare Phylogenetic Trees with Robinson-Foulds Distance
rf_distance <- RF.dist(COI_tree, CAD_tree)
print(paste("Robinson-Foulds distance between COI and CAD trees:", rf_distance))
# "Robinson-Foulds distance between COI and CAD trees: 26". The Robinson-Foulds (RF) distance between your COI and CAD trees is 26. This metric quantifies the difference between the two tree topologies. An RF distance of 30 indicates that there are 30 splits (or branches) that differ between the two trees.

# Define a function to assign colors to groups
assign_colors_phylo <- function(tree, n) {
  # Perform hierarchical clustering on the tree
  hc <- hclust(dist(cophenetic(tree)))
  # Cut the tree into n groups
  groups <- cutree(hc, k = n)
  # Assign a unique color to each group
  colors <- c("darkred", "darkblue", "darkgreen", "darkorange", "purple", "darkcyan")
  group_colors <- colors[groups]
  return(group_colors)
}

# Apply the function to the trees and plot them
n_groups <- 6 # Number of subfamilies at this family

# Assign colors to each tree
COI_colors <- assign_colors_phylo(COI_tree, n_groups)
CAD_colors <- assign_colors_phylo(CAD_tree, n_groups)

# Set up side-by-side plotting
par(mfrow = c(1, 2))

# Plot the COI and CAD tree with colors and longer branches
plot(COI_tree, main = "COI", cex = 0.9, tip.color = COI_colors, edge.width = 2, x.lim = c(0, 1.5))
plot(CAD_tree, main = "CAD", cex = 0.9, tip.color = CAD_colors, edge.width = 2, x.lim = c(0, 1.5))

# Add the main title centered above both plots with bold font
mtext("Phylogenetic Tree of COI and CAD", side = 3, outer = TRUE, line = -1, cex = 1.2, font = 2)
