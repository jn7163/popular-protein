# Popular Protein Tally Full Program v.1.0

####################### USER INPUT ######################

# Home Directory
home_directory <- "~/Documents/Ping Lab/Project Files/2015 Popular Proteins"

# Gene2PubMed File Location
gene2pubmed_location <- "Gene2PubMed file/20140811_gene2pubmed"

# PMID list location
pmid_location <- "PMID list/hliver or hepatic_pmid.txt"

# Species
taxonomy= c("10090")       #9606 for human, 10090 for mouse

# Annotation File location
annotation_location <- "~/Documents/Ping Lab/Project Files/2015 Popular Proteins/Annotation file/10090_annotations.csv"
##########################################################

### ALTERNATIVE: Combing all PMID files in the folder to get proteins for all orgams
#search_term_1 = read.table("PMID list/gut or intestine or intestinal_pmid.txt",header=0);
#search_term_2 = read.table("PMID list/heart or cardiac_pmid.txt",header=0);
#search_term_3 = read.table("PMID list/brain or cerebral_pmid.txt",header=0);
#search_term_4 = read.table("PMID list/kidney or kidneys or renal_pmid.txt",header=0);
#search_term_5 = read.table("PMID list/liver or hepatic_pmid.txt",header=0);
#search_term_6 = read.table("PMID list/lung or lungs or pulmonary_pmid.txt",header=0);
#search_term = rbind(search_term_1,search_term_2,search_term_3,search_term_4,search_term_5,search_term_6)
#search_term = unique(search_term)

###########################################################

# Loading the gene2pubmed file
setwd(home_directory)
gene2pubmed = read.table(gene2pubmed_location,header=1,sep = "\t");
summary(gene2pubmed)

# Subsetting the gene2pubmed file by taxonomy (9606 = homo sapiens, 10090 = mus musculus)
gene2pubmed_tax = gene2pubmed[ which(gene2pubmed$tax_id==taxonomy),]
dim(gene2pubmed_tax)
rm(gene2pubmed)
head(gene2pubmed_tax)
str(gene2pubmed_tax)

## Further subsetting the Gene2Pubmed file and retaining only the PMIDs that appear in the Keyword searches
search_term = read.table(pmid_location,header=0);
tally <- gene2pubmed_tax$PubMed_ID  %in% search_term[,1]
gene2pubmed_tax_term <- gene2pubmed_tax[tally,]

## Count the frequency of each gene in the subsetted Gene2Pubmed file
gene_count <- table(gene2pubmed_tax_term$GeneID) # Tabulate and count frequency
gene_count_table <- as.data.frame(gene_count) #Turn into data frame
sorted_gene_count_table <- gene_count_table[order(-gene_count_table[,2]),]
colnames(sorted_gene_count_table) <- c("gene","count")
rownames(sorted_gene_count_table) <- NULL

# Option to write out the first output (just frequency, no semantic distance or annotation)
#write.table(sorted_gene_count_table, file="First_output.txt",quote=FALSE, row.names=FALSE)

# Tallying the total number of gene-linked publications in the organism
total_linked_pub_count <- length(unique(gene2pubmed_tax$PubMed_ID))

# Tallying the total number of gene-linked publications in the organism 
total_linked_pub_count_term <- length(unique(gene2pubmed_tax_term$PubMed_ID))

# Output: how many times have each gene been linked in all of pubmed?
output <- paste("GeneID", "TotalPubCount", "TissuePubCount", "Semantic Distance", "Uniprot","GN","PN", sep = "\t")
write(output, file="Similarity_output.txt")

# Annotation file for Uniprot, GN, and PN
annot = read.csv(file = annotation_location);
dim(annot)
names(annot)
Gene2Uni = match(sorted_gene_count_table$gene, annot$GeneID)
annotUni = annot$Uniprot[Gene2Uni]
annotGN = annot$GN[Gene2Uni]
annotPN = annot$PN[Gene2Uni]
# The following is the number or probes without annotation:
sum(is.na(Gene2Uni))
# Should return 0.

# Loop through all genes in the first output table (for the particular search term)
for (c in 1:nrow(sorted_gene_count_table))
{
  print(paste("Now running Gene ", c, " of ", nrow(sorted_gene_count_table), ". ", round(c/nrow(sorted_gene_count_table)*100,2), "% done."))
  gene = sorted_gene_count_table$gene[c]
  # Getting the total number of papers for the gene (any search term)
  gene_count = (gene2pubmed_tax$GeneID==gene)
  
  # This is the total number of linked publication counts for the gene, in any search term
  linked_pub_count_all_tissues = length(gene2pubmed_tax$PubMed_ID[gene_count])
  # normalized distance between gene and search term (similarity should really read distance)
  similarity_numerator = max(log10(total_linked_pub_count_term),log10(linked_pub_count_all_tissues)) - log10(sorted_gene_count_table$count[c])
  similarity_denominator = log10(total_linked_pub_count) - min(log10(total_linked_pub_count_term),log10(linked_pub_count_all_tissues))
  output <- paste(gene, linked_pub_count_all_tissues, sorted_gene_count_table$count[c], round(similarity_numerator/similarity_denominator,3), annotUni[c], annotGN[c], annotPN[c], sep = "\t")
  write(output, file="Similarity_output.txt", append=T)
}



### Similarity matrix between any two genes; because I don't know how to do this efficiently, limit yourself to only the top handful of genes with most publications
depth <- 75

#output <- paste("GeneID", "TotalPubCount", "TissuePubCount", "Semantic Distance", "Uniprot","GN","PN", sep = "\t")
#write(output, file="Two_gene_cross_similarity_output.txt")

matching_matrix <- matrix(,depth,depth)
rownames(matching_matrix) <- 1:depth
colnames(matching_matrix) <- 1:depth

for (d in 1:depth)
{
  gene_1 = sorted_gene_count_table$gene[d]
  gene_count_1 = (gene2pubmed_tax_term$GeneID==gene_1) ####### NOTE THIS IS THE TOPIC SPECIFIC COUNT. USE gene2pubmed_tax for ALL TOPIC CLUSTERING
  linked_pub_1 = gene2pubmed_tax_term$PubMed_ID[gene_count_1] ##_term
  rownames(matching_matrix)[d] <- as.character(annotGN[d])
  for (e in 1:depth)
  {
    print(paste("Now running Comparison ", (d-1)*depth+e, " of ", depth^2, ". ", round((((d-1)*depth+e)/depth^2)*100,2), "% done."))
    gene_2 = sorted_gene_count_table$gene[e]
    gene_count_2 = (gene2pubmed_tax_term$GeneID==gene_2) #_term
    linked_pub_2 = gene2pubmed_tax_term$PubMed_ID[gene_count_2] #_term
    matching = match(linked_pub_2, linked_pub_1)
    matching_count = length(matching[!is.na(matching)])
    
    # normalized distance between two genes
    if (matching_count == 0){
      distance <- 1.1
    }
    else{
    distance_numerator = max(log10(length(linked_pub_1)),log10(length(linked_pub_2))) - log10(matching_count)
    distance_denominator = log10(total_linked_pub_count_term) - min(log10(length(linked_pub_1)),log10(length(linked_pub_2)))
    distance = distance_numerator/distance_denominator
    }
    colnames(matching_matrix)[e] <- as.character(annotGN[e])
    matching_matrix[d,e] <- 1-distance

   }
 
}
require("heatmap3")
require("RColorBrewer")
heatmap3(matching_matrix,col=brewer.pal(9,"YlOrRd"))

# Printing out the clustered genes in order
map <- heatmap3(matching_matrix,col=brewer.pal(9,"YlOrRd"))
for(i in 1:length(map$colInd)){
        print(rownames(matching_matrix)[map$colInd[i]])
}
