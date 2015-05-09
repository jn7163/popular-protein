
# Loading the gene2pubmed file
setwd("~/Documents/Ping Lab/Manuscripts/2015 JACC Popular Proteins/Data summary")
gene2pubmed = read.table("20140811_gene2pubmed",header=1,sep = "\t");
summary(gene2pubmed)

# Subsetting the gene2pubmed file by taxonomy (9606 = homo sapiens, 10090 = mus musculus)
taxonomy= c("9606")
gene2pubmed_tax = gene2pubmed[ which(gene2pubmed$tax_id==taxonomy),]
dim(gene2pubmed_tax)
rm(gene2pubmed)
head(gene2pubmed_tax)
str(gene2pubmed_tax)
# Loading the BD2KPubMed output file
# NOTE, please manually remove the gene name and protein name columns from the frequencies file prior to running.
tissue_list = read.table("heart_9606_gene_frequencies.txt",header=1,fill=1,sep = "\t");
summary(tissue_list)
head(tissue_list)

# Tallying the total number of gene-linked and heart-related publications (checking against sum of tissue_lists$Frequency)
# I couldn't figure out how to read the .csv file without error, so just download the PMID list from the Pubmed search for now.
#tissue_pubs = read.table("heart or cardiac PMID.txt", header=0)
#tally <- gene2pubmed_tax$PubMed_ID  %in% tissue_pubs[,1]
#total_tissue_linked_pub_count <- sum(tally == TRUE)

# Tallying the total number of gene-linked publications in the organism
total_linked_pub_count <- length(gene2pubmed_tax$PubMed_ID)

# Output: how many times have each gene been linked in all of pubmed?
output <- paste("GeneID", "TotalPubCount", "TissuePubCount", "Semantic Distance", "Uniprot","GN","PN", sep = "\t")
write(output, file="Similarity_output.txt")

# Annotation file for Uniprot, GN, and PN
annot = read.csv(file = "9606_annotations.csv");
dim(annot)
names(annot)
Gene2Uni = match(tissue_list$Gene, annot$GeneID)
annotUni = annot$Uniprot[Gene2Uni]
annotGN = annot$GN[Gene2Uni]
annotPN = annot$PN[Gene2Uni]
# The following is the number or probes without annotation:
sum(is.na(Gene2Uni))
# Should return 0.


for (c in 1:nrow(tissue_list))
{
  print(paste("Now running Gene ", c, " of ", nrow(tissue_list), ". ", round(c/nrow(tissue_list)*100,2), "% done."))
  gene = tissue_list$Gene[c]
  gene_count = (gene2pubmed_tax$GeneID==gene)
  linked_pub = gene2pubmed_tax$PubMed_ID[gene_count]
  linked_pub_count_all_tissues = length(linked_pub)
  # normalized distance
  similarity_numerator = max(log10(sum(tissue_list$Frequency)),log10(linked_pub_count_all_tissues)) - log10(tissue_list$Frequency[c])
  similarity_denominator = log10(total_linked_pub_count) - min(log10(sum(tissue_list$Frequency)),log10(linked_pub_count_all_tissues))
  similarity = similarity_numerator/similarity_denominator
  output <- paste(gene, linked_pub_count_all_tissues, tissue_list$Frequency[c], round(similarity,3), annotUni[c], annotGN[c], annotPN[c], sep = "\t")
  write(output, file="Similarity_output.txt", append=T)
}



### Similarity between two genes

