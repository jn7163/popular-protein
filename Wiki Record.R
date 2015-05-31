### GETTING WIKI ARTICLE RECORDS

# Home Directory
home_directory <- "~/Documents/Ping Lab/Project Files/2015 Gene Wiki"
gene_wiki_record_location <- "http://genewikiplus.org/data/genewiki-latest.txt"
gene_wiki_local_file <- "GeneWiki Record/genewiki-latest.txt"


counter_location <- "http://stats.grok.se/json/en/"
date <- "201204"

output_location <- paste("Wiki_Count",date,".txt",sep="")

### LOADING PACKAGES ###
require(dplyr)
require(jsonlite)

#########################

setwd(home_directory)

download.file(gene_wiki_record_location, destfile = gene_wiki_local_file)
genewiki_record <-read.table(gene_wiki_local_file, skip = 1, sep="\t", as.is=TRUE, fill = TRUE, header=TRUE)
genewiki_list <- distinct(select(genewiki_record, Gene.Name, Entrez.Id))

df <- data.frame(Gene.Name=character(),Entrez.Id=character(),Count=character(),stringsAsFactors=FALSE)
        
for (i in 1:nrow(genewiki_list)){
        print(paste("Now running Gene ", i, " of ", nrow(genewiki_list), ". ", round(i/nrow(genewiki_list)*100,2), "% done."))
        
        gene <- genewiki_list[i,1]
        gene <- gsub(" ","%20",gene)
        
        json_url <- paste(counter_location,date,"/",gene, sep="")
        jsonData <- fromJSON(json_url)
        jsonData2 <- as.data.frame(unlist(jsonData[[1]]))
        row <- cbind(genewiki_list[i,1],genewiki_list[i,2],sum(jsonData2))
        df <- rbind(df,row)
}

write.table(df,output_location,row.names=FALSE,sep="\t",append=TRUE,quote=FALSE)