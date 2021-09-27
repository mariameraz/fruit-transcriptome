#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                 DEG Analysis                                 #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

#Load libraries
library(tidyverse)
library(DESeq2)
library(ggpubr)
library(ggplot2)
library(ggvenn)

#Load count data
#Your path to the files here:
temp.path <- "./Count_data/"

#Load and merge the both count datasets in a single list
count_data <- list("pepper_count" = read.table(paste(temp.path, "pepper_count.txt", sep = ""), header=T),
                   "tomato_count" = read.table(paste(temp.path, "tomato_count.txt", sep = ""), header=T))

#Create metadata table
metadata <- data.frame(stage = rep(c("S1","S2","MG","B"), each=3),
                       sample = paste(rep(c("PRS1","PRS2","PRMG","PRB"), each=3), sep = "_") 
                       %>% paste(., rep(1:3, times=3), sep = "_")) 

rownames(metadata) <- metadata$sample
metadata$stage = factor(metadata$stage)

#::::::::::::::::::::::::::::::::::::::::::::#
#      DEG analysis for ALL the genes        #  
#::::::::::::::::::::::::::::::::::::::::::::#

names(count_data)
for (i in names(count_data)) {
  #DESEq analysis
      #DESEq analysis
  if (all(colnames(count_data[[1]]) == rownames(metadata)) == T) {
    ddseq <- DESeqDataSetFromMatrix(countData= count_data[[i]], colData= metadata, 
                                    design = ~ stage)
    keep = rowSums(counts(ddseq)) >= 5
    table(keep) 
    ddseq = ddseq[keep,]
    ddseq 
    
    dds <- DESeq(ddseq)
    
    res[[gsub("count","s1_s2",i)]] <- as.data.frame(lfcShrink(dds, contrast= c("stage","S1","S2")))
    res[[gsub("count","s2_mg",i)]] <- as.data.frame(lfcShrink(dds, contrast = c("stage", "S2","MG")))
    res[[gsub("count","mg_b",i)]] <- as.data.frame(lfcShrink(dds, contrast = c("stage","MG","B")))

    } else {
      stop()
      print("Metadata sample column and colnames of count data don't have the same order")
    }
    
}

res <- list("pepper_s1_s2" =pepper_s1_s2, "pepper_s2_mg"=pepper_s2_mg, "pepper_mg_b"=pepper_mg_b,
            "tomato_s1_s2" =tomato_s1_s2, "tomato_s2_mg"=tomato_s2_mg, "tomato_mg_b"=tomato_mg_b)

#How many significant (abs LFC > 1 and p adjusted < 0.1) DEG are in each contrast?
for (i in c("pepper","tomato")) {
  
  print(paste(paste("S1 vs S2 DEG in", i, sep = " "), nrow(res[[paste(i,"s1_s2", sep = "_")]] %>% 
                                                             dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))
  print(paste(paste("S2 vs MG DEG in", i, sep = " "), nrow(res[[paste(i,"s2_mg", sep = "_")]] %>%
                                                             dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))
  print(paste(paste("MG vs B DEG in", i, sep = " "), nrow(res[[paste(i,"mg_b", sep = "_")]] %>%
                                                            dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))

}

contras <- c("s1_s2","s2_mg","mg_b")
#Create MA plots
for (i in contras) {
  for (j in c("pepper","tomato")) {
    #Determine the colors
    res[[paste(j,i,sep = "_")]]<-
      res[[paste(j,i,sep = "_")]] %>% 
      na.omit() %>% 
      dplyr::mutate(color=ifelse(log2FoldChange > 1 & padj < 0.01, "red", 
                                 ifelse(log2FoldChange < -1 & padj < 0.01,"blue","black")))
    #Plot MA
    assign(
      paste(paste(j,i,sep = "_"),"ma", sep = "_"),  
      res[[paste(j,i,sep = "_")]] %>% 
        ggplot(aes(log2(baseMean), log2FoldChange, color=as.factor(color))) +
        geom_point(size=1.5) +
        theme(legend.position = "none") +
        labs(x="Log2 Expresión", y= "LFC") +
        ggtitle(str_to_upper(gsub("_"," vs ", i))) +
        scale_color_manual(name = "thershold", values = c("black" ="#929ea1","blue"= "#308285","red"= "#c4922b" ))
    )
  }
}

#Merge and plot
ggarrange(pepper_s1_s2_ma, pepper_s2_mg_ma, pepper_mg_b_ma,
          tomato_s1_s2_ma, tomato_s2_mg_ma, tomato_mg_b_ma,
          ncol = 3, nrow = 2, labels = c("A","B","C","D","F"))

#Intersection of DEG 
res_sig <- list()
for (i in contras) {
  for (j in c("pepper","tomato")) {
   res_sig[[paste(paste(j, i, sep = "_"), "sig", sep = "_")]] <- res[[paste(j,i, sep = "_")]] %>% 
            dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01))
     
  }
}

#Create logical table

for (i in c("pepper", "tomato")){
  
assign(
  paste("venn_ortho",i, sep="_"),
  tibble(count_data[[paste(i,"count", sep="_")]] %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(paste(i,"s1_s2", sep="_") = ortho %in% rownames(res_sig[[paste(paste(i, "s1_s2", sep = "_"), "sig", sep = "_")]]) %>%
  dplyr::mutate(paste(i,"s2_mg", sep="_")= ortho %in% rownames(res_sig[[paste(paste(i, "s2_mg", sep = "_"), "sig", sep = "_")]]) %>%
  dplyr::mutate(paste(i,"mg_b", sep="_")= ortho %in% rownames(res_sig[[paste(paste(i, "mg_b", sep = "_"), "sig", sep = "_")]])) 
  )

pepper_venn <- ggplot(venn_ortho_pepper) +
  geom_venn(aes(A = pepper_s1_s2, B = pepper_s2_mg, C= pepper_mg_b), 
            fill_color = c("#308285", "#c4922b", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 4, 
            text_size = 4, 
            text_color = "#1a1918", 
            set_names = c("S1 vs S2","S2 vs MG", "MG vs B")) +
  coord_fixed(xlim = c(-1.6,1.6), ylim = c(-1.8,1.8)) +
  theme_void() + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.title=element_text(size=15,hjust=0.5, vjust = -5)) 

tomato_venn <- ggplot(venn_ortho_tomato) +
  geom_venn(aes(A = tomato_s1_s2, B = tomato_s2_mg, C= tomato_mg_b), 
            fill_color = c("#308285", "#c4922b", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 4, 
            text_size = 4, 
            text_color = "#1a1918", 
            set_names = c("S1 vs S2","S2 vs MG", "MG vs B")) +
  coord_fixed(xlim = c(-1.6,1.6), ylim = c(-1.8,1.8)) +
  theme_void() +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.title=element_text(size=15,hjust=0.5, vjust = -5)) 

ggarrange(pepper_venn, tomato_venn)

#::::::::::::::::::::::::::::::::::::::::#
#    DEG analysis for 1:1 orthologous    #
#::::::::::::::::::::::::::::::::::::::::#

#Load count data

pepper_count_ort <- read.table(paste(temp.path, "pepper_count_ort.txt", sep = ""), header=T) 
tomato_count_ort <- read.table(paste(temp.path, "tomato_count_ort.txt", sep = ""), header=T)


#::::::::::::::::::::::::::::::::::::::::::::#
#      DEG analysis for ALL the genes        #  
#::::::::::::::::::::::::::::::::::::::::::::#

#Merge the both count datasets in a single list
count_data <- list("pepper_count" = pepper_count_ort, "tomato_count" = tomato_count_ort)
names(count_data)

res<- list()
for (i in names(count_data)) {
      #DESEq analysis
  if (all(colnames(count_data[[1]]) == rownames(metadata)) == T) {
  ddseq <- DESeqDataSetFromMatrix(countData= count_data[[i]], colData= metadata, 
                                  design = ~ stage)
  keep = rowSums(counts(ddseq)) >= 5
  table(keep) 
  ddseq = ddseq[keep,]
  ddseq 
  
  dds <- DESeq(ddseq)
  
  res[[gsub("count","s1_s2_ort",i)]] <- as.data.frame(lfcShrink(dds, contrast= c("stage","S1","S2")))
  res[[gsub("count","s2_mg_ort",i)]] <- as.data.frame(lfcShrink(dds, contrast = c("stage", "S2","MG")))
  res[[gsub("count","mg_b_ort",i)]] <- as.data.frame(lfcShrink(dds, contrast = c("stage","MG","B")))
  
    } else {
      stop()
      print("Metadata sample columns and colnames of count data don't have the same order")
    }
}


#How many significant (abs LFC > 1 and p adjusted < 0.1) DEG are in each contrast?
for (i in c("pepper","tomato")) {
  
  print(paste(paste("S1 vs S2 DEG in", i, sep = " "), nrow(res[[paste(i,"s1_s2", sep = "_")]] %>% 
                                                             dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))
  print(paste(paste("S2 vs MG DEG in", i, sep = " "), nrow(res[[paste(i,"s2_mg", sep = "_")]] %>%
                                                             dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))
  print(paste(paste("MG vs B DEG in", i, sep = " "), nrow(res[[paste(i,"mg_b", sep = "_")]] %>%
                                                            dplyr::filter(abs(log2FoldChange) > 1 & padj< 0.01)), sep = ": "))
  
}

contras <- c("s1_s2","s2_mg","mg_b")
#Create MA plots
for (i in contras) {
  for (j in c("pepper","tomato")) {
    #Determine the colors
    res[[paste(j,i,sep = "_")]]<-
      res[[paste(j,i,sep = "_")]] %>% 
      na.omit() %>% 
      dplyr::mutate(color=ifelse(log2FoldChange > 1 & padj < 0.01, "red", 
                                 ifelse(log2FoldChange < -1 & padj < 0.01,"blue","black")))
    #Plot MA
    assign(
      paste(paste(j,i,sep = "_"),"ma_ort", sep = "_"),  
      res[[paste(j,i,sep = "_")]] %>% 
        ggplot(aes(log2(baseMean), log2FoldChange, color=as.factor(color))) +
        geom_point(size=1.5) +
        theme(legend.position = "none") +
        labs(x="Log2 Expresión", y= "LFC") +
        ggtitle(str_to_upper(gsub("_"," vs ", i))) +
        scale_color_manual(name = "thershold", values = c("black" ="#929ea1","blue"= "#308285","red"= "#c4922b" ))
    )
  }
}

#Merge and plot
ggarrange(pepper_s1_s2_ma_ort, pepper_s2_mg_ma_ort, pepper_mg_b_ma_ort,
          tomato_s1_s2_ma_ort, tomato_s2_mg_ma_ort, tomato_mg_b_ma_ort,
          ncol = 3, nrow = 2, labels = c("A","B","C","D","F"))

#Intersection of DEG 
for (i in contras) {
  for (j in c("pepper","tomato")) {
    assign(paste(paste(paste(j, i, sep = "_"), "ort", sep = "_"), "sig", sep = "_"), res[[paste(j,i, sep = "_")]] %>% 
             dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01))
    
  }
}

res_sig <- list("pepper_s1_s2"=pepper_s1_s2_ort_sig, "pepper_s2_mg"=pepper_s2_mg_ort_sig, "pepper_mg_b"=pepper_mg_b_ort_sig,
                "tomato_s1_s2"=tomato_s1_s2_ort_sig, "tomato_s2_mg"=tomato_s2_mg_ort_sig, "tomato_mg_b"=tomato_mg_b_ort_sig)

#Create logical table
venn_ortho <- tibble(tomato_count_ort %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(tomato_s1_s2= ortho %in% rownames(tomato_s1_s2_ort_sig)) %>%
  dplyr::mutate(tomato_s2_mg= ortho %in% rownames(tomato_s2_mg_ort_sig)) %>%
  dplyr::mutate(tomato_mg_b= ortho %in% rownames(tomato_mg_b_ort_sig)) %>%
  dplyr::mutate(pepper_s1_s2= ortho %in% rownames(pepper_s1_s2_ort_sig)) %>%
  dplyr::mutate(pepper_s2_mg= ortho %in% rownames(pepper_s2_mg_ort_sig)) %>%
  dplyr::mutate(pepper_mg_b= ortho %in% rownames(pepper_mg_b_ort_sig)) %>%
  dplyr::mutate(pepper= ortho %in% c(rownames(pepper_s1_s2_ort_sig),
                                     rownames(pepper_s2_mg_ort_sig),
                                     rownames(pepper_mg_b_ort_sig))) %>%
  dplyr::mutate(tomato=ortho %in% c(rownames(tomato_s1_s2_ort_sig),
                                    rownames(tomato_s2_mg_ort_sig),
                                    rownames(tomato_mg_b_ort_sig)))


pepper_venn <- ggplot(venn_ortho) +
  geom_venn(aes(A = pepper_s1_s2, B = pepper_s2_mg, C= pepper_mg_b), 
            fill_color = c("#308285", "#c4922b", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 4, 
            text_size = 4, 
            text_color = "#1a1918", 
            set_names = c("S1 vs S2","S2 vs MG", "MG vs B")) +
  coord_fixed(xlim = c(-1.6,1.6), ylim = c(-1.8,1.8)) +
  theme_void() +
  labs(tag = "A") + 
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.title=element_text(size=15,hjust=0.5, vjust = -5)) 

tomato_venn <- ggplot(venn_ortho) +
  geom_venn(aes(A = tomato_s1_s2, B = tomato_s2_mg, C= tomato_mg_b), 
            fill_color = c("#308285", "#c4922b", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 4, 
            text_size = 4, 
            text_color = "#1a1918", 
            set_names = c("S1 vs S2","S2 vs MG", "MG vs B")) +
  coord_fixed(xlim = c(-1.6,1.6), ylim = c(-1.8,1.8)) +
  theme_void() +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.title=element_text(size=15,hjust=0.5, vjust = -5)) 

ggarrange(pepper_venn, tomato_venn)


s1_s2 <- ggplot(venn_ortho) +
  geom_venn(aes(A = pepper_s1_s2, B = tomato_s1_s2), 
            fill_color = c("#308285", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 5.5, 
            text_size = 5.5, 
            text_color = "#1a1918", 
            set_names = c("Chile","Tomate")) +
  coord_fixed(xlim = c(-1.57,1.55), ylim = c(-0.9,1.8)) +
  theme_void() +
  labs(title = "S1 vs S2", tag = "A") +
  theme(plot.tag = element_text(size = 15, face = "bold", vjust = -7, hjust = 1),
        plot.title=element_text(size=16,hjust=0.5, vjust = -10)) 

s2_mg <- ggplot(venn_ortho) +
  geom_venn(aes(A = pepper_s2_mg, B = tomato_s2_mg), 
            fill_color = c("#308285", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 5.5, 
            text_size = 5.5, 
            text_color = "#1a1918", 
            set_names = c("Chile","Tomate")) +
  coord_fixed(xlim = c(-1.57,1.55), ylim = c(-0.9,1.8)) +
  theme_void() +
  labs(title = "S2 vs MG", tag="B") +
  theme(plot.tag = element_text(size = 15, face = "bold", vjust = -7, hjust = 1),
        plot.title=element_text(size=16,hjust=0.5, vjust = -10)) 

mg_b <- ggplot(venn_ortho) +
  geom_venn(aes(A = pepper_mg_b, B = tomato_mg_b), 
            fill_color = c("#308285", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 5.5, 
            text_size = 5.5, 
            text_color = "#1a1918", 
            set_names = c("Chile","Tomate")) +
  coord_fixed(xlim = c(-1.57,1.55), ylim = c(-0.9,1.8)) +
  theme_void() +
  labs(title = "MG vs B", tag = "C") +
  theme(plot.tag = element_text(size = 15, face = "bold", vjust = -7, hjust = 1),
        plot.title=element_text(size=16,hjust=0.5, vjust = -10)) 

all <- ggplot(venn_ortho) +
  geom_venn(aes(A = pepper, B = tomato), 
            fill_color = c("#308285", "#a32314"), 
            stroke_color = "#302e2e", 
            stroke_size = 0.5,
            fill_alpha = 0.6, 
            set_name_size = 5.5, 
            text_size = 5.5, 
            text_color = "#1a1918", 
            set_names = c("Chile","Tomate")) +
  coord_fixed(xlim = c(-1.57,1.55), ylim = c(-0.9,1.8)) +
  labs(title = "Total", tag = "D") +
  theme_void() +
  theme(plot.tag = element_text(size = 15, face = "bold", vjust = -7, hjust = 1),
        plot.title=element_text(size=16,hjust=0.5, vjust = -10)) 

ggarrange(s1_s2, s2_mg, mg_b, all, ncol = 2, nrow = 2)


##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                         GO terms enrichment analysis                         #
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#Load libraries
library(biomaRt)
library(topGO)

#Gene annotation (Uniprot ID)
annotation <- read.table("./Annotation_data/annotation.txt")

#Single copy orthologs 
single_ortho <- read.table("./GOTerms_results/Single_copy_orthologs.txt") %>%
  dplyr::rename(tomato.id=lycopersicum.faa)

uniprot_ort <- merge(annotation, single_ortho, by = "tomato.id")
uniprot_ort_back <- uniprot_ort

uniprot_ort<- uniprot_ort %>% dplyr::select(ortho, uniprot.id)
head(uniprot_ort)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Working with genes DE in at least one stage shared in the two species.                      #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
shared_ort <- venn_ortho[venn_ortho$pepper == T & venn_ortho$tomato == T,]

list <- uniprot_ort %>% mutate(filt=ifelse(ortho %in% 
                                             shared_ort$ortho, "yes","no")) %>% filter(filt=="yes") %>%
  separate(uniprot.id, into = c('uniID', 'uniprot_ID'), sep = "\\|") %>%
  dplyr::select(ortho, uniprot_ID)

list$uniprot_ID<- gsub("_.*","",list$uniprot_ID)
list_back <- list

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
best_hits <-annotation %>%
  separate(uniprot.id, into = c(
    'uniID', 'uniprot_ID'), sep = "\\|") %>%
  dplyr::select(uniprot_ID)
best_hits=as.data.frame(best_hits)

GenesTotales<-gsub("_.*","",best_hits$uniprot_ID)
listMarts(host="plants.ensembl.org")
ensembl_plant <- useMart(host="https://plants.ensembl.org", 
                         biomart="plants_mart", 
                         port = 443,)

listDatasets(ensembl_plant)
db= useMart('plants_mart',dataset='athaliana_eg_gene', host="http://plants.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=GenesTotales,
              mart=db, useCache = FALSE)

#Build gene 2 GO annotation list
gene_2_GO=unstack(go_ids[,c(1,2)])

#Remove any candidate genes without GO annotation.
keep = list$uniprot_ID %in% go_ids[,2]
keep = which(keep==TRUE)
list=list$uniprot_ID[keep]
length(list) #426


list=factor(as.integer(GenesTotales %in% list))
names(list)= GenesTotales

GOdata = new('topGOdata', 
             ontology='BP', #Character string specifying the ontology of interest (BP, MF or CC)
             allGenes = list, #Named vector of type numeric or factor. Contains the genes identifiers. The genes listed in this object define the gene universe.
             annot = annFUN.gene2GO, #This function is used when the annotation are provided as a gene-to-GOs mapping.
             gene2GO = gene_2_GO)

GO_fisherW=runTest(GOdata, algorithm='weight01', statistic='fisher')

GO=usedGO(GOdata)
GO_res=GenTable(GOdata, weightFisher=GO_fisherW,
                orderBy='weightFisher', topNodes=length(GO)) 

#Top 30 of the most significant (padj > 0.01) enriched TopGO
GO_res_top <- GO_res %>%
  dplyr::filter(!stringr::str_detect(Term, "biological")) %>%
  dplyr::filter(as.numeric(weightFisher) < 0.01) %>%
  arrange(as.numeric(weightFisher)) %>%
  head(30) %>%
  dplyr::rename(Genes=Significant, p.adj =weightFisher) %>%
  dplyr::mutate(enrichment.factor = as.numeric(Genes)/as.numeric(Expected))

#Plot data
GO_plot <-     ggplot(GO_res_top) + 
  geom_point(aes(y=Term, x=as.numeric(enrichment.factor), color = as.numeric(p.adj), size = as.integer(Genes))) +
  theme_bw() +
  scale_color_gradientn(colors = c("#441650","#FFBB00"))+
  theme(axis.text.x = element_text(angle = 85,vjust = 0.5, size=10),
        axis.text.y = element_text(hjust = 1, size=10),
        axis.title.x = element_text(vjust = -2, size = 10),
        plot.title=element_text(size=11,face = 2,hjust=0, lineheight = 1), 
        plot.subtitle=element_text(size=9, face=1, hjust=0),
        legend.title = element_text(color = "darkslategrey",size = 9, face = "bold"), 
        plot.margin = unit(c(2, 2, 1,3.5), "cm")) +
  labs(x = "Radio de enrequecimiento", y = "", 
       color="valor padj", 
       size= "Número de genes")  
       # title="GO Enrichment analysis\n Pepper DEG1",
       # subtitle ="Biological process") + 
print(GO_plot)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Working with genes DE in at least one stage only in a specie.                               #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
specific_ort <- list("pepper" = venn_ortho[venn_ortho$pepper == T & venn_ortho$tomato == F,],
                     "tomato" = venn_ortho[venn_ortho$pepper == F & venn_ortho$tomato == T,])

for (i in c("pepper","tomato")) {
  assign(
  paste(i, "list",sep = "_"),
  uniprot_ort %>% mutate(filt=ifelse(ortho %in% 
                                             specific_ort[[i]]$ortho, "yes","no")) %>% 
  filter(filt=="yes") %>%
  separate(uniprot.id, into = c('uniID', 'uniprot_ID'), sep = "\\|") %>%
  dplyr::select(ortho, uniprot_ID) %>% 
  dplyr::mutate(uniprot_ID = gsub("_.*","",uniprot_ID)) %>%
  dplyr::filter(uniprot_ID %in% go_ids[,2])
  )
}

pepper_list <- factor(as.integer(GenesTotales %in% pepper_list$uniprot_ID))
tomato_list <- factor(as.integer(GenesTotales %in% tomato_list$uniprot_ID))
names(pepper_list) <- GenesTotales
names(tomato_list) <- GenesTotales

list_both <- list("pepper"=pepper_list, "tomato"=tomato_list)

for (i in c("pepper","tomato")) {
  assign(paste("GO_res",i, sep = "_"), 
         GenTable(new('topGOdata', 
                      ontology='BP',
                      allGenes = list_both[[i]],
                      annot = annFUN.gene2GO,
                      gene2GO = gene_2_GO), weightFisher=runTest(new('topGOdata', 
                                                                     ontology='BP',
                                                                     allGenes = list_both[[i]], 
                                                                     annot = annFUN.gene2GO, 
                                                                     gene2GO = gene_2_GO), algorithm='weight01', statistic='fisher'),
                  orderBy='weightFisher', topNodes=length(usedGO(new('topGOdata', 
                                                                     ontology='BP',
                                                                     allGenes = list_both[[i]], 
                                                                     annot = annFUN.gene2GO, 
                                                                     gene2GO = gene_2_GO))))
  )
}

#Top of the most significant (padj < 0.01) terms
GO_res_top_pepper <- GO_res_pepper %>%
  dplyr::filter(!stringr::str_detect(Term, "biological")) %>%
  dplyr::filter(as.numeric(weightFisher) < 0.01) %>%
  arrange(as.numeric(weightFisher)) %>%
  head(30) %>%
  dplyr::rename(Genes=Significant, p.adj =weightFisher) %>%
  dplyr::mutate(enrichment.factor = as.numeric(Genes)/as.numeric(Expected))

#Plot data
GO_plot_pepper <- ggplot(GO_res_top_pepper) + 
  geom_point(aes(y=Term, x=as.numeric(enrichment.factor), color = as.numeric(p.adj), size = as.integer(Genes))) +
  theme_bw() +
  scale_color_gradientn(colors = c("#441650","#FFBB00"))+
  theme(axis.text.x = element_text(angle = 85,vjust = 0.5, size=10),
        axis.text.y = element_text(hjust = 1, size=10),
        axis.title.x = element_text(vjust = -2, size = 10),
        plot.title=element_text(size=11,face = 2,hjust=0, lineheight = 1), 
        plot.subtitle=element_text(size=9, face=1, hjust=0),
        legend.title = element_text(color = "darkslategrey",size = 9, face = "bold"), 
        plot.margin = unit(c(2, 2, 1,3.5), "cm"),
        plot.tag = element_text(family = "bold")) +
  labs(x = "Radio de enrequecimiento", y = "", 
       color="valor padj", 
       size= "Número de genes",
       tag = "A")  

print(GO_plot_pepper)

#Top of the most significant (padj < 0.01) terms
GO_res_top_tomato <- GO_res_tomato %>%
  dplyr::filter(!stringr::str_detect(Term, "biological")) %>%
  dplyr::filter(as.numeric(weightFisher) < 0.01) %>%
  arrange(as.numeric(weightFisher)) %>%
  head(30) %>%
  dplyr::rename(Genes=Significant, p.adj =weightFisher) %>%
  dplyr::mutate(enrichment.factor = as.numeric(Genes)/as.numeric(Expected))

#Plot data
GO_plot_tomato <- ggplot(GO_res_top_tomato) + 
  geom_point(aes(y=Term, x=as.numeric(enrichment.factor), color = as.numeric(p.adj), size = as.integer(Genes))) +
  theme_bw() +
  scale_color_gradientn(colors = c("#441650","#FFBB00"))+
  theme(axis.text.x = element_text(angle = 85,vjust = 0.5, size=10),
        axis.text.y = element_text(hjust = 1, size=10),
        axis.title.x = element_text(vjust = -2, size = 10),
        plot.title=element_text(size=11,face = 2,hjust=0, lineheight = 1), 
        plot.subtitle=element_text(size=9, face=1, hjust=0),
        legend.title = element_text(color = "darkslategrey",size = 9, face = "bold"), 
        plot.margin = unit(c(2, 2, 1,3.5), "cm"),
        plot.tag = element_text(family = "bold")) +
  labs(x = "Radio de enrequecimiento", y = "", 
       color="valor padj", 
       size= "Número de genes",
       tag = "B")

print(GO_plot_tomato)

#Save session information
capture.output(sessionInfo(), file = "./R_sessionInfo.txt")


