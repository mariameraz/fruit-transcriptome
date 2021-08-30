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
temp.path <- "~/Data/"

pepper_count <- read.table(paste(temp.path, "pepper_count.txt", sep = ""), header=T)
tomato_count <- read.table(paste(temp.path, "tomato_count.txt", sep = ""), header=T)

#Create metadata table
metadata <- data.frame(stage = rep(c("S1","S2","MG","B"), each=3),
                       sample = paste(rep(c("PRS1","PRS2","PRMG","PRB"), each=3), sep = "_") 
                       %>% paste(., rep(1:3, times=3), sep = "_")) 

rownames(metadata) <- metadata$sample
metadata$stage = factor(metadata$stage)
all(colnames(pepper_count) == rownames(metadata))
all(colnames(tomato_count) == rownames(metadata))

#::::::::::::::::::::::::::::::::::::::::::::#
#      DEG analysis for ALL the genes        #  
#::::::::::::::::::::::::::::::::::::::::::::#

#Merge the both count datasets in a single list
count_data <- list("pepper_count" = pepper_count, "tomato_count" = tomato_count)
names(count_data)
for (i in names(count_data)) {
  #DESEq analysis
    all(colnames(count_data[[i]]) == rownames(metadata))
  
    ddseq <- DESeqDataSetFromMatrix(countData= count_data[[i]], colData= metadata, 
                                    design = ~ stage)
    keep = rowSums(counts(ddseq)) >= 5
    table(keep) 
    ddseq = ddseq[keep,]
    ddseq 
    
    dds <- DESeq(ddseq)
    
    assign(gsub("count","s1_s2",i) , as.data.frame(lfcShrink(dds, contrast= c("stage","S1","S2"))))
    assign(gsub("count","s2_mg",i), as.data.frame(lfcShrink(dds, contrast = c("stage", "S2","MG"))))
    assign(gsub("count","mg_b",i), as.data.frame(lfcShrink(dds, contrast = c("stage","MG","B"))))

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
for (i in contras) {
  for (j in c("pepper","tomato")) {
   assign(paste(paste(j, i, sep = "_"), "sig", sep = "_"), res[[paste(j,i, sep = "_")]] %>% 
            dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01))
     
  }
}

res_sig <- list("pepper_s1_s2"=pepper_s1_s2_sig, "pepper_s2_mg"=pepper_s2_mg_sig, "pepper_mg_b"=pepper_mg_b_sig,
                "tomato_s1_s2"=tomato_s1_s2_sig, "tomato_s2_mg"=tomato_s2_mg_sig, "tomato_mg_b"=tomato_mg_b_sig)
#Create logical table
venn_ortho_pepper <- tibble(pepper_count %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(pepper_s1_s2= ortho %in% rownames(pepper_s1_s2_sig)) %>%
  dplyr::mutate(pepper_s2_mg= ortho %in% rownames(pepper_s2_mg_sig)) %>%
  dplyr::mutate(pepper_mg_b= ortho %in% rownames(pepper_mg_b_sig))

venn_ortho_tomato <- tibble(tomato_count %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(tomato_s1_s2= ortho %in% rownames(tomato_s1_s2_sig)) %>%
  dplyr::mutate(tomato_s2_mg= ortho %in% rownames(tomato_s2_mg_sig)) %>%
  dplyr::mutate(tomato_mg_b= ortho %in% rownames(tomato_mg_b_sig))

names(res_sig)

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

all(colnames(pepper_count_ort) == rownames(metadata))
all(colnames(tomato_count_ort) == rownames(metadata))

#::::::::::::::::::::::::::::::::::::::::::::#
#      DEG analysis for ALL the genes        #  
#::::::::::::::::::::::::::::::::::::::::::::#

#Merge the both count datasets in a single list
count_data <- list("pepper_count" = pepper_count_ort, "tomato_count" = tomato_count_ort)
names(count_data)
for (i in names(count_data)) {
  #DESEq analysis
  all(colnames(count_data[[i]]) == rownames(metadata))
  
  ddseq <- DESeqDataSetFromMatrix(countData= count_data[[i]], colData= metadata, 
                                  design = ~ stage)
  keep = rowSums(counts(ddseq)) >= 5
  table(keep) 
  ddseq = ddseq[keep,]
  ddseq 
  
  dds <- DESeq(ddseq)
  
  assign(gsub("count","s1_s2",i) , as.data.frame(lfcShrink(dds, contrast= c("stage","S1","S2"))))
  assign(gsub("count","s2_mg",i), as.data.frame(lfcShrink(dds, contrast = c("stage", "S2","MG"))))
  assign(gsub("count","mg_b",i), as.data.frame(lfcShrink(dds, contrast = c("stage","MG","B"))))
  
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
for (i in contras) {
  for (j in c("pepper","tomato")) {
    assign(paste(paste(j, i, sep = "_"), "sig", sep = "_"), res[[paste(j,i, sep = "_")]] %>% 
             dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01))
    
  }
}

res_sig <- list("pepper_s1_s2"=pepper_s1_s2_sig, "pepper_s2_mg"=pepper_s2_mg_sig, "pepper_mg_b"=pepper_mg_b_sig,
                "tomato_s1_s2"=tomato_s1_s2_sig, "tomato_s2_mg"=tomato_s2_mg_sig, "tomato_mg_b"=tomato_mg_b_sig)
#Create logical table
venn_ortho_pepper <- tibble(pepper_count_ort %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(pepper_s1_s2= ortho %in% rownames(pepper_s1_s2_sig)) %>%
  dplyr::mutate(pepper_s2_mg= ortho %in% rownames(pepper_s2_mg_sig)) %>%
  dplyr::mutate(pepper_mg_b= ortho %in% rownames(pepper_mg_b_sig))

venn_ortho_tomato <- tibble(tomato_count_ort %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(tomato_s1_s2= ortho %in% rownames(tomato_s1_s2_sig)) %>%
  dplyr::mutate(tomato_s2_mg= ortho %in% rownames(tomato_s2_mg_sig)) %>%
  dplyr::mutate(tomato_mg_b= ortho %in% rownames(tomato_mg_b_sig))

names(res_sig)

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

venn_ortho <- tibble(pepper_count_ort %>% rownames_to_column(., "ortho")) %>%
  dplyr::select("ortho") %>%
  dplyr::mutate(tomato_s1_s2= ortho %in% rownames(tomato_s1_s2_sig)) %>%
  dplyr::mutate(tomato_s2_mg= ortho %in% rownames(tomato_s2_mg_sig)) %>%
  dplyr::mutate(tomato_mg_b= ortho %in% rownames(tomato_mg_b_sig)) %>%
  dplyr::mutate(pepper_s1_s2= ortho %in% rownames(pepper_s1_s2_sig)) %>%
  dplyr::mutate(pepper_s2_mg= ortho %in% rownames(pepper_s2_mg_sig)) %>%
  dplyr::mutate(pepper_mg_b= ortho %in% rownames(pepper_mg_b_sig)) %>%
  dplyr::mutate(pepper= ortho %in% c(rownames(pepper_s1_s2_sig),
                                     rownames(pepper_s2_mg_sig),
                                     rownames(pepper_mg_b_sig))) %>%
  dplyr::mutate(tomato=ortho %in% c(rownames(tomato_s1_s2_sig),
                                    rownames(tomato_s2_mg_sig),
                                    rownames(tomato_mg_b_sig)))

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
