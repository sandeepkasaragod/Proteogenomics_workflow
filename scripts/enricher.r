library(STRINGdb)
STRINGdb$methods()

#### Specify taxonomy ID from STRINGdb for species ####
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=400)

#### StringDB_Annotation ####
#annotations <- string_db$get_annotations(proteins_vector)

#### Fetch the list of proteins or genes for STRINGdb search ####
setwd("C:/Users/Sandeep/Desktop/snakemake/")
proteins <- read.csv("Enrichment_analysis_genes.txt", sep = "\t", header = F)
colnames(proteins)

#### StringDB_enrichment ####
enrich_proteins <- string_db$get_enrichment(proteins)
#enrich_proteins$Modification <- rep("STRINGdb_enrichment",each=dim(enrich_proteins)[1])
write.table(enrich_proteins,file = "DEX_Proteins_Gene_Symbol_STRINGdb_enrichment.txt", sep = "\t")

#### plot PPI network ####
string_db$plot_network(enrich_proteins)

####
DEX_enrich_proteins <- read.csv("DEX_Proteins_Gene_Symbol_STRINGdb_enrichment.txt", sep = '\t', header = T)

DEX_CC <- DEX_enrich_proteins[DEX_enrich_proteins$category == "Component", ]
DEX_MF <- DEX_enrich_proteins[DEX_enrich_proteins$category == "Function", ]
DEX_BP <- DEX_enrich_proteins[DEX_enrich_proteins$category == "Process", ]


#### Bubble plot ####
library(viridis)
library(ggrepel)

data_DEX <- rbind(DEX_CC,DEX_MF,DEX_BP)

ggplot(data_DEX, aes(term,-log(p_value),label=description))+
  geom_point(aes(color = category, size = log(number_of_genes), alpha=0.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "darkslateblue", "blueviolet", "chartreuse", "violetred1"),#, "navyblue", "hotpink4", "salmon", "mediumaquamarine"),
                     labels=c("CC","MF","BP")) +
  scale_size(range = c(0.5, 12))+
  facet_wrap(~category,  ncol = 3, scales = "fixed") +
  geom_text_repel(data=subset(data_DEX, -log(p_value) > 10),aes(term,-log(p_value),label=description), size=4)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

summary <- string_db$get_summary(proteins_vector)
all <- string_db$load()
plot(all)


#### Dotplot ####
CC <- enrichment[1:15,]
dotplot <- ggplot(String_enrich, aes(x=category, y=term, fill=category)) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45,  hjust = 1), plot.title = element_text(size=10))

plot(dotplot)


