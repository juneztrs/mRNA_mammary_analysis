library("tidyverse")

counts <- read.csv("Junez_R_proyecto/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")

sampleinfo <- read.csv("Junez_R_proyecto/GSE60450_filtered_metadata.csv")


colnames(counts)[1] <- "gene_id"
colnames(sampleinfo)[1] <- "sample_id"

view(counts)
view(sampleinfo)

seqdata <- pivot_longer(counts, cols = starts_with("GSM"),
                        names_to = "Sample", values_to = "Count")

seqdata <- pivot_longer(counts, cols =GSM1480291:GSM1480302
, names_to = "Sample", values_to = "Count" )

seqdata <- pivot_longer(counts, cols = -c ("gene_id", "gene_symbol"),
                        names_to = "Sample", values_to = "Count")


view(seqdata)

allinfo <- full_join(seqdata, sampleinfo, by = c("Sample"="sample_id"))
                     
                     
view(allinfo)

ggplot(data = allinfo, mapping = aes(x =Sample, y = Count)) +
  geom_boxplot()

ggplot(data = allinfo, mapping = aes(x =Sample, y = log2(Count))) +
  geom_boxplot()


ggplot(data = allinfo, mapping = aes(x =Sample, y = log2(Count + 1))) +
  geom_boxplot()
 


ggplot(data = allinfo, mapping = aes(x =Sample, y = log2(Count))) +
  geom_violin()


ggplot(data = allinfo, mapping = aes(x =Sample, y = log2(Count +1), colour = Sample)) +
  geom_violin()

ggplot(data = allinfo, mapping = aes(x =Sample, weight = log2(Count +1),
                                     fill= Sample)) +
  geom_bar()
  
data("women")

view(women) 

ggplot(data = women, mapping = aes(x =height, y = weight))+
  geom_point()


ggplot(data = women, mapping = aes(x =height, y = weight))+
  geom_path()







ggplot(data = women, mapping = aes(x = height))+
  geom_density()

ggplot(data = allinfo, mapping = aes(x = log2(Count +1), fill = Sample)) +
  geom_density()


ggplot(data = allinfo, mapping = aes(x = log2(Count +1), color = Sample)) +
  geom_density()

pdf("My_first_plot.pdf")
ggplot(data = allinfo, mapping = aes(x = log2(Count +1), color = Sample)) +
  geom_density()
dev.off()

view(sampleinfo)


colnames(sampleinfo)

?mutate

view(women)

women <- mutate(women, BMI = weight / (height**2) * 703)















allinfo <- mutate(allinfo, Group = case_when(
  
  str_detect(characteristics, "basal.*virgin") ~ "bvirg",
  str_detect(characteristics, "basal.*preg") ~ "bpreg",
  str_detect(characteristics, "basal.*lact") ~ "blact",
  str_detect(characteristics, "luminal.*virgin") ~ "lvirg",
  str_detect(characteristics, "luminal.*preg") ~ "lpreg",
  str_detect(characteristics, "luminal.*lact") ~ "llact",
  ))


view(allinfo)
 
my_genes <- allinfo %>%
  group_by(gene_symbol) %>%
  summarise(Total_count = sum(Count)) %>%
  arrange(desc(Total_count)) %>%
  head(n = 10) %>%
  pull(gene_symbol)

my_genes

my_genes_count <- filter(allinfo, gene_symbol %in% my_genes)

view(my_genes_count)

ggplot(my_genes_count, mapping = aes(x = Group, y = log2(Count+1), 
                                     colour = Group))+
  geom_point()+
  facet_wrap( ~ gene_symbol)


pdf("Mammary_Gland_RNA-seq_data.pdf")
ggplot(my_genes_count, mapping = aes(x = Group, y = log2(Count+1), 
                                     colour = Group))+
  geom_point()+
  facet_wrap( ~ gene_symbol)+
  labs(x = "Cell type and stage", y = "Count", title = "Mammary gland RNA-seq data")+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

 