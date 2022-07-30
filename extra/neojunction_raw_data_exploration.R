library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggsci)

events <- read_delim("~/Downloads/merge_graphs_intron_retention_C2.confirmed.txt.gz", delim = "\t")

# intron retention explanation from spladder:
# features intron_retention:
# e1_cov: mean coverage of the left flanking exon (in genomic coordinates, ignoring strand)
# e2_cov: mean coverage of the retained intron
# e3_cov: mean coverage of the right flanking exon (in genomic coordinates, ignoring strand)
# e1e3_conf: number of spliced alignments spanning the intron
# e2_cov_region: fraction of positions in the intron that have a coverage > 0
# to count the number of spliced alignments that confirm the connection of exon segments e1 and e3 in an exon skip, the corresponding feature name would be e1e3_conf




events %>% filter(strand == '+') %>% select(ends_with(":e1_cov")) -> tmp
names(tmp) <- str_replace_all(names(tmp), c("ohuman_ML-" = "", "Aligned.sortedByCoord.unique.out:e1_cov" =""))
tmp %>% pivot_longer(everything()) -> e1_cov_pos
summary(e1_cov_pos$value)

events %>% filter(strand == '+') %>% select(ends_with("e2_cov")) -> tmp
names(tmp) <- str_replace_all(names(tmp), c("ohuman_ML-" = "", "Aligned.sortedByCoord.unique.out:e2_cov" =""))
tmp %>% pivot_longer(everything()) -> e2_cov_pos
summary(e2_cov_pos$value)

events %>% filter(strand == '+') %>% select(ends_with("e3_cov")) -> tmp
names(tmp) <- str_replace_all(names(tmp), c("ohuman_ML-" = "", "Aligned.sortedByCoord.unique.out:e3_cov" =""))
tmp %>% pivot_longer(everything()) -> e3_cov_pos
summary(e3_cov_pos$value)

events %>% filter(strand == '+') %>% select(ends_with(":e1e3_conf")) -> tmp
names(tmp) <- str_replace_all(names(tmp), c("ohuman_ML-" = "", "Aligned.sortedByCoord.unique.out:e1e3_conf" =""))
tmp %>% pivot_longer(everything()) -> e1e3_conf_pos
summary(e1e3_conf_pos$value)

events %>% filter(strand == '+') %>% select(ends_with(":psi")) -> tmp
names(tmp) <- str_replace_all(names(tmp), c("ohuman_ML-" = "", "Aligned.sortedByCoord.unique.out:psi" =""))
tmp %>% pivot_longer(everything()) -> psi_pos
summary(psi_pos$value)

a <- ggplot(e1_cov_pos, aes(x=name,y = value, fill = name)) + 
  geom_boxplot() + 
  ylim(0, summary(e1_cov_pos$value)['Median'] + summary(e1_cov_pos$value)['3rd Qu.']) + 
  theme_pubclean() +
  #labs(fill = "samples") +
  theme(legend.position = "none") +
  xlab("sample") + 
  ylab('exon1 coverage') + 
  ggtitle("Pre-retention exon coverage (+ strand, exon1, trunc axes)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_npg(palette = c("nrc"))

b<-ggplot(e2_cov_pos, aes(x=name,y = value, fill = name)) + 
  geom_boxplot() + 
  ylim(0, summary(e2_cov_pos$value)['Median'] + summary(e2_cov_pos$value)['3rd Qu.']) + 
  theme_pubclean() + 
  #labs(fill = "samples") +
  xlab("sample") + 
  theme(legend.position = "none") +
  ylab('exon2 coverage') + 
  ggtitle("intron retention coverage (+ strand, exon2, trunc axes)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_npg(palette = c("nrc"))

c<-ggplot(e3_cov_pos, aes(x=name,y = value, fill = name)) + 
  geom_boxplot() + 
  ylim(0, summary(e3_cov_pos$value)['Median'] + summary(e3_cov_pos$value)['3rd Qu.']) + 
  theme_pubclean() +
  #labs(fill = "samples") +
  theme(legend.position = "none") +
  xlab("sample") + 
  ylab('exon3 coverage') + 
  ggtitle("Post-retention exon coverage (+ strand, exon3, trunc axes)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_npg(palette = c("nrc"))



d<-ggplot(e1e3_conf_pos, aes(x=name,y = value, fill = name)) + 
  geom_boxplot() + 
  ylim(0, summary(e1e3_conf_pos$value)['Median'] + summary(e1e3_conf_pos$value)['3rd Qu.']) + 
  theme_pubclean() +
  #labs(fill = "samples") +
  theme(legend.position = "none") +
  xlab("sample") + 
  ylab('exon1-exon3 spliced alignments') + 
  ggtitle("Number of splice alignments connecting exon1-exon3 (+ strand, trunc axes)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_npg(palette = c("nrc"))


ggarrange(a, b, c, d, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = T)

