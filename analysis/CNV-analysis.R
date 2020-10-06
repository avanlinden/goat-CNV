### Wild and domestic goat Copy Number Variation
### 5-5-20

library(tidyverse)

### import domestic CNV data ========

CNV_col_names <- c("Chromosome", "Start", "End", "CN", "Type")

domCNVs <- read_tsv("ERR340337_vs_REF.sort.bam_CNVs", col_names = CNV_col_names)

domCNVs %>% 
  filter(CN != 0) %>% 
  filter(Type == 'loss') %>% 
  count(Chromosome)

wildCNVs <- read_tsv("ERR340341_vs_REF.sort.bam_CNVs", col_names = CNV_col_names)

wildCNVs %>% 
  arrange(desc(CN)) %>% 
  select(chrName, Start, End, CN, Type, Goat) %>% 
  ggplot(aes(x = CN)) +
  geom_dotplot(binwidth = 5) +
  xlim(0, 650) +
  labs(title = "Wild CN")+
  scale_y_continuous(NULL, breaks = NULL) +
  theme_classic() 


domCNVs %>% 
  arrange(desc(CN)) %>% 
  select(chrName, Start, End, CN, Type, Goat) %>% 
  ggplot(aes(x = CN)) +
    geom_dotplot(binwidth = 5) +
    xlim(0, 650) +
  labs(title = "Domestic CN") +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_classic() 

wildCNVs %>% 
  inner_join(domCNVs, by = c('Chromosome', 'Start', 'End'), suffix = c(".wild", ".dom")) %>% 
  left_join(chrLength, by = 'Chromosome') %>% 
  filter(CN.wild != CN.dom)
  

wildCNVs <- wildCNVs %>% 
  left_join(chrLength, by = 'Chromosome') %>% 
  mutate(Goat = "wild")

domCNVs <- domCNVs %>% 
  left_join(chrLength, by = 'Chromosome') %>% 
  mutate(Goat = 'domestic')

allCNVs <- bind_rows(wildCNVs, domCNVs)

allCNVs %>% 
  arrange(desc(CN)) %>% 
  ggplot(aes(x = CN)) +
    geom_density(aes(fill = Goat))

allCNVs %>% 
  select(-chrLen) %>% 
  mutate(chr = chrNum) %>% 
  filter(CN < 160) %>% 
  filter(chr %in% 17:31) %>% 
  ggplot(aes(x = Start, y = CN)) +
    geom_point(aes(colour = Goat, shape = Type), alpha = 7/10, size = 1, position = position_jitter(height = 1), stroke = 0.75) +
    facet_wrap(vars(chr), labeller = label_both, ncol = 4) +
    scale_x_continuous(name = "CN start position on chromosome", breaks = c(0, 1e8), limits = c(0, 151032000)) +
    ylim(0, 160) + 
    scale_shape_manual(values = c(19, 1)) + 
    theme_minimal() 

### import chromosome names =========

chrNames <- read_tsv("chr_len_reformat.txt", col_names = c("chrNum", "Chromosome", "chrLen"))

chrLength <- chrNames %>% 
  mutate(chr = "chr") %>% 
  mutate(chrName = str_c(chr, chrNum)) %>% 
  mutate(chrName = if_else(chrName == 'chr30', 'chrX', chrName)) %>% 
  mutate(chrName = if_else(chrName == 'chr31', 'chrY', chrName)) %>% 
  select(-chr)


wildCNVs %>% 
  arrange(desc(CN))