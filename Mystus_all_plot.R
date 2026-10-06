setwd("E:/RProject/RSCU_Mystus")
library(RSCUcaller)
library(cubar)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(patchwork)

KaKs <- read_excel("E:/Manuscript/2025_Mystus_Singaringan/RSCU_mystus.xlsx",
                   sheet = "KaKs")
KaKs_long <- KaKs %>%
  pivot_longer(
    cols = everything(),
    names_to = "Gene",
    values_to = "KaKs"
  )

head(KaKs_long)
gene_order 
KaKs_summary <- KaKs_long %>%
  group_by(Gene) %>%
  summarise(
    n = n(),
    Mean = mean(KaKs),
    SD = sd(KaKs),
    SE = SD/sqrt(n),
    Median = median(KaKs),
    Min = min(KaKs),
    Max = max(KaKs)
  )

KaKs_summary
KaKs_plot = ggplot(KaKs_long,
                   aes(Gene, KaKs, fill = Gene)) +
  
  geom_boxplot(width = 0.7,
               alpha = 0.85,
               outlier.shape = NA) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.7) +
  
  scale_fill_brewer(palette = "Paired") +
  
  labs(
    x = "Mitochondrial gene",
    y = expression(K[a]/K[s])
  ) + 
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    )
  )
KaKs_plot


library(Biostrings)
cds <- readDNAStringSet("E:/RProject/RSCU_Mystus/each_gene_Bagrid.fasta")
cds_full <- readDNAStringSet("E:/RProject/RSCU_Mystus/Bagrid_noSTOP_extraction.fasta")
mito_code <- get_codon_table(gcid = '2')
head(mito_code)
bagrid_mt_qc <- check_cds(
  cds,
  codon_table = mito_code,
  check_stop = FALSE,
  rm_stop = FALSE,
  check_len = FALSE,
  start_codons = c('ATG', 'ATA', 'ATT', "GTG"))

bagrid_mt_qc 
mt_cf <- count_codons(bagrid_mt_qc)
enc = get_enc(mt_cf, codon_table = mito_code)
enc

bagrid_mtfull_qc <- check_cds(
  cds_full,
  codon_table = mito_code,
  check_stop = FALSE,
  rm_stop = FALSE,
  check_len = FALSE,
  start_codons = c('ATG', 'ATA', 'ATT', "GTG"))
bagrid_mtfull_qc 
mtfull_cf <- count_codons(bagrid_mtfull_qc)
enc_full = get_enc(mtfull_cf, codon_table = mito_code)
enc_full


write.table(enc_full,
            file="G:/ENC_full_bagridae.txt",
            row.names=T, sep="\t")

#calculate CAI of each gene
cf_all <- count_codons(cds)
rscu <- est_rscu(mt_cf, codon_table = mito_code)
head(rscu)
cai <- get_cai(cf_all, rscu)
head(cai)
hist(cai, main = "Distribution of CAI values")
write.table(cai,
            file="G:/CAI_gene_bagridae.txt",
            row.names=T, sep="\t")

tai = get_tai(mt_cf, w, w_format = "cubar")
write.table(tai,
            file="G:/tAI_gene_bagridae.txt",
            row.names=T, sep="\t")


cf_full <- count_codons(cds_full)
rscu_full <- est_rscu(cf_full, codon_table = mito_code)
head(rscu_full)
cai_full <- get_cai(cf_full, rscu_full)
head(cai_full)
hist(cai_full, main = "Distribution of CAI values", xlab = "Codon Adaptation Index")




write.table(cai_full,
            file="G:/CAIfull_bagridae.txt",
            row.names=T, sep="\t")

write.table(tai,
            file="G:/tAIfull_bagridae.txt",
            row.names=T, sep="\t")



#calculate GC12
# Use seqinr's built-in GC bias functions on your sequences
# GC1, GC2, and GC3
library(Biostrings)
library(dplyr)

# Function to calculate nucleotide composition
calc_gc_stats <- function(seq){
  
  seq <- toupper(as.character(seq))
  
  # split into codons
  codons <- substring(seq,
                      seq(1, nchar(seq)-2, 3),
                      seq(3, nchar(seq), 3))
  
  # remove ambiguous codons
  codons <- codons[grepl("^[ATGC]{3}$", codons)]
  
  # codon positions
  pos1 <- substring(codons,1,1)
  pos2 <- substring(codons,2,2)
  pos3 <- substring(codons,3,3)
  
  GC1 <- mean(pos1 %in% c("G","C"))
  GC2 <- mean(pos2 %in% c("G","C"))
  GC3 <- mean(pos3 %in% c("G","C"))
  GC12 <- (GC1 + GC2)/2
  
  A3 <- sum(pos3=="A")
  T3 <- sum(pos3=="T")
  G3 <- sum(pos3=="G")
  C3 <- sum(pos3=="C")
  
  A3_AT3 <- A3/(A3+T3)
  G3_GC3 <- G3/(G3+C3)
  
  data.frame(
    GC1 = GC1,
    GC2 = GC2,
    GC3 = GC3,
    GC12 = GC12,
    A3_AT3 = A3_AT3,
    G3_GC3 = G3_GC3
  )
}

#Apply to every CDS:

gc_stats <- bind_rows(
  lapply(cds, calc_gc_stats),
  .id = "Gene"
)

gc_stats
write.table(gc_stats,
            file="G:/gc_stats_gene_bagridae.txt",
            row.names=T, sep="\t")

library(readxl)
library(GGally)
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

#------------------------------------------------------------
# Read Excel file
#------------------------------------------------------------
dat <- read_excel("E:/Manuscript/2025_Mystus_Singaringan/RSCU_mystus.xlsx",
                  sheet = "Sheet2")

# Variables to plot
df <- dat %>%
  select(ENC, CAI, tAI, GC1, GC2, GC3, GC12, A3_AT3, G3_GC3)

#------------------------------------------------------------
# Upper panel: correlation coefficient + significance
#------------------------------------------------------------
cor_panel <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ct <- cor.test(x, y)
  
  r <- unname(ct$estimate)
  p <- ct$p.value
  
  stars <- ifelse(p < 0.001,"***",
                  ifelse(p < 0.01,"**",
                         ifelse(p < 0.05,"*","")))
  
  # Background colour
  fill <- colorRampPalette(c("#B39DDB","white","#39FF14"))(100)[
    round(rescale(r,c(1,100)))
  ]
  
  ggplot()+
    annotate("rect",
             xmin=-Inf,xmax=Inf,
             ymin=-Inf,ymax=Inf,
             fill=fill)+
    annotate("text",
             x=0.5,
             y=0.5,
             label=paste0(sprintf("%.3f",r),stars),
             size=7,
             fontface="bold")+
    theme_void()
}

#------------------------------------------------------------
# Lower panel: scatter + loess
#------------------------------------------------------------
lower_panel <- function(data,mapping,...){
  
  ggplot(data=data,mapping=mapping)+
    geom_point(size=1.6,
               shape=21,
               fill="white",
               colour="black")+
    geom_smooth(method="loess",
                se=FALSE,
                colour="red",
                linewidth=0.7)+
    theme_bw()
}

#------------------------------------------------------------
# Diagonal: histogram + density
#------------------------------------------------------------
diag_panel <- function(data,mapping,...){
  
  ggplot(data=data,mapping=mapping)+
    geom_histogram(aes(y=after_stat(density)),
                   bins=8,
                   fill="#18DDE3",
                   colour="black")+
    geom_density(colour="red",
                 linewidth=0.8)+
    theme_bw()
}


library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

ENC_GC3_plot = ggplot(dat, aes(GC3, ENC)) +
  
  # Wright expected curve
  geom_line(data = enc_curve,
            aes(GC3, ENC),
            colour = "darkorange",
            linewidth = 1,
            linetype = "dashed",
            inherit.aes = FALSE) +
  
  # points
  geom_point(
    aes(fill = gene),
    shape = 21,
    colour = "black",
    size = 3.5,
    stroke = 0.4
  ) +
  
  scale_fill_brewer(palette = "Set3") +
  
  labs(
    x = "GC3s",
    y = "ENC",
    fill = "Gene"
  ) +
  
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(20,61)
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.title = element_text(face="bold", size=15),
    axis.text = element_text(colour="black", size=12),
    legend.position="right",
    legend.title=element_text(face="bold"),
    panel.border=element_rect(fill=NA, colour="black")
  ) +
  
  ggtitle("A")


#---------------------------------------------------------
# Read data
#---------------------------------------------------------
ENC_data <- read_excel("ENC_data_Bagridae.xlsx")

# Make factors
ENC_data <- ENC_data %>%
  mutate(
    genus = factor(genus),
    gene = factor(gene,
                  levels = c("ATP8","ATP6","COX1","COX2","COX3",
                             "CYTB","ND1","ND2","ND3",
                             "ND4","ND4L","ND5","ND6"))
  )
library(dplyr)

delta_genus <- ENC_data %>%
  group_by(genus) %>%
  summarise(
    n = n(),
    Mean = mean(delta_ENC, na.rm = TRUE),
    SD = sd(delta_ENC, na.rm = TRUE),
    SE = SD / sqrt(n),
    Min = min(delta_ENC, na.rm = TRUE),
    Max = max(delta_ENC, na.rm = TRUE)
  ) %>%
  arrange(desc(Mean))

delta_genus

gene_order <- c("ATP6","ATP8","COX1","COX2","COX3",
                "CYTB","ND1","ND2","ND3",
                "ND4","ND4L","ND5","ND6")

delta_gene <- ENC_data %>%
  mutate(gene = factor(gene)) %>%
  group_by(gene) %>%
  summarise(
    n = n(),
    Mean = mean(delta_ENC, na.rm = TRUE),
    SD = sd(delta_ENC, na.rm = TRUE),
    SE = SD / sqrt(n),
    Min = min(delta_ENC, na.rm = TRUE),
    Max = max(delta_ENC, na.rm = TRUE)
  ) %>%
  arrange(gene)

delta_gene
#---------------------------------------------------------
# Figure A : by Genus
#---------------------------------------------------------
pB <- ggplot(ENC_data,
             aes(x = genus,
                 y = delta_ENC,
                 fill = genus)) +
  
  geom_boxplot(width = 0.65,
               outlier.shape = NA,
               alpha = 0.8) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.7) +
  
  labs(x = "Genus",
       y = "ENC ratio") +
  
  scale_fill_brewer(palette = "Set2") +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               face = "italic"),
    plot.title = element_text(face = "bold")
  ) +
  
  ggtitle("B")
pB
#---------------------------------------------------------
# Figure B : by Gene
#---------------------------------------------------------
pC <- ggplot(ENC_data,
             aes(x = gene,
                 y = delta_ENC,
                 fill = gene)) +
  
  geom_boxplot(width = 0.65,
               outlier.shape = NA,
               alpha = 0.8) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.7) +
  
  labs(x = "Mitochondrial gene",
       y = "ENC ratio") +
  
  scale_fill_brewer(palette = "Paired") +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               face = "italic"),
    plot.title = element_text(face = "bold")
  ) +
  
  ggtitle("C")
pC
#---------------------------------------------------------
# Combine panels
#---------------------------------------------------------
ENC_final_fig6 <- ggarrange(
  ENC_GC3_plot,
  ggarrange(
    pB,
    pC,
    ncol = 1
  ),
  ncol = 2,
  widths = c(1.2, 1)   # A is 40% wider than B+C
)
ENC_final_fig6


#------------------------------------------------------------
# Read Excel file
#------------------------------------------------------------
dat <- read_excel("E:/Manuscript/2025_Mystus_Singaringan/RSCU_mystus.xlsx",
                  sheet = "Sheet2")

# Variables to plot
df <- dat %>%
  select(ENC, CAI, tAI, GC1, GC2, GC3, GC12, A3_AT3, G3_GC3)

#------------------------------------------------------------
# Upper panel: correlation coefficient + significance
#------------------------------------------------------------
cor_panel <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ct <- cor.test(x, y)
  
  r <- unname(ct$estimate)
  p <- ct$p.value
  
  stars <- ifelse(p < 0.001,"***",
                  ifelse(p < 0.01,"**",
                         ifelse(p < 0.05,"*","")))
  
  # Background colour
  fill <- colorRampPalette(c("#B39DDB","white","#39FF14"))(100)[
    round(rescale(r,c(1,100)))
  ]
  
  ggplot()+
    annotate("rect",
             xmin=-Inf,xmax=Inf,
             ymin=-Inf,ymax=Inf,
             fill=fill)+
    annotate("text",
             x=0.5,
             y=0.5,
             label=paste0(sprintf("%.3f",r),stars),
             size=7,
             fontface="bold")+
    theme_void()
}

#------------------------------------------------------------
# Lower panel: scatter + loess
#------------------------------------------------------------
lower_panel <- function(data,mapping,...){
  
  ggplot(data=data,mapping=mapping)+
    geom_point(size=1.6,
               shape=21,
               fill="white",
               colour="black")+
    geom_smooth(method="loess",
                se=FALSE,
                colour="red",
                linewidth=0.7)+
    theme_bw()
}

#------------------------------------------------------------
# Diagonal: histogram + density
#------------------------------------------------------------
diag_panel <- function(data,mapping,...){
  
  ggplot(data=data,mapping=mapping)+
    geom_histogram(aes(y=after_stat(density)),
                   bins=8,
                   fill="#18DDE3",
                   colour="black")+
    geom_density(colour="red",
                 linewidth=0.8)+
    theme_bw()
}

#------------------------------------------------------------
# Create plot
#------------------------------------------------------------
library(GGally)
fig7 <- ggpairs(
  df,
  upper=list(continuous=cor_panel),
  lower=list(continuous=lower_panel),
  diag=list(continuous=diag_panel),
  columnLabels=names(df),
  progress=FALSE
)+
  theme_bw(base_size=14)+
  theme(
    strip.background=element_rect(fill="grey90"),
    strip.text=element_text(size=16,face="bold"),
    panel.grid=element_blank()
  )

fig7