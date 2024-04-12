# Soil microbiota and herbivory drive the assembly of plant-associated microbial communities through different mechanisms

**Antonino Malacrinò, Alison E. Bennett**

## Abstract

Plant-associated microbial communities are key to shaping many aspects of plant biology. In this study, we tested whether soil microbial communities and herbivory influence the bacterial community of tomato plants and whether their influence in different plant compartments is driven by microbial spillover between compartments or whether plants are involved in mediating this effect. We grew our plants in soils hosting three different microbial communities and covered (or not) the soil surface to prevent (or allow) passive microbial spillover between compartments, and we exposed them (or not) to herbivory by Manduca sexta. Here we show that the soil-driven effect on aboveground compartments is consistently detected regardless of soil coverage, whereas soil cover influences the herbivore-driven effect on belowground microbiota. Together, our results suggest that the soil microbiota influences aboveground plant and insect microbial communities via changes in plant metabolism and physiology or by sharing microorganisms via xylem sap. In contrast, herbivores influence the belowground plant microbiota via a combination of microbial spillover and changes in plant metabolism. These results demonstrate the important role of plants in linking aboveground and belowground microbiota, and can foster further research on soil microbiota manipulation for sustainable pest management.

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject number `PRJNA910821`.

Our pipeline included:
* nf-core/ampliseq v2.3.2 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* MAFFT [https://academic.oup.com/nar/article/30/14/3059/2904316](https://academic.oup.com/nar/article/30/14/3059/2904316)
* FastTree [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)
* R v4.3.2 [https://www.R-project.org/](https://www.R-project.org/)

AM is supported by the Italian Ministry of University and Research (MUR) through the PRIN 2022 PNRR program (project P2022KY74N, financed by the European Union - NextGenerationEU).

# Data processing

Run ampliseq

```bash
nextflow run nf-core/ampliseq -r 2.3.2 -profile singularity \
--input $INDIR \
--FW_primer "GTGYCAGCMGCCGCGGTAA" \
--RV_primer "GGACTACNVGGGTWTCTAAT" \
--outdir $OUTDIR \
--illumina_novaseq \
--extension "/*_{1,2}.fastq.gz" \
--trunclenf 0 \
--trunclenr 0 \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--ignore_empty_input_files \
--email "antonino.malacrino@gmail.com" \
--max_cpus 16 \
--max_memory '128.GB'

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```


# Data analysis

## Load libraries

```r
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library("DESeq2")
library("vegan")
library("ape")
library("plyr") 
library("ggrepel")
library("dplyr")
library("broom")
library("picante")
library("Rmisc")
library("emmeans")
library("car")
library("lme4")
library("tibble")
library("data.table")
library("limma")
library("microbiome")
library("ggvenn")
library("patchwork")
library("RColorBrewer")
library("RVAideMemoire")
library("data.table")
library("MuMIn")
```

## Load and clean microbiome data

```r
asv.table <- as.matrix(read.table('data/ASV_table.tsv', sep = "\t", header = T, row.names = 1))
tax.table <- as.matrix(read.table('data/ASV_tax_species.tsv', sep = "\t", header = T, row.names = 1)[,c(2:8)])
mapping <- read.table('data/metadata.txt', sep = "\t", header = T, row.names = 1)
tree <- read.tree('data/tree.tre')

ps.16S <- phyloseq(otu_table(asv.table, taxa_are_rows = T), 
                         sample_data(mapping), 
                         tax_table(tax.table), 
                         phy_tree(tree))

ps.16S <- subset_taxa(ps.16S, Class !="Chloroplast")
ps.16S <- subset_taxa(ps.16S, Order !="Chloroplast")
ps.16S <- subset_taxa(ps.16S, Family !="Mitochondria")

ps.16S <- subset_samples(ps.16S, Compartment != "Control") 
ps.16S <- subset_samples(ps.16S, Soil != "No_treatment") 
ps.16S <- filter_taxa(ps.16S, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.16S.L <- subset_samples(ps.16S, Compartment == "Leaves") 
ps.16S.R <- subset_samples(ps.16S, Compartment == "Roots") 
ps.16S.S <- subset_samples(ps.16S, Compartment == "Soil") 
ps.16S.H <- subset_samples(ps.16S, Compartment == "Herbivore")
```

## Phylogenetic diversity

Code for: Fig. 2A, Tab. S1, and Figs. S1-S3.

### Calculate Phylogenetic Diversity.

```r
calc.div <- function(ps){
  div <- microbiome::alpha(ps, index = "diversity_shannon")
  otus <- as.data.frame(t(otu_table(ps)))
  tree <- phy_tree(ps)
  div.pd <- pd(otus, tree, include.root = FALSE)
  div.2 <- cbind(sample_data(ps), div)
  div.2 <- cbind(div.2, div.pd)
  return(div.2)
}

div.df <- calc.div(ps.16S)
```

### Model, run post-hoc contrasts, and estimate variance explained for **rhizosphere** samples

```r
model.div <- function(x){
  model.qv <- lmer(PD ~ Soil * Herbivore * Treatment * (1|Block), data = x)
  anova.qv <- Anova(model.qv)
  anova.qv <- data.table::setDT(anova.qv, keep.rownames = TRUE)[]
  r2lmer <- data.frame(Factor = c("Soil", "Herbivore", "Treatment"),
                       R2 = c(r.squaredGLMM(lmer(PD ~ Soil + (1|Block), data = x))[1],
                             r.squaredGLMM(lmer(PD ~ Herbivore + (1|Block), data = x))[1],
                             r.squaredGLMM(lmer(PD ~ Treatment + (1|Block), data = x))[1]))
  model.sum <- merge(anova.qv, r2lmer, by.x = "rn", by.y = "Factor")
  return(model.sum)
}

div.df.S <- div.df[which(div.df$Compartment == "Soil"),]
model <- lmer(PD ~ Soil * Herbivore * Treatment * (1|Block), data = div.df.S)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)
m1 <- emmeans(model, "Treatment", by = c("Soil", "Herbivore"))
pairs(m1)
model.div(div.df.S)
```

### Figure S1

```r
plot1 <- ggplot(div.df.S, aes(x = Herbivore, y = PD, fill = Herbivore)) +
                      facet_grid(Soil ~ Treatment) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")

plot2 <- ggplot(div.df.S, aes(x = Treatment, y = PD, fill = Treatment)) +
                      facet_grid(Soil ~ Herbivore) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#b3b3b3", "#8da0cb"),
                                      breaks = c("Covered", "Uncovered"),
                                      labels=c("Covered", "Uncovered"))

px <- ggpubr::ggarrange(plot1, plot2, ncol = 2, nrow = 1, align = "hv", common.legend = F)
```
<img width="1079" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/c0aac4f3-b6df-4ab7-a307-99942965bef3">

### Model, run post-hoc contrasts, and estimate variance explained for **root** samples

```r
div.df.R <- div.df[which(div.df$Compartment == "Roots"),]
model <- lmer(PD ~ Soil * Herbivore * Treatment * (1|Block), data = div.df.R)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)
m1 <- emmeans(model, "Treatment", by = c("Soil", "Herbivore"))
pairs(m1)
model.div(div.df.R)
```

### Figure S2

```r
plot1 <- ggplot(div.df.R, aes(x = Herbivore, y = PD, fill = Herbivore)) +
                      facet_grid(Soil ~ Treatment) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")

plot2 <- ggplot(div.df.R, aes(x = Treatment, y = PD, fill = Treatment)) +
                      facet_grid(Soil ~ Herbivore) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#b3b3b3", "#8da0cb"),
                                      breaks = c("Covered", "Uncovered"),
                                      labels=c("Covered", "Uncovered"),)

px <- ggpubr::ggarrange(plot1, plot2, ncol = 2, nrow = 1, align = "hv", common.legend = F)
```
<img width="1079" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/0f91a5a9-48ed-44be-b0c8-ce2f8f3acd89">

### Model, run post-hoc contrasts, and estimate variance explained for **leaf** samples

```r
div.df.L <- div.df[which(div.df$Compartment == "Leaves"),]
model <- lmer(PD ~ Soil * Herbivore * Treatment * (1|Block), data = div.df.L)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)
m1 <- emmeans(model, "Treatment", by = c("Soil", "Herbivore"))
pairs(m1)
model.div(div.df.L)
```

### Figure S3

```r
plot1 <- ggplot(div.df.L, aes(x = Herbivore, y = PD, fill = Herbivore)) +
                      facet_grid(Soil ~ Treatment) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")

plot2 <- ggplot(div.df.L, aes(x = Treatment, y = PD, fill = Treatment)) +
                      facet_grid(Soil ~ Herbivore) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#b3b3b3", "#8da0cb"),
                                      breaks = c("Covered", "Uncovered"),
                                      labels=c("Covered", "Uncovered"),)

px <- ggpubr::ggarrange(plot1, plot2, ncol = 2, nrow = 1, align = "hv", common.legend = F)
```
<img width="1079" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/a2c78f0a-ee15-43ad-adbf-8c5baf2be906">

### Model, run post-hoc contrasts, and estimate variance explained for **herbivore** samples

```r
model.div.H <- function(x){
  model.qv <- lmer(PD ~ Soil * Treatment * (1|Block), data = x)
  anova.qv <- Anova(model.qv)
  anova.qv <- data.table::setDT(anova.qv, keep.rownames = TRUE)[]
  r2lmer <- data.frame(Factor = c("Soil", "Treatment"),
                       R2 = c(r.squaredGLMM(lmer(PD ~ Soil + (1|Block), data = x))[1],
                             r.squaredGLMM(lmer(PD ~ Treatment + (1|Block), data = x))[1]))
  model.sum <- merge(anova.qv, r2lmer, by.x = "rn", by.y = "Factor")
  return(model.sum)
}

div.df.H <- div.df[which(div.df$Compartment == "Herbivore"),]
model <- lmer(PD ~ Soil * Treatment * (1|Block), data = div.df.H)
Anova(model)
model.div.H(div.df.H)
m1 <- emmeans(model, "Soil")
pairs(m1)
```

### Figure 2A

```r
plot1 <- ggplot(div.df.H, aes(x = Soil, y = PD, fill = Soil)) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "Phylogenetic diversity") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#d73027", "#fee08b", "#1a9850"),
                        breaks = c("Agricultural", "Margin", "Prairie"),
                        labels=c("Agricultural", "Margin", "Prairie"),
                                      guide = "none")
```
<img width="421" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/b9c643a8-afb1-448b-9bf5-dc9cfdec9296">

## Multivariate analyses

Code for: Tab 1, Fig. 2B, Tabs. S2-S3, Figs. S4-S6.

### PERMANOVA and posthoc contrasts for rhizosphere samples

```r
permanova.compartment <- function(psobject, compartment){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  sampledf <- data.frame(sample_data(ps))
  dist.mat <- phyloseq::distance(ps, method = "unifrac")
  perm <- how(nperm = 999)
  set.seed(100)
  pmv <- adonis2(dist.mat ~ Soil * Herbivore * Treatment, data = sampledf, permutations = perm)
  return(pmv)
}

posthoc.soilXherbivoryXTreatment <- function(psobject, compartment){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  sampledf <- data.frame(sample_data(ps))
  dist.mat <- phyloseq::distance(ps, method = "unifrac")
  pairwise.perm.manova(dist.mat, paste0(sampledf$Soil, sampledf$Herbivore, sampledf$Treatment), nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)
}

permanova.compartment(ps.16S, "Soil")
posthoc.soilXherbivoryXTreatment(ps.16S, "Soil")
```

### PERMANOVA and posthoc contrasts for root samples

```r
permanova.compartment(ps.16S, "Roots")
posthoc.soilXherbivoryXTreatment(ps.16S, "Roots")
```

### PERMANOVA and posthoc contrasts for leaf samples

```r
permanova.compartment(ps.16S, "Leaves")
posthoc.soilXherbivoryXTreatment(ps.16S, "Leaves")
```

### PERMANOVA for herbivore samples

```r
permanova.compartment.H <- function(psobject, compartment){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  sampledf <- data.frame(sample_data(ps))
  dist.mat <- phyloseq::distance(ps, method = "unifrac")
  perm <- how(nperm = 999)
  set.seed(100)
  pmv <- adonis2(dist.mat ~ Soil * Treatment, data = sampledf, permutations = perm)
  return(pmv)
}

permanova.compartment.H(ps.16S, "Herbivore")
```

### NMDS according to soil inoculum

Figure S4.

```r
nmds.soil <- function(psobject, compartment, label){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = "unifrac", formula = ~ 1)
  cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
    stat_ellipse(mapping = aes(fill = Soil),
                 alpha = 0.4,
                 geom = "polygon",
                 show.legend=T) +
    theme_bw(base_size = 15) +
    geom_point(mapping = aes(color = Soil), size = 2.5) +
    scale_colour_manual(values=c("#d73027", "#fee08b", "#1a9850"),
                        breaks = c("Agricultural", "Margin", "Prairie"),
                        labels=c("Agricultural", "Margin", "Prairie"),
                        guide = "none") +
    scale_fill_manual(name = "Legend",
                      values=c("#d73027", "#fee08b", "#1a9850"),
                      breaks = c("Agricultural", "Margin", "Prairie"),
                      labels=c("Agricultural", "Margin", "Prairie"),
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 12, family="Helvetica"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          panel.grid = element_blank()) +
    ggtitle(paste0(label))
  return(cap_plot)}

a <- nmds.soil(ps.16S, "Soil", "Rhizosphere soil")
b <- nmds.soil(ps.16S, "Roots", "Roots")
c <- nmds.soil(ps.16S, "Leaves", "Leaves")

px <- ggpubr::ggarrange(a, b, c, ncol = 3, nrow = 1, align = "hv", common.legend = T)
```
<img width="1290" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/d01712e5-3156-41db-9425-f460aff92aca">

Figure 2B.

```r
nmds.soil <- function(psobject, compartment){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = "unifrac", formula = ~ 1)
  cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
    stat_ellipse(mapping = aes(fill = Soil),
                 alpha = 0.4,
                 geom = "polygon",
                 show.legend=T) +
    theme_bw(base_size = 15) +
    geom_point(mapping = aes(color = Soil), size = 2.5) +
    scale_colour_manual(values=c("#d73027", "#fee08b", "#1a9850"),
                        breaks = c("Agricultural", "Margin", "Prairie"),
                        labels=c("Agricultural", "Margin", "Prairie"),
                        guide = "none") +
    scale_fill_manual(name = "Legend",
                      values=c("#d73027", "#fee08b", "#1a9850"),
                      breaks = c("Agricultural", "Margin", "Prairie"),
                      labels=c("Agricultural", "Margin", "Prairie"),
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.title=element_blank(),
          legend.position="none",
          legend.text = element_text(size = 12, family="Helvetica"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          panel.grid = element_blank())
  return(cap_plot)}

d <- nmds.soil(ps.16S, "Herbivore")
```
<img width="424" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/ab75cc01-00f0-4b74-ac0a-a93b58619012">

### NMDS according to herbivory

```r
nmds.herb <- function(psobject, compartment, label){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = "unifrac", formula = ~ 1)
  cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
    stat_ellipse(mapping = aes(fill = Herbivore),
                 alpha = 0.4,
                 geom = "polygon",
                 show.legend=T) +
    theme_bw(base_size = 15) +
    geom_point(mapping = aes(color = Herbivore), size = 2.5) +
    scale_colour_manual(values=c("#ff7f00", "#377eb8"),
                        breaks = c("Present", "Absent"),
                        labels=c("Present", "Absent"),
                        guide = "none") +
    scale_fill_manual(name = "Legend",
                      values=c("#ff7f00", "#377eb8"),
                      breaks = c("Present", "Absent"),
                      labels=c("Present", "Absent"),
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 12, family="Helvetica"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          panel.grid = element_blank()) +
    ggtitle(paste0(label))
  return(cap_plot)}

a <- nmds.herb(ps.16S, "Soil", "Rhizosphere soil")
b <- nmds.herb(ps.16S, "Roots", "Roots")
c <- nmds.herb(ps.16S, "Leaves", "Leaves")

px <- ggpubr::ggarrange(a, b, c, ncol = 3, nrow = 1, align = "hv", common.legend = T)
```
<img width="1290" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/833ccbfb-8eeb-403e-982a-6c39a66005b8">

### NMDS according to treatment

```r
nmds.treatment <- function(psobject, compartment, label){
  ks <- sample_data(psobject)[["Compartment"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = "unifrac", formula = ~ 1)
  cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
    stat_ellipse(mapping = aes(fill = Treatment),
                 alpha = 0.4,
                 geom = "polygon",
                 show.legend=T) +
    theme_bw(base_size = 15) +
    geom_point(mapping = aes(color = Treatment), size = 2.5) +
    scale_colour_manual(values=c("#b3b3b3", "#8da0cb"),
                        breaks = c("Covered", "Uncovered"),
                        labels=c("Covered", "Uncovered"),
                        guide = "none") +
    scale_fill_manual(name = "Legend",
                      values=c("#b3b3b3", "#8da0cb"),
                      breaks = c("Covered", "Uncovered"),
                      labels=c("Covered", "Uncovered"),
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 12, family="Helvetica"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          panel.grid = element_blank()) +
    ggtitle(paste0(label))
  return(cap_plot)}

a <- nmds.treatment(ps.16S, "Soil", "Rhizosphere soil")
b <- nmds.treatment(ps.16S, "Roots", "Roots")
c <- nmds.treatment(ps.16S, "Leaves", "Leaves")

px <- ggpubr::ggarrange(a, b, c, ncol = 3, nrow = 1, align = "hv", common.legend = T)
```
<img width="1290" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/c4d87281-3927-4c94-9c88-5c160bb5ff8c">

## MNTD

Code for: Fig. 2C, Tab. 4, and Figs. S7-S10.

```r
mntd.calc <- function(ps){
  comm <- as.data.frame(t(otu_table(ps)))
  phy <- phy_tree(ps)
  phy.dist <- cophenetic(phy)
  comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
  ks <- sample_data(ps)
  bnti <- cbind(ks, comm.sesmpd)
  return(bnti)
}

mntd.df <- mntd.calc(ps.16S)
mntd.df.H <- mntd.df[which(mntd.df$Compartment == "Herbivore"),]
mntd.df.L <- mntd.df[which(mntd.df$Compartment == "Leaves"),]
mntd.df.R <- mntd.df[which(mntd.df$Compartment == "Roots"),]
mntd.df.S <- mntd.df[which(mntd.df$Compartment == "Soil"),]
```

### Model, post-hoc, and plot for rhizosphere samples.

Tab. S4 and Fig. S7.

```r
model <- lmer(mntd.obs ~ Soil * Herbivore * Treatment * (1|Block), data = mntd.df.S)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)

b_plot <- ggplot(mntd.df.S, aes(x = Herbivore, y = mntd.obs, fill = Herbivore)) +
                      facet_grid(Treatment ~ Soil) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")
```
<img width="1205" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/84b45969-9b8e-4f6b-ba94-d8525132674c">

### Model, post-hoc, and plot for root samples.

Tab. S4 and Fig. S8.

```r
model <- lmer(mntd.obs ~ Soil * Herbivore * Treatment * (1|Block), data = mntd.df.R)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)

b_plot <- ggplot(mntd.df.R, aes(x = Herbivore, y = mntd.obs, fill = Herbivore)) +
                      facet_grid(Treatment ~ Soil) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")

```
<img width="1205" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/2c576e03-7c1c-4a09-99bd-40ad79c30a95">

Fig. S10.

```r
m1 <- emmeans(model, "Soil")
pairs(m1)
b_plot1 <- ggplot(mntd.df.R, aes(x = Soil, y = mntd.obs, fill = Soil)) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#d73027", "#fee08b", "#1a9850"),
                        breaks = c("Agricultural", "Margin", "Prairie"),
                        labels=c("Agricultural", "Margin", "Prairie"),
                                      guide = "none") 

m1 <- emmeans(model, "Treatment")
pairs(m1)

b_plot2 <- ggplot(mntd.df.R, aes(x = Treatment, y = mntd.obs, fill = Treatment)) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#b3b3b3", "#8da0cb"),
                                      breaks = c("Covered", "Uncovered"),
                                      labels=c("Covered", "Uncovered"),
                                      guide = "none") 
px <- ggpubr::ggarrange(b_plot1, b_plot2, ncol = 2, nrow = 1, align = "hv", common.legend = T)
```
<img width="1205" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/9de5f41a-3cdb-43c4-9c18-1618d3ec9495">

### Model, post-hoc, and plot for leaf samples.

Tab. S4 and Fig. S9.

```r
model <- lmer(mntd.obs ~ Soil * Herbivore * Treatment * (1|Block), data = mntd.df.L)
Anova(model)
m1 <- emmeans(model, "Herbivore", by = c("Soil", "Treatment"))
pairs(m1)

b_plot <- ggplot(mntd.df.L, aes(x = Herbivore, y = mntd.obs, fill = Herbivore)) +
                      facet_grid(Treatment ~ Soil) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#ff7f00", "#377eb8"),
                                      breaks = c("Present", "Absent"),
                                      labels=c("Present", "Absent"),
                                      guide = "none")
```
<img width="1205" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/a1256d4a-3faf-4e8f-b516-a2bab6a3ed97">

### Model, post-hoc, and plot for herbivore samples.

Tab. S4 and Fig. 2C

```r
model <- lmer(mntd.obs ~ Soil * Treatment * (1|Block), data = mntd.df.H)
Anova(model)
m1 <- emmeans(model, "Soil")
pairs(m1)

b_plot <- ggplot(mntd.df.H, aes(x = Soil, y = mntd.obs, fill = Soil)) +
                      theme_bw(base_size = 15) +
                      geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
                      labs(y = "MNTD") +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_text(color="black"),
                            axis.text.y = element_text(color="black"),
                            axis.title.y = element_text(color="black"),
                            panel.grid = element_blank(),
                            strip.background = element_blank(),
                            legend.position = "none") +
                    scale_fill_manual(values=c("#d73027", "#fee08b", "#1a9850"),
                        breaks = c("Agricultural", "Margin", "Prairie"),
                        labels=c("Agricultural", "Margin", "Prairie"),
                                      guide = "none") 
```
<img width="476" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/8a609d6f-bbcb-4835-af1f-f6fb0c3dae09">

## Upset plot

Code for: Fig. 3 and Tab. S5.

```r
library('ComplexUpset')
library('metagMisc')

upsetvec <- function(ps){
  ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
  glom <- microbiome::transform(ps, "compositional")
  dat <- psmelt(glom)
  dat2 <- dat %>% group_by(OTU) %>% dplyr::summarize(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) %>% filter(cs > 0)
  return(dat2$OTU)
}

genupset.df <- function(x){
  ps.16S <- filter_taxa(x, function (x) {sum(x > 0) > 1}, prune=TRUE)
  
  ks <- sample_data(ps.16S)[["Compartment"]] %in% "Leaves"
  ps.16S.L <- prune_samples(samples = ks, ps.16S)
  ks <- sample_data(ps.16S)[["Compartment"]] %in% "Roots"
  ps.16S.R <- prune_samples(samples = ks, ps.16S)
  ks <- sample_data(ps.16S)[["Compartment"]] %in% "Soil"
  ps.16S.S <- prune_samples(samples = ks, ps.16S)
  ks <- sample_data(ps.16S)[["Compartment"]] %in% "Herbivore"
  ps.16S.H <- prune_samples(samples = ks, ps.16S)
  
  upset.H <- upsetvec(ps.16S.H)
  upset.L <- upsetvec(ps.16S.L)
  upset.R <- upsetvec(ps.16S.R)
  upset.S <- upsetvec(ps.16S.S)
  
  rn <- Reduce(union, list(upset.H, upset.L, upset.R, upset.S))
  dat <- data.frame(herbivore = as.integer(rn %in% upset.H),
                    leaves = as.integer(rn %in% upset.L),
                    roots = as.integer(rn %in% upset.R),
                    rhizosphere = as.integer(rn %in% upset.S))
  row.names(dat) <- rn
  compartments <- colnames(dat)
  return(dat)
}


compartments <- factor(levels = c("herbivore", "leaves", "roots", "rhizosphere"))
ps.covered <- subset_samples(ps.16S, Treatment == "Covered") 
upset.df.covered <- genupset.df(ps.covered)
upset.covered <- upset(upset.df.covered, compartments, name=' ', width_ratio=0.1, sort_intersections=FALSE, sort_sets=FALSE,
    intersections=list("herbivore", "leaves", "roots", "rhizosphere",
                       c("herbivore", "leaves"), c("herbivore", "roots"), c("herbivore", "rhizosphere"), c("leaves", "roots"), c("leaves", "rhizosphere"), c("roots", "rhizosphere"),
                       c("herbivore", "leaves", "roots"), c("herbivore", "leaves", "rhizosphere"), c("leaves", "roots", "rhizosphere"), c("herbivore", "roots", "rhizosphere"),
                       c("herbivore", "leaves", "roots", "rhizosphere")
                       ),
    base_annotations = list('Intersection size'=(intersection_size()+
                                                   ylim(c(0, 4000))+ 
                                                   theme(panel.grid = element_blank(),
                                                         axis.ticks = element_line(color="black"),
                                                         axis.line = element_line(colour = "black"))+
                                                   ylab('# of ASVs - covered'))),
    set_sizes=F)


ps.uncovered <- subset_samples(ps.16S, Treatment == "Uncovered") 
upset.df.uncovered <- genupset.df(ps.uncovered)
upset.uncovered <- upset(upset.df.uncovered, compartments, name=' ', width_ratio=0.1, sort_intersections=FALSE, sort_sets=FALSE,
    intersections=list("herbivore", "leaves", "roots", "rhizosphere",
                       c("herbivore", "leaves"), c("herbivore", "roots"), c("herbivore", "rhizosphere"), c("leaves", "roots"), c("leaves", "rhizosphere"), c("roots", "rhizosphere"),
                       c("herbivore", "leaves", "roots"), c("herbivore", "leaves", "rhizosphere"), c("leaves", "roots", "rhizosphere"), c("herbivore", "roots", "rhizosphere"),
                       c("herbivore", "leaves", "roots", "rhizosphere")
                       ),
    base_annotations = list('Intersection size'=(intersection_size()+
                                                   ylim(c(0, 4000))+ 
                                                   theme(panel.grid = element_blank(),
                                                         axis.ticks = element_line(color="black"),
                                                         axis.line = element_line(colour = "black"))+
                                                   ylab('# of ASVs - uncovered'))),
    set_sizes=F)

px <- ggpubr::ggarrange(upset.covered, upset.uncovered, ncol = 1, nrow = 2, align = "hv", common.legend = T)
```
<img width="964" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/1239b98f-1b6e-46a8-91c1-9ae21365eae6">


```r
compare.upset.df <- function(df1, df2){
  intersections_data.1 = upset_data(df1,  c("herbivore", "leaves", "roots", "rhizosphere"))
  intersection_sizes.1 = unique(intersections_data.1$with_sizes[, c('intersection', 'exclusive_intersection_size')])
  intersections_data.2 = upset_data(df2,  c("herbivore", "leaves", "roots", "rhizosphere"))
  intersection_sizes.2 = unique(intersections_data.2$with_sizes[, c('intersection', 'exclusive_intersection_size')])
  compare <- merge(intersection_sizes.1, intersection_sizes.2, by = "intersection")
  colnames(compare) <- c("intersection", "g1", "g2")
  df <- compare %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      chi_sq = chisq.test(c(g1, g2))$statistic,
      p_val = chisq.test(c(g1, g2))$p.value,
      sig = ifelse(p_val<0.05, "*", " "))
  return(df)
}

compare.upset.df(upset.df.covered, upset.df.uncovered)
```

## Differentially abundant taxa

Code for: Fig. 4, Fig. 5, Tab. S6, and Fig. S11.

```r
cal.diff.taxa.cover <- function(object){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$Treatment))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", "Covered", "Uncovered"), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  return(tx)
}

cal.diff.taxa.herb <- function(object){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$Herbivore))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", "Present", "Absent"), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  return(tx)
}

plot.diff.taxa <- function(ps, title, contrast, labels){
  df.diff <- if(contrast == "cover"){cal.diff.taxa.cover(ps)} else if(contrast == "herbivory"){cal.diff.taxa.herb(ps)}
  df.diff$diffexpressed <- "no changes"
  df.diff$diffexpressed[df.diff$log2FoldChange > 0.1 & df.diff$padj < 0.05] <- if(contrast == "cover"){"Covered"} else if(contrast == "herbivory"){"Present"}
  df.diff$diffexpressed[df.diff$log2FoldChange < -0.1 & df.diff$padj < 0.05] <- if(contrast == "cover"){"Uncovered"} else if(contrast == "herbivory"){"Absent"}
  df.diff$Genus <- ifelse(df.diff$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobia(*)", df.diff$Genus)
  df.diff$Genus[df.diff$Genus==""]<-NA
  df.diff$Genus <- ifelse(is.na(df.diff$Genus) == T, df.diff$Family, df.diff$Genus)
  plot <- ggplot(data=df.diff) +
        theme_bw(base_size = 12) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed)) +
        geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = if(labels == T){ifelse(diffexpressed != "no changes", Genus,"")} else if(labels == F){""}), position = "dodge", size = 3, max.overlaps = Inf, force = 10, box.padding = 0.2, min.segment.length = 0, seed = 42) +
        scale_color_manual(values=c("blue", "red", "black")) +
        geom_vline(xintercept=0, col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = if(contrast == "cover"){"#f7fbff"} else if(contrast == "herbivory"){"#ffffcc"}, 
                                              colour = if(contrast == "cover"){"#f7fbff"} else if(contrast == "herbivory"){"#ffffcc"}, 
                                              size = 0.5, linetype = "solid"),
              legend.position="none",
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        ggtitle(paste0(title, if(contrast == "cover"){" - covered vs uncovered"} else if(contrast == "herbivory"){" - herbivory vs control"})) +
        xlab(expression(paste(Log[2], " Fold Changes"))) +
        ylab(expression(paste(-Log[10], " P"))) +
        scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c", "#000000"), labels = if(contrast == "cover"){c("Treatment", "Covered", "Uncovered")} else if(contrast == "herbivory"){c("Herbivore", "Present", "Absent")},
                           breaks = if(contrast == "cover"){c("Treatment", "Covered", "Uncovered")} else if(contrast == "herbivory"){c("Herbivore", "Present", "Absent")})
  return(plot)
}

plot1 <- plot.diff.taxa(ps.16S.S, "Rhizosphere",  "cover", T)
plot2 <- plot.diff.taxa(ps.16S.R, "Roots", "cover", F)
plot3 <- plot.diff.taxa(ps.16S.L, "Leaves", "cover", T)
plot4 <- plot.diff.taxa(ps.16S.H, "Herbivore", "cover", T)
plot5 <- plot.diff.taxa(ps.16S.S, "Rhizosphere", "herbivory", T)
plot6 <- plot.diff.taxa(ps.16S.R, "Roots",  "herbivory", T)
plot7 <- plot.diff.taxa(ps.16S.L, "Leaves", "herbivory", T)


px <- ggpubr::ggarrange(plot1, plot5, plot2, plot6, plot3, plot7, plot4, ncol = 2, nrow = 4, align = "hv", common.legend = T)
```
<img width="826" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/1c17dd0e-c009-4893-b683-fc261788991b">


```r
cal.diff.taxa.soil <- function(object, g1, g2){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$Soil))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", g1, g2), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  return(tx)
}


plot.diff.taxa.soil <- function(ps, title, g1, g2, labels){
  df.diff <- cal.diff.taxa.soil(ps, g1, g2)
  df.diff$diffexpressed <- "no changes"
  df.diff$diffexpressed[df.diff$log2FoldChange > 0.1 & df.diff$padj < 0.05] <- paste0(g1)
  df.diff$diffexpressed[df.diff$log2FoldChange < -0.1 & df.diff$padj < 0.05] <- paste0(g2)
  df.diff$Genus <- ifelse(df.diff$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobia(*)", df.diff$Genus)
  df.diff$Genus[df.diff$Genus==""]<-NA
  df.diff$Genus <- ifelse(is.na(df.diff$Genus) == T, df.diff$Family, df.diff$Genus)
  plot <- ggplot(data=df.diff) +
        theme_bw(base_size = 12) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed)) +
        geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = if(labels == T){ifelse(diffexpressed != "no changes", Genus,"")} else if(labels == F){""}), position = "dodge", size = 3, max.overlaps = Inf, force = 10, box.padding = 0.2, min.segment.length = 0, seed = 42) +
        scale_color_manual(values=c("blue", "red", "black")) +
        geom_vline(xintercept=0, col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", 
                                              colour = "white", 
                                              size = 0.5, linetype = "solid"),
              legend.position="none",
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        ggtitle(paste0(title, " - ", g1, " vs ", g2)) +
        xlab(expression(paste(Log[2], " Fold Changes"))) +
        ylab(expression(paste(-Log[10], " P"))) +
        scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c", "#000000"), labels = c("Soil", g1, g2), breaks =c("Soil", g1, g2))
  return(plot)
}

plot1 <- plot.diff.taxa.soil(ps.16S.S, "Rhizosphere", "Agricultural", "Margin", T)
plot2 <- plot.diff.taxa.soil(ps.16S.S, "Rhizosphere", "Agricultural", "Prairie", T)
plot3 <- plot.diff.taxa.soil(ps.16S.S, "Rhizosphere", "Margin", "Prairie", T)

plot4 <- plot.diff.taxa.soil(ps.16S.R, "Roots", "Agricultural", "Margin", F)
plot5 <- plot.diff.taxa.soil(ps.16S.R, "Roots", "Agricultural", "Prairie", F)
plot6 <- plot.diff.taxa.soil(ps.16S.R, "Roots", "Margin", "Prairie", F)

plot7 <- plot.diff.taxa.soil(ps.16S.L, "Leaves", "Agricultural", "Margin", T)
plot8 <- plot.diff.taxa.soil(ps.16S.L, "Leaves", "Agricultural", "Prairie", T)
plot9 <- plot.diff.taxa.soil(ps.16S.L, "Leaves", "Margin", "Prairie", T)

plot10 <- plot.diff.taxa.soil(ps.16S.H, "Herbivore", "Agricultural", "Margin", T)
plot11 <- plot.diff.taxa.soil(ps.16S.H, "Herbivore", "Agricultural", "Prairie", T)
plot12 <- plot.diff.taxa.soil(ps.16S.H, "Herbivore", "Margin", "Prairie", T)

px <- ggpubr::ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, ncol = 3, nrow = 4, align = "hv", common.legend = F)
```
<img width="1133" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/1a2e1bcb-77d9-4d19-9921-2ce73ff9d053">


```r
table.diff.taxa <- function(ps, g1, g2){
  df.diff <- cal.diff.taxa.soil(ps, g1, g2)
  df.diff <- df.diff[which(df.diff$padj < 0.05),]
  df.diff$diffexpressed[df.diff$log2FoldChange > 0.1 & df.diff$padj < 0.05] <- paste0(g1)
  df.diff$diffexpressed[df.diff$log2FoldChange < -0.1 & df.diff$padj < 0.05] <- paste0(g2)
  return(df.diff)
}

table.diff.taxa(ps.16S.R, "Agricultural", "Margin") %>% write.table("diffasv_AM.txt", sep = "\t")
table.diff.taxa(ps.16S.R, "Agricultural", "Prairie") %>% write.table("diffasv_AP.txt", sep = "\t")
table.diff.taxa(ps.16S.R, "Margin", "Prairie") %>% write.table("diffasv_MP.txt", sep = "\t")
```


```r
library("ggvenn")
a <- table.diff.taxa(ps.16S.R, "Agricultural", "Margin")$rn
b <- table.diff.taxa(ps.16S.R, "Agricultural", "Prairie")$rn
c <- table.diff.taxa(ps.16S.R, "Margin", "Prairie")$rn

list.venn <- list(AvsM = a,
          AvsP = b,
          MvsP = c)

ggvenn(list.venn) 
```
<img width="420" alt="image" src="https://github.com/amalacrino/malacrino_and_bennett_CommsBio/assets/21124426/7c5fe4d7-3f74-4763-b21f-9a0540ca02b4">

## SEM

Code for: Fig. 6.

```r
library("tidyverse")
library("piecewiseSEM")

sem.dataset <- read.table('data/biomass.txt', sep = "\t", header = T)

mds.calc <- function(ps, y){
  sampledf <- data.frame(sample_data(ps.16S.L)) %>% tibble::rownames_to_column("sample")
  dist.mat <- phyloseq::distance(ps.16S.L, method = "unifrac")
  cap_ord <- ordinate(physeq = ps.16S.L, method = "NMDS", distance = dist.mat, formula = ~ 1)$points %>% as.data.frame %>% tibble::rownames_to_column("sample")
  cap_ord <- merge(sampledf, cap_ord, by = "sample")[,c(3,9)]
  colnames(cap_ord) <- c("plantID", paste0("MDS1_", y))
  return(cap_ord)
}

nmds.H <- mds.calc(ps.16S.H, "herbivore")
nmds.L <- mds.calc(ps.16S.L, "leaves")
nmds.R <- mds.calc(ps.16S.R, "roots")
nmds.S <- mds.calc(ps.16S.S, "soil")
df_list <- list(nmds.H, nmds.L, nmds.R, nmds.S)
nmds.all <- df_list %>% reduce(full_join, by='plantID')
sem.dataset <- merge(sem.dataset, nmds.all, by = "plantID")
sem.dataset <- as_tibble(sem.dataset)
sem.dataset$soil <- as.numeric(factor(sem.dataset$soil))
sem.dataset$herbivore <- as.numeric(factor(sem.dataset$herbivore))
sem.dataset$cover <- as.numeric(factor(sem.dataset$treatment))

modelList <- psem(
  lmer(shoot_biomass ~ cover + soil + herbivore + MDS1_leaves + MDS1_soil + MDS1_roots + (1|Block), sem.dataset),
  lmer(root_biomass ~ cover + soil + herbivore + MDS1_roots + MDS1_soil + (1|Block), sem.dataset),
  lmer(insect_biomass ~  cover + soil + MDS1_leaves + MDS1_soil + MDS1_roots + MDS1_herbivore + (1|Block), sem.dataset),
  lmer(MDS1_soil ~ cover + soil + herbivore + (1|Block), sem.dataset),
  lmer(MDS1_roots ~ MDS1_soil + MDS1_leaves + cover + soil + herbivore + (1|Block), sem.dataset),
  lmer(MDS1_leaves ~ MDS1_soil + MDS1_roots + cover + soil + herbivore + (1|Block), sem.dataset),
  lmer(MDS1_herbivore ~ MDS1_soil + MDS1_roots + MDS1_leaves + cover + soil + (1|Block), sem.dataset)
)

a <- coefs(modelList)
plot(modelList, alpha = 0.05)
```
