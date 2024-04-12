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



```r

```



```r

```



```r

```



```r

```



```r

```



```r

```



```r

```



```r

```






