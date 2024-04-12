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


## Data analysis

### Load libraries

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

### Load and clean microbiome data

```r
NTHREADS = 8

'%ni%' <- Negate('%in%')

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))))
}

asv.table <- as.matrix(read.table('data/ASV_table.tsv', sep = "\t", header = T, row.names = 1))
tax.table <- as.matrix(read.table('data/ASV_tax_species.tsv', sep = "\t", header = T, row.names = 1)[,c(2:8)])
mapping <- read.table('data/metadata.txt', sep = "\t", header = T, row.names = 1)
tree <- read.tree('data/tree.tre')

ps.16S <- phyloseq(otu_table(asv.table, taxa_are_rows = T), 
                         sample_data(mapping), 
                         tax_table(tax.table), 
                         phy_tree(tree))

ps.16S.inoculum <- subset_samples(ps.16S, Treatment == "Inoculum")

ps.16S <- subset_samples(ps.16S, Treatment != "Herbivore_artif_food") 
ps.16S <- subset_samples(ps.16S, Treatment != "Inoculum") 
ps.16S <- subset_samples(ps.16S, Treatment != "Sterile_soil") 

ps.16S <- subset_taxa(ps.16S, Class !="Chloroplast")
ps.16S <- subset_taxa(ps.16S, Order !="Chloroplast")
ps.16S <- subset_taxa(ps.16S, Family !="Mitochondria")

remove.cont <- function(ps){
  sample_data(ps)$is.neg <- sample_data(ps)$Treatment == "DNA_extraction"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.05)
  cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
  cont.remove <- row.names(cont.remove)
  allTaxa = taxa_names(ps)
  allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
  ps <-  prune_taxa(allTaxa, ps)
  temp <- sample_data(ps)[["Treatment"]] %ni% "DNA_extraction"
  ps <- prune_samples(samples = temp, ps)
  return(ps)
}

wrench.norm.ps <- function(ps){
  count_tab <- as.matrix((data.frame(otu_table(ps))))
  group <- sample_data(ps)$Compartment
  W <- wrench(count_tab, condition=group)
  norm_factors <- W$nf
  head(norm_factors)
  norm_counts <- sweep(count_tab, 2, norm_factors, FUN = '/')
  norm_counts_trim <- data.frame(t(data.frame(norm_counts)))                                                  
  norm_counts_trim[] <- lapply(norm_counts_trim, function(x) DescTools::Winsorize(x, probs = c(0, 0.97), type = 1))
  norm_counts_trim <- data.frame(t(norm_counts_trim))
  norm_counts_trim[norm_counts_trim == 0] <- 1
  norm_counts_trim[norm_counts_trim < 0] <- 1
  norm_counts_trim <- log2(norm_counts_trim)
  colnames(norm_counts_trim) <- gsub("\\.", "-", colnames(norm_counts_trim))
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(norm_counts_trim, taxa_are_rows =  TRUE)
  return(ps_norm)
}

ps.16S <- subset_samples(ps.16S, Compartment != "Control") 
ps.16S <- subset_samples(ps.16S, Soil != "No_treatment") 
ps.16S <- filter_taxa(ps.16S, function (x) {sum(x > 0) > 1}, prune=TRUE)
```

