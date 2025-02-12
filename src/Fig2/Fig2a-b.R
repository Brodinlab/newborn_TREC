library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(EnhancedVolcano)
library(rstatix)
library(tibble)

data <- read_csv("data/TREC_KREC_combine_Baby_4.csv")
subtype_names <- colnames(data)[23:36]
protein_names <- colnames(data)[37:124]

data %>%
    filter(!is.na(Neutrophils)) %>%
    dim() # check the number of samples we have both data

# Fig 2a
res_subtype <- lapply(subtype_names, function(name) {
    cor.test(as.formula(paste0("~ log10TREC + ", name)), data, use = "complete.obs", method = "spearman")
})
df_subtype <- data.frame(
    r = sapply(res_subtype, function(x) x$estimate),
    p = sapply(res_subtype, function(x) x$p.value),
    subtype = subtype_names
) %>% adjust_pvalue(method = "BH")
EnhancedVolcano(df_subtype,
    x = "r", y = "p.adj",
    lab = df_subtype$subtype, xlim = c(-0.75, 0.75),
    pCutoff = 0.05, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "TREC and subpopulation frequencies (CyTOF)",
    subtitle = ""
) + xlab("Spearman correlation")
ggsave("figures/Fig2/Fig2a.pdf", height = 8, width = 8)

# Fig 2b
na_count <- data %>%
    select(all_of(protein_names)) %>%
    is.na() %>%
    colSums()
hist(na_count / dim(data)[1])
protein_names_qc <- protein_names[na_count / dim(data)[1] <= 0.9]
res_protein <- lapply(protein_names_qc, function(name) {
    cor.test(as.formula(paste0("~ log10TREC + `", name, "`")), data, use = "complete.obs", method = "spearman")
})
df_protein <- data.frame(
    r = sapply(res_protein, function(x) x$estimate),
    p = sapply(res_protein, function(x) x$p.value),
    name = protein_names_qc
) %>% adjust_pvalue(method = "BH")

EnhancedVolcano(df_protein,
    x = "r", y = "p.adj",
    lab = df_protein$name, xlim = c(-0.75, 0.75),
    pCutoff = 0.05, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "TREC and protein profiles (Olink)",
    subtitle = "",
    arrowheads = FALSE
) + xlab("Spearman correlation")
ggsave("figures/Fig2/Fig2b.pdf", width = 8, height = 8)
