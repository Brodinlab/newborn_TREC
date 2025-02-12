library(dplyr)
library(ggplot2)
library(readr)
library(rstatix)
library(wesanderson)
library(ggpubr)

my_pal <- wes_palette("Moonrise3")
theme_set(
    theme_pubr() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_blank()
        )
)

data_parents_raw <- read_csv("data/TREC_KREC_parents.csv")
data_parents <- data_parents_raw %>%
    filter(type == "Mother") %>%
    rename(
        SNP = `GENOTYPE SNP rs2204985`,
        log10TREC = `LOG10 sjTREC`
    ) %>%
    select(child_id, SNP, log10TREC)

raw <- read_csv("data/TREC_KREC_combine_Baby_4.csv")
data <- raw %>%
    filter(!is.na(age_bin)) %>%
    mutate(age_bin = factor(age_bin,
        levels = c(
            "CB", "0-3 Days", "4-7 Days", "8-30 Days", "31-57 Days",
            "58-78 Days", "79-87 Days", "88-100 Days", "101-120 Days",
            "121-150 Days", "151-200 Days", "201-400 Days", ">400 Days"
        )
    )) %>%
    mutate(group_SNP = if_else(group == "Preterm", paste0(group, "_", SNP_merge), group))

# suppl fig 2a
data %>%
    filter(baby_age > 100) %>%
    ggplot(aes(x = group, y = log10TREC, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    stat_compare_means(comparisons = list(c("Preterm", "Term")), method = "wilcox.test") +
    scale_fill_manual(values = my_pal) +
    ggtitle("baby_age > 100")
ggsave("figures/suppl/suppl fig2a.pdf", width = 8, height = 8)

# suppl fig 2b
data %>%
    filter(age_bin %in% c("0-3 Days", "4-7 Days")) %>%
    ggplot(aes(x = group_SNP, y = log10TREC, fill = group_SNP)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    ylim(2, 4.25) +
    stat_compare_means(comparisons = list(c("Preterm_AA", "Preterm_GA/GG"), c("Preterm_GA/GG", "Term"), c("Preterm_AA", "Term")), method = "wilcox.test") +
    scale_fill_manual(values = my_pal)
ggsave("figures/suppl/suppl fig2b.pdf", width = 8, height = 8)

# suppl fig 2c
data_parents %>%
    ggplot(aes(x = SNP, y = log10TREC, fill = SNP)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    stat_compare_means(comparisons = list(c("AA", "GA"), c("GA", "GG"), c("AA", "GG")), method = "wilcox.test") +
    scale_fill_manual(values = my_pal)
ggsave("figures/suppl/suppl fig2c.pdf", width = 8, height = 8)