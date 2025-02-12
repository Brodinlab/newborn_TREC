library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)

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
    mutate(group_SNP = if_else(group == "Preterm", paste0(group, "_", SNP), group))

df <- left_join(
    data %>%
        filter(week == 1) %>%
        select(child_id, log10TREC, SNP, group) %>%
        rename(w1_TREC = log10TREC),
    data %>%
        filter(week == 12) %>%
        select(child_id, log10TREC) %>%
        rename(w12_TREC = log10TREC),
    by = "child_id"
) %>%
    left_join(
        data_parents %>% select(-SNP) %>% rename(mother_TREC = log10TREC),
        by = "child_id"
    )

# fig 1f
ggplot(df, aes(x = mother_TREC, y = w1_TREC)) +
    geom_point() +
    stat_cor(
        label.x = 2,
        label.x.npc = "right",
        method = "spearman"
    ) +
    annotate(
        geom = "text", x = 2, y = Inf, vjust = 2, hjust = 0, label = "Spearman correlation"
    ) +
    theme_pubr()
ggsave("figures/Fig1/Fig1f.pdf", width = 8, height = 8)

# fig 1g
ggplot(df, aes(x = mother_TREC, y = w12_TREC)) +
    geom_point() +
    stat_cor(
        label.x = 2,
        label.x.npc = "right",
        method = "spearman"
    ) +
    annotate(
        geom = "text", x = 2, y = Inf, vjust = 2, hjust = 0, label = "Spearman correlation"
    ) +
    theme_pubr()
ggsave("figures/Fig1/Fig1g.pdf", width = 8, height = 8)

# fig 1h
ggplot(df, aes(x = w1_TREC, y = w12_TREC)) +
    geom_point() +
    stat_cor(
        label.x = 3,
        label.x.npc = "right",
        method = "spearman"
    ) +
    annotate(
        geom = "text", x = 3, y = Inf, vjust = 2, hjust = 0, label = "Spearman correlation"
    ) +
    theme_pubr()
ggsave("figures/Fig1/Fig1h.pdf", width = 8, height = 8)
