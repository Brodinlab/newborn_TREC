library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(wesanderson)

raw <- read_csv("data/TREC_KREC_combine_Baby_4.csv")
data <- raw %>%
    filter(!is.na(age_bin)) %>%
    mutate(age_bin = factor(age_bin,
        levels = c(
            "CB", "0-3 Days", "4-14 Days", "15-30 Days", "31-57 Days",
            "58-78 Days", "79-87 Days", "88-100 Days", "101-120 Days",
            "121-150 Days", "151-200 Days", "201-400 Days", ">400 Days"
        )
    ))

# data %>% group_by(age_bin) %>% tally()
# data %>% select(child_id, mode_delivery, group) %>% distinct() %>% group_by(mode_delivery, group) %>% tally()

my_pal <- wes_palette("Moonrise3")
theme_set(
    theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_blank()
        )
)

data %>%
    arrange(child_id, age_bin) %>%
    ggplot(aes(x = age_bin, y = log10TREC)) +
    geom_boxplot(outlier.shape = NA, fill = "grey") +
    geom_line(aes(group = child_id), position = position_jitter(width = 0.2, seed = 123), alpha = 0.2) +
    geom_point(aes(group = child_id), position = position_jitter(width = 0.2, seed = 123)) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("figures/Fig1/1b.pdf", width = 8, height = 8)

data %>%
    filter(baby_age <= 7 & age_bin != "CB") %>% # first week, no cord blood included. 
    ggplot(aes(x = group, y = log10TREC, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    ylim(2, 4.25) +
    stat_compare_means(comparisons = list(c("Preterm", "Term"))) +
    scale_fill_manual(values = my_pal)

ggsave("figures/Fig1/1c.pdf", width = 8, height = 8)

tmp <- data %>%
    select(child_id, age_bin, log10TREC) %>%
    filter(age_bin %in% c(
        "CB", "0-3 Days", "88-100 Days", "101-120 Days",
        "121-150 Days", "151-200 Days", "201-400 Days", ">400 Days"
    )) %>%
    pivot_wider(
        names_from = age_bin,
        values_from = log10TREC,
        values_fill = NA,
        values_fn = first
    ) %>%
    mutate(
        birth = if_else(is.na(`0-3 Days`), CB, `0-3 Days`),
        later = coalesce(`88-100 Days`, `101-120 Days`, `121-150 Days`, `151-200 Days`, `201-400 Days`, `>400 Days`)
    ) %>%
    mutate(delta = later - birth) %>%
    left_join(data %>% select(child_id, SNP) %>% distinct(), by = "child_id")

stat.test.birth <- tmp %>%
    filter(!is.na(birth)) %>%
    wilcox_test(birth ~ SNP, comparisons = list(c("AA", "GA"), c("AA", "GG"), c("GA", "GG"))) %>%
    adjust_pvalue(method = "fdr") %>%
    add_xy_position(x = "SNP")

tmp %>%
    filter(!is.na(birth)) %>%
    ggplot(aes(x = SNP, y = birth)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    stat_pvalue_manual(stat.test.birth, label = "p.adj") +
    theme_pubr() +
    ggtitle("birth") +
    labs(subtitle = paste0(
        "AA: n=", sum(tmp$SNP == "AA" & !is.na(tmp$birth)),
        ", GA: n=", sum(tmp$SNP == "GA" & !is.na(tmp$birth)),
        ", GG: n=", sum(tmp$SNP == "GG" & !is.na(tmp$birth))
    ))
ggsave("figures/Fig1/1d.pdf", width = 8, height = 8)


stat.test.later <- tmp %>%
    filter(!is.na(later)) %>%
    wilcox_test(later ~ SNP, comparisons = list(c("AA", "GA"), c("AA", "GG"), c("GA", "GG"))) %>%
    adjust_pvalue(method = "fdr") %>%
    add_xy_position(x = "SNP")

tmp %>%
    filter(!is.na(later)) %>%
    ggplot(aes(x = SNP, y = later)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    stat_pvalue_manual(stat.test.later, label = "p.adj") +
    theme_pubr() +
    ggtitle("later") +
    labs(subtitle = paste0(
        "AA: n=", sum(tmp$SNP == "AA" & !is.na(tmp$later)),
        ", GA: n=", sum(tmp$SNP == "GA" & !is.na(tmp$later)), 
        ", GG: n=", sum(tmp$SNP == "GG" & !is.na(tmp$later))
    ))
ggsave("figures/Fig1/1e.pdf", width = 8, height = 8)

stat.test.delta <- tmp %>%
    filter(!is.na(delta)) %>%
    wilcox_test(delta ~ SNP, comparisons = list(c("AA", "GA"), c("AA", "GG"), c("GA", "GG"))) %>%
    adjust_pvalue(method = "fdr") %>%
    add_xy_position(x = "SNP")

tmp %>%
    filter(!is.na(delta)) %>%
    ggplot(aes(x = SNP, y = delta)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    stat_pvalue_manual(stat.test.delta, label = "p.adj") +
    theme_pubr() +
    ggtitle("delta") +
    labs(subtitle = paste0(
        "AA: n=", sum(tmp$SNP == "AA" & !is.na(tmp$delta)),
        ", GA: n=", sum(tmp$SNP == "GA" & !is.na(tmp$delta)),
        ", GG: n=", sum(tmp$SNP == "GG" & !is.na(tmp$delta))
    ))

ggsave("figures/Fig1/1f.pdf", width = 8, height = 8)

# SNP-group sample sizes for the legend (unique subjects)
snp_levels <- c("AA", "GA", "GG")
snp_n <- data %>%
    filter(baby_age != 0, !is.na(SNP)) %>% # remove CB
    distinct(child_id, SNP) %>%
    count(SNP, name = "n")
snp_n_vec <- snp_n$n[match(snp_levels, snp_n$SNP)]
snp_n_vec[is.na(snp_n_vec)] <- 0
snp_labs <- setNames(paste0(snp_levels, "(n=", snp_n_vec, ")"), snp_levels)

data %>%
    filter(baby_age != 0, !is.na(SNP)) %>%
    ggplot(aes(x = log10(baby_age), y = log10TREC, color = SNP)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(
        values = setNames(my_pal[seq_along(snp_levels)], snp_levels),
        breaks = snp_levels,
        labels = snp_labs
    ) +
    theme_pubr()
ggsave("figures/Fig1/1g.pdf", width = 8, height = 8)
