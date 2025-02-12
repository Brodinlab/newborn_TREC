library(dplyr)
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
            "CB", "0-3 Days", "4-7 Days", "8-30 Days", "31-57 Days",
            "58-78 Days", "79-87 Days", "88-100 Days", "101-120 Days",
            "121-150 Days", "151-200 Days", "201-400 Days", ">400 Days"
        )
    ))

my_pal <- wes_palette("Moonrise3")
theme_set(
    theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_blank()
        )
)

# Fig 1b
ggplot(data, aes(x = age_bin, y = log10TREC)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig1/Fig1b.pdf", width = 8, height = 8)

# Fig 1c
data %>%
    filter(age_bin %in% c("0-3 Days", "4-7 Days")) %>%
    ggplot(aes(x = group, y = log10TREC, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    ylim(2, 4.25) +
    stat_compare_means(comparisons = list(c("Preterm", "Term"))) +
    scale_fill_manual(values = my_pal)
ggsave("figures/Fig1/Fig1c.pdf", width = 8, height = 8)

# Fig 1d
data %>%
    filter(age_bin != "CB") %>% # remove CB
    ggplot(aes(x = log10(baby_age), y = log10TREC, color = SNP)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = my_pal)
ggsave("figures/Fig1/Fig1d.pdf", width = 8, height = 8)
