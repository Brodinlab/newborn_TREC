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
    theme_pubr() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_blank()
        )
)

# Fig 3a
ggplot(data, aes(x = age_bin, y = log10KREC)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/Fig3a.pdf", width = 8, height = 8)

# Fig 3b
ggplot(data, aes(x = age_bin, y = `LOG10 CjInt`)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/Fig3b.pdf", width = 8, height = 8)

# Fig 3c
ggplot(data, aes(x = age_bin, y = `B cell Div`)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/Fig3c.pdf", width = 8, height = 8)
