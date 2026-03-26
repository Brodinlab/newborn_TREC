library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(wesanderson)

raw <- read_csv("data/TREC_KREC_combine_Baby_4.csv")
data <- raw %>%
    mutate(
        age_bin_B = case_when(
            baby_age == 0 ~ "CB",
            baby_age >= 1 & baby_age <= 3 ~ "0-3 Days",
            baby_age >= 4 & baby_age <= 7 ~ "4-7 Days",
            baby_age >= 8 & baby_age <= 30 ~ "8-30 Days",
            baby_age >= 31 & baby_age <= 57 ~ "31-57 Days",
            baby_age >= 58 & baby_age <= 78 ~ "58-78 Days",
            baby_age >= 79 & baby_age <= 87 ~ "79-87 Days",
            baby_age >= 88 & baby_age <= 100 ~ "88-100 Days",
            baby_age >= 101 & baby_age <= 120 ~ "101-120 Days",
            baby_age >= 121 & baby_age <= 150 ~ "121-150 Days",
            baby_age >= 151 & baby_age <= 200 ~ "151-200 Days",
            baby_age >= 201 & baby_age <= 400 ~ "201-400 Days",
            baby_age > 400 ~ ">400 Days",
            TRUE ~ NA_character_
        ),
        age_bin_B = factor(
            age_bin_B,
            levels = c(
                "CB",
                "0-3 Days",
                "4-7 Days",
                "8-30 Days",
                "31-57 Days",
                "58-78 Days",
                "79-87 Days",
                "88-100 Days",
                "101-120 Days",
                "121-150 Days",
                "151-200 Days",
                "201-400 Days",
                ">400 Days"
            )
        )
    ) %>%
    filter(!is.na(age_bin_B))

my_pal <- wes_palette("Moonrise3")
theme_set(
    theme_pubr() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_blank()
        )
)

# Fig 3a
ggplot(data, aes(x = age_bin_B, y = log10KREC)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/3a.pdf", width = 8, height = 8)

# Fig 3b
ggplot(data, aes(x = age_bin_B, y = `LOG10 CjInt`)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/3b.pdf", width = 8, height = 8)

# Fig 3c
ggplot(data, aes(x = age_bin_B, y = `B cell Div`)) +
    geom_boxplot(outlier.shape = NA, fill = my_pal[[1]]) +
    geom_jitter()
ggsave("figures/Fig3/3c.pdf", width = 8, height = 8)
