library(dplyr)
library(ggpubr)
library(rstatix)

data <- readxl::read_excel("data/sjTREC_150kcells_21112024.xlsx") %>%
    mutate(log10age = log10(Age))

# Fig 4c
data %>%
    ggplot(aes(x = SNP, y = sjTREC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.5) +
    theme_pubr() +
    ylab("log10 TREC") +
    stat_compare_means(
        comparisons = list(
            c("AA", "GA"),
            c("GA", "GG"),
            c("AA", "GG")
        ),
        method = "wilcox.test"
    )
ggsave("figures/Fig4/Fig4c.pdf", width = 8, height = 6)

# Fig 4d
ggplot(data, aes(x = log10age, y = sjTREC)) +
    geom_point() +
    stat_cor(method = "spearman", label.x = 0.5) +
    theme_pubr() +
    ylab("log10 TREC")

ggsave("figures/Fig4/Fig4d.pdf", width = 8, height = 6)
