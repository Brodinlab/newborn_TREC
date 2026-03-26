library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(lme4)
library(rstatix)

data <- read_csv("data/TREC_KREC_combine_Baby_4.csv") %>%
    filter(type == "Baby", baby_age > 0) %>% # remove cord blood
    mutate(log10age = log10(baby_age))

m <- glm(log10KREC ~ log10age + SNP_merge + sex + mode_delivery + group, data = data)

# Function to perform bootstrapped GLM
bootstrap_glm <- function(data, formula, family = gaussian(), n_bootstrap = 100, sample_fraction = 0.9) {
    # Initialize storage for bootstrap results
    m <- glm(formula, data = data, family = family)
    bootstrap_coefs <- matrix(NA,
        nrow = n_bootstrap,
        ncol = length(coef(m))
    )
    colnames(bootstrap_coefs) <- names(coef(m))

    # Set random seed for reproducibility
    set.seed(123)

    # Perform bootstrap
    for (i in 1:n_bootstrap) {
        # Sample 90% of data
        bootstrap_sample <- data[sample(nrow(data),
            size = floor(sample_fraction * nrow(data)),
            replace = FALSE
        ), ]

        # Fit GLM
        model <- glm(formula, data = bootstrap_sample, family = family)

        # Store coefficients
        bootstrap_coefs[i, ] <- coef(model)
    }

    # Compute summary statistics
    bootstrap_summary <- data.frame(
        coefficient = colMeans(bootstrap_coefs),
        std_error = apply(bootstrap_coefs, 2, sd),
        lower_ci = apply(bootstrap_coefs, 2, quantile, probs = 0.025),
        upper_ci = apply(bootstrap_coefs, 2, quantile, probs = 0.975)
    )

    return(list(
        bootstrap_summary = bootstrap_summary,
        all_bootstraps = bootstrap_coefs
    ))
}

result <- bootstrap_glm(data,
                        formula = log10KREC ~ log10age + SNP_merge + sex + mode_delivery + group,
                        family = gaussian(),
                        n_bootstrap = 500,
                        sample_fraction = 0.9)

result$bootstrap_summary %>%
    rownames_to_column("var_name") %>%
    arrange(desc(coefficient)) %>%
    mutate(
        var_name = factor(var_name,
            levels = unique(var_name)
        )
    ) %>%
    filter(var_name != "(Intercept)") %>%
    ggplot(aes(x = var_name, y = coefficient)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
        width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
        title = "Bootstrap Coefficient Estimates",
        y = "Coefficient Value",
        x = "Predictor"
    ) +
    coord_flip() +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_blank()
    )

rstatix::Anova(m) # add p values to the plot according to the ANOVA table

ggsave("figures/Fig3/3d.pdf", width = 8, height = 8)
