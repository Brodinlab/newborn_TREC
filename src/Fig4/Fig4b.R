library(dplyr)
library(ggplot2)
library(lme4)
library(rstatix)
library(tibble)

data <- readxl::read_excel("data/sjTREC_150kcells_21112024.xlsx") %>%
    mutate(
        log10age = log10(Age),
        SNP_merged = if_else(SNP == "AA", SNP, "GA/GG")
    )

m <- glm(sjTREC ~ log10age + SNP_merged + Gender, data = data)

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
                        formula = sjTREC ~ log10age + SNP_merged + Gender,
                        family = gaussian(),
                        n_bootstrap = 500,
                        sample_fraction = 0.9)

# Get p-values from Anova test
anova_res <- rstatix::Anova(m)
p_values <- anova_res$`Pr(>Chisq)`
names(p_values) <- c("log10age", "SNP_mergedGA/GG", "GenderMale")

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
    ylim(-0.2, 0.5) +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_blank()
    ) +
    # Add p-values as text annotations
    geom_text(
        aes(
            x = var_name,
            y = upper_ci + 0.02,
            label = sprintf("p = %.3f", p_values[as.character(var_name)])
        ),
        hjust = 0
    )

ggsave("figures/Fig4/4b.pdf", width = 8, height = 8)
