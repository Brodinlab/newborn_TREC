library(dplyr)
library(tibble)
library(readr)
library(lme4)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(patchwork) 

df <- read_csv("data/TREC_KREC_combine_Baby_4.csv")

df$SNP <- factor(df$SNP)
df$SNP_merge <- factor(df$SNP_merge)
df$group <- factor(df$group)
df$sex <- factor(df$sex)
df$steroids <- factor(df$steroids)
df$mode_delivery<- factor(df$mode_delivery)
df$age_bin<- factor(df$age_bin)
df$type<- factor(df$type)
df$log10age <- log10(df$baby_age)
df$age_bin <- factor(df$age_bin, levels = c("CB", "0-3 Days", "4-14 Days", "15-30 Days", 
                                            "31-57 Days", "58-78 Days", "79-87 Days", 
                                            "88-100 Days", "101-120 Days", "121-150 Days", 
                                            "151-200 Days", "201-400 Days", ">400 Days"))

#Supp 1b
colnames(df)
mass_spectrometry <- c("B.cells","Basophils","CD4T","CD8T","gdT","Monocytes","Neutrophils","NK.cells","pDC" ,"mDC","Tregs","Eosinophils","Platelets","MSC")
olink <- c(
  "IL8", "VEGFA", "CD8A", "MCP-3", "GDNF", "CDCP1",
  "CD244", "IL7", "OPG", "LAP TGF-beta-1", "uPA", "IL6",
  "IL-17C", "MCP-1", "IL-17A", "CXCL11", "AXIN1", "TRAIL",
  "CXCL9", "CST5", "OSM", "CXCL1", "CCL4", "CD6",
  "SCF", "IL18", "TGF-alpha", "MCP-4", "CCL11", "TNFSF14",
  "FGF-23", "IL-10RA", "MMP-1", "LIF-R", "FGF-21", "CCL19",
  "IL-15RA", "IL-10RB", "IL-22 RA1", "IL-18R1", "PD-L1", "CXCL5",
  "TRANCE", "HGF", "IL-12B", "MMP-10", "IL10", "TNF",
  "CCL23", "CD5", "CCL3", "Flt3L", "CXCL6", "CXCL10",
  "4E-BP1", "SIRT2", "CCL28", "DNER", "EN-RAGE", "CD40",
  "IFN-gamma", "FGF-19", "MCP-2", "CASP-8", "CCL25", "CX3CL1",
  "TNFRSF9", "NT-3", "TWEAK", "CCL20", "ST1A1", "STAMBP",
  "ADA", "TNFB", "CSF-1", "SLAMF1", "FGF-5", "Beta-NGF",
  "STAMPB", "IL-2RB", "ARTN", "IL-20RA", "IL-1 alpha", "TSLP",
  "IL13", "IL4", "LIF", "IL5"
)

all_vars <- c(mass_spectrometry, olink)

df$mass_spectrometry <- rowSums(!is.na(df[, intersect(names(df), mass_spectrometry)])) > 0
df$olink <- rowSums(!is.na(df[, intersect(names(df), olink)])) > 0

df$weekcb <- as.character(df$week)
df$weekcb[df$age_bin == "CB"] <- "CB"

semaines_numeriques <- sort(as.numeric(na.omit(as.character(df$weekcb[df$weekcb != "CB"]))))
niveaux <- unique(c("CB", as.character(semaines_numeriques)))
df$weekcb <- factor(df$weekcb, levels = niveaux)

df <- df %>%
  mutate(
    assay = case_when(
      !is.na(log10TREC) & olink == TRUE & mass_spectrometry == TRUE ~ "TREC_olink_MS",
      !is.na(log10TREC) & olink == TRUE & mass_spectrometry == FALSE ~ "TREC_olink",
      !is.na(log10TREC) ~ "TREC",
      TRUE ~ NA_character_
    )
  )

niveaux <- c("TREC_olink_MS","TREC_olink","TREC")
df$assay <- factor(df$assay, levels = niveaux)
summary(df$assay)

df_plot <- df %>% filter(!is.na(weekcb), !is.na(child_id), !is.na(assay))

df_plot <- df_plot %>%
  mutate(
    weekcb = droplevels(weekcb),
    child_id = droplevels(factor(child_id))
  )

p <-ggplot(df_plot, aes(x = weekcb, y = child_id, color = assay, group = child_id)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +  
  geom_point(aes(shape = as.factor(group)), size = 2.5) + 
  geom_point(size = 1) +                      
  theme_minimal() +
  labs(x = "Children Age(week)",
       y = "Children",
       color = "Analysis perfomed",
       shape = "group"
  ) + 
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),     
    axis.ticks.y = element_blank(), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    
  )

ggsave(
  filename = "figures/supplementary/Supp1b.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 2a
df_filtered <- df %>%
  filter(baby_age > 100, !is.na(log10TREC), !is.na(group))

p <- ggplot(df_filtered, aes(x = group, y = log10TREC)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +  
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = max(df_filtered$log10TREC, na.rm = TRUE) * 1.05
  ) +
  theme_minimal() +
  labs(
    x = "Group",
    y = "Log10TREC",
    
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

ggsave(
  filename = "figures/supplementary/Supp2a.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 2b

df_Preterm <- subset(df, group == "Preterm")
df_Preterm_clean <- df_Preterm %>% 
  filter(
    !is.na(log10TREC),
    !is.na(SNP),
    baby_age != 0
  )
df_Preterm_clean$log10age <- log10(df_Preterm_clean$baby_age)

m <- glm(log10TREC ~ log10age + SNP , data = df_Preterm_clean)
model_pvals <- summary(m)$coefficients[, "Pr(>|t|)"]
names(model_pvals)


bootstrap_glm <- function(data, formula, family = gaussian(), n_bootstrap = 100, sample_fraction = 0.9) {
  set.seed(123)  # pour reproductibilité
  
  n <- nrow(data)
  sample_size <- floor(sample_fraction * n)
  
  # Ajuster un modèle une fois pour connaître tous les noms de coefficients
  initial_model <- glm(formula, data = data, family = family)
  coef_names <- names(coef(initial_model))
  
  # Initialiser la matrice pour stocker les coefficients bootstrap
  bootstrap_coefs <- matrix(NA, nrow = n_bootstrap, ncol = length(coef_names))
  colnames(bootstrap_coefs) <- coef_names
  
  for (i in 1:n_bootstrap) {
    sample_indices <- sample(1:n, size = sample_size, replace = TRUE)
    sample_data <- data[sample_indices, ]
    
    model <- try(glm(formula, data = sample_data, family = family), silent = TRUE)
    
    if (inherits(model, "try-error")) {
      next  # sauter cette itération si le modèle échoue
    }
    
    model_coefs <- coef(model)
    
    # S'assurer que les bons coefficients vont aux bons endroits
    matched_coefs <- rep(NA, length(coef_names))
    names(matched_coefs) <- coef_names
    matched_coefs[names(model_coefs)] <- model_coefs
    
    bootstrap_coefs[i, ] <- matched_coefs
  }
  
  # Calculer les statistiques bootstrap
  bootstrap_summary <- data.frame(
    coefficient = colMeans(bootstrap_coefs, na.rm = TRUE),
    std_error = apply(bootstrap_coefs, 2, sd, na.rm = TRUE),
    lower = apply(bootstrap_coefs, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(bootstrap_coefs, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  bootstrap_summary <- tibble::rownames_to_column(bootstrap_summary, "var_name")
  
  return(list(
    bootstrap_summary = bootstrap_summary,
    all_bootstraps = bootstrap_coefs
  ))
}

result <- bootstrap_glm(df_Preterm_clean,
                        formula = log10TREC ~ log10age + SNP,
                        family = gaussian(),
                        n_bootstrap = 500)

summary_df <- result$bootstrap_summary %>%
  mutate(
    classical_pval = model_pvals[var_name],
    signif = ifelse(classical_pval < 0.05, "Significatif", "NS"),
    label = sprintf("p = %.4g", classical_pval)  # ✅ p-val exacte, notation intelligente
  ) %>%
  filter(var_name != "(Intercept)")


p <- ggplot(summary_df, aes(x = var_name, y = coefficient, color = signif)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(aes(label = label), hjust = -0.1, nudge_y = -0.25, size = 3.5)+
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("NS" = "black", "Significatif" = "red")) +
  theme_minimal() +
  
  coord_flip()

ggsave(
  filename = "figures/supplementary/Supp2b.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 2c
df_term <- subset(df, group == "Term")
df_term_clean <- df_term %>% 
  filter(
    !is.na(log10TREC),
    !is.na(SNP),
    baby_age != 0
  )
df_term_clean$log10age <- log10(df_term_clean$baby_age)

m <- glm(log10TREC ~ log10age + SNP , data = df_term_clean)
model_pvals <- summary(m)$coefficients[, "Pr(>|t|)"]
names(model_pvals)


result <- bootstrap_glm(df_term_clean,
                        formula = log10TREC ~ log10age + SNP,
                        family = gaussian(),
                        n_bootstrap = 500)

summary_df <- result$bootstrap_summary %>%
  mutate(
    classical_pval = model_pvals[var_name],
    signif = ifelse(classical_pval < 0.05, "Significatif", "NS"),
    label = sprintf("p = %.4g", classical_pval)  # ✅ p-val exacte, notation intelligente
  ) %>%
  filter(var_name != "(Intercept)")

p <- ggplot(summary_df, aes(x = var_name, y = coefficient, color = signif)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(aes(label = label), hjust = -0.1, nudge_y = -0.25, size = 3.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("NS" = "black", "Significatif" = "red")) +
  theme_minimal() +
  coord_flip()

ggsave(
    filename = "figures/supplementary/Supp2c.pdf",
    plot = p,
    width = 9,
    height = 9,
    units = "in"
  )
#Supp 3a
p2 <- ggplot(df %>% filter(!is.na(age_bin), !is.na(TRANCE)) %>% droplevels(),
             aes(x = age_bin, y = TRANCE)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Preterm" = "orange", "Term" = "purple")) +
  labs(title = "Scatter plot of RANKL depending on baby age and SNP",
       x = "Baby age",
       y = "RANKL") +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3a.pdf",
  plot = p2,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3b
p2 <- ggplot(df %>% filter(!is.na(age_bin), !is.na(TNFB)) %>% droplevels(),
             aes(x = age_bin, y = TNFB)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Preterm" = "orange", "Term" = "purple")) +
  labs(title = "Scatter plot of LTa depending on baby age and SNP",
       x = "Baby age",
       y = "LTa") +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3b.pdf",
  plot = p2,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3c
cor_test <- cor.test(df$log10TREC, df$TRANCE, method = "spearman", exact = FALSE)
rho <- round(cor_test$estimate, 3)
pval <- cor_test$p.value
p_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 3)))
label_text <- paste0("Spearman rho = ", rho, "\n", p_label)

p1 <- ggplot(df, aes(x = log10TREC, y = TRANCE, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") + # droite globale
  scale_color_manual(values = c("Preterm" = "orange", "Term" = "purple")) +
  annotate("text",
    x = min(df$log10TREC, na.rm = TRUE),
    y = max(df$TRANCE, na.rm = TRUE),
    label = label_text,
    hjust = 0, size = 5, fontface = "bold"
  ) +
  labs(
    title = "Correlation of TRANCE, vs log10TREC",
    x = "log10TREC",
    y = "TRANCE,",
  ) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(
  filename = "figures/supplementary/Supp3c.pdf",
  plot = p1,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3d
p1 <- ggplot(df_Preterm, aes(x = log10TREC, y = TRANCE)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = min(df_term$log10TREC, na.rm = TRUE),
           label.y = max(df$TRANCE, na.rm = TRUE), size = 5) +
  
  labs(title = "Scatter plot of RANKL vs log10sjTREC, preterm",
       x = "log10TREC",
       y = "RANKL") +
  theme_minimal()

p2 <- ggplot(df_term, aes(x = log10TREC, y = TRANCE)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(
    method = "spearman", label.x = min(df_term$log10TREC, na.rm = TRUE),
    label.y = max(df$TRANCE, na.rm = TRUE), size = 5
  ) +
  labs(
    title = "Scatter plot of RANKL vs log10sjTREC, term",
    x = "log10TREC",
    y = "RANKL"
  ) +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3d.pdf",
  plot = p1 / p2,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3e
cor_test <- cor.test(df$CD4T, df$TRANCE, method = "spearman", exact = FALSE)
rho <- round(cor_test$estimate, 3)
pval <- cor_test$p.value
p_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 3)))
label_text <- paste0("Spearman rho = ", rho, "\n", p_label)
p <- ggplot(df, aes(x = CD4T, y = TRANCE, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("Preterm" = "orange", "Term" = "purple")) +
  labs(
    title = "Correlation of TRANCE, vs CD4T cells",
    x = "CD4T",
    y = "TRANCE",
  ) +
  annotate("text",
    x = min(df$CD4T, na.rm = TRUE),
    y = max(df$TRANCE, na.rm = TRUE),
    label = label_text,
    hjust = 0, size = 5, fontface = "bold"
  ) +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3e.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3f
cor_test <- cor.test(df$CD4T, df$TNFB, method = "spearman", exact = FALSE)
rho <- round(cor_test$estimate, 3)
pval <- cor_test$p.value
p_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 3)))
label_text <- paste0("Spearman rho = ", rho, "\n", p_label)
p <- ggplot(df, aes(x = CD4T, y = TNFB, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("Preterm" = "orange", "Term" = "purple")) +
  labs(
    title = "Correlation of Lta vs CD4T cells",
    x = "CD4T",
    y = "LTa",
  ) +
  annotate("text",
    x = min(df$CD4T, na.rm = TRUE),
    y = max(df$TNFB, na.rm = TRUE),
    label = label_text,
    hjust = 0, size = 5, fontface = "bold"
  ) +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3f.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3g
p <- ggplot(df, aes(x = log10TREC, y = TRANCE)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = min(df_term$log10TREC, na.rm = TRUE),
           label.y = max(df$TRANCE, na.rm = TRUE), size = 5) +
  
  labs(title = "Scatter plot of RANKL vs log10sjTREC",
       x = "log10TREC",
       y = "RANKL") +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3g.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
#Supp 3h
p <- ggplot(df, aes(x = CD4T, y = TRANCE)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", label.x = min(df_term$log10TREC, na.rm = TRUE),
           label.y = max(df$TRANCE, na.rm = TRUE), size = 5) +
  
  labs(title = "Scatter plot of RANKL vs CD4T",
       x = "CD4T",
       y = "RANKL") +
  theme_minimal()
ggsave(
  filename = "figures/supplementary/Supp3h.pdf",
  plot = p,
  width = 9,
  height = 9,
  units = "in"
)
