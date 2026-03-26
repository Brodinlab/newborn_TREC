library(dplyr)
library(readxl)
library(ggpubr)
library(rstatix)

RANKL_TEC <- read_excel("data/RANKL_TEC.xlsx")

p <- ggplot(RANKL_TEC, aes(x = cTEC, y = gMFI)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_pubr() +
  ylab("RANKL gMFI")
ggsave(
  filename = "figures/Fig4/4e.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
p <- ggplot(RANKL_TEC, aes(x = mTEC, y = gMFI)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_pubr() +
  ylab("RANKL gMFI")
ggsave(
  filename = "figures/Fig4/4f.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
