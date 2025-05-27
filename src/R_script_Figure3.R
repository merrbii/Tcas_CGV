# ====================================================================
# PoolSeq FET Analysis:
# Tribolium castaneum RNAi and 17-DMAG normal versus reduced-eye lines
# Author: Mohammed Errbii 2025
# ====================================================================

# -----------------------
# Load required libraries
# -----------------------
library(tidyverse)
library(cowplot)
library(patchwork)

# ------------------------------
# Load and prepare PoolSeq data
# ------------------------------
# Load Fisher Exact Test (FET) results
dt1 <- read_delim(
  "results/RNe_RRe_DNe_DRe.covFiltered.edited.fet",
  col_names = FALSE
)

# Assign column names
colnames(dt1) <- c(
  "CHR", "CENTER", "N_SNP", "COV", "MEAN_COV",
  "RNe_RRe", "RNe_DNe", "RNe_DRe", "RRe_DNe", "RRe_DRe", "DNe_DRe"
)

# Keep relevant columns and preview
dt <- dt1[, c("CHR", "CENTER", "RNe_RRe", "DNe_DRe")]
head(dt)

# -------------------------------
# Apply FDR correction to p-values
# -------------------------------
dt$pvals1 <- 10^(-dt$RNe_RRe)
dt$pvals2 <- 10^(-dt$DNe_DRe)

# FDR-adjusted and transformed
dt$FDRpvals1 <- -log10(p.adjust(dt$pvals1, method = "fdr"))
dt$FDRpvals2 <- -log10(p.adjust(dt$pvals2, method = "fdr"))

# Keep only linkage groups (LG)
dt <- dt %>% filter(str_detect(CHR, "LG"))

# ------------------------
# Reshape data for plotting
# ------------------------
fstmelted <- reshape2::melt(
  dt,
  id = c("CHR", "CENTER", "pvals1", "pvals2", "RNe_RRe", "DNe_DRe")
)

# Define chromosome order
fstmelted$CHR <- factor(fstmelted$CHR, levels = paste0("LG", c(2:10, "X")))

# Label definitions
labs <- c("FDRpvals1" = "RNAi line", "FDRpvals2" = "17-DMAG line")

# Subset for FDR plots
sub <- fstmelted %>%
  filter(variable %in% c("FDRpvals1", "FDRpvals2")) %>%
  mutate(CHR = factor(CHR, levels = paste0("LG", c(2:10, "X"))))

nCHR <- length(unique(sub$CHR))

# -------------------------------
# Plot A: Manhattan-style plot
# -------------------------------
A <- ggplot(sub, aes(CENTER / 1e6, value)) +
  geom_point(aes(color = CHR), alpha = 0.55, size = 0.5) +
  scale_color_manual(values = rep(c("grey", "#183059"), nCHR)) +
  scale_y_continuous(
    name = expression(paste("-log"["10"], " (FDR-corrected ", italic("p"), ")"))
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 30, 2)
  ) +
  xlab("Chromosome") +
  theme_bw() +
  theme(
    panel.spacing.x=unit(0,"lines"),
    panel.spacing.y=unit(.5,"lines"),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 7, family = "Arial"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 7, family = "Arial"),
    strip.background = element_rect(fill = 'white', colour = 'black'),
    panel.spacing = unit(0.1, "cm")
  ) +
  facet_grid(
    variable ~ CHR,
    scales = "free_x",
    space = "free_x",
    switch = "x",
    labeller = labeller(variable = labs)
  )


# -------------------------------
# Plot B: Correlation scatter plot (LG10 only)
# -------------------------------
sub2 <- dt %>% filter(CHR == "LG10")

B <- ggplot(sub2, aes(FDRpvals1, FDRpvals2)) +
  geom_point(
    color = ifelse(sub2$FDRpvals1 >= 20 & sub2$FDRpvals2 >= 20, "#fd7f20", "grey"),
    size = 0.5,
    alpha = ifelse(sub2$FDRpvals1 >= 20 & sub2$FDRpvals2 >= 20, 0.8, 0.3)
  ) +
  scale_x_continuous(name = expression(paste("RNAi: -log"["10"], " (FDR-corrected ", italic("p"), ")"))) +
  scale_y_continuous(name = expression(paste("17-DMAG: -log"["10"], " (FDR-corrected ", italic("p"), ")"))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 7, family = "Arial")
  )

# -------------------------------
# Plot C: Zoomed Manhattan (LG10 only)
# -------------------------------
sub3 <- sub %>% filter(CHR == "LG10")

C <- ggplot(sub3, aes(CENTER / 1e6, value)) +
  geom_point(
    color = ifelse(
      (sub3$variable == "FDRpvals1" | sub3$variable == "FDRpvals2") &
        sub3$value >= 20 & sub3$CENTER >= 3871149 & sub3$CENTER <= 3906982,
      "#fd7f20", "grey"
    ),
    size = 0.5,
    alpha = ifelse(
      (sub3$variable == "FDRpvals1" | sub3$variable == "FDRpvals2") &
        sub3$value >= 20 & sub3$CENTER >= 3871149 & sub3$CENTER <= 3906982,
      0.8, 0.3
    )
  ) +
  coord_cartesian(xlim = c(3.8, 3.951)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(3.8, 3.95, 0.05)) +
  scale_y_continuous(name = expression(paste("-log"["10"], " (FDR-corrected ", italic("p"), ")"))) +
  xlab("Position (Mb)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x=unit(.3,"lines"),
    panel.spacing.y=unit(.5,"lines"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 7, family = "Arial"),
    axis.title = element_text(size = 7, family = "Arial"),
    strip.text = element_text(size = 7, family = "Arial"),
    strip.background = element_rect(fill = 'white', colour = 'black'),
    strip.text.x.bottom = element_blank(),
    panel.spacing = unit(0.1, "cm")
  ) +
  facet_grid(
    variable ~ CHR,
    scales = "free_x",
    space = "free_x",
    switch = "x",
    labeller = labeller(variable = labs)
  )


# -------------------------------
# Combine A, B, and C into final figure
# -------------------------------
BC <- plot_grid(B, plot_spacer(), C, nrow = 1, rel_widths = c(0.35, 0.03, 0.65))
final_plot <- plot_grid(A, plot_spacer(), BC, nrow = 3, rel_heights = c(0.5, 0.05, 0.5))

# Save final figure
ggsave(
  filename = "results/RNe_RRe_DNe_DRe_perSNP_fdr_correctedP_Figure.Arial2.png",
  plot = final_plot,
  width = 10,
  height = 6.5,
  dpi = 500
)

# -------------------------------
# Extract candidate region (FDR ≥ 20 in both lines)
# -------------------------------
candidate <- sub2 %>% filter(FDRpvals1 >= 20 & FDRpvals2 >= 20)

# Save candidate SNPs
write_delim(
  candidate,
  "candidates/candidateSNPs_logFDRpval20.txt",
  delim = "\t",
  col_names = TRUE
)

# Print candidate region boundaries
cat("Candidate region:\n")
cat(paste0("LG10:", min(candidate$CENTER), "-", max(candidate$CENTER)), "\n")

# Print candidate region ±32kb to obtain a 100kb candidate region
flank <- 32083.5
cat("Candidate region ±32kb:\n")
cat(paste0("LG10:", min(candidate$CENTER) - flank, "-", max(candidate$CENTER) + flank), "\n")