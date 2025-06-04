library(dplyr)
library(ggplot2)
library(future.apply)
library(tidyr)
library(purrr)
library(broom)
library(patchwork)
library(scales)
plan(multisession)


### directories
snp_dir <- "/home/berdreyminn/Desktop/ultima/posteriors_snp/"
wgs_dir <- "/home/berdreyminn/Desktop/ultima/posteriors_wgs/"

### get folders with our files and ID extraction
snp_files <- list.files(snp_dir, pattern = "_snp_posterior$", full.names = TRUE)
wgs_files <- list.files(wgs_dir, pattern = "_wgs_posterior$", full.names = TRUE)

#
get_id <- function(path) sub("(_snp_posterior|_wgs_posterior)$", "", basename(path))
snp_ids <- sapply(snp_files, get_id)
wgs_ids <- sapply(wgs_files, get_id)
common_ids <- intersect(snp_ids, wgs_ids)

### calculate allele freq from three compatible genotypes
get_freq <- function(a1, a2, probs) {
  g1 <- paste0(a1, a1)
  g2 <- paste0(a1, a2)
  g3 <- paste0(a2, a2)
  p1 <- as.numeric(probs[g1])
  p2 <- as.numeric(probs[g2])
  p3 <- as.numeric(probs[g3])
  if (any(is.na(c(p1, p2, p3)))) return(NA)
  total <- p1 + p2 + p3
  return((p1 + 0.5 * p2) * total)
}


### summaries
bias_list <- list()
valid_alleles <- c("A", "C", "G", "T")

### Calculate directional mean difference
bias_list <- future_lapply(common_ids, function(id) {
  cat("processing", id, "
")
  snp_file <- file.path(snp_dir, paste0(id, "_snp_posterior"))
  wgs_file <- file.path(wgs_dir, paste0(id, "_wgs_posterior"))
  
  snp <- tryCatch(read.table(snp_file, header = TRUE, sep = "	"), error = function(e) NULL)
  wgs <- tryCatch(read.table(wgs_file, header = TRUE, sep = "	"), error = function(e) NULL)
  if (is.null(snp) || is.null(wgs)) return(NULL)
  
  merged <- merge(snp, wgs, by = "Variant_ID", suffixes = c("_snp", "_wgs"))
  
  merged <- merged %>%
    mutate(
      ref_snp = toupper(trimws(Ref_snp)),
      alt_snp = toupper(trimws(Alt_snp)),
      ref_wgs = toupper(trimws(Ref_wgs)),
      alt_wgs = toupper(trimws(Alt_wgs))
    ) %>%
    filter(
      ref_snp %in% valid_alleles & alt_snp %in% valid_alleles &
        ref_wgs %in% valid_alleles & alt_wgs %in% valid_alleles
    ) %>%
    rowwise() %>%
    mutate(
      snp_pair = paste(sort(c(ref_snp, alt_snp)), collapse = ""),
      wgs_pair = paste(sort(c(ref_wgs, alt_wgs)), collapse = ""),
      keep = snp_pair == wgs_pair
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    rowwise() %>%
    mutate(
      allele1 = sort(c(ref_snp, alt_snp))[1],
      allele2 = sort(c(ref_snp, alt_snp))[2],
      allele_pair = paste0(allele1, "-", allele2),
      freq1_snp = get_freq(allele1, allele2,
                           c(AA = AA_snp, AC = AC_snp, AG = AG_snp, AT = AT_snp,
                             CC = CC_snp, CG = CG_snp, CT = CT_snp,
                             GG = GG_snp, GT = GT_snp, TT = TT_snp)),
      freq1_wgs = get_freq(allele1, allele2,
                           c(AA = AA_wgs, AC = AC_wgs, AG = AG_wgs, AT = AT_wgs,
                             CC = CC_wgs, CG = CG_wgs, CT = CT_wgs,
                             GG = GG_wgs, GT = GT_wgs, TT = TT_wgs)),
      freq_diff = freq1_snp - freq1_wgs
    ) %>%
    ungroup() %>%
    filter(!is.na(freq_diff))
  
  summary <- merged %>%
    group_by(allele_pair) %>%
    summarise(
      Bias = mean(freq_diff, na.rm = TRUE),
      n = n()
    ) %>%
    filter(n >= 1) %>%
    mutate(Sample = id) %>%
    ungroup()
  
  return(summary)
})

### combining summaries
bias_long <- do.call(bind_rows, bias_list)

### exclude two samples with not enough data (optional)
bias_long <- bias_long %>% filter(!Sample %in% c("AL2921", "OL4117"))

### Define your desired sample order
sample_order <- c(
  "ROM-3102240001", "ROM-15038",
  "AL2541", "AL2744", "AL2657", "MCZ-BOM-11183",
  "LS0121", "LS0122", "LS0123", "LS0124", "LS0125", "LS0126",
  "ANH5", "AC02",
  "KT0033", "KT0039", "KT0041", "KT0043", "KT0052", "KT0056",
  "ACAD17371"
)

### generate combinations using the specified order
all_combinations <- expand.grid(
  Sample = sample_order,  # Use your ordered samples here
  allele_pair = unique(bias_long$allele_pair),
  stringsAsFactors = FALSE
)

### join
bias_long_complete <- left_join(all_combinations, bias_long, by = c("Sample", "allele_pair")) %>%
  mutate(
    n_label = ifelse(is.na(n), "NA", as.character(n)),
    Sample = factor(Sample, levels = sample_order)  # Enforce the order
  )

### plotting the heatmap (unchanged except for sample order)
ggplot(bias_long_complete, aes(x = Sample, y = allele_pair, fill = Bias)) +
  geom_tile(color = "white", size = 0.4) +
  geom_text(aes(label = n_label), color = "black", angle = 90, size = 3, vjust = 0.5) +
  scale_fill_gradient2(
    low = "yellow", mid = "white", high = "purple", midpoint = 0,
    na.value = "gray90",
    name = "Frequency difference \nfor the first allele\n(SNP - WGS)"
  ) +
  labs(
    title = "",
    x = "Sample", y = "Allele Pair"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


##SAVE optional
ggsave("/home/berdreyminn/Desktop/try4.png", width = 10, height = 6, dpi = 500, bg='white')

##########################################stats
###############################################

########################################################
# A. Weighted absolute bias per sample (normalized) (based on directional)
sample_weighted_abs_bias <- bias_long_complete %>%
  filter(!is.na(Bias) & !is.na(n)) %>%
  group_by(Sample) %>%
  mutate(
    total_n = sum(n, na.rm = TRUE),
    weight = n / total_n,
    weighted_abs_bias = abs(Bias) * weight
  ) %>%
  summarise(
    WeightedAbsBias = sum(weighted_abs_bias, na.rm = TRUE),
    TotalSites = unique(total_n)
  ) %>%
  ungroup()

########################################################
# Weighted bootstrap function
weighted_bootstrap_p <- function(bias, weights, R = 10000, center = 0) {
  probs <- weights / sum(weights, na.rm = TRUE)
  boot_means <- replicate(R, {
    sampled <- sample(bias, size = length(bias), replace = TRUE, prob = probs)
    mean(sampled, na.rm = TRUE)
  })
  observed <- sum(bias * probs, na.rm = TRUE)
  p_value <- mean(abs(boot_means - center) >= abs(observed - center))
  return(list(p = p_value, observed = observed))
}


########################################################
# B. Statistical test per sample (Wilcoxon on raw Bias)
bias_tests_weighted <- bias_long_complete %>%
  filter(!is.na(Bias) & !is.na(n)) %>%
  group_by(Sample) %>%
  summarise(
    test = list(weighted_bootstrap_p(Bias, n)),
    p.value = test[[1]]$p,
    statistic = NA_real_,  # placeholder; bootstrap doesn't produce a "statistic"
    method = "Weighted bootstrap"
  ) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),
    Significant_FDR = ifelse(p.adj < 0.05, "Yes", "No"),
    sig_label = ifelse(Significant_FDR == "Yes", "*", "")
  )

########################################################
# C. Plot data: merge sample-level weighted bias with test results
bias_plot_data_weighted <- left_join(sample_weighted_abs_bias, bias_tests_weighted, by = "Sample")

########################################################
# D. Plot: per-sample weighted absolute bias (based on directional)
gg1 <- ggplot(bias_plot_data_weighted, aes(x = Sample, y = WeightedAbsBias)) +
  geom_col(fill = "#7B3294", width = 0.7) +
  geom_text(aes(label = sig_label), vjust = -0.2, size = 5.5, color = "black", fontface = "bold") +
  geom_text(aes(label = comma(TotalSites)), vjust = -1.6, size = 3.2, color = "black", fontface = "bold") +
  labs(
    x = "Specimen",
    y = "Specimen Bias Based on Directional Differences",
    title = "Weighted SNP vs WGS Bias per Sample",
    subtitle = "* = FDR-significant (0.05) | Numbers above bars = total sites"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

gg1

########################################################
# E. Weighted mean directional bias per allele pair
allele_weighted_mean_bias <- bias_long_complete %>%
  filter(!is.na(Bias) & !is.na(n)) %>%
  group_by(allele_pair) %>%
  mutate(
    total_n = sum(n, na.rm = TRUE),
    weight = n / total_n,
    weighted_bias = Bias * weight
  ) %>%
  summarise(
    WeightedMeanBias = sum(weighted_bias, na.rm = TRUE),
    TotalSites = unique(total_n)
  ) %>%
  ungroup()

########################################################
# F. Wilcoxon test per allele pair (on unweighted Bias)
allele_tests_weighted <- bias_long_complete %>%
  filter(!is.na(Bias) & !is.na(n)) %>%
  group_by(allele_pair) %>%
  summarise(
    test = list(weighted_bootstrap_p(Bias, n)),
    p.value = test[[1]]$p,
    statistic = NA_real_,
    method = "Weighted bootstrap"
  ) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),
    Significant_FDR = ifelse(p.adj < 0.05, "Yes", "No"),
    sig_label = ifelse(Significant_FDR == "Yes", "*", "")
  )


########################################################
# G. Merge for allele-pair plot
allele_plot_data <- left_join(allele_weighted_mean_bias, allele_tests_weighted, by = "allele_pair")

########################################################
# H. Plot: per-allele pair weighted directional bias
gg1N <- ggplot(allele_plot_data, aes(x = allele_pair, y = WeightedMeanBias)) +
  geom_col(fill = "#7B3294", width = 0.7) +
  geom_text(aes(label = sig_label), vjust = -0.2, size = 6, color = "black", fontface = "bold") +
  geom_text(aes(label = comma(TotalSites)), vjust = -1.6, size = 3.2, color = "black", fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    x = "Allele Pair",
    y = "Weighted Mean Directional Bias",
    title = "Directional Bias of SNP vs WGS by Allele Pair",
    subtitle = "* = FDR-significant (0.05) | Numbers above bars = total sites"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

gg1N

########################################################
# I. Final combined plot
plotstats <- gg1 + gg1N
plotstats

ggsave("/home/berdreyminn/Desktop/biasstats.png", plotstats,
       width = 14, height = 8, dpi = 500, bg = "white")








