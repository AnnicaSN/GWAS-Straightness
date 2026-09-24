library(qqman)
library(qvalue)
library(data.table)
library(dplyr)
library(tidyr)


result <- fread("BLUPCurveRatioGWAS_savar_old_only_biallelic_maf01_missing10_ordered_+_phenoBLUPCurveRatio_stats.txt")


#results variables
trait <- NULL
lambda <- NULL
cutoff05 <- NULL
cutoff01 <- NULL
cutoff001 <- NULL
fdr05 <- NULL
fdr01 <- NULL
fdr001 <- NULL

#loop through GWAS results for each trait
for(i in unique(result$Trait)){
  #Estimate genomic inflation (lambda)
  data <- qchisq(result$p[which(result$Trait==i)],1,lower.tail=F)
  data <- sort(data)
  ppoi <- ppoints(data) #Generates the sequence of probability points
  ppoi <- sort(qchisq(ppoi,df=1,lower.tail=F))
  s <- summary(lm(data ~ 0 + ppoi))$coeff
  
  #add trait and lambda to results variables
  trait <- c(trait,i)
  lambda <- c(lambda,s[1,1])
  
  lambda_summary <- list(
    median = median(lambda, na.rm = TRUE),
    q1     = quantile(lambda, 0.25, na.rm = TRUE),
    q3     = quantile(lambda, 0.75, na.rm = TRUE),
    IQR    = IQR(lambda, na.rm = TRUE),
    min    = min(lambda, na.rm = TRUE),
    max    = max(lambda, na.rm = TRUE)
  )
  lambda_summary
  
  cat(sprintf("λ across loops — median = %.3f, IQR = %.3f–%.3f (range %.3f–%.3f)\n",
              lambda_summary$median, lambda_summary$q1, lambda_summary$q3,
              lambda_summary$min, lambda_summary$max))
  
  #estimate empirical cutoffs at alpha=0.05,0.01,0.001; then add to results variables
  cutoff05 <- c(cutoff05,quantile(-log10(result$p[which(result$Trait==i&result$p<0.05)]),0.95))
  cutoff01 <- c(cutoff01,quantile(-log10(result$p[which(result$Trait==i&result$p<0.05)]),0.99))
  cutoff001 <- c(cutoff001,quantile(-log10(result$p[which(result$Trait==i&result$p<0.05)]),0.999))
  
  #FDR at alpha=0.05,0.01,0.001
  qval05 <- qvalue(result$p[which(result$Trait==i)],fdr.level=0.05)
  qval01 <- qvalue(result$p[which(result$Trait==i)],fdr.level=0.01)
  qval001 <- qvalue(result$p[which(result$Trait==i)],fdr.level=0.001)
  
  #Add FDR cutoffs at alpha=0.05,0.01,0.001 to results variables; if no SNPs pass FDR, save cutoff as NA
  if(length(which(qval05$significant==T))>0){
    fdr05 <- c(fdr05,-log10(max(qval05$pvalues[which(qval05$significant==T)])))
  }else{fdr05 <- c(fdr05,NA)}
  if(length(which(qval01$significant==T))>0){
    fdr01 <- c(fdr01,-log10(max(qval01$pvalues[which(qval01$significant==T)])))
  }else{fdr01 <- c(fdr01,NA)}
  if(length(which(qval001$significant==T))>0){
    fdr001 <- c(fdr001,-log10(max(qval001$pvalues[which(qval001$significant==T)])))
  }else{fdr001 <- c(fdr001,NA)}
}

traits <- sort(unique(result$Trait))

# Storage
lambda_slope_all  <- setNames(rep(NA_real_, length(traits)), traits)
lambda_slope_trim <- setNames(rep(NA_real_, length(traits)), traits)
lambda_gc_all     <- setNames(rep(NA_real_, length(traits)), traits)   # median-χ² method
lambda_gc_trim    <- setNames(rep(NA_real_, length(traits)), traits)

for (tr in traits) {
  p <- result$p[result$Trait == tr]
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  n <- length(p)
  if (n < 50) next  # skip tiny sets safely
  
  # Chi-square stats and expected under null
  obs <- sort(qchisq(p, df = 1, lower.tail = FALSE))
  exp <- sort(qchisq(ppoints(n), df = 1, lower.tail = FALSE))
  
  # ---- λ (QQ-slope) using ALL points
  lambda_slope_all[tr] <- coef(lm(obs ~ 0 + exp))[1]
  
  # ---- λ (median-χ² / "genomic control")
  lambda_gc_all[tr] <- median(obs) / qchisq(0.5, df = 1)  # 0.456...
  
  # ---- Trim top 1% (drop the largest 1% of obs / smallest 1% p)
  k <- floor(0.99 * n)
  if (k >= 10) {
    lambda_slope_trim[tr] <- coef(lm(obs[1:k] ~ 0 + exp[1:k]))[1]
    lambda_gc_trim[tr]    <- median(obs[1:k]) / qchisq(0.5, df = 1)
  }
}

# Summaries to report
summ <- function(x) c(median = median(x, na.rm=TRUE),
                      q1     = quantile(x, .25, na.rm=TRUE),
                      q3     = quantile(x, .75, na.rm=TRUE),
                      IQR    = IQR(x, na.rm=TRUE))

s_all   <- summ(lambda_slope_all)
s_trim  <- summ(lambda_slope_trim)
gc_all  <- summ(lambda_gc_all)
gc_trim <- summ(lambda_gc_trim)

cat(sprintf("\nλ (QQ-slope) — ALL:  median=%.3f, IQR=%.3f–%.3f\n",
            s_all['median'], s_all['q1'], s_all['q3']))
cat(sprintf("λ (QQ-slope) — TRIM: median=%.3f, IQR=%.3f–%.3f\n",
            s_trim['median'], s_trim['q1'], s_trim['q3']))

cat(sprintf("\nλ (median-χ²) — ALL:  median=%.3f, IQR=%.3f–%.3f\n",
            gc_all['median'], gc_all['q1'], gc_all['q3']))
cat(sprintf("λ (median-χ²) — TRIM: median=%.3f, IQR=%.3f–%.3f\n\n",
            gc_trim['median'], gc_trim['q1'], gc_trim['q3']))

## ===== Genomic inflation across all loops (and trimmed 1%) =====

# helper to compute lambda (QQ-slope) from a vector of p-values
compute_lambda <- function(p, trim_frac = 0) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  n <- length(p)
  if (n < 50) return(NA_real_)  # skip tiny sets safely
  
  obs <- sort(qchisq(p, df = 1, lower.tail = FALSE))
  exp <- sort(qchisq(ppoints(n), df = 1, lower.tail = FALSE))
  
  if (trim_frac > 0) {
    k <- floor((1 - trim_frac) * n)        # keep lower (1 - trim_frac)
    if (k < 10) return(NA_real_)
    obs <- obs[1:k]
    exp <- exp[1:k]
  }
  as.numeric(coef(lm(obs ~ 0 + exp))[1])   # QQ-slope estimate of lambda
}

traits <- sort(unique(result$Trait))

# per-loop lambdas
lambda_all  <- sapply(traits, function(tr) compute_lambda(result$p[result$Trait == tr], trim_frac = 0))
lambda_trim <- sapply(traits, function(tr) compute_lambda(result$p[result$Trait == tr], trim_frac = 0.01))  # drop top 1%

# tidy table (optional to inspect/save)
lambda_table <- data.frame(
  Trait            = traits,
  lambda_all       = as.numeric(lambda_all),
  lambda_trim_top1 = as.numeric(lambda_trim),
  stringsAsFactors = FALSE
)

# summaries to report
med_all <- median(lambda_all, na.rm = TRUE)
q1_all  <- quantile(lambda_all, 0.25, na.rm = TRUE)
q3_all  <- quantile(lambda_all, 0.75, na.rm = TRUE)
rng_all <- range(lambda_all, na.rm = TRUE)

med_trim <- median(lambda_trim, na.rm = TRUE)
q1_trim  <- quantile(lambda_trim, 0.25, na.rm = TRUE)
q3_trim  <- quantile(lambda_trim, 0.75, na.rm = TRUE)
rng_trim <- range(lambda_trim, na.rm = TRUE)

cat(sprintf("λ across loops — median = %.3f, IQR = %.3f–%.3f (range %.3f–%.3f)\n",
            med_all, q1_all, q3_all, rng_all[1], rng_all[2]))
cat(sprintf("λ (trim top 1%%) — median = %.3f, IQR = %.3f–%.3f (range %.3f–%.3f)\n",
            med_trim, q1_trim, q3_trim, rng_trim[1], rng_trim[2]))

# optional: write per-loop values for supplement
# write.csv(lambda_table, "lambda_per_loop.csv", row.names = FALSE)



fdrList <- matrix(NA, nrow = 15, ncol = 100)

for (i in 1:100) {
  pheno_name <- paste("Pheno", i, sep = "")
  
  # Filter rows where the Trait column matches the current phenotype
  filtered_data <- result[result$Trait == pheno_name, ] %>% drop_na(Pos)
  
  #fdr05s <- qvalue(p=filtered_data$p,fdr.level=0.05)
  fdrsignificant_indices <- which(-log10(filtered_data$p) > cutoff001[i]) #which(fdr05$pvalues < 0.05)       
  # Extract the markers corresponding to those indices
  fdrmarkers <- filtered_data$Marker[fdrsignificant_indices]
  filtered_data$p[fdrsignificant_indices]
  
  # Fill the i-th column with the markers
  fdrList[, i] <- c(fdrmarkers, rep(NA, nrow(fdrList) - length(fdrmarkers)))
}  

# Flatten the matrix into a single vector
flat_matrix <- as.vector(fdrList)

# Count the frequency of each unique value
fdrvalue_counts <- table(flat_matrix)

# Sort the frequencies in descending order and pick the top 30 most frequent values
fdrtop_30_valuesNames <- names(sort(fdrvalue_counts, decreasing = TRUE)[1:30])
fdrtop_30_values <- sort(fdrvalue_counts, decreasing = TRUE)[1:30]

# Sort your data
sorted_fdr <- sort(fdrvalue_counts, decreasing = TRUE)

# Turn it into a data frame (a table)
fdr_table <- data.frame(
  Value = names(sorted_fdr),
  Count = as.numeric(sorted_fdr)
)

# Save it as a CSV file (Excel can open this)
write.csv(fdr_table, "sorted_fdr_table.csv", row.names = FALSE)


##FRD RESULTS SMA_10431774_1058 67 SMA_10434783_5058 63  SMA_7066810_174 63  SMA_36867_13907 60  SMA_277070_7618 58

# --- SETTINGS YOU CAN TUNE ---
y_max <- 60  # common y-limit for all QQ panels

library(dplyr)
library(tidyr)

## Global mention counts & palette (blue -> yellow; 0 stays black)
mention_counts_tbl <- table(na.omit(as.vector(fdrList)))
max_count <- if (length(mention_counts_tbl)) max(mention_counts_tbl) else 0L
pal <- colorRampPalette(c("#0000FF", "#FF0000"))

count_to_col_vec <- function(cnt_vec) {
  cols <- rep("#000000", length(cnt_vec))
  if (max_count > 0L) {
    grad_cols <- pal(max_count)
    nz <- which(cnt_vec > 0L)
    cols[nz] <- grad_cols[pmin(cnt_vec[nz], max_count)]
  }
  cols
}

draw_mini_colorbar <- function() {
  if (max_count <= 0L) return(invisible(NULL))
  usr <- par("usr")
  x1 <- usr[1] + 0.05 * (usr[2] - usr[1])
  x2 <- usr[1] + 0.08 * (usr[2] - usr[1])
  y1 <- usr[4] - 0.35 * (usr[4] - usr[3])
  y2 <- usr[4] - 0.08 * (usr[4] - usr[3])
  
  grad_cols <- pal(max_count)
  n_steps   <- max_count
  yy <- seq(y1, y2, length.out = n_steps + 1)
  for (j in 1:n_steps) rect(x1, yy[j], x2, yy[j + 1], col = grad_cols[j], border = NA)
  rect(x1, y1, x2, y2, border = "black")
  
  # Ticks strictly 1..max_count
  tick_vals <- if (max_count == 1L) 1L else {
    tv <- pretty(c(1, max_count), n = 5)
    tv <- unique(pmax(1, pmin(max_count, round(tv))))
    tv[tv >= 1 & tv <= max_count]
  }
  y_ticks <- y1 + (tick_vals - 1) / max(1, (max_count - 1)) * (y2 - y1)
  text(x2 + 0.01 * (usr[2] - usr[1]), y_ticks, labels = tick_vals, cex = 0.7, adj = 0)
  
  text(x1, y2 + 0.02 * (usr[4] - usr[3]), "Mentions", adj = 0, cex = 0.8)
  # Optional 0 swatch
  rect(x1, y2 + 0.04 * (usr[4] - usr[3]),
       x1 + 0.02 * (usr[2] - usr[1]),
       y2 + 0.07 * (usr[4] - usr[3]),
       col = "black", border = "black")
  text(x1 + 0.024 * (usr[2] - usr[1]), y2 + 0.055 * (usr[4] - usr[3]), "0", adj = 0, cex = 0.7)
}

# ---- Plot 3x3 with shared axis labels ----
op <- par(no.readonly = TRUE)
on.exit(par(op))

par(mfrow = c(3, 3),
    mar   = c(1.8, 1.8, 1.6, 0.6),   # inner margins per panel (tight)
    oma   = c(5.5, 5.5, 2, 1.2),     # OUTER margins (space for shared labels)
    mgp   = c(1.6, 0.6, 0))          # axis title/label/line spacing a bit tighter

for (i in 1:9) {
  pheno_name <- paste0("Pheno", i)
  
  filtered_data <- result %>%
    filter(Trait == pheno_name) %>%
    drop_na(Pos, p)
  
  if (nrow(filtered_data) == 0) { plot.new(); title(pheno_name); next }
  
  mentions <- ifelse(filtered_data$Marker %in% names(mention_counts_tbl),
                     as.integer(mention_counts_tbl[filtered_data$Marker]), 0L)
  cols <- count_to_col_vec(mentions)
  
  ord  <- order(filtered_data$p, decreasing = FALSE)
  pvec <- filtered_data$p[ord]
  cols <- cols[ord]
  
  col_markers_i <- na.omit(fdrList[, i])
  smallLine <- NA_real_
  if (length(col_markers_i)) {
    p_in_col <- filtered_data$p[filtered_data$Marker %in% col_markers_i]
    if (length(p_in_col)) smallLine <- min(-log10(p_in_col))
  }
  
  # IMPORTANT: blank per-panel labels; axes still draw ticks
  qq(pvec, main = pheno_name, ylim = c(0, y_max), col = cols, pch = 20,
     xlab = "", ylab = "")
  
  if (!is.na(smallLine)) abline(h = smallLine, lwd = 2, col = "#CC79A7")
  
  # Draw the compact legend only once (top-left panel: i == 1)
  if (i == 1) draw_mini_colorbar()
}

# Shared axis labels in the OUTER margin (visible once)
mtext(expression(Expected~ -log[10](p)), side = 1, line = 3.6, outer = TRUE, cex = 1.05)
mtext(expression(Observed~ -log[10](p)), side = 2, line = 3.6, outer = TRUE, cex = 1.05)


png(width=2000,height=2000,res=300)

dev.copy2pdf(
  file = "Figure2n.pdf",
  width = 500 / 25.4,
  height = 300 / 25.4
)




plotqq <- result %>% drop_na(Pos)

snp1 <- plotqq$p[(which(plotqq$Marker%in%fdrtop_30_valuesNames))]
snp3 <- plotqq$p[(which(plotqq$Marker%in%fdrList))]
count <- length(which(snp3%in%snp1))
if (count > 0) {
  snps <- snp3[-which(snp3%in%snp1)]
} else {
  snps <- snp3
}

snp2 <- plotqq$p[-c((which(plotqq$Marker%in%fdrtop_30_valuesNames)), (which(plotqq$Marker%in%fdrList)))]
length(snp2)
length(plotqq$p) - length(snp1) - length(snps)

smallLine <- min(-log10(snp3))


qq(main="All loops", c(snps, snp1, snp2),col=c(rep("green",length(snps)), rep("red",length(snp1)) ,rep("black", length(snp2)))[order(c(snps,snp1, snp2), decreasing = FALSE)]) 


length(unique(na.omit(as.vector(fdrList))))











## --- assumes `result` is already loaded and has columns: Trait, p ---

compute_lambda <- function(p, trim_frac = 0) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  n <- length(p); if (n < 50) return(NA_real_)
  obs <- sort(qchisq(p, df = 1, lower.tail = FALSE))
  exp <- sort(qchisq(ppoints(n), df = 1, lower.tail = FALSE))
  if (trim_frac > 0) {
    k <- floor((1 - trim_frac) * n)
    if (k < 10) return(NA_real_)
    obs <- obs[1:k]; exp <- exp[1:k]
  }
  as.numeric(coef(lm(obs ~ 0 + exp))[1])  # QQ-slope λ
}

traits <- sort(unique(result$Trait))
lambda_all  <- sapply(traits, function(tr) compute_lambda(result$p[result$Trait==tr], 0))
lambda_trim <- sapply(traits, function(tr) compute_lambda(result$p[result$Trait==tr], 0.01))  # drop top 1%

# Summaries for annotation
med  <- median(lambda_all,  na.rm = TRUE)
q1   <- quantile(lambda_all, .25, na.rm = TRUE)
q3   <- quantile(lambda_all, .75, na.rm = TRUE)

cat(sprintf("λ across loops — median = %.3f, IQR = %.3f–%.3f\n", med, q1, q3))




# ---- Histogram + density for λ ----
op <- par(no.readonly = TRUE); on.exit(par(op))
par(mar=c(4,4,1.5,1))
hist(lambda_all, breaks = 20, col = "grey85", border = "grey50",
     xlab = expression(lambda), main = "Distribution of genomic inflation (λ) across loops")
rug(lambda_all, col = "grey40")
dens <- density(lambda_all, na.rm = TRUE)
lines(dens, lwd = 2)

abline(v = 1, col = "darkgreen", lwd = 2, lty = 2)          # ideal calibration
abline(v = med, col = "firebrick", lwd = 2)                 # median
abline(v = q1,  col = "firebrick", lwd = 1, lty = 3)        # IQR
abline(v = q3,  col = "firebrick", lwd = 1, lty = 3)
legend("topright",
       legend = c("Density", "Median", "IQR", "λ = 1"),
       lwd    = c(2,2,1,2), lty=c(1,1,3,2),
       col    = c("black","firebrick","firebrick","darkgreen"),
       bty = "n")



# ---- Optional: side-by-side violin/boxplot comparing λ vs trimmed λ ----
# install.packages("ggplot2") if needed
library(ggplot2)
df <- data.frame(
  lambda = c(as.numeric(lambda_all), as.numeric(lambda_trim)),
  type   = rep(c("All tests","Trim top 1%"), each = length(lambda_all))
)

ggplot(df, aes(type, lambda)) +
  geom_violin(fill = "grey85", color = "grey40", width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, color = "darkgreen") +
  labs(x = NULL, y = expression(lambda),
       title = "Genomic inflation per loop: all tests vs. trimmed (top 1% removed)") +
  theme_classic(base_size = 12)







# Packages
library(dplyr)
library(tidyr)
library(qvalue)
library(ggplot2)
library(forcats)

# 1) Count SNPs passing FDR in each loop
traits <- sort(unique(result$Trait))

fdr_counts <- do.call(rbind, lapply(traits, function(tr) {
  p <- result$p[result$Trait == tr]
  p <- p[is.finite(p) & !is.na(p)]
  if (length(p) < 10) {
    return(data.frame(Trait = tr, FDR_0.05 = NA_integer_, FDR_0.01 = NA_integer_, FDR_0.001 = NA_integer_))
  }
  c05  <- sum(qvalue(p, fdr.level = 0.05)$significant,  na.rm = TRUE)
  c01  <- sum(qvalue(p, fdr.level = 0.01)$significant,  na.rm = TRUE)
  c001 <- sum(qvalue(p, fdr.level = 0.001)$significant, na.rm = TRUE)
  data.frame(Trait = tr, FDR_0.05 = c05, FDR_0.01 = c01, FDR_0.001 = c001)
}))

# Optional: extract numeric loop index if your traits are "Pheno1..100"
fdr_counts$loop <- suppressWarnings(as.integer(gsub("\\D", "", fdr_counts$Trait)))

# Save (optional)
# write.csv(fdr_counts, "fdr_counts_per_loop.csv", row.names = FALSE)

# 2) Long format for plotting
long <- fdr_counts %>%
  pivot_longer(cols = starts_with("FDR_"), names_to = "alpha", values_to = "count") %>%
  mutate(alpha = factor(alpha, levels = c("FDR_0.05","FDR_0.01","FDR_0.001")))

# 3A) Bar plot, loops sorted by count within each alpha
ggplot(
  long %>% group_by(alpha) %>% arrange(desc(count), .by_group = TRUE) %>%
    mutate(Trait_sorted = fct_reorder(Trait, count, .desc = TRUE)),
  aes(Trait_sorted, count, fill = alpha)
) +
  geom_col(width = 0.8, show.legend = FALSE) +
  facet_wrap(~ alpha, ncol = 1, scales = "free_y") +
  geom_hline(
    data = long %>% group_by(alpha) %>% summarise(med = median(count, na.rm = TRUE)),
    aes(yintercept = med), linetype = 2, color = "grey30"
  ) +
  labs(x = NULL, y = "SNPs passing FDR per loop",
       title = "Discovery counts per bootstrap loop",
       subtitle = "Dashed line = median within each FDR level") +
  coord_flip() +
  theme_classic(base_size = 12)

# 3B) (Alternative) counts along loop index (keeps Pheno order on x)
# ggplot(long, aes(loop, count, color = alpha)) +
#   geom_point() + geom_line(alpha = 0.5) +
#   labs(x = "Loop (Pheno index)", y = "SNPs passing FDR",
#        title = "FDR discovery counts across bootstrap loops") +
#   theme_classic(base_size = 12)








flat_matrix <- as.vector(fdrList)

fdrvalue_counts <- table(flat_matrix)

sorted_fdr <- sort(
  fdrvalue_counts,
  decreasing = TRUE
)

fdr_table <- data.frame(
  Value = names(sorted_fdr),
  Count = as.numeric(sorted_fdr)
)


recurrence_distribution <-
  as.data.frame(
    table(as.numeric(fdrvalue_counts))
  )

names(recurrence_distribution) <-
  c("Loops", "SNPs")

recurrence_distribution$Loops <-
  as.numeric(as.character(
    recurrence_distribution$Loops
  ))

recurrence_distribution$SNPs <-
  as.numeric(
    recurrence_distribution$SNPs
  )

ggplot(
  recurrence_distribution,
  aes(Loops, SNPs)
) +
  geom_col(
    fill = "grey70",
    width = 0.88
  ) +
  geom_vline(
    xintercept = c(10, 20, 30),
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    title =
      "Frequency spectrum of recurrent SNPs across 100 bootstrap loops",
    subtitle =
      "Dashed lines at 10, 20, 30 loops",
    x =
      "Number of loops a SNP was rediscovered in",
    y =
      "Number of SNPs"
  ) +
  theme_classic(base_size = 12)


library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# FIGURE 5: Pairwise LD among the 17 highest-recurrence SNPs
# ============================================================

# Read PLINK pairwise LD output
# Change the filename if yours is called something else
ld <- fread(
  "recurrent_all_LD.ld",
  data.table = FALSE,
  check.names = FALSE
)

# 17 highest-recurrence candidate SNPs
top17 <- c(
  "SMA_10425984_38060",
  "SMA_10429260_10473",
  "SMA_106974_4066",
  "SMA_176455_9193",
  "SMA_18020_38016",
  "SMA_185331_32429",
  "SMA_427734_5763",
  "SMA_500288_3414",
  "SMA_8183_33268",
  "SMA_90177_8949",
  "SMA_953813_140",
  "SMA_10428241_11763",
  "SMA_140292_6115",
  "SMA_31083_12506",
  "SMA_318388_997",
  "SMA_467548_3209",
  "SMA_97450_20455"
)

# Keep only pairwise LD involving the 17 candidates
ld17 <- ld %>%
  filter(
    SNP_A %in% top17,
    SNP_B %in% top17
  )

# Create all possible SNP pairs
all_pairs <- expand.grid(
  SNP_A = top17,
  SNP_B = top17,
  stringsAsFactors = FALSE
)

# Add directly reported LD values
ld_full <- all_pairs %>%
  left_join(
    ld17 %>%
      select(SNP_A, SNP_B, R2),
    by = c("SNP_A", "SNP_B")
  )

# Add values reported in the opposite orientation
ld_reverse <- ld17 %>%
  transmute(
    SNP_A = SNP_B,
    SNP_B = SNP_A,
    R2_reverse = R2
  )

ld_full <- ld_full %>%
  left_join(
    ld_reverse,
    by = c("SNP_A", "SNP_B")
  ) %>%
  mutate(
    R2 = coalesce(R2, R2_reverse)
  ) %>%
  select(-R2_reverse)

# IMPORTANT:
# This reproduces the old plot appearance by displaying
# unreported pairs as 0.
# If PLINK omitted pairs because they were outside the 1 Mb
# window, they are technically "not calculated", not truly zero.
ld_full$R2[is.na(ld_full$R2)] <- 0

# Add numeric positions so we can keep only one triangle
ld_full <- ld_full %>%
  mutate(
    x_index = match(SNP_A, top17),
    y_index = match(SNP_B, top17)
  ) %>%
  filter(x_index < y_index)

# Set plotting order
ld_full$SNP_A <- factor(
  ld_full$SNP_A,
  levels = top17
)

ld_full$SNP_B <- factor(
  ld_full$SNP_B,
  levels = rev(top17)
)

# Create Figure 5
p_ld <- ggplot(
  ld_full,
  aes(
    x = SNP_A,
    y = SNP_B,
    fill = R2
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.2f", R2)),
    size = 3.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 1),
    name = expression(r^2)
  ) +
  labs(
    title = "Pairwise LD among the 17 highest-recurrence candidate SNPs",
    x = NULL,
    y = NULL
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.title = element_blank(),
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 9
    )
  )

# Show plot
p_ld

# Save journal-quality vector PDF
ggsave(
  "Figure5_pairwise_LD.pdf",
  plot = p_ld,
  device = cairo_pdf,
  width = 180,
  height = 150,
  units = "mm"
)

# Print strongest LD pairs
print(
  ld17 %>%
    arrange(desc(R2)) %>%
    select(SNP_A, SNP_B, R2)
)

# Check if any pair exceeds the high-LD threshold
print(
  ld17 %>%
    filter(R2 > 0.8) %>%
    arrange(desc(R2))
)








# ============================================================
# SUPPLEMENTARY FIGURES S1 AND S2
# Schematic illustrations only — not based on measured data
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)


# ============================================================
# FIGURE S1
# Sävar site schematic:
# 9 rows x 15 positions, arranged as five groups of
# three adjacent clonal ramets
# ============================================================

site_schematic <- expand_grid(
  Row = 1:9,
  Column = 1:15
) %>%
  mutate(
    Genotype_group = ceiling(Column / 3),
    Clone = ((Column - 1) %% 3) + 1
  )


# Lines connecting each trio
trio_lines <- site_schematic %>%
  group_by(Row, Genotype_group) %>%
  summarise(
    x_start = min(Column),
    x_end   = max(Column),
    .groups = "drop"
  )


p_S1 <- ggplot() +
  
  # horizontal lines joining the three clonal ramets
  geom_segment(
    data = trio_lines,
    aes(
      x = x_start,
      xend = x_end,
      y = Row,
      yend = Row
    ),
    colour = "grey75",
    linewidth = 0.5
  ) +
  
  # tree symbols
  geom_point(
    data = site_schematic,
    aes(
      x = Column,
      y = Row
    ),
    shape = 17,
    size = 4.2,
    colour = "darkgreen"
  ) +
  
  # clone numbers 1, 2, 3
  geom_text(
    data = site_schematic,
    aes(
      x = Column,
      y = Row + 0.20,
      label = Clone
    ),
    size = 3.2,
    colour = "grey25"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 4, 7, 10, 13),
    limits = c(0.5, 15.5),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = 1:9,
    limits = c(0.5, 9.6),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Sävar site schematic (9 × 15) – three adjacent clones per genotype",
    x = "Columns",
    y = "Rows"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    ),
    
    axis.title = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    axis.line = element_line(
      linewidth = 0.7
    ),
    
    axis.ticks = element_line(
      linewidth = 0.6
    ),
    
    plot.margin = margin(
      10, 10, 10, 10
    )
  )


# Show S1
p_S1




# ============================================================
# FIGURE S2
# Illustrative random selection:
# one of three clonal ramets selected for each genotype
# in each GWAS loop
#
# This intentionally shows only a small illustrative subset
# rather than all 100 loops / all genotypes.
# ============================================================

set.seed(12345)

n_example_loops <- 12
n_example_genotypes <- 17


selection_example <- expand_grid(
  Loop = 1:n_example_loops,
  Genotype = 1:n_example_genotypes
) %>%
  mutate(
    Clone = sample(
      1:3,
      size = n(),
      replace = TRUE
    ),
    
    Clone_factor = factor(
      Clone,
      levels = c(1, 2, 3),
      labels = c(
        "Clone 1",
        "Clone 2",
        "Clone 3"
      )
    ),
    
    Text_colour = ifelse(
      Clone == 1,
      "white",
      "dark"
    )
  )


# Reverse loop factor so Loop 1 appears at the top
selection_example <- selection_example %>%
  mutate(
    Loop_label = factor(
      paste("Loop", Loop),
      levels = paste(
        "Loop",
        n_example_loops:1
      )
    )
  )


p_S2 <- ggplot(
  selection_example,
  aes(
    x = Genotype,
    y = Loop_label,
    fill = Clone_factor
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.45,
    width = 1,
    height = 1
  ) +
  
  geom_text(
    aes(
      label = Clone,
      colour = Text_colour
    ),
    size = 3.6
  ) +
  
  scale_fill_manual(
    values = c(
      "Clone 1" = "#1B7F3A",
      "Clone 2" = "#7FC97F",
      "Clone 3" = "#D7EED2"
    ),
    name = NULL
  ) +
  
  scale_colour_manual(
    values = c(
      "white" = "white",
      "dark" = "#145A22"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 17),
    limits = c(0.5, 17.5),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Random selection per loop: one clone (1–3) chosen per trio",
    subtitle = "Illustrative subset of the ramet-subsampling procedure",
    x = "Example genotype trios",
    y = "GWAS loops"
  ) +
  
  guides(
    fill = guide_legend(
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      colour = "grey35"
    ),
    
    axis.title = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.width = grid::unit(
      0.7,
      "cm"
    ),
    
    legend.key.height = grid::unit(
      0.45,
      "cm"
    ),
    
    plot.margin = margin(
      10, 10, 10, 10
    )
  )


# Show S2
p_S2


