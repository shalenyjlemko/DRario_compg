library(ggplot2)
#/../data/Info_KaKs/all_kaks.tsv
kaks <- read.table("/Users/raulduran/Documents/M2_GENIOMHE/comparative_genomics/DRario_compg/data/Info_KaKS/all_kaks.tsv",
                   header = TRUE, sep = "\t", check.names = TRUE,
                   comment.char = "")

# See column names
names(kaks)
unique(kaks$Method)
kaks_ma <- subset(kaks, Method == "MA")
kaks_ma <- subset(kaks_ma,
                  Ks > 0 & !is.na(Ka.Ks) & Ka.Ks >= 0 & Ka.Ks < 5)
summary(kaks_ma$Ka.Ks)
hist(kaks_ma$Ka.Ks,
     breaks = 30,
     main   = "Ka/Ks distribution – Danio rerio paralogs",
     xlab   = "Ka/Ks",
     ylab   = "Number of pairs")

abline(v = 1, col = "red", lty = 2)  # neutral line
plot(density(kaks_ma$Ka.Ks),
     main = "Density of Ka/Ks ratios – Danio rerio paralogs",
     xlab = "Ka/Ks")

abline(v = 1, col = "red", lty = 2)

hist(kaks_ma$Ks,
     breaks = 40,  # more or fewer bins if you like
     main   = "Ks distribution – Danio rerio duplicate pairs",
     xlab   = "Synonymous divergence (Ks, substitutions/site)",
     ylab   = "Number of duplicate pairs")
abline(v = 2, col = "red", lty = 2)

k6 <- subset(kaks_ma,
                  Ks > 2)
hist(k6$Ks,
     breaks = 40,  # more or fewer bins if you like
     main   = "Ks distribution – Danio rerio duplicate pairs",
     xlab   = "Synonymous divergence (Ks, substitutions/site)",
     ylab   = "Number of duplicate pairs")
abline(v = 2, col = "red", lty = 2)

k7 <- subset(kaks_ma,
             Ks < 2)
hist(k7$Ks,
     breaks = 40,  # more or fewer bins if you like
     main   = "Ks distribution – Danio rerio duplicate pairs",
     xlab   = "Synonymous divergence (Ks, substitutions/site)",
     ylab   = "Number of duplicate pairs")
abline(v = 2, col = "red", lty = 2)

k4 <- subset(kaks,
             Ka > 0 & Ka < 2 &
               Ks > 0 & Ks < 5)

plot(k4$Ks, k4$Ka,
     pch = 16, cex = 0.4,
     xlab = "Ks",
     ylab = "Ka",
     main = "Ka vs Ks – Duplicate pairs")
abline(a = 0, b = 1, col = "red", lty = 2)  # Ka = Ks line

k5 <- subset(kaks, Divergence.Time > 0 & Divergence.Time < 500)  # adjust cutoff

hist(k5$Divergence.Time,
     breaks = 40,
     main   = "Distribution of divergence times",
     xlab   = "Divergence time (units from KaKs_Calculator)",
     ylab   = "Number of duplicate pairs")

# GC.1.2.3. looks like: "0.501066(0.541578:0.420043:0.541578)"
# we’ll extract overall GC and positions if needed
library(tidyr)
library(dplyr)

gc_df <- kaks %>%
  mutate(GC_all = as.numeric(sub("\\(.*", "", GC.1.2.3.)))  # before first "("

k_gc <- subset(gc_df,
               Ks > 0 & Ks < 5 & !is.na(GC_all))

plot(k_gc$Ks, k_gc$GC_all,
     pch = 16, cex = 0.4,
     xlab = "Ks",
     ylab = "GC content (overall)",
     main = "GC vs Ks")

abline(h = mean(k_gc$GC_all, na.rm = TRUE), col = "red", lty = 2)

k6 <- subset(kaks, Ks > 0 & Ks < 5)

plot(k6$Ks, k6$Substitutions,
     pch = 16, cex = 0.4,
     xlab = "Ks",
     ylab = "Total substitutions",
     main = "Substitutions vs Ks – saturation check")

library(ggplot2)

ggplot(kaks_ma, aes(x = Ks)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(
    title = "Ks distribution – Danio rerio duplicate pairs",
    x     = "Synonymous divergence (Ks)",
    y     = "Number of duplicate pairs"
  ) +
  theme_minimal()


k2 <- subset(kaks,
             Ks > 0 & Ks < 5 &
               !is.na(Ka.Ks) & Ka.Ks >= 0 & Ka.Ks < 3)

ggplot(k2, aes(x = Ks, y = Ka.Ks)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Ka/Ks vs Ks",
       x = "Ks",
       y = "Ka/Ks") +
  theme_minimal()








