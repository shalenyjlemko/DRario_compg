library(ggplot2)

kaks <- read.table("/../data/Info_KaKs/all_kaks.tsv",
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
