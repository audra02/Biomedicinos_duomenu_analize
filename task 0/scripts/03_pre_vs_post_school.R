# scripts/03_pre_vs_post_school.R
# 3) Kaip skiriasi nusiskundimai tarp ikimokyklinio ir mokyklinio amžiaus vaikų?
# Pagrindinė analizė: 2 grupės (0–6) vs (7–18)
# Papildoma analizė: 4 grupės (0–6, 7–11, 12–14, 15–18)
# Papildomai: amžius × lytis

library(ggplot2)

dir.create("plots", showWarnings = FALSE)

# Atrenkame vaikus iki 18 imtinai (0–18)
idx <- dat$age <= 18
d_minor <- dat[idx, , drop = FALSE]
tm_minor <- tm[idx, , drop = FALSE]

cat("Vaikų (0–18) įrašų skaičius:", sum(idx), "\n")

# Amžiaus grupės (pagrindinė analizė)
# ikimokyklinis: 0–6
# mokyklinis: 7–18
d_minor$age_group <- ifelse(d_minor$age <= 6, "ikimokyklinis (0–6)", "mokyklinis (7–18)")

cat("\nPagrindinės 2 grupės (N):\n")
print(table(d_minor$age_group))

# Term matrix padalinimas (suporavimas išlaikytas)
tm_pre  <- tm_minor[d_minor$age_group == "ikimokyklinis (0–6)", , drop = FALSE]
tm_post <- tm_minor[d_minor$age_group == "mokyklinis (7–18)", , drop = FALSE]

n_pre  <- nrow(tm_pre)
n_post <- nrow(tm_post)

# Dažniai ir procentai kiekvienoje grupėje
freq_pre  <- colSums(tm_pre)
freq_post <- colSums(tm_post)

pct_pre  <- 100 * freq_pre  / n_pre
pct_post <- 100 * freq_post / n_post

rez_age <- data.frame(
  simptomas = colnames(tm_minor),
  pre_proc = round(as.numeric(pct_pre), 2),
  post_proc = round(as.numeric(pct_post), 2),
  diff_post_minus_pre = round(as.numeric(pct_post - pct_pre), 2),
  stringsAsFactors = FALSE
)

# TOP 15
rez_sorted <- rez_age[order(abs(rez_age$diff_post_minus_pre), decreasing = TRUE), ]
top_diff <- head(rez_sorted, 15)

cat("\nTOP 15 simptomų su didžiausiu procentiniu skirtumu (mokyklinis - ikimokyklinis):\n")
print(top_diff)

# Grafikas: TOP 10 simptomų
overall_freq <- freq_pre + freq_post
top_symptoms <- names(sort(overall_freq, decreasing = TRUE))[1:10]

plot_data <- data.frame(
  simptomas = rep(top_symptoms, 2),
  procentai = c(pct_pre[top_symptoms], pct_post[top_symptoms]),
  grupe = rep(c("ikimokyklinis (0–6)", "mokyklinis (7–18)"),
              each = length(top_symptoms))
)

plot_data$simptomas <- factor(plot_data$simptomas, levels = rev(top_symptoms))

g_top10 <- ggplot(plot_data, aes(x = simptomas, y = procentai, fill = grupe)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Simptomų palyginimas pagal amžiaus grupę",
    subtitle = "Rodoma procentinė dalis kiekvienoje grupėje (TOP 10 simptomų)",
    x = NULL,
    y = "Procentai (%)",
    fill = "Amžiaus grupė"
  ) +
  theme_minimal()

print(g_top10)

ggsave("plots/03_amziaus_palyginimas_top10.png",
       plot = g_top10,
       width = 10,
       height = 6,
       dpi = 300)

# Chi-square test TOP 10 simptomų (2 grupės)
rez_test <- data.frame(simptomas = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (s in top_symptoms) {
  x <- tm_minor[, s] > 0
  tab <- table(d_minor$age_group, x)
  tst <- suppressWarnings(chisq.test(tab))
  rez_test <- rbind(rez_test, data.frame(simptomas = s, p_value = round(tst$p.value, 6)))
}

cat("\nChi-square test (TOP 10 simptomų, 2 grupės):\n")
print(rez_test)

# Papildoma analizė: amžius × lytis (2 grupės)
cat("\n--- Papildoma analizė: amžius × lytis (2 grupės) ---\n")

keep_sex <- d_minor$sex %in% c("male", "female")
d_sex <- d_minor[keep_sex, , drop = FALSE]
tm_sex <- tm_minor[keep_sex, , drop = FALSE]

cat("Įrašų su aiškia lytimi (male/female):", nrow(d_sex), "\n")
cat("Grupių dydžiai (amžius × lytis):\n")
print(table(d_sex$age_group, d_sex$sex))

top6_sym <- head(top_diff$simptomas, 6)

rez_age_sex <- data.frame()

for (s in top6_sym) {
  x <- tm_sex[, s] > 0
  
  for (ag in unique(d_sex$age_group)) {
    for (sx in c("male", "female")) {
      
      idx_g <- (d_sex$age_group == ag) & (d_sex$sex == sx)
      n_g <- sum(idx_g)
      p_g <- ifelse(n_g > 0, 100 * sum(x[idx_g]) / n_g, NA)
      
      rez_age_sex <- rbind(rez_age_sex, data.frame(
        simptomas = s,
        age_group = ag,
        sex = sx,
        procentai = round(p_g, 2),
        N = n_g,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\nTOP 6 simptomai: procentai pagal amžių ir lytį (2 grupės):\n")
print(rez_age_sex)

rez_age_sex$simptomas <- factor(rez_age_sex$simptomas, levels = rev(top6_sym))

g_sex <- ggplot(rez_age_sex, aes(x = simptomas, y = procentai, fill = age_group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ sex) +
  labs(
    title = "Amžiaus grupių skirtumai pagal lytį",
    subtitle = "Rodoma procentinė dalis (TOP 6 simptomai su didžiausiu amžiaus skirtumu)",
    x = NULL,
    y = "Procentai (%)",
    fill = "Amžiaus grupė"
  ) +
  theme_minimal()

print(g_sex)

ggsave("plots/03_amzius_x_lytis_top6.png",
       plot = g_sex,
       width = 11,
       height = 6,
       dpi = 300)

# Papildoma analizė: 4 grupės (ikimokyklinukai, pradinukai, progimnazistai, gimnazistai)
cat("\n--- Papildoma analizė: 4 amžiaus grupės ---\n")

d_minor$age_group4 <- cut(
  d_minor$age,
  breaks = c(-Inf, 6, 11, 14, 18),
  labels = c("ikimokyklinukai (0–6)", "pradinukai (7–11)", "progimnazistai (12–14)", "gimnazistai (15–18)"),
  right = TRUE
)

cat("4 grupių dydžiai (N):\n")
print(table(d_minor$age_group4, useNA = "ifany"))

keep_agegrp <- !is.na(d_minor$age_group4)
d4 <- d_minor[keep_agegrp, , drop = FALSE]
tm4 <- tm_minor[keep_agegrp, , drop = FALSE]

age_levels <- levels(d4$age_group4)

freq_by_group <- sapply(age_levels, function(g) colSums(tm4[d4$age_group4 == g, , drop = FALSE]))
n_by_group <- as.integer(table(d4$age_group4)[age_levels])

pct_by_group <- sweep(freq_by_group, 2, n_by_group, FUN = "/") * 100
pct_by_group <- round(pct_by_group, 2)

overall4 <- rowSums(freq_by_group)
top_symptoms4 <- names(sort(overall4, decreasing = TRUE))[1:10]

plot_data4 <- do.call(rbind, lapply(age_levels, function(g) {
  data.frame(
    simptomas = top_symptoms4,
    procentai = pct_by_group[top_symptoms4, g],
    grupe = g,
    stringsAsFactors = FALSE
  )
}))

plot_data4$simptomas <- factor(plot_data4$simptomas, levels = rev(top_symptoms4))
plot_data4$grupe <- factor(plot_data4$grupe, levels = age_levels)

g4 <- ggplot(plot_data4, aes(x = simptomas, y = procentai, fill = grupe)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Simptomų palyginimas pagal 4 amžiaus grupes",
    subtitle = "Rodoma procentinė dalis kiekvienoje grupėje (TOP 10 simptomų)",
    x = NULL,
    y = "Procentai (%)",
    fill = "Amžiaus grupė"
  ) +
  theme_minimal()

print(g4)

ggsave("plots/03_4grupes_top10.png",
       plot = g4,
       width = 12,
       height = 6,
       dpi = 300)