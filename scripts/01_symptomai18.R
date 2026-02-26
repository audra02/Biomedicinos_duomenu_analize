# scripts/01_symptomai18.R

# 2 atrenkam nepilnamečius 
idx <- dat$age < 18

# pasidarom d18 anksti, kad viršuje esantys check'ai veiktų
d18 <- dat[idx, , drop = FALSE]

names(d18)
ls()
names(dat)

head(d18$code, 20)
length(unique(d18$code))

cat("Nepilnamečių įrašų skaičius:", sum(idx), "\n")

# 3) paimam tik nepilnamečių term matrix dalį
tm18 <- tm[idx, , drop = FALSE]

# 4) kiek kartų pasitaiko kiekvienas simptomas 
freq <- colSums(tm18)

# 5) surikiuojam mažėjančiai ir paimam TOP
freq <- sort(freq, decreasing = TRUE)

print(head(freq, 10))

# 6) lentelė su procentais ( 15)
topN <- 15
top <- head(freq, topN)

rez_symptomai18 <- data.frame(
  simptomas = names(top),
  kiekis = as.integer(top),
  procentai = round(100 * as.integer(top) / sum(idx), 2)
)

cat("\nTOP 15 simptomų \n")
print(rez_symptomai18)

## grafikai 
topN <- 15
top <- head(freq, topN)

par(mar = c(5, 18, 4, 2))
dotchart(
  x = as.integer(top),
  labels = names(top),
  pch = 16,
  main = "Pavojingiausi simptomai nepilnamečiams ",
  xlab = "Įrašų skaičius"
)

## ggploto
topN <- 15
top <- head(freq, topN)

df <- data.frame(
  simptomas = names(top),
  kiekis = as.integer(top)
)

df$procentai <- round(100 * df$kiekis / sum(idx), 2)

# kad būtų tvarkingas rikiavimas grafike
df$simptomas <- factor(df$simptomas, levels = df$simptomas[order(df$kiekis)])

ggplot(df, aes(x = simptomas, y = kiekis)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = paste0(procentai, "%")),
            hjust = -0.1, size = 3.5) +
  labs(
    title = "Simptomai nepilnamečiams (age < 18)",
    subtitle = paste0("N = ", sum(idx), " įrašų"),
    x = NULL,
    y = "Įrašų skaičius"
  ) +
  theme_minimal() +
  expand_limits(y = max(df$kiekis) * 1.12)

#Grupiu palyginimas
# Nepilnamečiai berniukai vs mergaitės
# Ar skiriasi TOP 10 simptomų dažniai

# (idx, d18, tm18 jau turim iš viršaus)
# tik male ir female
keep_sex <- d18$sex %in% c("male", "female")
d18_sex <- d18[keep_sex, , drop = FALSE]
tm18_sex <- tm18[keep_sex, , drop = FALSE]

nrow(d18_sex)
print(table(d18_sex$sex))

# TOP 10 simptomų tarp nepilnamečių (male+female)
freq_sex <- sort(colSums(tm18_sex), decreasing = TRUE)

topN_test <- 10
top_symptoms_test <- names(head(freq_sex, topN_test))

print(top_symptoms_test)

# Rezultatų lentel
rez_sex_symptoms <- data.frame(
  simptomas = character(),
  male_kiekis = integer(),
  male_proc = numeric(),
  female_kiekis = integer(),
  female_proc = numeric(),
  p_value = numeric(),
  reiksminga_05 = character(),
  stringsAsFactors = FALSE
)

n_male <- sum(d18_sex$sex == "male")
n_female <- sum(d18_sex$sex == "female")

# Ciklas per simptomus
for (s in top_symptoms_test) {
  
  # ar simptomas yra (TRUE/FALSE)
  x <- tm18_sex[, s] > 0
  
  # 2x2 lentelė: lytis x simptomas (yra/nera)
  tab <- table(d18_sex$sex, x)
  
  cat("Simptomas:", s, "\n")
  print(tab)
  
  # kiek turi simptomą
  male_yes <- sum(x[d18_sex$sex == "male"])
  female_yes <- sum(x[d18_sex$sex == "female"])
  
  # procentai grupėse
  male_proc <- 100 * male_yes / n_male
  female_proc <- 100 * female_yes / n_female
  
  # chi-square expected reikšmės 
  chi_tmp <- suppressWarnings(chisq.test(tab))
  exp_counts <- chi_tmp$expected
  cat("Expected reikšmės (chi-square):\n")
  print(round(exp_counts, 2))
  
  # Kadangi expected didelės -> naudojam chi-square testą
  tst <- chisq.test(tab)
  
  reiksminga <- ifelse(tst$p.value < 0.05, "taip", "ne")
  
  # įrašom rezultatą į lentelę
  rez_sex_symptoms <- rbind(
    rez_sex_symptoms,
    data.frame(
      simptomas = s,
      male_kiekis = as.integer(male_yes),
      male_proc = round(male_proc, 2),
      female_kiekis = as.integer(female_yes),
      female_proc = round(female_proc, 2),
      p_value = round(tst$p.value, 6),
      reiksminga_05 = reiksminga,
      stringsAsFactors = FALSE
    )
  )
}

# Surikiuojam pagal p reikšme
if (nrow(rez_sex_symptoms) > 0) {
  rez_sex_symptoms <- rez_sex_symptoms[order(rez_sex_symptoms$p_value), ]
}

cat("Rezultatai: TOP 10 simptomų palyginimas (chi-square)\n")
cat("Statistiškai reikšminga, jei p < 0.05\n")
print(rez_sex_symptoms)