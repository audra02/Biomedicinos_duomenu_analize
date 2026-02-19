# scripts/01_symptomai18.R

library(ggplot2)

names(d18)
ls()
names(dat)

head(d18$code, 20)
length(unique(d18$code))

       
# 2) atrenkam nepilnamečius (tas pats indeksas taikomas ir dat, ir tm eilutėms)
idx <- dat$age < 18
       
cat("Nepilnamečių įrašų skaičius:", sum(idx), "\n")
       
# 3) paimam tik nepilnamečių term matrix dalį
tm18 <- tm[idx, , drop = FALSE]
       
# 4) kiek kartų pasitaiko kiekvienas simptomas (stulpelių sumos)
 freq <- colSums(tm18)
       
# 5) surikiuojam mažėjančiai ir paimam TOP
    freq <- sort(freq, decreasing = TRUE)
       
cat("\nTOP 10 simptomų:\n")
print(head(freq, 10))
       
# 6) lentelė su procentais (TOP 15)
topN <- 15
top <- head(freq, topN)
       
       rez_symptomai18 <- data.frame(
         simptomas = names(top),
         kiekis = as.integer(top),
         procentai = round(100 * as.integer(top) / sum(idx), 2)
       )
       
cat("\nTOP 15 simptomų (su procentais nuo visų nepilnamečių įrašų):\n")
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