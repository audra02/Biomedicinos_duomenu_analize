dat <- readRDS("data/dat.rds")

minor <- dat[dat$age < 18, ]
unique(minor$outcome)

minor$severity_score <- 0

minor$severity_score <- minor$severity_score +
  10 * grepl("death", minor$outcome) +
  6 * grepl("life threatening", minor$outcome) +
  5 * grepl("required intervention", minor$outcome) +
  4 * grepl("hospitalization", minor$outcome) +
  4 * grepl("disability", minor$outcome) +
  4 * grepl("congenital anomaly", minor$outcome) +
  2 * grepl("visited emergency room", minor$outcome) +
  1 * grepl("visited a health care provider", minor$outcome) +
  1 * grepl("injury", minor$outcome) +
  1 * grepl("allergic reaction", minor$outcome) +
  1 * grepl("other serious", minor$outcome)

library(dplyr)

product_risk <- minor %>%
  group_by(product) %>%
  summarise(
    cases = n(),
    total_severity = sum(severity_score),
    avg_severity = mean(severity_score)
  ) %>%
  arrange(desc(total_severity))

head(product_risk, 5)

library(ggplot2)
top5 <- head(product_risk, 5)
grafikas <- ggplot(top5, aes(x = reorder(product, total_severity),
                 y = total_severity)) +
  geom_col() +
  geom_text(aes(label = paste0("Cases: ", cases,
                               "\nAvg: ", round(avg_severity, 2))),
            hjust = -0.1,
            size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "Top 5 pavojingiausi produktai nepilnamečiams",
    subtitle = "Bendras sunkumo balas, atvejų skaičius ir vidutinis ligos sunkumas",
    x = "Produktas",
    y = "Bendras ligos sunkumo balas"
  ) +
  theme_minimal()

ggsave("plots/top5_pavojingiausi.png",
       plot = grafikas,
       width = 10,
       height = 6,
       dpi = 300)