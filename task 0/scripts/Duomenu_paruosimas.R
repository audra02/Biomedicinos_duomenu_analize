## scripts/Task_01.R
# ------------------------------------------------------------
# DUOMENŲ PARUOŠIMO SKRIPTAS
#
# Šio failo paskirtis:
#  1) Užkrauti pradinį CSV failą (HFCS duomenis).
#  2) Išvalyti ir sutvarkyti duomenis analizei.
#  3) Sukurti term matrix (tm) objektą iš simptomų stulpelio.
#  4) Išsaugoti du paruoštus failus:
#        - data/dat.rds  (pagrindinė lentelė)
#        - data/tm.rds   (simptomų 0/1 matrica)
#
# SVARBU:
#  - Analizės skriptai turi prasidėti nuo:
#        dat <- readRDS("data/dat.rds")
#        tm  <- readRDS("data/tm.rds")
#  - Visi komandos nariai generuoja tuos pačius dat.rds ir tm.rds failus.
#
# Komandos struktūra:
#   scripts/01_symptomai18.R
#   scripts/02_top_produktai.R
#   scripts/03_pre_vs_post_school.R
#
# Bendras skriptas scripts/Task0.R tik sujungia visas analizes per source().
# ------------------------------------------------------------

Sys.setlocale("LC_ALL", "C")
dat <- read.csv("data/HFCS-Quarterly-20250930--CSV_PRODUCT-BASED.csv")

# Perpavadiname stulpelius
colnames(dat) <- c(
  "date_report",  # Date_FDA_First_Received_Case
  "id",           # Report_ID
  "date_event",   # Date_Event
  "suspect",      # Product_Type (suspect/concomitant)
  "product",      # Product
  "code",         # Product_Code
  "category",     # Description (industry description)
  "age",          # Patient_Age
  "age_units",    # Age_Units
  "sex",          # Sex
  "symptoms",     # Case_MedDRA_Preferred_Terms
  "outcome"       # Case_Outcome
)

# Konvertuojame tekstinius laukus į mažąsias raides
# (tai padeda išvengti problemų su skirtingu raidžių registru)
dat$suspect  <- tolower(as.character(dat$suspect))
dat$product  <- tolower(as.character(dat$product))
dat$category <- tolower(as.character(dat$category))
dat$sex      <- tolower(as.character(dat$sex))
dat$symptoms <- tolower(as.character(dat$symptoms))
dat$outcome  <- tolower(as.character(dat$outcome))

# Pašaliname pasikartojančius įrašus pagal unikalų ID
dat <- dat[!duplicated(dat$id), ]

# Sutvarkome lyties kintamąjį:
# nežinomos ar neapibrėžtos reikšmės suvedamos į "unknown"
dat$sex[dat$sex %in% c("", "not reported", "unknown",
                       "unspecified (or another gender identity)")] <- "unknown"

dat$age <- suppressWarnings(as.numeric(as.character(dat$age)))

# Skirtingi vienetai konvertuojami į metus
dat$age[dat$age_units == "week(s)"]   <- dat$age[dat$age_units == "week(s)"] / 52
dat$age[dat$age_units == "month(s)"]  <- dat$age[dat$age_units == "month(s)"] / 12
dat$age[dat$age_units == "day(s)"]    <- dat$age[dat$age_units == "day(s)"] / 365
dat$age[dat$age_units == "decade(s)"] <- dat$age[dat$age_units == "decade(s)"] * 10

# Pašaliname įrašus be amžiaus arba su nerealiais amžiais
dat <- dat[!is.na(dat$age), ]
dat <- dat[dat$age < 600, ]

# Amžiaus vienetų stulpelis nebenaudojamas
dat$age_units <- NULL

# Simptomų stulpelis suskaidomas pagal kablelį
# Kiekvienas simptomas tampa atskiru stulpeliu (0/1)
syms <- strsplit(dat$symptoms, ", ")

# Sukuriame matricą:
#   - eilutė = vienas atvejis
#   - stulpelis = konkretus simptomas
tm <- matrix(0, 
             nrow = length(syms), 
             ncol = length(unique(unlist(syms))))

colnames(tm) <- unique(unlist(syms))

# Užpildome matricą: jei simptomas yra – žymime 1
for(i in 1:length(syms)) {
  tm[i, syms[[i]]] <- 1
}

# Pašaliname simptomų tekstinį stulpelį iš pagrindinės lentelės
dat$symptoms <- NULL
dat$symptoms <- NULL

# Išsaugome paruoštus objektus, kuriuos naudos analizės skriptai
saveRDS(dat, "data/dat.rds")
saveRDS(tm,  "data/tm.rds")
