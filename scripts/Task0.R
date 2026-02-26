# scripts/Task0.R
# ------------------------------------------------------------
# PAGRINDINIS KOMANDOS ANALIZĖS FAILAS (master script)
#
# Šio failo tikslas:
#  1) Užkrauti paruoštus duomenis iš:
#       - data/dat.rds  (pagrindinė lentelė su įrašais)
#       - data/tm.rds   (term matrix su simptomais, 0/1)
#  2) Užtikrinti, kad dat ir tm eilutės būtų teisingai suporuotos
#  3) Paleisti kiekvieno komandos nario analizės skriptus per source().
#
# Pastaba:
#  - Šis failas neatlieka analizės pats — jis tik sujungia visus darbus.
#  - Kiekvienas analizės skriptas privalo naudoti jau įkeltus objektus:
#       dat ir tm
# ------------------------------------------------------------

# 1) Užkrauname paruoštus duomenis
dat <- readRDS("data/dat.rds")
tm  <- readRDS("data/tm.rds")

# 2) Jei dat ir tm nesutampa eilučių skaičius, stabdom darbą
#    (tai reiškia, kad suporavimas tarp lentelių sugadintas)
stopifnot(nrow(dat) == nrow(tm))

# 3) Paleidžiame komandos analizes

# 1. Kokiais simptomais dažniausiai skundžiasi nepilnamečiai?
source("scripts/01_symptomai18.R")

# 2. Top 5 pavojingiausi produktai kuriuos reikėtų riboti nepilnamečiams?
source("scripts/02_top_produktai.R")

# 3. Kaip skiriasi nusiskundimai tarp iki mokylinio ir po mokyklinio amžiaus vaikų?
# source("scripts/03_pre_vs_post_school.R")   # įjungti, kai bus paruoštas
