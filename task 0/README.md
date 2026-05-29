# Biomedicinos_duomenu_analize

Šis repo yra darbas su FDA duomenimis.

## Failų logika

### Duomenų paruošimas: scripts/Duomenu_paruosimas.R

scripts/Duomenu_paruosimas.R yra **duomenų paruošimo skriptas**, paremtas dėstytojo pateiktu šablonu.
Jo paskirtis – vienodai paruošti duomenis visiems komandos nariams, kad analizės būtų daromos su tais pačiais objektais.

Šis skriptas sugeneruoja ir išsaugo paruoštus objektus į data/:

* data/dat.rds – pagrindinė lentelė
* data/tm.rds – simptomų term-matrix

### Pagrindinis analizės failas: scripts/Task0.R

scripts/Task0.R yra **pagrindinis (main) analizės failas**,
įkelia paruoštus objektus, sugeneruotus Duomenu_paruosimas.R:

dat <- readRDS("data/dat.rds")
tm  <- readRDS("data/tm.rds")

Kiekvienas komandos narys rašė savo atskirą analizės skriptą, o Duomenu_paruosimas.R juos sujungė į vieną bendrą darbą.

Naudojami analizės skriptai:

* scripts/01_symptomai18.R
* scripts/02_top_produktai.R
* scripts/03_pre_vs_post_school.R

## Duomenys

Šiuo metu **duomenys jau yra įkelti į GitHub**:

* HFCS-Quarterly-20250930--CSV_PRODUCT-BASED.csv
* dat.rds
* tm.rds
