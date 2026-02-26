# Biomedicinos_duomenu_analize


Šis repo yra darbas su FDA duomenimis.

## Failų logika


### 1) Duomenų paruošimas: scripts/Task_01.R
scripts/Task_01.R yra **duomenų paruošimo skriptas**, paremtas dėstytojo pateiktu šablonu.
Jo paskirtis vienodai paruošti duomenis visiems komandos nariams, kad analizės būtų daromos su tais pačiais objektais.

Šis skriptas sugeneruoja ir išsaugo paruoštus objektus į data/:
- data/dat.rds – pagrindinė lentelė 
- data/tm.rds – simptomų term-matrix 

### 2) Pagrindinis analizės failas: scripts/Task0.R
scripts/Task0.R yra **pagrindinis (main) analizės failas**,
Įkelia paruoštus objektus, sugeneruotus Task_01.R:
   ```r
   dat <- readRDS("data/dat.rds")
   tm  <- readRDS("data/tm.rds")

Kiekvienas komandos narys raš2ė savo atskirą analizės skriptą, o Task0.R juos sujungė į vieną bendrą darbą.
Naudojami analizės skriptai:
scripts/01_symptomai18.R
scripts/02_top_produktai.R
scripts/03_pre_vs_school.R

## Duomenys
Šiuo metu **duomenys jau yra įkelti į GitHub **:
- HFCS-Quarterly-20250930--CSV_PRODUCT-BASED.csv
- dat.rds
- tm.rds


**jei nereikės pačių duomenų repozitorijoje**, tada duomenis laikysime tik lokaliai tai įkelsimė duomenys į .gitignore, ar tiesiog pratrinsim. 
```gitignore
data/*.csv
data/*.rds

