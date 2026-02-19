# Biomedicinos_duomenu_analize


Šis repo yra darbas su FDA duomenimis.

## Failų logika

### `scripts/Task_01.R` – duomenų paruošimas (dėstytojo šablonas)
`Task_01.R` yra **duomenų paruošimo skriptas**, paremtas dėstytojo duotu kodu.  
Jis skirtas tik tam, kad visi turėtų vienodai paruoštus objektus analizėms.

## Duomenys (`data/`) ir GitHub
Šiuo metu **duomenys jau yra įkelti į GitHub (netyčia xd)**:
- `data/HFCS-Quarterly-20250930--CSV_PRODUCT-BASED.csv`
- `data/dat.rds`
- `data/tm.rds`

**jei nereikės pačių duomenų repozitorijoje**, tada duomenis laikysime tik lokaliai tai įkelsimė duomenys į .gitignore, ar tiesiog pratrinsim. 
```gitignore
data/*.csv
data/*.rds

Paleidus `Task_01.R` sugeneruojami lokaliai:
- `data/dat.rds`
- `data/tm.rds`
## `Task0.R` (main analizės paleidėjas)

`scripts/Task0.R` yra **pagrindinis (main) analizės failas**, kuris,
 įkelia paruoštus objektus (sugeneruotus iš `Task_01.R`):
```r
dat <- readRDS("data/dat.rds")
tm  <- readRDS("data/tm.rds")

paleidžia kiekvieno komandos nario analizės dalį per source(...).
Idėja tokia: kiekvienas  rašo savo atskirą skriptą
scripts/01_symptomai18.R
scripts/02_top_produktai.R
scripts/03_pre_vs_school.R
O Task0.R juos sujungia

RStudio Console:
```r
source("scripts/Task_01.R")
