## scripts/task_01.R
# Tai cia tikriausiai bus komandinis failas, kur bus duomenų paruošimas, ikelsime destytojo koda. 

# Kaip supratu kiekvienam reikes atsiusti CSV (duomenu) failą  ir kiekvienas laikys ta csv faila pas save lokaliai, 
# nes duomenu tikrai i githuba nelabai ir reikia kelti 
# Šis skriptas sugeneruojs du failus lokaliai:
#    - data/dat.rds
#    - data/tm.rds
#  Analizės skriptai turi prasidėti nuo:
#      dat <- readRDS("data/dat.rds")
#      tm  <- readRDS("data/tm.rds")


# ir manau, nezinau kaip cia geriau padaryti, pvz kiekvienas daro savo analizės failą 
#     scripts/01_symptomai18.R
#     scripts/02_top_produktai.R
#     scripts/03_pre_vs_school.R
#  o bendras scripts/task0.R tik sujungia viską:





#######

# duomenų paruošimas komandai

#--- init
Sys.setlocale("LC_ALL", "C")

#--- read ----------------------------------------------------------------------

dat <- read.csv("data/HFCS-Quarterly-20250930--CSV_PRODUCT-BASED.csv")

#--- clean ---------------------------------------------------------------------

#rename columns 
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

# to lower (apsaugom nuo NA)
dat$suspect  <- tolower(as.character(dat$suspect))
dat$product  <- tolower(as.character(dat$product))
dat$category <- tolower(as.character(dat$category))
dat$sex      <- tolower(as.character(dat$sex))
dat$symptoms <- tolower(as.character(dat$symptoms))
dat$outcome  <- tolower(as.character(dat$outcome))

# remove duplicated ids
dat <- dat[!duplicated(dat$id), ]

# fix sex
dat$sex[dat$sex %in% c("", "not reported", "unknown",
                       "unspecified (or another gender identity)")] <- "unknown"

# fix age units -> years
dat$age <- suppressWarnings(as.numeric(as.character(dat$age)))

dat$age[dat$age_units == "week(s)"]   <- dat$age[dat$age_units == "week(s)"] / 52
dat$age[dat$age_units == "month(s)"]  <- dat$age[dat$age_units == "month(s)"] / 12
dat$age[dat$age_units == "day(s)"]    <- dat$age[dat$age_units == "day(s)"] / 365
dat$age[dat$age_units == "decade(s)"] <- dat$age[dat$age_units == "decade(s)"] * 10

dat <- dat[!is.na(dat$age), ]   # remove entries with unknown age
dat <- dat[dat$age < 600, ]     # remove extreme age values (just one sample)

dat$age_units <- NULL

#--- create term matrix --------------------------------------------------------
# Patikimiau splitinti: kablelis + bet koks tarpas
syms <- strsplit(dat$symptoms, ", ")
tm   <- matrix(0, nrow = length(syms), ncol = length(unique(unlist(syms))))
colnames(tm) <- unique(unlist(syms))

for(i in 1:length(syms)) {
  tm[i, syms[[i]]] <- 1
}

dat$symptoms <- NULL
dat$symptoms <- NULL

#--- save ----------------------------------------------------------------------

saveRDS(dat, "data/dat.rds")
saveRDS(tm,  "data/tm.rds")

