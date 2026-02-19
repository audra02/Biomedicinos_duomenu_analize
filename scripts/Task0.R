#Task0.R

# analizės paleidimas-sujungimas musu main scriptas
# to failo tikslas isikleti paruostas duomnis iriskves to to kiekvieno is musu dali per source()

#pvz ikelimas duomenu
dat<- readRDS("data/dat.rds")
tm <- readRDS("data/tm.rds") 
# ir po to per source() kiekvieno dalis pvz
source("scripts/01_simptomai.R")
source("scripts/02_top_produktai.R")