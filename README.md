# Biomedicinos_duomenu_analize

Šioje repozitorijoje pateikiamos dvi atskiros ir tarpusavyje nesusijusios užduotys: `task0` ir `task4`.

## task0

`task0` dalyje atliekama FDA duomenų analizė.

Pirmiausia duomenys paruošiami bendru duomenų paruošimo skriptu, kad visi analizės žingsniai naudotų tuos pačius objektus. Po to pagrindinis analizės failas įkelia paruoštus duomenis ir paleidžia atskirus analizės skriptus.

Šioje dalyje nagrinėjami simptomai, produktai ir skirtumai tarp pasirinktų laikotarpių. Rezultatai naudojami aprašomajai analizei, grafikams ir išvadoms parengti.

## task4

`task4` dalyje kuriamas epigenetinis laikrodis, t. y. modelis, kuris pagal pilno kraujo DNR metilinimo duomenis prognozuoja individo chronologinį amžių.

Prieš modeliavimą atliekamas duomenų paruošimas: pritaikoma mėginių kokybės kontrolė, paliekami tik visiems duomenų rinkiniams bendri CpG, pašalinami mėginiai be amžiaus informacijos ir paruošiama amžiaus transformacija modeliavimui.

Prieš paleidžiant epigenetinių laikrodžių skriptus, reikia paleisti du paruošiamuosius skriptus iš `task4/scripts/` aplanko:

1. `prepare_data.R`

   Šis skriptas paruošia pradinius metilinimo duomenis modeliavimui. Jis pritaiko kokybės kontrolės filtrus, palieka tik visiems duomenų rinkiniams bendrus CpG, pašalina mėginius be amžiaus informacijos ir išsaugo išvalytas kohortas.

2. `02_select_age_related_cpgs.R`

   Šis skriptas naudoja išvalytas kohortas ir atrenka CpG, kurie labiausiai susiję su amžiumi. Šie požymiai vėliau naudojami epigenetinių laikrodžių modeliams kurti.

Tik atlikus šiuos paruošimo žingsnius galima paleisti atskirus epigenetinių laikrodžių kūrimo skriptus:

- `Justinas.R`
- `Audra.R`
- `Daniel.R`

Modeliai kuriami naudojant regresijos ar mašininio mokymosi metodus.
