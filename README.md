# Biomedicinos_duomenu_analize

Šioje repozitorijoje pateikiamos dvi atskiros ir tarpusavyje nesusijusios užduotys: `task0` ir `task4`.

## task0

`task0` dalyje atliekama FDA duomenų analizė.

Pirmiausia duomenys paruošiami bendru duomenų paruošimo skriptu, kad visi analizės žingsniai naudotų tuos pačius objektus. Po to pagrindinis analizės failas įkelia paruoštus duomenis ir paleidžia atskirus analizės skriptus.

Šioje dalyje nagrinėjami simptomai, produktai ir skirtumai tarp pasirinktų laikotarpių. Rezultatai naudojami aprašomajai analizei, grafikams ir išvadoms parengti.

## task4

`task4` dalyje kuriamas epigenetinis laikrodis, t. y. modelis, kuris pagal pilno kraujo DNR metilinimo duomenis prognozuoja individo chronologinį amžių.

Šioje užduotyje naudojami keli DNR metilinimo duomenų rinkiniai. Prieš modeliavimą atliekamas duomenų paruošimas, kokybės kontrolė, pasirenkami bendri CpG visiems rinkiniams ir paruošiami duomenys amžiaus prognozavimo modeliams.

Svarbi užduoties sąlyga – negalima naudoti jau egzistuojančių epigenetinių laikrodžių. Todėl kiekvienas modelis kuriamas savarankiškai, naudojant regresijos ar mašininio mokymosi metodus.

`task4` dalyje yra epigenetinių laikrodžių kūrimo skriptai:

- `Justinas.R`
- `Audra.R`
- `Daniel.R`

Galutiniame palyginime bus vertinama, kuris iš sukurtų epigenetinių laikrodžių geriausiai prognozuoja amžių. Modeliai bus lyginami pagal prognozavimo paklaidą ir kitus tikslumo įverčius.
