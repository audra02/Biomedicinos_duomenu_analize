# Biomedicinos duomenų analizė

Šioje repozitorijoje pateikiamos dvi atskiros užduotys:

* **task0** – FDA duomenų analizė.
* **task4** – epigenetinio laikrodžio kūrimas naudojant DNR metilinimo duomenis.

---

# Task 0

Šioje užduotyje atliekama FDA nepageidaujamų reiškinių duomenų analizė.

Duomenys pirmiausia paruošiami bendru duomenų paruošimo skriptu, kad visi analizės žingsniai naudotų vienodus duomenų objektus. Toliau vykdoma simptomų, produktų ir skirtingų laikotarpių analizė.

Rezultatai naudojami:

* aprašomajai statistinei analizei;
* grafikų sudarymui;
* rezultatų interpretacijai;
* išvadų formulavimui.

---

# Task 4

Šioje užduotyje kuriamas epigenetinis laikrodis – modelis, prognozuojantis individo chronologinį amžių pagal pilno kraujo DNR metilinimo duomenis.

Analizei buvo naudojamos septynios nepriklausomos kohortos:

* Fraga
* Johansson
* Kalyakulina
* Kurushima
* Mahdi
* Quinn
* Xu

Po kokybės kontrolės galutiniame duomenų rinkinyje liko:

* 2235 mėginiai;
* 369543 bendri CpG lokusai.

---

# Duomenų paruošimas

Prieš kuriant modelius būtina paleisti du paruošiamuosius skriptus:

```r
source("prepare_data.R")
source("02_select_age_related_cpgs.R")
```

## prepare_data.R

Šis skriptas:

* atlieka mėginių kokybės kontrolę;
* pašalina nekokybiškus mėginius;
* palieka tik visoms kohortoms bendrus CpG;
* pašalina mėginius be amžiaus informacijos;
* sukuria išvalytų kohortų failus tolimesnei analizei.

## 02_select_age_related_cpgs.R

Šis skriptas:

* apskaičiuoja koreliacijas tarp CpG ir amžiaus;
* atrenka su amžiumi labiausiai susijusius CpG;
* sudaro CpG reitingus tolimesniam modeliavimui.

Tik atlikus šiuos žingsnius galima kurti epigenetinius laikrodžius.

---

# Epigenetinių laikrodžių modeliai

Po duomenų paruošimo paleidžiami modelių kūrimo skriptai:

```r
source("Justinas.R")
source("Audra.R")
```

Buvo išbandyti keli regresijos ir mašininio mokymosi metodai:

* PCA + Linear Regression
* PCA + kNN Regression
* Random Forest
* M-value Ridge Regression
* Signature Kernel

Modeliai buvo vertinami naudojant:

* MAE (Mean Absolute Error);
* RMSE (Root Mean Squared Error);
* Pearson koreliaciją.

Vertinimui naudotas 80 % mokymo ir 20 % testinis duomenų padalinimas kiekvienoje kohortoje.

---

# Rezultatai

Geriausią rezultatą pasiekė modelis:

**m_value_ridge_calibrated**

Modelio rezultatai testinėje aibėje:

| Metrika     | Reikšmė |
| ----------- | ------: |
| MAE         |   2.603 |
| RMSE        |   3.869 |
| Koreliacija |   0.986 |

Tai reiškia, kad modelis vidutiniškai klysta maždaug 2,6 metų prognozuodamas individo amžių pagal DNR metilinimo duomenis.

Papildomai nustatyta, kad Quinn kohortoje (0–3 metų amžiaus tiriamieji) modeliai pasiekė itin aukštą tikslumą. M-value Ridge modelio vidutinė absoliuti paklaida šioje kohortoje siekė vos 0,155 metų.

---

# Projekto struktūra

```text
task4/
├── data/
├── scripts/
│   ├── prepare_data.R
│   ├── 02_select_age_related_cpgs.R
│   ├── Justinas.R
│   └── Audra.R
├── results/
│   ├── cleaned_cohorts/
│   ├── age_cpg_selection/
│   └── model_comparison/
```

Visi modelių rezultatai, grafikai ir tarpinių analizės žingsnių failai išsaugomi `results/` kataloge.
