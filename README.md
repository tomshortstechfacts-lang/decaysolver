# decaysolver

Bibliothèque C++20 qui calcule l'évolution temporelle d'une chaîne de décroissance radioactive
(équations de Bateman) par plusieurs méthodes indépendantes, et qui documente la vérification de
ses résultats. Destinée aux ingénieurs qui vieillissent des inventaires radiologiques et veulent
savoir ce que vaut le chiffre obtenu.

[![CI](https://github.com/tomshortstechfacts-lang/decaysolver/actions/workflows/ci.yml/badge.svg)](https://github.com/tomshortstechfacts-lang/decaysolver/actions/workflows/ci.yml)
[![Licence CeCILL-C](https://img.shields.io/badge/licence-CeCILL--C-blue.svg)](LICENSE)

> **État du projet : lot 0, socle.** Le cœur numérique n'est pas encore livré. Ce README décrit
> le périmètre visé et marque explicitement ce qui existe (✅) et ce qui est prévu (⬜).

| Lot | Contenu | État |
|---|---|---|
| 0 | Structure, CMake strict, conversions d'unités, provenance, CI | ✅ |
| 1 | Bibliothèque de nucléides, Bateman analytique avec cas dégénérés, Euler explicite/implicite, Crank–Nicolson, RK4, mode inventaire et ligne de commande | ✅ |
| 2 | Oracle multiprécision, ordres de convergence mesurés, cas dégénérés, rapport de vérification | ⬜ |
| 3 | Exponentielle de matrice, non-régression, sanitizers, couverture, bibliothèque complète | ⬜ |
| 4 | Évaluation croisée avec `radioactivedecay`, DOI | ⬜ |

## 1. Le problème mathématique

Pour une chaîne de $n$ nucléides, la population $N_i(t)$ du nucléide $i$ obéit à

$$
\frac{\mathrm{d}N_i}{\mathrm{d}t} = \sum_{j} b_{ji}\,\lambda_j\,N_j \;-\; \lambda_i\,N_i ,
\qquad \lambda_i = \frac{\ln 2}{T_{1/2,i}} ,
$$

où $b_{ji}$ est le rapport d'embranchement de $j$ vers $i$ ($\sum_i b_{ji} = 1$).
Sous forme matricielle, $\dot{\mathbf N} = A\,\mathbf N$, $\mathbf N(0) = \mathbf N_0$, dont la
solution exacte est $\mathbf N(t) = e^{At}\,\mathbf N_0$. L'activité est $A_i = \lambda_i N_i$.

Propriétés de $A$ exploitées et testées : $A_{ii} = -\lambda_i \le 0$ et $A_{ij} \ge 0$ pour
$i \ne j$ (matrice de Metzler), donc $e^{At}$ est positive ; si aucune voie ne quitte le système,
les colonnes de $A$ somment à zéro et le nombre total d'atomes est conservé.

## 2. Schémas prévus

| Méthode | Ordre | Quand l'utiliser | État |
|---|---|---|---|
| Formule analytique de Bateman, par différences divisées de l'exponentielle | exacte | référence ; les $\lambda$ proches ou égaux sont traités par série de Taylor (McCurdy, Ng & Parlett) au lieu de la somme fermée, qui perd ses chiffres par annulation | ✅ |
| Exponentielle de matrice (SciPy, Padé + scaling-and-squaring) | exacte | référence indépendante, côté Python | ⬜ lot 3 |
| Euler implicite | 1 | problèmes raides ; L-stable, positivité garantie (A de Metzler) | ✅ |
| Crank–Nicolson | 2 | non raide seulement : A-stable mais **non L-stable**, produit des concentrations négatives sur les modes raides | ✅ |
| Euler explicite | 1 | pédagogique ; stable si $h\lambda_{max} \le 2$ | ✅ |
| RK4 explicite | 4 | problèmes non raides, mesure d'ordre ; stable si $h\lambda_{max} \le 2{,}785$ | ✅ |

Les schémas implicites se résolvent par substitution avant : dans l'ordre topologique, $A$ est
triangulaire inférieure.

## 3. Hypothèses

- Décroissance spontanée seule : coefficients constants, système linéaire.
- **Mode inventaire, convention des filles.** Par défaut (`input-only`), la sortie ne liste que
  les nucléides de l'entrée : les filles hors liste sont calculées mais masquées, comme dans les
  inventaires industriels qui ne listent pas les filles à vie courte en équilibre (Y-90 sous
  Sr-90, Ba-137m sous Cs-137). Avec `all`, elles apparaissent et l'activité totale change. Les deux
  résultats sont exacts ; ils ne répondent pas à la même question.
- Les totaux α et β/γ classent chaque nucléide par son **mode de décroissance principal**
  (Pu-241, β⁻ à 99,998 %, compte en β/γ).
- Un nucléide stable n'a pas d'activité : il n'apparaît ni en entrée ni en sortie du mode
  inventaire.
- Rapports d'embranchement et demi-vies tirés d'une bibliothèque évaluée, avec leur provenance.
- Convention d'année : **année julienne, 365,25 j = 31 557 600 s**.
- Arithmétique IEEE-754 `double` stricte : le projet **refuse `-ffast-math`** et désactive la
  contraction FMA. Ces options changent l'ordre et l'arrondi des opérations, et rendent caduc tout
  raisonnement sur l'erreur d'arrondi (voir `cmake/FloatingPoint.cmake`).

## 4. Limitations

- Pas de flux neutronique, pas de capture, pas de fission : ce n'est pas un code d'évolution
  sous irradiation.
- Pas de transport, pas de dosimétrie, pas de spectres d'émission.
- Pas de conversion activité → masse, pas de classement réglementaire des déchets.
- La méthode CRAM, référence des codes industriels pour les systèmes très raides, n'est **pas**
  implémentée.
- Pas de parallélisme.

## 5. Domaine de validité

Premiers chiffres, à consolider dans le rapport de vérification (lot 2) :

- **Solution analytique.** Cas sains (D4, $\lambda = (1,2,3)$) : erreur relative $\le 10^{-14}$.
  Cas quasi dégénérés (D2, écarts $10^{-7}$ ; D3, écarts $10^{-11}$) : $\le 10^{-12}$ grâce aux
  différences divisées ; la somme fermée classique donne respectivement $2\cdot10^{-3}$ d'erreur et
  **exactement 0**. Chaîne du Ra-226 (constantes étalées sur 14 ordres de grandeur) : semi-groupe
  vérifié à $10^{-10}$.
- **Conservation du nombre d'atomes** : exacte à $10^{-14}$ sur des données à rapports
  d'embranchement exacts ; sur ICRP-107, limitée par les arrondis des rapports (Bi-210 :
  $1 + 1{,}3\cdot10^{-6}$), et par les voies de fission spontanée non suivies ($\le 1{,}4\cdot10^{-6}$).
- **Intégrateurs** : ordres observés 1, 1, 2, 4 à $\pm 0{,}1$ sur le problème non raide
  $\lambda = (1,2,3,0)$, $T = 4$. En régime raide, seul Euler implicite est utilisable.
- **Données** : les incertitudes des demi-vies évaluées (non fournies par ICRP-107) dépassent
  l'erreur numérique ; voir `data/PROVENANCE.md`.

## 6. Vérification

Le vocabulaire suit le guide ASN n°28 (*Qualification des outils de calcul scientifique utilisés
dans la démonstration de sûreté nucléaire*, 2017), pris comme **référentiel méthodologique
d'inspiration** : ce guide porte formellement sur les outils de la démonstration de sûreté, et un
solveur de décroissance n'entre pas dans son champ d'application. On en retient la distinction
**vérification** (l'outil fait ce que l'on a voulu qu'il fasse : réalisation informatique et
numérique correcte) / **validation** (l'outil représente correctement les phénomènes physiques
dans son domaine de validation).

| Niveau | Contenu | Terme ASN | État |
|---|---|---|---|
| T1 | Unitaires : conversions, parsing, noms de nucléides, validation de la bibliothèque ($\sum b = 1$, filles présentes, absence de cycle) | vérification | ✅ |
| T2 | Solutions analytiques : formules fermées à un et deux corps, Sr-90/Y-90, équilibre séculaire ; oracle mpmath sur le jeu D | vérification (cas de validation analytiques) | ✅ partiel (oracle complet : lot 2) |
| T3 | Ordres de convergence observés vs théoriques | vérification | ✅ contrôle à deux raffinements (protocole complet : lot 2) |
| T4 | Invariants : $N(0) = N_0$ bit à bit, positivité, conservation, semi-groupe $\Phi(t_1+t_2)=\Phi(t_2)\circ\Phi(t_1)$ | vérification | ✅ |
| T5 | Cas dégénérés D1–D3, formule naïve en `[known-limitation]`, raideur (D5) : L-stabilité, concentrations négatives de Crank–Nicolson, divergence d'Euler explicite | vérification | ✅ |
| T6 | Non-régression, sorties figées avec tolérances | vérification | ⬜ |
| T7 | Évaluation croisée avec un outil de référence (`radioactivedecay`) | validation | ⬜ |

## 7. Installation

Prérequis : CMake ≥ 3.21, Ninja, un compilateur C++20 (GCC ≥ 13, Clang ≥ 17, MSVC 2022).
Catch2 est récupéré automatiquement (tag figé) s'il n'est pas installé.

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Sous Windows avec MSVC, lancer ces commandes depuis un *x64 Native Tools Command Prompt* (ou
après `vcvars64.bat`).

## 8. Exemple : vieillir un inventaire

Un inventaire est un fichier `nuclide;valeur`, virgule ou point décimal, graphies usuelles des
noms acceptées (`Cs137`, `Cs-137`, `137Cs`, `Ag108m`). Ici trois nucléides en proportions
normalisées, vieillis de six ans :

```
Cs137;1,10E-02
Pu241;1.21e-1
Am241;3.93e-4
```

```bash
./build/dev/apps/decaysolver age --input spectre.csv --kind fraction --age 6a
```

```
# decaysolver: 0.1.0
# git_sha: 7bb3639a5fc1 (capturé à la configuration CMake)
# compiler: MSVC 19.44.35217.0
# build_type: Debug
# generated_utc: 2026-09-05T14:16:35Z
# data: source: ICRP Publication 107, Nuclear Decay Data for Dosimetric Calculations (2008),
# data: […]  (toute la provenance des données est recopiée)
# input: spectre.csv vieilli de 6a
# input_kind: fraction
# age_s: 189345600 (année julienne = 31 557 600 s)
# daughters: input-only
# total_activity: 0.101535
# alpha_activity: 0.00139493
# beta_gamma_activity: 0.10014
# beta_gamma_over_alpha: 71.7884
nuclide;activity;fraction;primary_mode
Cs-137;9.583412e-03;9.438554e-02;beta-
Pu-241;9.055640e-02;8.918760e-01;beta-
Am-241;1.394930e-03;1.373845e-02;alpha
```

L'Am-241 a été alimenté par la décroissance β⁻ du Pu-241 ; avec `--daughters all`, les filles
hors liste (par exemple Ba-137m sous Cs-137) apparaissent et modifient toutes les fractions.
Options : `--kind bq|fraction`, `--daughters input-only|all`, `--library`, `--output`.
Durées : `6a`, `30j`, `12h`, `90min`, `0s` (année julienne).

## 9. Données nucléaires

Aucune donnée n'est codée en dur. `data/nuclides_icrp107.csv` contient 139 nucléides (les 35 de la
liste standard de déclaration des déchets et la fermeture complète de leurs descendants), extraits
d'**ICRP-107** par un script versionné, avec demi-vies dans leur valeur et unité d'origine.
Source, chemin d'extraction, licence, SHA-256, voies rétablies et limitations (pas d'incertitudes,
divergences entre bibliothèques) : [`data/PROVENANCE.md`](data/PROVENANCE.md).

## 10. Références

- H. Bateman, *Solution of a system of differential equations occurring in the theory of
  radio-active transformations*, Proc. Cambridge Phil. Soc. 15 (1910) 423–427.
- C. Moler, C. Van Loan, *Nineteen Dubious Ways to Compute the Exponential of a Matrix,
  Twenty-Five Years Later*, SIAM Review 45(1) (2003) 3–49. doi:10.1137/S00361445024180
- E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II: Stiff and
  Differential-Algebraic Problems*, Springer, 2e éd., 1996. doi:10.1007/978-3-642-05221-7
- ASN, *Guide n°28 : Qualification des outils de calcul scientifique utilisés dans la
  démonstration de sûreté nucléaire — première barrière*, 2017.

## 11. Citer, licence, contribuer

Citation : voir [`CITATION.cff`](CITATION.cff).

Licence : **CeCILL-C** (voir [`LICENSE`](LICENSE), version française dans
[`LICENSE.fr`](LICENSE.fr)). Choix motivé par la cohérence avec l'écosystème français du logiciel
scientifique (CEA, CNRS, Inria) et par son régime de type LGPL : réutilisable dans un logiciel
propriétaire, modifications de la bibliothèque à publier.

Contributions : toute modification qui change un résultat numérique doit venir avec la mise à
jour des tests de non-régression et une entrée `Numerics` dans le [CHANGELOG](CHANGELOG.md).
