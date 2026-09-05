# Changelog

Format : [Keep a Changelog 1.1.0](https://keepachangelog.com/fr/1.1.0/).
Versionnement : [SemVer 2.0.0](https://semver.org/lang/fr/), avec une règle supplémentaire propre
aux codes de calcul : **tout changement de résultat numérique impose au minimum un incrément
MINOR**, et une entrée sous la rubrique `Numerics` décrivant le changement et son ampleur.

## [Unreleased]

### Added
- Vérification (lot 2) : oracle mpmath indépendant (`verification/scripts/oracle_common.py`,
  formule fermée à 50 chiffres, superposition des chemins) pour la chaîne (1,2,3,0) et la chaîne
  réelle Ra-226 → Pb-206 à 30 j, 1 a, 100 a ; test V1 confrontant la solution analytique à ces
  oracles ; exécutable `decaysolver_convergence` (V2) mesurant les ordres des quatre schémas sur
  15 raffinements, normes L∞ et L2, critère d'acceptation ±0,1 sur ≥ 4 raffinements au-dessus du
  plancher d'arrondi, code de retour exploité par ctest et par un job CI dédié ; figures log-log
  et tableaux générés par `verification/scripts/plot_convergence.py` ; rapport de vérification.
- Mode inventaire (`decaysolver/inventory.hpp`) : lecture `nuclide;valeur` (virgule ou point
  décimal, graphies usuelles), vieillissement par la solution analytique, renormalisation,
  totaux α / β-γ par mode principal, convention des filles `input-only` (défaut) ou `all`,
  sortie CSV avec provenance du code, des données et du cas.
- Ligne de commande `decaysolver age --input … --age 6a [--kind bq|fraction]
  [--daughters input-only|all] [--library …] [--output …]`.
- Système de décroissance (`decaysolver/decay_system.hpp`) : fermeture par filiation, ordre
  topologique (A triangulaire inférieure), termes de production, second membre.
- Solution analytique de Bateman (`decaysolver/bateman.hpp`) sur graphe branché quelconque par
  superposition des chemins ; somme de Bateman calculée comme différence divisée de l'exponentielle,
  cas dégénérés et quasi dégénérés traités par série de Taylor autour des groupes de nœuds proches
  (McCurdy, Ng & Parlett 1984), reste de série majoré rigoureusement. Formule naïve conservée à
  part pour documenter la pathologie d'annulation.
- Intégrateurs à pas fixe (`decaysolver/integrator.hpp`) : Euler explicite, Euler implicite,
  Crank–Nicolson, RK4 ; schémas implicites par substitution avant.
- Oracle mpmath `verification/scripts/oracle_degenerate.py` pour les cas D1–D4.
- Tests T2 (formules fermées, Sr-90/Y-90, équilibre séculaire), T3 (ordres observés des quatre
  schémas à ±0,1), T4 (N(0)=N₀ bit à bit, positivité, conservation, semi-groupe sur la chaîne du
  Ra-226), T5 (D1–D3 corrigés, D2–D3 naïfs en `[known-limitation]`, raideur : L-stabilité
  d'Euler implicite, concentration négative de Crank–Nicolson, divergence d'Euler explicite).
- Bibliothèque de nucléides (`decaysolver/nuclide.hpp`, `decaysolver/nuclide_library.hpp`) :
  noms canoniques (graphies usuelles acceptées), modes de décroissance, voies avec rapports
  d'embranchement, chargement CSV validé (Σ b = 1 à 5e-4, filles présentes, absence de cycle).
- Données `data/nuclides_icrp107.csv` : 139 nucléides ICRP-107 extraits par script versionné,
  provenance et SHA-256 dans `data/PROVENANCE.md`, deux voies β⁻ rétablies (At-219, At-217).
- Socle du projet : bibliothèque `DecaySolver::decaysolver`, exécutable `decaysolver`, CMake ≥ 3.21
  avec presets `dev`, `release`, `asan`, `ubsan`, `coverage`.
- Avertissements stricts (GCC, Clang, MSVC), `-Werror` en CI, refus explicite de `-ffast-math`
  et de la contraction FMA (`cmake/FloatingPoint.cmake`).
- Conversions d'unités (`decaysolver/units.hpp`) : demi-vie ↔ constante de décroissance,
  année julienne 365,25 j, avec tests unitaires (T1).
- En-tête de provenance des sorties (`decaysolver/provenance.hpp`) : version, SHA git,
  compilateur, type de build, horodatage UTC.
- Intégration continue GitHub Actions : GCC 14 et Clang 18 (Ubuntu 24.04), MSVC 2022, Debug et
  Release, vérification `clang-format` bloquante.

### Numerics
- Première version du cœur numérique. Conservation du nombre d'atomes sur données ICRP-107 exacte
  au niveau des arrondis des rapports d'embranchement (jusqu'à 1,3e-6 sur Bi-210), pas de
  l'arithmétique.
