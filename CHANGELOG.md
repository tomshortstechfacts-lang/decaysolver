# Changelog

Format : [Keep a Changelog 1.1.0](https://keepachangelog.com/fr/1.1.0/).
Versionnement : [SemVer 2.0.0](https://semver.org/lang/fr/), avec une règle supplémentaire propre
aux codes de calcul : **tout changement de résultat numérique impose au minimum un incrément
MINOR**, et une entrée sous la rubrique `Numerics` décrivant le changement et son ampleur.

## [Unreleased]

### Added
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
- Aucun résultat numérique produit à ce stade.
