# ADR 0003 — Cœur numérique en C++ volontairement simple, vérification en Python

**Statut :** accepté, 2026-09-05.

## Contexte

Le solveur doit être un outil de calcul scientifique lisible et défendable ligne par ligne, et sa
vérification doit reposer sur des références indépendantes du code vérifié. Le projet est
développé avec une assistance IA ; la maîtrise du résultat, pas son origine, est le critère.

## Décision

1. **C++20, style simple** : structures à champs publics, fonctions libres, `std::vector`,
   `std::string`, exceptions typées, boucles explicites. Aucun template dans l'API publique, pas
   de métaprogrammation, aucune bibliothèque à templates (Eigen écarté). Chaque idiome moderne
   (`constexpr`, `[[nodiscard]]`, `std::span`, `std::numbers`) est commenté une fois.
2. **Toute la vérification en Python** : oracle mpmath (formule fermée, 50 chiffres), figures,
   tableaux, exponentielle de matrice SciPy, évaluation croisée `radioactivedecay`. Les scripts
   sont versionnés ; aucune valeur de référence n'est committée sans le script qui la produit.
3. **Un seul point de contact** entre les deux mondes : des fichiers CSV avec en-tête de
   provenance, lus par les tests C++ et par les scripts Python.

## Alternatives écartées

- **Tout en C++** (oracle en `boost::multiprecision`) : plus de C++ à défendre, moins de recul :
  l'oracle partagerait le langage et une partie du code du solveur.
- **Tout en Python** : ne répond pas à l'exigence du domaine (codes de calcul en C++/Fortran) et
  n'exerce pas l'arithmétique `double` avec les options de compilation réelles.
- **C++ « moderne expert »** (concepts, ranges, expression templates) : plus court, moins lisible
  pour un physicien, et impossible à défendre honnêtement par l'auteur dans le délai.

## Conséquences

- Le C++ ressemble à du code de physique, pas à du code de bibliothèque : c'est voulu.
- La vérification est reproductible avec un interpréteur Python et deux paquets (mpmath, SciPy).
- L'assistance IA à l'écriture du code est déclarée ; la revue humaine porte sur chaque tolérance,
  chaque hypothèse et chaque limitation écrite.
