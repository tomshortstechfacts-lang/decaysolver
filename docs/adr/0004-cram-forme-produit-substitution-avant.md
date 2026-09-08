# ADR 0004 — CRAM en forme produit, résolutions par substitution avant

**Statut :** accepté, 2026-09-08.

## Contexte

La version 0.1 déclarait CRAM comme limitation : « référence des codes industriels pour les
systèmes très raides, non implémentée ». La mesure qui la justifie existait déjà : l'exponentielle
de matrice Padé + *scaling-and-squaring* (SciPy) perd quatre chiffres sur la chaîne du Ra-226 à
100 ans, parce que ‖At‖ ≈ 10¹³ impose une quarantaine d'élévations au carré, chacune doublant
l'erreur relative (`verification/report/expm_comparison.md`). La solution analytique par
différences divisées n'a pas ce défaut, mais elle énumère les chemins du graphe, ce qui la
condamne sur un graphe large ou cyclique : c'est la voie « exacte » du solveur, pas celle des
codes d'évolution.

## Décision

1. **Approximation rationnelle de Tchebychev de e^x sur ]−∞, 0]** (Pusa 2010, 2011, 2016), ordres
   16 et 48, coefficients publiés (Pusa 2016), relus depuis la transcription publique d'OpenMC et
   **vérifiés** : r_k(x) contre e^x sur ]−10¹², 0] à 1,6·10⁻¹⁵ (k = 16) et 3,2·10⁻¹⁵ (k = 48)
   absolu, deux fois (C++ `test_cram.cpp`, Python `cram_python.py` contre mpmath).
2. **Forme produit (IPF)** : y ← y + 2 Re[α_i (At − θ_i I)⁻¹ y], k/2 fois, puis α₀ y. Elle évite
   la forme en fractions partielles, dont les résidus atteignent 10⁵ et s'annulent entre eux.
3. **Résolution par substitution avant** en arithmétique complexe : dans l'ordre topologique,
   At − θ_i I est triangulaire inférieure. Aucune factorisation, aucune dépendance. Ce choix
   restreint CRAM aux graphes acycliques, ce que la bibliothèque garantit au chargement.
4. **Ordre 48 par défaut** dans l'API : le surcoût (24 substitutions au lieu de 8) est nul à
   l'échelle des problèmes traités, et l'ordre 48 ramène l'erreur *relative* au niveau de
   l'arrondi jusqu'à des populations de 10⁻²⁶ ; l'ordre 16 plafonne en relatif dès 10⁻¹³
   (plancher absolu α₀ ≈ 2·10⁻¹⁶). Le mode inventaire garde **Bateman par défaut** : positivité
   garantie et N(0) restitué bit à bit, ce que CRAM ne promet pas.

## Alternatives écartées

- **Padé + scaling-and-squaring en C++** : reproduirait la perte de quatre chiffres mesurée ; la
  raison d'être de CRAM est précisément d'y échapper.
- **Fractions partielles** (somme des α_i (At − θ_i I)⁻¹ N₀) : même précision théorique, mais
  annulation entre termes de module 10⁵ ; Pusa 2016 recommande la forme produit.
- **Factorisation LU générale** (systèmes avec cycles) : soixante lignes de plus, sans cas d'usage
  dans le périmètre (décroissance spontanée seule). Reportée, notée comme extension.
- **Sous-pas** (`substeps` de Serpent/OpenMC) : utiles quand ‖At‖ dépasse la plage où les
  coefficients sont valables pour les termes sources ; sans terme source, inutiles ici.

## Conséquences

- Trois voies dans le code pour le même résultat : analytique (chemins), CRAM (algèbre linéaire),
  intégrateurs (pas de temps) ; elles s'accordent à 10⁻¹⁴ absolu sur la liste des 35 nucléides
  et sa fermeture (139 nucléides, 10 ans).
- Propriétés déclarées et testées : précision **absolue** (pas relative), positivité **non
  garantie**, N(0) à l'arrondi près. Le test `[cram][T4]` vérifie la borne, pas le signe.
- Nouvelle option `--method bateman|cram16|cram48` du mode inventaire, tracée dans l'en-tête de
  provenance (`# method:`). Aucun résultat existant ne change (défaut inchangé).
- L'extension naturelle est le **terme source** (production sous flux ou par un procédé), où CRAM
  s'étend par matrice augmentée ; c'est là que la forme produit et les sous-pas prennent leur sens.
