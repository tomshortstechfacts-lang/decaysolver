# ADR 0002 — Solution analytique par différences divisées de l'exponentielle

**Statut :** accepté, 2026-09-05.

## Contexte

La formule fermée de Bateman est une somme alternée dont les termes se compensent quand deux
constantes de décroissance sont proches : en `double`, on perd environ log₁₀(λ/|Δλ|) chiffres, et
le résultat vaut exactement 0 pour Δλ/λ ≈ 10⁻¹¹ (cas D3). Le cas exactement dégénéré (Δλ = 0)
divise par zéro. Or les bibliothèques évaluées produisent réellement des constantes égales après
arrondi, et le cahier des charges impose de détecter et traiter ces cas.

## Décision

Écrire la somme de Bateman comme la différence divisée exp[−λ₀, …, −λ_k] de x ↦ e^{xt}, et la
calculer par la table classique, sauf entre nœuds proches ((x_{i+1} − x_i)·t ≤ 0,1) où l'on
utilise la série de Taylor autour du centre du groupe (McCurdy, Ng & Parlett, Math. Comp. 1984),
avec un reste majoré rigoureusement. Les nœuds confondus donnent les dérivées, soit exactement les
termes tᵏe^{−λt}/k! de la littérature : le cas dégénéré n'est pas un cas particulier, c'est la
limite continue de la formule.

La formule fermée naïve est conservée à part, pour documenter la pathologie par un test
`[known-limitation]`, jamais utilisée pour un calcul.

## Alternatives écartées

- **Perturber les constantes égales** (λ → λ(1+ε)) : introduit une erreur de modèle arbitraire et
  ne traite pas les constantes proches.
- **Forme confluente au cas par cas** (multiplicité 2, puis 3…) : combinatoire, et n'apporte rien
  pour le quasi-dégénéré.
- **Exponentielle de matrice en `double`** : perd quatre chiffres sur la chaîne du Ra-226 à 100 ans
  (41 élévations au carré pour ‖At‖ ≈ 10¹³, voir `verification/report/expm_comparison.md`).
  Excellente comme troisième voie de contrôle, mauvaise comme méthode de production sur un
  problème raide. C'est précisément la raison d'être de CRAM dans les codes industriels.
- **Multiprécision partout** : coût et dépendance pour un problème qui a une solution en `double`.

## Conséquences

- Précision : ≤ 10⁻¹² relatif sur D1–D3, 4·10⁻¹³ sur la chaîne du Ra-226 (15 nucléides, 14 ordres
  de grandeur de constantes) contre l'oracle multiprécision.
- Un seuil de regroupement (0,1) et un garde-fou sur la longueur de série (200 termes), tous deux
  dans `BatemanOptions`, documentés avec leur justification.
- Le produit λ₀⋯λ_{k−1} d'un chemin peut descendre vers 10⁻¹⁰⁰ sur les séries d'actinides ;
  dans la plage du `double`, mais limite connue et déclarée.
