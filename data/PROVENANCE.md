# Provenance des données nucléaires

Aucune donnée nucléaire n'est codée en dur dans decaysolver : tout vient de ce dossier, et tout
fichier de données y est accompagné du script qui l'a produit.

## `nuclides_icrp107.csv`

| | |
|---|---|
| **Source** | ICRP Publication 107, *Nuclear Decay Data for Dosimetric Calculations*, Ann. ICRP 38(3), 2008. Copyright © 2008 A. Endo et K.F. Eckerman. doi:10.1016/j.icrp.2008.10.004 |
| **Chemin d'extraction** | paquet Python [`radioactivedecay`](https://github.com/radioactivedecay/radioactivedecay) 0.6.1 (MIT), jeu `icrp107_ame2020_nubase2020`, qui redistribue ICRP-107 |
| **Script** | [`scripts/extract_icrp107.py`](scripts/extract_icrp107.py), déterministe (tri par nom) |
| **Date d'extraction** | 2026-09-05 (UTC) |
| **Licence des données** | [`LICENSE.ICRP-07`](LICENSE.ICRP-07) : usage, copie et redistribution autorisés à des fins éducatives, de recherche et non lucratives, à condition de joindre la notice. Ce dépôt relève de cet usage. |
| **Contenu** | 35 nucléides de base (liste standard de déclaration des déchets) + fermeture complète de leurs descendants : 139 nucléides, 171 voies, 35 stables |
| **SHA-256** | `14dc5f9c41a371a9fc851c2fe7e0b3902ba7f307752936f877558d79fc093254` (fichier [`nuclides_icrp107.sha256`](nuclides_icrp107.sha256)) |
| **Corrections** | 2 voies β⁻ rétablies (At-219, At-217), voir ci-dessous |

Vérification de l'empreinte :

```bash
sha256sum -c data/nuclides_icrp107.sha256
```

### Format

Une ligne par voie de décroissance, séparateur `;`, en-tête `#` recopié dans la provenance des
sorties du solveur :

```
nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction
Cs-137;30.1671;y;beta-;Ba-137m;0.94399
Cs-137;30.1671;y;beta-;Ba-137;0.056005
Ba-137;stable;;stable;;
```

- `half_life_value` et `half_life_unit` sont la **valeur et l'unité d'origine** (s, m, h, d, y).
  La conversion en secondes est faite par decaysolver avec **1 y = 365,25 j** (année julienne).
  `radioactivedecay` utilise 365,2422 j : écart relatif 2·10⁻⁵ sur les demi-vies exprimées en
  années, à retrancher avant d'interpréter une évaluation croisée.
- `mode` ∈ {`alpha`, `beta-`, `beta+/EC`, `IT`, `SF`, `stable`}. ICRP-107 regroupe β⁺ et
  capture électronique.
- `SF` (fission spontanée) : voie conservée avec `daughter` vide ; les produits de fission ne
  sont pas suivis. Les rapports concernés sont ≤ 1,4·10⁻⁶ sur ce jeu.

### Voies rétablies par le script

Le jeu redistribué par `radioactivedecay` ne contient, pour **At-219** et **At-217**, que la voie α
(0,97 et 0,99988), alors qu'ICRP-107 leur attribue aussi une voie β⁻. Sans correction, la somme
des rapports vaut 0,97 et 0,99988, et le chargeur de decaysolver refuse le fichier pour At-219.
Le script rétablit ces deux voies :

| Nucléide | Voie rétablie | Rapport | Fille | Contrôle indépendant (ENSDF, API Livechart AIEA, 2026-09-05) |
|---|---|---|---|---|
| At-219 | β⁻ | 1 − 0,97 = 0,03 | Rn-219 | β⁻ 6,4 % (évaluation actuelle) |
| At-217 | β⁻ | 1 − 0,99988 = 1,2·10⁻⁴ | Rn-217 | β⁻ 0,007 % |

Le rapport retenu est le complément à 1 de la valeur ICRP-107, pour rester cohérent avec la
bibliothèque déclarée ; la valeur ENSDF, différente, illustre la dispersion des évaluations sur ces
nucléides très mineurs de la chaîne du Np-237. Rn-217 (0,54 ms, α → Po-213) entre dans la fermeture
du fait de cette correction. Ces deux voies pèsent moins de 10⁻⁶ sur tout inventaire réaliste.

**Leçon retenue :** une voie manquante de 1,2·10⁻⁴ est invisible pour un contrôle de somme à
5·10⁻⁴ ; seule une comparaison à une seconde source la révèle. Le script écrit dans l'en-tête du
fichier tout écart de somme supérieur à 10⁻⁵ (aujourd'hui : Fr-223, +6·10⁻⁵, arrondi ICRP-107).

### Limitations connues de ce jeu

1. **Pas d'incertitudes.** ICRP-107 ne publie pas les incertitudes sur les demi-vies. Or, pour
   plusieurs nucléides, l'incertitude des données évaluées dépasse largement l'erreur numérique du
   solveur : c'est elle qui borne la précision réelle d'un inventaire vieilli. Prévu : colonne
   d'incertitude alimentée depuis DDEP/LNHB pour les 35 nucléides de base.
2. **Les évaluations divergent entre bibliothèques.** Exemple : Cs-137 vaut 30,1671 a dans
   ICRP-107, 30,05 a dans DDEP (2006), 30,08 a dans ENSDF (2020) ; Sn-121m vaut 43,9 a ici contre
   55 a dans les évaluations antérieures à 2000. Un écart de quelques 10⁻³ entre deux outils de
   vieillissement est le plus souvent un écart de bibliothèque, pas de numérique.
3. **Le mode inventaire suit par défaut les filles qui figurent dans la liste d'entrée**
   (convention `input-only`) ; le jeu contient pourtant toutes les filles, jusqu'aux stables, pour
   que le mode `all` soit exact.
