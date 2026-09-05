# ADR 0001 — Catch2 v3 comme cadre de test

**Statut :** accepté, 2026-09-05.

## Contexte

Un solveur numérique se teste avec des tolérances, jamais avec `==`. Le cadre de test doit rendre
la tolérance explicite dans l'assertion, la justifier lisiblement, et produire un message d'échec
qui montre les deux valeurs et l'écart.

## Décision

Catch2 v3 (tag figé v3.8.1 via FetchContent, ou installation système via `find_package`), avec
les matchers flottants `WithinRel`, `WithinAbs`, `WithinULP`, composables (`||`) pour traiter les
zéros numériques. Découverte des tests par `catch_discover_tests`, un cas ctest par `TEST_CASE`.

## Alternatives écartées

- **GoogleTest** : `EXPECT_NEAR` n'a qu'une tolérance absolue ; les tolérances relatives, seules
  pertinentes ici, demandent du code maison.
- **doctest** : plus léger, mais sans matchers flottants intégrés.
- **Assertions maison** : réinventer ce que Catch2 fait bien, sans le rapport d'échec.

## Conséquences

- Chaque tolérance est écrite dans le test avec sa justification en commentaire (ulp, arrondis
  de la source de données, plancher d'arrondi).
- Les noms de `TEST_CASE` restent en ASCII : ctest les transmet en filtre et les consoles Windows
  mutilent les accents. Les `SECTION` peuvent être accentuées.
- Catch2 est traité comme en-tête système : nos avertissements stricts (`-Werror` en CI) portent
  sur notre code, pas sur la dépendance.
