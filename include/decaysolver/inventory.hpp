#pragma once

/// @file inventory.hpp
/// Mode inventaire : vieillissement d'une liste de radionucléides avec leurs activités.
///
/// Entrée : un fichier `nuclide;valeur` (activités en Bq, ou proportions normalisées). Traitement :
/// N_i = A_i / λ_i, évolution par la solution analytique de Bateman, retour en activités
/// A_i = λ_i N_i, renormalisation. Sortie : activités, fractions et agrégats α / β-γ.
///
/// Convention des filles, hypothèse à déclarer :
///  - `input_only` : la sortie ne contient que les nucléides de l'entrée. Les filles hors liste
///    sont calculées mais non affichées. C'est la pratique des inventaires industriels, qui ne
///    listent pas les filles à vie courte en équilibre (Y-90 sous Sr-90, Ba-137m sous Cs-137) ;
///  - `all` : tous les nucléides radioactifs produits sont listés.
/// Les deux résultats diffèrent : avec `all`, l'activité totale compte les filles en équilibre, ce
/// qui change toutes les fractions.
///
/// Les nucléides stables ont une activité nulle par construction : ils n'apparaissent jamais dans
/// un inventaire d'activités, en entrée comme en sortie.

#include <decaysolver/nuclide_library.hpp>

#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace decaysolver {

enum class DaughterPolicy { input_only, all };
[[nodiscard]] DaughterPolicy daughter_policy_from_string(std::string_view text);
[[nodiscard]] std::string_view to_string(DaughterPolicy policy);

/// Nature des valeurs d'entrée. Déclarée par l'utilisateur, jamais devinée.
enum class ValueKind { activity_bq, fraction };
[[nodiscard]] ValueKind value_kind_from_string(std::string_view text);
[[nodiscard]] std::string_view to_string(ValueKind kind);

struct InventoryEntry {
    std::string nuclide; ///< forme canonique
    double value;        ///< Bq ou proportion, selon `ValueKind`
};

struct Inventory {
    ValueKind kind;
    std::vector<InventoryEntry> entries;
};

/// Lit un inventaire : lignes `nuclide;valeur` (séparateur `;`, `,` ou tabulation), virgule ou
/// point décimal, notation scientifique, lignes vides et lignes `#` ignorées, une ligne d'en-tête
/// textuelle optionnelle. Les noms sont canonisés, les doublons et les valeurs négatives rejetés.
/// @throws DataError sur toute ligne illisible.
[[nodiscard]] Inventory read_inventory(std::istream& in, ValueKind kind);

/// Durée : nombre suivi d'une unité `s`, `min`, `h`, `j` (jour), `a` (année julienne 365,25 j) ;
/// un `+` initial est accepté ("+6a"). @throws std::invalid_argument sinon.
[[nodiscard]] double parse_duration_s(std::string_view text);

struct AgedEntry {
    std::string nuclide;
    double activity; ///< dans l'unité de l'entrée (Bq ou relative)
    double fraction; ///< activité / total des nucléides listés
    DecayMode primary_mode;
};

struct AgedInventory {
    double age_s;
    DaughterPolicy policy;
    ValueKind kind;
    std::vector<AgedEntry> entries;
    double total;      ///< somme des activités listées
    double alpha;      ///< somme des émetteurs α (mode principal)
    double beta_gamma; ///< somme des autres émetteurs
};

/// Vieillit `inventory` de `age_s` secondes. Un nucléide stable ou absent de la bibliothèque en
/// entrée est une erreur (DataError / std::out_of_range).
[[nodiscard]] AgedInventory age_inventory(const NuclideLibrary& library, const Inventory& inventory,
                                          double age_s, DaughterPolicy policy);

/// Écrit le résultat en CSV avec l'en-tête de provenance (code, données, cas de calcul).
void write_aged_inventory(std::ostream& out, const AgedInventory& aged,
                          const NuclideLibrary& library, std::string_view input_description);

} // namespace decaysolver
