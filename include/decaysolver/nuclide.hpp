#pragma once

/// @file nuclide.hpp
/// Description d'un nucléide : nom canonique, demi-vie, voies de décroissance.
///
/// Style volontairement simple : des structures à champs publics et des fonctions libres.
/// Les invariants (rapports d'embranchement, filles existantes, absence de cycle) sont vérifiés
/// au chargement de la bibliothèque (nuclide_library.hpp), pas ici.

#include <string>
#include <string_view>
#include <vector>

namespace decaysolver {

/// Modes de décroissance retenus. `beta_plus_ec` regroupe β⁺ et capture électronique, qui
/// mènent à la même fille ; ICRP-107 les regroupe de la même façon. `unlisted` désigne une voie
/// dont le jeu de données ne donne pas la fille (rapport = complément à 1 des voies connues) :
/// comme la fission spontanée, les noyaux quittent le système. Elle n'apparaît que dans la
/// bibliothèque complète, pour sept nucléides exotiques (voir data/PROVENANCE.md).
enum class DecayMode {
    alpha,
    beta_minus,
    beta_plus_ec,
    isomeric_transition,
    spontaneous_fission,
    unlisted,
    stable
};

/// Lecture du mode depuis le fichier de données : "alpha", "beta-", "beta+/EC", "IT", "SF",
/// "unlisted", "stable".
/// @throws std::invalid_argument pour tout autre texte.
[[nodiscard]] DecayMode decay_mode_from_string(std::string_view text);

/// Opération inverse, avec le même vocabulaire.
[[nodiscard]] std::string_view to_string(DecayMode mode);

/// Une voie de décroissance.
struct DecayBranch {
    DecayMode mode;
    /// Fille dans sa forme canonique ; vide pour la fission spontanée et les voies non
    /// répertoriées, dont les produits ne sont pas suivis (les noyaux quittent alors le système :
    /// la conservation ne s'applique plus).
    std::string daughter;
    double branching_fraction;
};

struct Nuclide {
    std::string name;                  ///< forme canonique, par exemple "Cs-137" ou "Ba-137m"
    double half_life_s;                ///< +inf pour un nucléide stable
    std::vector<DecayBranch> branches; ///< vide pour un nucléide stable

    [[nodiscard]] bool is_stable() const;
    /// λ = ln 2 / T½ ; 0 pour un nucléide stable.
    [[nodiscard]] double decay_constant_per_s() const;
    /// Mode de la voie de plus grand rapport d'embranchement ; `stable` si aucune voie.
    /// Sert à classer un nucléide en émetteur α ou β/γ pour les agrégats du mode inventaire.
    [[nodiscard]] DecayMode primary_mode() const;
};

/// Normalise une graphie usuelle vers la forme canonique "Xx-NNN", "Xx-NNNm" ou "Xx-NNNn"
/// (second isomère, notation ICRP-107 : "Bi-212n").
///
/// Graphies acceptées : "Cs-137", "Cs137", "cs137", "CS-137", "137Cs", "Ag-108m", "Ag108m".
/// La forme à nombre de masse en tête ("137Cs") n'accepte pas d'indicateur d'isomère, pour ne
/// pas confondre "m" avec la première lettre d'un symbole (Mg, Mn, Mo...).
/// Le symbole n'est pas vérifié contre la table périodique : c'est la recherche dans la
/// bibliothèque qui rejette un nucléide inconnu.
/// @throws std::invalid_argument si la graphie n'est pas reconnue.
[[nodiscard]] std::string canonical_name(std::string_view raw);

} // namespace decaysolver
