#pragma once

/// @file bateman.hpp
/// Solution analytique des équations de Bateman sur un système de décroissance quelconque
/// (graphe acyclique avec embranchements).
///
/// Principe. Pour un chemin de nucléides j₀ → j₁ → … → j_k, de constantes λ₀…λ_k et de rapports
/// d'embranchement b₁…b_k, la population apportée au dernier nucléide par N₀ atomes de j₀ vaut
/// (Bateman, 1910) :
///
///     N_k(t) = N₀ · (b₁⋯b_k) · (λ₀⋯λ_{k−1}) · Σ_m exp(−λ_m t) / Π_{l≠m} (λ_l − λ_m).
///
/// La somme est exactement la **différence divisée** de la fonction x ↦ exp(x t) aux nœuds
/// x_m = −λ_m : exp[x₀,…,x_k]. Cette lecture apporte deux choses :
///  - le cas dégénéré (λ_p = λ_q) est la limite continue : les nœuds confondus donnent des
///    dérivées, et exp[x,…,x] (k+1 fois) = tᵏ e^{xt} / k!, ce qui redonne les termes en
///    tᵏ e^{−λt}/k! de la littérature ;
///  - le cas quasi dégénéré (λ_p ≈ λ_q), où la somme de Bateman perd ses chiffres par annulation
///    catastrophique, se traite par développement de Taylor autour du centre du groupe de nœuds
///    proches (McCurdy, Ng & Parlett, 1980), qui converge d'autant plus vite que les nœuds sont
///    proches.
///
/// Le système général se résout par superposition : on énumère tous les chemins issus de chaque
/// nucléide de population initiale non nulle, et chaque préfixe de chemin contribue au nucléide
/// où il s'arrête.
///
/// Références :
///  - H. Bateman, Proc. Cambridge Phil. Soc. 15 (1910) 423–427.
///  - A.C. McCurdy, K.C. Ng, B.N. Parlett, "Accurate computation of divided differences of the
///    exponential function", Math. Comp. 43 (1984) 501–528. doi:10.1090/S0025-5718-1984-0758202-4
///
/// Complexité : nombre de chemins × k² pour la table des différences divisées ; négligeable pour
/// les chaînes de décroissance réelles (k ≤ 15, quelques centaines de chemins au plus).

#include <decaysolver/decay_system.hpp>

#include <cstddef>
#include <vector>

namespace decaysolver {

struct BatemanOptions {
    /// Deux nœuds consécutifs (triés) x_i < x_{i+1} sont regroupés si (x_{i+1} − x_i)·t ≤ seuil.
    /// Justification : la différence (e^{x_{i+1}t} − e^{x_i t})/(x_{i+1} − x_i) perd environ
    /// log10(1/((x_{i+1} − x_i) t)) chiffres par annulation quand cette quantité est < 1.
    /// Avec 0,1 on perd au plus un chiffre hors des groupes, et dans un groupe la série de Taylor
    /// a des termes en (0,1)ⁿ/n! : précision machine en une quinzaine de termes.
    double cluster_threshold = 0.1;
    /// Garde-fou sur la longueur de la série de Taylor.
    std::size_t max_taylor_terms = 200;
};

/// Différence divisée exp[x₀,…,x_k] de x ↦ exp(x·t) aux nœuds donnés (ordre quelconque, nœuds
/// répétés admis). Fonction de base, exposée pour les tests de vérification.
[[nodiscard]] double divided_difference_exp(std::vector<double> nodes, double t,
                                            const BatemanOptions& options = {});

/// Même quantité par la formule fermée de Bateman, sans traitement de la dégénérescence.
/// **Sert uniquement à documenter la pathologie d'annulation** (tests `[known-limitation]`,
/// rapport de vérification) ; ne pas utiliser pour un calcul.
/// @throws std::domain_error si deux nœuds sont exactement égaux (division par zéro).
[[nodiscard]] double divided_difference_exp_naive(const std::vector<double>& nodes, double t);

/// Populations N(t) de tous les nucléides du système, pour les populations initiales `n0`
/// (même indexation que `system.names()`) et le temps `t_s` en secondes.
/// @throws std::invalid_argument si `n0.size() != system.size()` ou si `t_s < 0`.
[[nodiscard]] std::vector<double> solve_bateman(const DecaySystem& system,
                                                const std::vector<double>& n0, double t_s,
                                                const BatemanOptions& options = {});

} // namespace decaysolver
