#pragma once

/// @file cram.hpp
/// Méthode CRAM (*Chebyshev Rational Approximation Method*) pour N(t) = e^{At} N₀.
///
/// Principe. Les valeurs propres de A sont les −λ_i ≤ 0 : il suffit d'approcher e^x sur l'axe réel
/// négatif. La meilleure approximation rationnelle (au sens de Tchebychev, erreur uniforme
/// minimale) de e^x sur ]−∞, 0] par une fraction de degré k a une erreur qui décroît comme
/// 9,29^{−k} : 10⁻¹⁵ dès k = 16. Ses pôles θ_i sont deux à deux conjugués ; en forme produit
/// (*incomplete partial fractions*, Pusa 2016) :
///
///     r_k(x) = α₀ · Π_{i=1}^{k/2} ( 1 + 2 Re[ α_i / (x − θ_i) ] ),
///
/// et pour la matrice, chaque facteur est une résolution linéaire complexe :
///
///     y ← y + 2 Re[ α_i (A t − θ_i I)⁻¹ y ],   puis  N(t) = α₀ y.
///
/// Dans l'ordre topologique, A t − θ_i I est triangulaire inférieure : substitution avant en
/// arithmétique complexe, k/2 fois. Le coût ne dépend **pas** de la raideur : pas de mise à
/// l'échelle ni d'élévations au carré (contrairement à Padé + *scaling-and-squaring*, qui perd un
/// chiffre par facteur 2 de ‖At‖, soit quatre chiffres sur la chaîne du Ra-226 à 100 ans).
///
/// Propriétés à connaître avant d'utiliser :
///  - l'erreur de l'approximation est bornée en **absolu** sur tout l'axe négatif (max ≈ 1,6·10⁻¹⁵
///    pour k = 16, 3,2·10⁻¹⁵ pour k = 48, atteint près de x = 0, relativement à ‖N₀‖), avec un
///    plancher égal à α₀ : une population très petite devant la population initiale n'est pas
///    résolue *relativement*. Pour k = 16 (α₀ ≈ 2·10⁻¹⁶) l'erreur relative sur e^x vaut 10⁻⁷ à
///    x = −20 et 10⁻³ à x = −30 ; pour k = 48 (α₀ ≈ 2·10⁻⁴⁷) elle reste ≤ 10⁻¹⁴ jusqu'à x ≈ −60,
///    soit des populations de 10⁻²⁶ (mesures : `verification/report/cram_comparison.md`) ;
///  - conséquence : **la positivité n'est pas garantie** (une population vraie de 10⁻²⁰ peut sortir
///    négative de 10⁻¹⁶), contrairement à la solution analytique et à l'Euler implicite ;
///  - r_k(0) = 1 à l'arrondi près, pas exactement : N(0) n'est pas restitué bit à bit ;
///  - la méthode ne suppose pas le graphe acyclique : c'est la résolution linéaire qui l'exige
///    ici (substitution avant). Un système avec cycles (activation sous flux) demanderait une
///    factorisation LU générale, non fournie.
///
/// Références :
///  - M. Pusa, J. Leppänen, "Computing the Matrix Exponential in Burnup Calculations", Nucl. Sci.
///    Eng. 164 (2010) 140–150.
///  - M. Pusa, "Rational Approximations to the Matrix Exponential in Burnup Calculations", Nucl.
///    Sci. Eng. 169 (2011) 155–167.
///  - M. Pusa, "Higher-Order Chebyshev Rational Approximation Method and Application to Burnup
///    Equations", Nucl. Sci. Eng. 182 (2016) 297–318 (forme IPF, ordres 16 et 48, coefficients).

#include <decaysolver/decay_system.hpp>

#include <complex>
#include <string_view>
#include <vector>

namespace decaysolver {

/// Ordre k de l'approximation rationnelle (k pôles, k/2 stockés).
enum class CramOrder { order_16, order_48 };

[[nodiscard]] std::string_view to_string(CramOrder order);

struct CramCoefficients {
    double alpha0;                           ///< facteur α₀ (≈ limite de r_k en −∞)
    std::vector<std::complex<double>> alpha; ///< α_i, i = 1..k/2
    std::vector<std::complex<double>> theta; ///< θ_i, i = 1..k/2, partie imaginaire > 0
};

/// Table de coefficients (Pusa 2016), voir `src/cram_coefficients.cpp` pour la provenance.
[[nodiscard]] const CramCoefficients& cram_coefficients(CramOrder order);

/// r_k(x) ≈ e^x pour x ≤ 0, évalué en forme produit. Exposé pour les tests : c'est le contrôle
/// le plus direct des coefficients (T2, erreur absolue ≤ 10⁻¹⁴ sur ]−10¹², 0]).
[[nodiscard]] double cram_scalar(double x, CramOrder order);

/// Populations N(t) pour les populations initiales `n0` (indexation de `system.names()`) et le
/// temps `t_s` en secondes, par CRAM d'ordre `order`.
/// @throws std::invalid_argument si `n0.size() != system.size()` ou si `t_s < 0`.
[[nodiscard]] std::vector<double> solve_cram(const DecaySystem& system,
                                             const std::vector<double>& n0, double t_s,
                                             CramOrder order = CramOrder::order_48);

} // namespace decaysolver
