#pragma once

/// @file integrator.hpp
/// Intégration temporelle de dN/dt = A·N à pas fixe, par quatre schémas classiques.
///
/// | Schéma           | Ordre | Stabilité                         | Positivité |
/// |------------------|-------|-----------------------------------|------------|
/// | Euler explicite  | 1     | h·λ_max ≤ 2                       | h·λ_max ≤ 1 |
/// | RK4 explicite    | 4     | h·λ_max ≤ 2,785                   | non garantie |
/// | Euler implicite  | 1     | inconditionnelle, A- et L-stable  | garantie (A de Metzler) |
/// | Crank–Nicolson   | 2     | inconditionnelle, A-stable, **non L-stable** | non : R(z) → −1
/// quand z → −∞, les modes raides oscillent et N peut devenir négatif |
///
/// A est triangulaire inférieure dans l'ordre topologique du système : les schémas implicites se
/// résolvent par substitution avant, nucléide par nucléide, sans système linéaire général.
///
/// Référence : E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II*, Springer,
/// 1996, chap. IV.3 (fonctions de stabilité) et IV.4 (L-stabilité).

#include <decaysolver/decay_system.hpp>

#include <cstddef>
#include <string_view>
#include <vector>

namespace decaysolver {

enum class Scheme { euler_explicit, euler_implicit, rk4, crank_nicolson };

[[nodiscard]] std::string_view to_string(Scheme scheme);
/// Noms acceptés : "euler-explicit", "euler-implicit", "rk4", "crank-nicolson".
/// @throws std::invalid_argument sinon.
[[nodiscard]] Scheme scheme_from_string(std::string_view text);

/// Ordre de convergence théorique du schéma (1, 1, 4, 2).
[[nodiscard]] int theoretical_order(Scheme scheme);

/// Un pas de longueur `h_s` depuis `populations`.
[[nodiscard]] std::vector<double>
step(const DecaySystem& system, const std::vector<double>& populations, double h_s, Scheme scheme);

/// Intègre de 0 à `t_end_s` en `n_steps` pas égaux (h = t_end_s / n_steps : le nombre de pas est
/// un entier par construction, ce qui évite tout arrondi sur T/h).
/// @throws std::invalid_argument si n0 a la mauvaise taille, si t_end_s < 0 ou si n_steps == 0.
[[nodiscard]] std::vector<double> integrate(const DecaySystem& system,
                                            const std::vector<double>& n0, double t_end_s,
                                            std::size_t n_steps, Scheme scheme);

} // namespace decaysolver
