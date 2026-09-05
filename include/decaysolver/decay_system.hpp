#pragma once

/// @file decay_system.hpp
/// Système d'équations dN/dt = A·N construit depuis la bibliothèque pour un ensemble de
/// nucléides de départ et tous leurs descendants.
///
/// Les nucléides sont numérotés dans un **ordre topologique** (tout parent avant ses filles).
/// Dans cet ordre, A est triangulaire inférieure : A_ii = −λ_i, A_ij = b_ji λ_j si j est parent
/// de i, 0 sinon. Deux conséquences exploitées par les solveurs :
///  - les schémas implicites se résolvent par substitution avant, sans système linéaire général ;
///  - les valeurs propres sont exactement les −λ_i.

#include <decaysolver/nuclide_library.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace decaysolver {

/// Un terme de production : le nucléide `parent` alimente le nucléide courant avec le taux
/// `branching_fraction · λ_parent`.
struct Production {
    std::size_t parent;
    double branching_fraction;
};

class DecaySystem {
public:
    /// Fermeture des `seeds` par filiation, dans la bibliothèque `library`.
    /// @throws std::out_of_range si un nucléide de départ est absent de la bibliothèque.
    [[nodiscard]] static DecaySystem build(const NuclideLibrary& library,
                                           const std::vector<std::string>& seeds);

    [[nodiscard]] std::size_t size() const { return names_.size(); }
    /// Noms canoniques dans l'ordre topologique.
    [[nodiscard]] const std::vector<std::string>& names() const { return names_; }
    /// @throws std::out_of_range si le nucléide n'est pas dans le système.
    [[nodiscard]] std::size_t index_of(std::string_view name) const;
    [[nodiscard]] bool contains(std::string_view name) const;

    /// λ_i en s⁻¹ (0 pour un nucléide stable).
    [[nodiscard]] const std::vector<double>& decay_constants_per_s() const { return lambdas_; }
    /// Pour chaque nucléide i, ses termes de production (parents j avec b_ji).
    [[nodiscard]] const std::vector<std::vector<Production>>& productions() const {
        return productions_;
    }
    /// Pour chaque nucléide j, les indices de ses filles avec leur rapport d'embranchement.
    [[nodiscard]] const std::vector<std::vector<Production>>& daughters() const {
        return daughters_;
    }
    /// Le nucléide `index` de la bibliothèque d'origine (demi-vie, mode principal…).
    [[nodiscard]] const Nuclide& nuclide(std::size_t index) const { return *nuclides_[index]; }

    /// Second membre : (A·N)_i = −λ_i N_i + Σ_j b_ji λ_j N_j.
    /// @throws std::invalid_argument si `populations.size() != size()`.
    [[nodiscard]] std::vector<double> rate(const std::vector<double>& populations) const;

private:
    std::vector<std::string> names_;
    std::vector<const Nuclide*> nuclides_;
    std::vector<double> lambdas_;
    std::vector<std::vector<Production>> productions_;
    std::vector<std::vector<Production>> daughters_;
};

} // namespace decaysolver
