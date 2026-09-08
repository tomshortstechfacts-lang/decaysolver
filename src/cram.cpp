#include <decaysolver/cram.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace decaysolver {

namespace detail {
// Tables définies dans cram_coefficients.cpp (fichier de données, généré depuis la source citée).
extern const CramCoefficients cram16;
extern const CramCoefficients cram48;
} // namespace detail

std::string_view to_string(CramOrder order) {
    switch (order) {
    case CramOrder::order_16:
        return "cram16";
    case CramOrder::order_48:
        return "cram48";
    }
    throw std::invalid_argument("ordre CRAM hors énumération");
}

const CramCoefficients& cram_coefficients(CramOrder order) {
    switch (order) {
    case CramOrder::order_16:
        return detail::cram16;
    case CramOrder::order_48:
        return detail::cram48;
    }
    throw std::invalid_argument("ordre CRAM hors énumération");
}

double cram_scalar(double x, CramOrder order) {
    const CramCoefficients& c = cram_coefficients(order);
    double y = 1.0;
    for (std::size_t i = 0; i < c.alpha.size(); ++i) {
        // Le facteur 2 Re(·) restitue la paire de pôles conjugués (θ_i, conj θ_i).
        y += 2.0 * std::real(c.alpha[i] * y / (x - c.theta[i]));
    }
    return c.alpha0 * y;
}

std::vector<double> solve_cram(const DecaySystem& system, const std::vector<double>& n0, double t_s,
                               CramOrder order) {
    if (n0.size() != system.size()) {
        throw std::invalid_argument("solve_cram : n0 de taille " + std::to_string(n0.size()) +
                                    ", attendu " + std::to_string(system.size()));
    }
    if (std::isnan(t_s) || t_s < 0.0) {
        throw std::invalid_argument("solve_cram : temps négatif ou NaN");
    }
    const CramCoefficients& c = cram_coefficients(order);
    const std::vector<double>& lambdas = system.decay_constants_per_s();
    const std::size_t n = system.size();

    std::vector<double> y = n0;
    std::vector<std::complex<double>> x(n);
    for (std::size_t i = 0; i < c.alpha.size(); ++i) {
        // (A t − θ_i I) x = y, par substitution avant : dans l'ordre topologique, la ligne r ne
        // fait intervenir que les parents de r, déjà calculés.
        //   (−λ_r t − θ_i) x_r + Σ_j b_jr λ_j t x_j = y_r.
        // Le dénominateur n'est jamais nul : Im θ_i ≠ 0.
        for (std::size_t r = 0; r < n; ++r) {
            std::complex<double> rhs = y[r];
            for (const Production& p : system.productions()[r]) {
                rhs -= p.branching_fraction * lambdas[p.parent] * t_s * x[p.parent];
            }
            x[r] = rhs / (-lambdas[r] * t_s - c.theta[i]);
        }
        for (std::size_t r = 0; r < n; ++r) {
            y[r] += 2.0 * std::real(c.alpha[i] * x[r]);
        }
    }
    for (double& value : y) {
        value *= c.alpha0;
    }
    return y;
}

} // namespace decaysolver
