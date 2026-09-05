#include <decaysolver/integrator.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace decaysolver {

namespace {

struct SchemeName {
    Scheme scheme;
    std::string_view text;
    int order;
};

constexpr SchemeName scheme_names[] = {
    {Scheme::euler_explicit, "euler-explicit", 1},
    {Scheme::euler_implicit, "euler-implicit", 1},
    {Scheme::rk4, "rk4", 4},
    {Scheme::crank_nicolson, "crank-nicolson", 2},
};

const SchemeName& entry_of(Scheme scheme) {
    for (const SchemeName& entry : scheme_names) {
        if (entry.scheme == scheme) {
            return entry;
        }
    }
    throw std::invalid_argument("schéma hors énumération");
}

// y + h·r, terme à terme.
std::vector<double> axpy(const std::vector<double>& y, double h, const std::vector<double>& r) {
    std::vector<double> out(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        out[i] = y[i] + h * r[i];
    }
    return out;
}

std::vector<double> euler_explicit_step(const DecaySystem& system, const std::vector<double>& n,
                                        double h) {
    return axpy(n, h, system.rate(n));
}

// Runge–Kutta classique d'ordre 4 : quatre évaluations du second membre.
std::vector<double> rk4_step(const DecaySystem& system, const std::vector<double>& n, double h) {
    const std::vector<double> k1 = system.rate(n);
    const std::vector<double> k2 = system.rate(axpy(n, h / 2.0, k1));
    const std::vector<double> k3 = system.rate(axpy(n, h / 2.0, k2));
    const std::vector<double> k4 = system.rate(axpy(n, h, k3));
    std::vector<double> out(n.size());
    for (std::size_t i = 0; i < n.size(); ++i) {
        out[i] = n[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return out;
}

// Schémas implicites à un pas, forme générale θ :
//   (1 + θ h λ_i) N_i⁺ = (1 − (1−θ) h λ_i) N_i + h Σ_j b_ji λ_j [ (1−θ) N_j + θ N_j⁺ ].
// θ = 1 : Euler implicite ; θ = 1/2 : Crank–Nicolson. Comme A est triangulaire inférieure
// (parents avant filles), N_j⁺ est déjà connu quand on calcule N_i⁺ : substitution avant.
std::vector<double> theta_step(const DecaySystem& system, const std::vector<double>& n, double h,
                               double theta) {
    const std::vector<double>& lambdas = system.decay_constants_per_s();
    std::vector<double> out(n.size());
    for (std::size_t i = 0; i < n.size(); ++i) {
        double production_old = 0.0;
        double production_new = 0.0;
        for (const Production& p : system.productions()[i]) {
            const double rate = p.branching_fraction * lambdas[p.parent];
            production_old += rate * n[p.parent];
            production_new += rate * out[p.parent];
        }
        const double numerator = (1.0 - (1.0 - theta) * h * lambdas[i]) * n[i] +
                                 h * ((1.0 - theta) * production_old + theta * production_new);
        out[i] = numerator / (1.0 + theta * h * lambdas[i]);
    }
    return out;
}

} // namespace

std::string_view to_string(Scheme scheme) {
    return entry_of(scheme).text;
}

Scheme scheme_from_string(std::string_view text) {
    for (const SchemeName& entry : scheme_names) {
        if (entry.text == text) {
            return entry.scheme;
        }
    }
    throw std::invalid_argument("schéma inconnu : '" + std::string(text) + "'");
}

int theoretical_order(Scheme scheme) {
    return entry_of(scheme).order;
}

std::vector<double> step(const DecaySystem& system, const std::vector<double>& populations,
                         double h_s, Scheme scheme) {
    switch (scheme) {
    case Scheme::euler_explicit:
        return euler_explicit_step(system, populations, h_s);
    case Scheme::euler_implicit:
        return theta_step(system, populations, h_s, 1.0);
    case Scheme::rk4:
        return rk4_step(system, populations, h_s);
    case Scheme::crank_nicolson:
        return theta_step(system, populations, h_s, 0.5);
    }
    throw std::invalid_argument("schéma hors énumération");
}

std::vector<double> integrate(const DecaySystem& system, const std::vector<double>& n0,
                              double t_end_s, std::size_t n_steps, Scheme scheme) {
    if (n0.size() != system.size()) {
        throw std::invalid_argument("integrate : n0 de taille " + std::to_string(n0.size()) +
                                    ", attendu " + std::to_string(system.size()));
    }
    if (std::isnan(t_end_s) || t_end_s < 0.0) {
        throw std::invalid_argument("integrate : temps négatif ou NaN");
    }
    if (n_steps == 0) {
        throw std::invalid_argument("integrate : au moins un pas");
    }
    const double h = t_end_s / static_cast<double>(n_steps);
    std::vector<double> n = n0;
    for (std::size_t k = 0; k < n_steps; ++k) {
        n = step(system, n, h, scheme);
    }
    return n;
}

} // namespace decaysolver
