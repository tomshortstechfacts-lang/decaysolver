#include <decaysolver/bateman.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace decaysolver {

namespace {

// exp[x_first, …, x_last] pour un groupe de nœuds proches, par série de Taylor autour du centre c
// du groupe (McCurdy, Ng & Parlett). Avec δ_m = x_m − c et k = nombre de nœuds − 1 :
//
//     exp[x_0..x_k] = e^{ct} · Σ_{n≥0} t^{n+k} / (n+k)! · h_n(δ_0, …, δ_k)
//
// où h_n est le polynôme symétrique complet de degré n (somme de tous les monômes de degré n).
// h_n se calcule nœud par nœud : en ajoutant δ, h_n ← h_n + δ · h_{n−1}.
double taylor_cluster(const std::vector<double>& nodes, std::size_t first, std::size_t last,
                      double t, const BatemanOptions& options) {
    const std::size_t count = last - first + 1;
    const std::size_t k = count - 1;
    double center = 0.0;
    for (std::size_t m = first; m <= last; ++m) {
        center += nodes[m];
    }
    center /= static_cast<double>(count);

    // h_n pour n = 0..max_terms, construits nœud par nœud. h_0 = 1 est poussé avant le
    // redimensionnement : le vecteur n'est jamais vide, ce que GCC vérifie (-Wnull-dereference).
    if (options.max_taylor_terms == 0) {
        throw std::invalid_argument("série de Taylor : au moins un terme");
    }
    std::vector<double> h;
    h.reserve(options.max_taylor_terms + 1);
    h.push_back(1.0);
    h.resize(options.max_taylor_terms + 1, 0.0);
    double spread = 0.0; // S = Σ|δ_m|, pour la majoration du reste
    for (std::size_t m = first; m <= last; ++m) {
        const double delta = nodes[m] - center;
        spread += std::abs(delta);
        for (std::size_t n = 1; n < h.size(); ++n) {
            h[n] += delta * h[n - 1];
        }
    }

    // Arrêt de la série par majoration du reste, pas sur la petitesse d'un terme : pour un groupe
    // symétrique (δ = ±d), h_n est nul pour tout n impair, et un test « terme négligeable »
    // s'arrêterait sur ce zéro avant convergence. Comme |h_n| ≤ Sⁿ, le reste après le terme n est
    // majoré par une série géométrique de raison r = t·S/(n+k+2) dès que r < 1 :
    //     |reste_n| ≤ M_{n+1} / (1 − r),   avec M_n = t^{n+k}/(n+k)! · Sⁿ.
    double coefficient = 1.0; // t^{n+k}/(n+k)!, mis à jour de proche en proche
    for (std::size_t p = 1; p <= k; ++p) {
        coefficient *= t / static_cast<double>(p);
    }
    double majorant = coefficient; // M_n
    double sum = 0.0;
    for (std::size_t n = 0; n < h.size(); ++n) {
        sum += coefficient * h[n];
        const double next_index = static_cast<double>(n + k + 1);
        coefficient *= t / next_index;
        majorant *= t * spread / next_index;
        const double ratio = t * spread / (next_index + 1.0);
        if (ratio < 0.5 &&
            majorant / (1.0 - ratio) <= std::numeric_limits<double>::epsilon() * std::abs(sum)) {
            break;
        }
    }
    return std::exp(center * t) * sum;
}

} // namespace

double divided_difference_exp(std::vector<double> nodes, double t, const BatemanOptions& options) {
    if (nodes.empty()) {
        throw std::invalid_argument("différence divisée : aucun nœud");
    }
    std::sort(nodes.begin(), nodes.end());
    const std::size_t k = nodes.size() - 1;

    // Groupes de nœuds proches : coupure entre i et i+1 si (x_{i+1} − x_i)·t > seuil.
    std::vector<std::size_t> cluster_of(nodes.size(), 0);
    for (std::size_t i = 1; i <= k; ++i) {
        const bool same_cluster = (nodes[i] - nodes[i - 1]) * t <= options.cluster_threshold;
        cluster_of[i] = same_cluster ? cluster_of[i - 1] : cluster_of[i - 1] + 1;
    }

    // Table des différences divisées D[i][j] = exp[x_i..x_j], remplie par longueur croissante.
    // Dans un même groupe : série de Taylor. Entre groupes : récurrence classique
    // exp[x_i..x_j] = (exp[x_{i+1}..x_j] − exp[x_i..x_{j−1}]) / (x_j − x_i).
    std::vector<std::vector<double>> table(nodes.size(), std::vector<double>(nodes.size(), 0.0));
    for (std::size_t length = 0; length <= k; ++length) {
        for (std::size_t i = 0; i + length <= k; ++i) {
            const std::size_t j = i + length;
            if (cluster_of[i] == cluster_of[j]) {
                table[i][j] = taylor_cluster(nodes, i, j, t, options);
            } else {
                table[i][j] = (table[i + 1][j] - table[i][j - 1]) / (nodes[j] - nodes[i]);
            }
        }
    }
    return table[0][k];
}

double divided_difference_exp_naive(const std::vector<double>& nodes, double t) {
    double sum = 0.0;
    for (std::size_t m = 0; m < nodes.size(); ++m) {
        double denominator = 1.0;
        for (std::size_t l = 0; l < nodes.size(); ++l) {
            if (l != m) {
                const double difference = nodes[m] - nodes[l];
                if (!(std::abs(difference) > 0.0)) { // égalité exacte : le cas interdit
                    throw std::domain_error("formule de Bateman naïve : deux constantes égales");
                }
                denominator *= difference;
            }
        }
        sum += std::exp(nodes[m] * t) / denominator;
    }
    return sum;
}

std::vector<double> solve_bateman(const DecaySystem& system, const std::vector<double>& n0,
                                  double t_s, const BatemanOptions& options) {
    if (n0.size() != system.size()) {
        throw std::invalid_argument("solve_bateman : n0 de taille " + std::to_string(n0.size()) +
                                    ", attendu " + std::to_string(system.size()));
    }
    if (std::isnan(t_s) || t_s < 0.0) {
        throw std::invalid_argument("solve_bateman : temps négatif ou NaN");
    }
    const std::vector<double>& lambdas = system.decay_constants_per_s();
    std::vector<double> result(system.size(), 0.0);

    // Un chemin en cours d'exploration : nucléide courant, prochaine fille à visiter, et les
    // facteurs accumulés (produit des rapports d'embranchement, produit des λ des nucléides
    // traversés sauf le dernier).
    struct Step {
        std::size_t nuclide;
        std::size_t next_daughter;
        double branching_product;
        double lambda_product;
    };

    for (std::size_t seed = 0; seed < system.size(); ++seed) {
        if (!(std::abs(n0[seed]) > 0.0)) {
            continue; // population initiale nulle : aucune contribution
        }
        std::vector<Step> path{{seed, 0, 1.0, 1.0}};
        std::vector<double> nodes{-lambdas[seed]};
        // Le préfixe réduit au nucléide de départ contribue lui aussi.
        result[seed] += n0[seed] * divided_difference_exp(nodes, t_s, options);

        while (!path.empty()) {
            Step& step = path.back();
            const std::vector<Production>& daughters = system.daughters()[step.nuclide];
            if (step.next_daughter >= daughters.size()) {
                path.pop_back();
                nodes.pop_back();
                continue;
            }
            const Production& link = daughters[step.next_daughter];
            ++step.next_daughter;
            const Step next{link.parent, 0, step.branching_product * link.branching_fraction,
                            step.lambda_product * lambdas[step.nuclide]};
            nodes.push_back(-lambdas[link.parent]);
            result[link.parent] += n0[seed] * next.branching_product * next.lambda_product *
                                   divided_difference_exp(nodes, t_s, options);
            path.push_back(next);
        }
    }
    return result;
}

} // namespace decaysolver
