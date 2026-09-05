// V2 — mesure des ordres de convergence des intégrateurs (T3, CDC §4.3).
//
// Problème non raide λ = (1, 2, 3, 0) s⁻¹, N(0) = (1, 0, 0, 0), T = 4 s. Référence : oracle mpmath
// (oracle_lambda123.csv, 50 chiffres). Pas h_k = T / 2^k, k = 2..16 ; erreurs à T en normes L∞ et
// L2 sur le vecteur des populations ; ordre observé p_obs = log2(E_k / E_{k+1}).
//
// Critère d'acceptation : pour chaque schéma, au moins quatre p_obs consécutifs dans
// [p − 0,1 ; p + 0,1], parmi les raffinements dont l'erreur reste au-dessus du plancher d'arrondi
// (E_{k+1} ≥ 1e-12). Au-delà, l'arrondi domine la troncature et p_obs s'effondre : c'est attendu,
// mesuré, et rapporté, pas caché.
//
// Sortie : CSV sur la sortie standard (ou dans le fichier donné en argument), lignes de synthèse
// préfixées par `#`. Code de retour 0 si tous les schémas satisfont le critère, 1 sinon.

#include <decaysolver/integrator.hpp>

#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "chain_fixture.hpp"

namespace {

constexpr double t_end_s = 4.0;
constexpr int k_min = 2;
constexpr int k_max = 16;
constexpr double roundoff_floor = 1e-12;
constexpr double order_tolerance = 0.1;
constexpr int required_consecutive = 4;

std::map<std::string, double> read_oracle(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("oracle introuvable : " + path);
    }
    std::map<std::string, double> values;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line.front() == '#' || line.rfind("case;", 0) == 0) {
            continue;
        }
        std::vector<std::string> fields;
        std::stringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ';')) {
            fields.push_back(field);
        }
        if (fields.size() != 4 || fields[0] != "lambda123") {
            continue;
        }
        double value = 0.0;
        const auto result =
            std::from_chars(fields[3].data(), fields[3].data() + fields[3].size(), value);
        if (result.ec != std::errc{}) {
            throw std::runtime_error("valeur d'oracle illisible : " + fields[3]);
        }
        values[fields[2]] = value;
    }
    if (values.size() != 4) {
        throw std::runtime_error("oracle incomplet : 4 nucléides attendus");
    }
    return values;
}

struct Row {
    int k;
    std::size_t n_steps;
    double h;
    double error_inf;
    double error_2;
};

std::vector<Row> measure(const decaysolver::DecaySystem& system, const std::vector<double>& n0,
                         const std::vector<double>& reference, decaysolver::Scheme scheme) {
    std::vector<Row> rows;
    for (int k = k_min; k <= k_max; ++k) {
        const std::size_t n_steps = std::size_t{1} << k;
        const std::vector<double> n = decaysolver::integrate(system, n0, t_end_s, n_steps, scheme);
        double error_inf = 0.0;
        double error_2 = 0.0;
        for (std::size_t i = 0; i < n.size(); ++i) {
            const double difference = std::abs(n[i] - reference[i]);
            error_inf = std::max(error_inf, difference);
            error_2 += difference * difference;
        }
        rows.push_back(
            {k, n_steps, t_end_s / static_cast<double>(n_steps), error_inf, std::sqrt(error_2)});
    }
    return rows;
}

// Plus longue suite de p_obs consécutifs acceptables au-dessus du plancher ; rend sa longueur et
// le k de départ.
std::pair<int, int> longest_valid_run(const std::vector<Row>& rows, int theoretical,
                                      bool use_inf_norm) {
    int best = 0;
    int best_start = 0;
    int current = 0;
    int current_start = 0;
    for (std::size_t i = 0; i + 1 < rows.size(); ++i) {
        const double coarse = use_inf_norm ? rows[i].error_inf : rows[i].error_2;
        const double fine = use_inf_norm ? rows[i + 1].error_inf : rows[i + 1].error_2;
        const double observed = std::log2(coarse / fine);
        const bool valid = fine >= roundoff_floor &&
                           std::abs(observed - static_cast<double>(theoretical)) <= order_tolerance;
        if (valid) {
            if (current == 0) {
                current_start = rows[i].k;
            }
            ++current;
            if (current > best) {
                best = current;
                best_start = current_start;
            }
        } else {
            current = 0;
        }
    }
    return {best, best_start};
}

} // namespace

int main(int argc, char** argv) {
    try {
        const std::string oracle_path = argc > 1 ? argv[1] : DECAYSOLVER_CONVERGENCE_ORACLE;
        std::ofstream file;
        if (argc > 2) {
            file.open(argv[2]);
            if (!file) {
                throw std::runtime_error(std::string("impossible d'écrire : ") + argv[2]);
            }
        }
        std::ostream& out = argc > 2 ? file : std::cout;

        const std::map<std::string, double> oracle = read_oracle(oracle_path);
        const ChainFixture chain = make_linear_chain({1.0, 2.0, 3.0, 0.0});
        std::vector<double> reference;
        for (const std::string& name : chain.system.names()) {
            reference.push_back(oracle.at(name));
        }
        const std::vector<double> n0{1.0, 0.0, 0.0, 0.0};

        out << "# decaysolver_convergence : lambda = (1,2,3,0) s^-1, T = 4 s, N0 = (1,0,0,0)\n";
        out << "# oracle : " << oracle_path << "\n";
        out << "# plancher d'arrondi : " << roundoff_floor << " ; tolerance sur l'ordre : +-"
            << order_tolerance << " ; raffinements consecutifs requis : " << required_consecutive
            << "\n";
        out << "scheme;order;k;n_steps;h;E_inf;E_2;p_inf;p_2\n";
        out << std::scientific << std::setprecision(6);

        bool all_ok = true;
        for (const decaysolver::Scheme scheme :
             {decaysolver::Scheme::euler_explicit, decaysolver::Scheme::euler_implicit,
              decaysolver::Scheme::crank_nicolson, decaysolver::Scheme::rk4}) {
            const int order = decaysolver::theoretical_order(scheme);
            const std::vector<Row> rows = measure(chain.system, n0, reference, scheme);
            for (std::size_t i = 0; i < rows.size(); ++i) {
                out << decaysolver::to_string(scheme) << ';' << order << ';' << rows[i].k << ';'
                    << rows[i].n_steps << ';' << rows[i].h << ';' << rows[i].error_inf << ';'
                    << rows[i].error_2 << ';';
                if (i + 1 < rows.size()) {
                    out << std::fixed << std::setprecision(3)
                        << std::log2(rows[i].error_inf / rows[i + 1].error_inf) << ';'
                        << std::log2(rows[i].error_2 / rows[i + 1].error_2) << std::scientific
                        << std::setprecision(6);
                } else {
                    out << ';';
                }
                out << '\n';
            }
            const auto [run_inf, start_inf] = longest_valid_run(rows, order, true);
            const auto [run_2, start_2] = longest_valid_run(rows, order, false);
            const bool ok = run_inf >= required_consecutive && run_2 >= required_consecutive;
            all_ok = all_ok && ok;
            out << "# " << decaysolver::to_string(scheme) << " : ordre theorique " << order
                << ", suite valide L_inf = " << run_inf << " (depuis k=" << start_inf
                << "), L2 = " << run_2 << " (depuis k=" << start_2 << ") -> "
                << (ok ? "OK" : "ECHEC") << '\n';
        }
        out << "# resultat global : " << (all_ok ? "OK" : "ECHEC") << '\n';
        return all_ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "erreur : " << error.what() << '\n';
        return 2;
    }
}
