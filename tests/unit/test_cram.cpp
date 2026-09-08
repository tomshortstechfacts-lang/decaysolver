// CRAM : contrôle des coefficients (T2, scalaire), accord avec l'oracle multiprécision (V1) et
// avec la solution analytique (T2), invariants (T4) et propriétés déclarées (positivité non
// garantie, précision absolue et non relative).
//
// Tolérances. L'approximation rationnelle a une erreur *absolue* uniforme sur ]−∞, 0] : ≈ 10⁻¹⁵
// pour k = 16, quelques 10⁻¹⁵ pour k = 48 (mesurées par verification/scripts/cram_python.py). Les
// comparaisons matricielles portent donc sur l'écart absolu rapporté à ‖N₀‖ (ici 1), à 10⁻¹⁴ ;
// l'écart relatif n'est exigé que pour k = 48 et des populations ≥ 10⁻³⁰, à 10⁻¹³.

#include <decaysolver/bateman.hpp>
#include <decaysolver/cram.hpp>
#include <decaysolver/decay_system.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "chain_fixture.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using decaysolver::cram_scalar;
using decaysolver::CramOrder;
using decaysolver::DecaySystem;
using decaysolver::NuclideLibrary;
using decaysolver::solve_bateman;
using decaysolver::solve_cram;
namespace units = decaysolver::units;

namespace {

const NuclideLibrary& icrp_library() {
    static const NuclideLibrary library =
        NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    return library;
}

struct OracleRow {
    std::string case_name;
    double t_s;
    std::string nuclide;
    double value;
};

// Même lecteur que test_oracle_cases.cpp (fichiers `case;t_s;nuclide;value`).
std::vector<OracleRow> read_oracle(const std::filesystem::path& path) {
    std::ifstream file(path);
    REQUIRE(file.good());
    std::vector<OracleRow> rows;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line.front() == '#' || line.rfind("case;", 0) == 0) {
            continue;
        }
        std::stringstream stream(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(stream, field, ';')) {
            fields.push_back(field);
        }
        REQUIRE(fields.size() == 4);
        OracleRow row{fields[0], 0.0, fields[2], 0.0};
        std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), row.t_s);
        std::from_chars(fields[3].data(), fields[3].data() + fields[3].size(), row.value);
        rows.push_back(row);
    }
    REQUIRE_FALSE(rows.empty());
    return rows;
}

const std::filesystem::path verification_dir =
    std::filesystem::path(DECAYSOLVER_DATA_DIR) / ".." / "verification";

const CramOrder both_orders[] = {CramOrder::order_16, CramOrder::order_48};

} // namespace

TEST_CASE("cram: coefficients", "[cram][T1]") {
    const auto& c16 = decaysolver::cram_coefficients(CramOrder::order_16);
    const auto& c48 = decaysolver::cram_coefficients(CramOrder::order_48);
    REQUIRE(c16.alpha.size() == 8);
    REQUIRE(c16.theta.size() == 8);
    REQUIRE(c48.alpha.size() == 24);
    REQUIRE(c48.theta.size() == 24);
    // Pôles dans le demi-plan supérieur : le dénominateur x − θ_i ne s'annule jamais sur l'axe
    // réel.
    for (const auto& theta : c16.theta) {
        REQUIRE(theta.imag() > 0.0);
    }
    for (const auto& theta : c48.theta) {
        REQUIRE(theta.imag() > 0.0);
    }
    REQUIRE(decaysolver::to_string(CramOrder::order_16) == "cram16");
    REQUIRE(decaysolver::to_string(CramOrder::order_48) == "cram48");
}

TEST_CASE("cram: r_k(x) reproduit exp(x) sur l'axe negatif", "[cram][T2]") {
    SECTION("erreur absolue uniforme, x de -1e-8 a -1e12") {
        for (const CramOrder order : both_orders) {
            double worst = 0.0;
            for (int k = -32; k <= 48; ++k) { // x = −10^{k/4}
                const double x = -std::pow(10.0, k / 4.0);
                worst = std::max(worst, std::abs(cram_scalar(x, order) - std::exp(x)));
            }
            INFO(decaysolver::to_string(order) << " : erreur absolue max " << worst);
            REQUIRE(worst <= 1e-14);
        }
    }
    SECTION("r_k(0) = 1 a l'arrondi pres, pas exactement") {
        REQUIRE_THAT(cram_scalar(0.0, CramOrder::order_16), WithinAbs(1.0, 1e-15));
        REQUIRE_THAT(cram_scalar(0.0, CramOrder::order_48), WithinAbs(1.0, 1e-15));
    }
    SECTION("erreur relative : bonne pour k = 48 jusqu'a x = -60, pas pour k = 16") {
        for (int k = 1; k <= 60; ++k) {
            const double x = -static_cast<double>(k);
            REQUIRE_THAT(cram_scalar(x, CramOrder::order_48), WithinRel(std::exp(x), 1e-13));
        }
        // Limitation déclarée : à x = −30, e^x ≈ 10⁻¹³ et l'erreur absolue de l'ordre 16 (≈ 10⁻¹⁶)
        // pèse 10⁻³ en relatif. Ce test documente le fait, il ne le corrige pas.
        const double x = -30.0;
        const double relative = std::abs(cram_scalar(x, CramOrder::order_16) / std::exp(x) - 1.0);
        REQUIRE(relative > 1e-6);
        REQUIRE(relative < 1e-2);
    }
    SECTION("limite en -infini : alpha0, pas zero") {
        REQUIRE_THAT(cram_scalar(-1e300, CramOrder::order_16),
                     WithinAbs(2.124853710495224e-16, 1e-30));
        REQUIRE_THAT(cram_scalar(-1e300, CramOrder::order_48),
                     WithinAbs(2.258038182743983e-47, 1e-60));
    }
}

TEST_CASE("cram: chaine (1,2,3,0) a T = 4 s vs oracle", "[cram][T2][verification][V1]") {
    const std::vector<OracleRow> oracle =
        read_oracle(verification_dir / "V2_order_of_accuracy" / "oracle_lambda123.csv");
    const ChainFixture chain = make_linear_chain({1.0, 2.0, 3.0, 0.0});
    for (const CramOrder order : both_orders) {
        const std::vector<double> n = solve_cram(chain.system, {1.0, 0.0, 0.0, 0.0}, 4.0, order);
        for (const OracleRow& row : oracle) {
            INFO(decaysolver::to_string(order) << " " << row.nuclide);
            REQUIRE_THAT(n[chain.system.index_of(row.nuclide)], WithinAbs(row.value, 1e-14));
            if (order == CramOrder::order_48) {
                REQUIRE_THAT(n[chain.system.index_of(row.nuclide)], WithinRel(row.value, 1e-13));
            }
        }
    }
}

TEST_CASE("cram: chaine Ra-226 a 30 j, 1 a, 100 a vs oracle",
          "[cram][T2][verification][V1][data]") {
    const std::vector<OracleRow> oracle =
        read_oracle(verification_dir / "V1_analytic_bateman" / "oracle_ra226.csv");
    const DecaySystem system = DecaySystem::build(icrp_library(), {"Ra-226"});
    std::vector<double> n0(system.size(), 0.0);
    n0[system.index_of("Ra-226")] = 1.0;

    for (const CramOrder order : both_orders) {
        std::map<std::string, std::vector<double>> solutions;
        double worst_absolute = 0.0;
        double worst_relative = 0.0;
        for (const OracleRow& row : oracle) {
            if (!solutions.contains(row.case_name)) {
                solutions[row.case_name] = solve_cram(system, n0, row.t_s, order);
            }
            const double computed = solutions[row.case_name][system.index_of(row.nuclide)];
            INFO(decaysolver::to_string(order) << " " << row.case_name << " " << row.nuclide);
            // Précision absolue : c'est ce que CRAM garantit, quel que soit l'ordre.
            REQUIRE_THAT(computed, WithinAbs(row.value, 1e-14));
            worst_absolute = std::max(worst_absolute, std::abs(computed - row.value));
            if (row.value >= 1e-30) {
                const double relative = std::abs(computed / row.value - 1.0);
                worst_relative = std::max(worst_relative, relative);
                if (order == CramOrder::order_48) {
                    REQUIRE_THAT(computed, WithinRel(row.value, 1e-13));
                }
            }
        }
        // Visible avec `decaysolver_unit_tests -s "[cram][V1]"` ; repris dans le rapport.
        SUCCEED(decaysolver::to_string(order) << " : écart absolu max " << worst_absolute
                                              << ", relatif max (N ≥ 1e-30) " << worst_relative);
    }
}

TEST_CASE("cram: accord avec la solution analytique sur les 35 nucleides", "[cram][T2][data]") {
    // Un atome de chacun des 35 nucléides de la liste standard, 10 ans : le graphe complet
    // (embranchements, filles communes), comparé à la voie analytique (chemins + différences
    // divisées). Deux méthodes sans rien en commun, hormis la matrice.
    // Même liste que tests/regression/test_regression.cpp (liste standard de déclaration).
    const std::vector<std::string> seeds = {
        "Be-10",   "C-14",    "Cl-36",  "Ca-41",  "Mn-54",  "Fe-55",  "Co-60",  "Ni-59",  "Ni-63",
        "Zn-65",   "Se-79",   "Sr-90",  "Zr-93",  "Nb-94",  "Mo-93",  "Tc-99",  "Pd-107", "Ag-108m",
        "Ag-110m", "Sn-121m", "Sn-126", "Sb-125", "I-129",  "Cs-134", "Cs-135", "Cs-137", "Sm-151",
        "U-235",   "U-238",   "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Am-241", "Cm-244",
    };
    const DecaySystem system = DecaySystem::build(icrp_library(), seeds);
    std::vector<double> n0(system.size(), 0.0);
    for (const std::string& seed : seeds) {
        n0[system.index_of(seed)] = 1.0;
    }
    const double t = 10.0 * units::seconds_per_year;
    const std::vector<double> reference = solve_bateman(system, n0, t);
    for (const CramOrder order : both_orders) {
        const std::vector<double> n = solve_cram(system, n0, t, order);
        double worst = 0.0;
        for (std::size_t i = 0; i < system.size(); ++i) {
            INFO(decaysolver::to_string(order) << " " << system.names()[i]);
            REQUIRE_THAT(n[i], WithinAbs(reference[i], 1e-14));
            worst = std::max(worst, std::abs(n[i] - reference[i]));
        }
        SUCCEED(decaysolver::to_string(order) << " : écart absolu max vs Bateman " << worst
                                              << " sur " << system.size() << " nucléides");
    }
}

TEST_CASE("cram: invariants", "[cram][T4]") {
    SECTION("t = 0 : N(0) = N0 a l'arrondi pres") {
        const ChainFixture chain = make_linear_chain({1.0, 2.0, 0.0});
        for (const CramOrder order : both_orders) {
            const std::vector<double> n = solve_cram(chain.system, {1.0, 0.5, 0.0}, 0.0, order);
            REQUIRE_THAT(n[0], WithinRel(1.0, 1e-15));
            REQUIRE_THAT(n[1], WithinRel(0.5, 1e-15));
            REQUIRE_THAT(n[2], WithinAbs(0.0, 1e-16));
        }
    }
    SECTION("conservation du nombre d'atomes, rapports exacts (Sr-90 -> Y-90 -> Zr-90)") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Sr-90"});
        std::vector<double> n0(system.size(), 0.0);
        n0[system.index_of("Sr-90")] = 1.0;
        for (const CramOrder order : both_orders) {
            for (const double years : {1.0, 30.0, 300.0, 3000.0}) {
                const std::vector<double> n =
                    solve_cram(system, n0, years * units::seconds_per_year, order);
                double total = 0.0;
                for (const double value : n) {
                    total += value;
                }
                INFO(decaysolver::to_string(order) << " " << years << " a");
                REQUIRE_THAT(total, WithinAbs(1.0, 1e-14));
            }
        }
    }
    SECTION("semi-groupe : Phi(1100 a) = Phi(1000 a) o Phi(100 a), Ra-226") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Ra-226"});
        std::vector<double> n0(system.size(), 0.0);
        n0[system.index_of("Ra-226")] = 1.0;
        for (const CramOrder order : both_orders) {
            const std::vector<double> direct =
                solve_cram(system, n0, 1100.0 * units::seconds_per_year, order);
            const std::vector<double> composed =
                solve_cram(system, solve_cram(system, n0, 100.0 * units::seconds_per_year, order),
                           1000.0 * units::seconds_per_year, order);
            for (std::size_t i = 0; i < system.size(); ++i) {
                REQUIRE_THAT(direct[i], WithinAbs(composed[i], 1e-14));
            }
        }
    }
    SECTION("positivite : non garantie, mais bornee par la precision absolue") {
        // Sr-90 à 3000 ans : population vraie ≈ 10⁻³¹ ; CRAM peut rendre un nombre négatif de
        // l'ordre de son erreur absolue. On vérifie la borne, pas le signe : limitation déclarée.
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Sr-90"});
        std::vector<double> n0(system.size(), 0.0);
        n0[system.index_of("Sr-90")] = 1.0;
        for (const CramOrder order : both_orders) {
            const std::vector<double> n =
                solve_cram(system, n0, 3000.0 * units::seconds_per_year, order);
            for (const double value : n) {
                REQUIRE(value >= -1e-14);
            }
        }
    }
}

TEST_CASE("cram: raideur extreme sans mise a l'echelle", "[cram][T5]") {
    // Constantes étalées sur 19 ordres de grandeur, t tel que λ_max t = 10²⁰ : Padé + scaling-and-
    // squaring devrait élever 67 fois au carré. CRAM résout k/2 systèmes triangulaires, point.
    // Référence : la solution analytique (constantes distinctes : fiable, voir V1).
    const ChainFixture chain = make_linear_chain({1e4, 1.0, 1e-5, 1e-10, 1e-15, 0.0});
    const double t = 1e16;
    const std::vector<double> reference = solve_bateman(chain.system, {1.0, 0, 0, 0, 0, 0}, t);
    for (const CramOrder order : both_orders) {
        const std::vector<double> n = solve_cram(chain.system, {1.0, 0, 0, 0, 0, 0}, t, order);
        for (std::size_t i = 0; i < n.size(); ++i) {
            INFO(decaysolver::to_string(order) << " A-" << i + 1);
            REQUIRE_THAT(n[i], WithinAbs(reference[i], 1e-14));
        }
        // La seule population non négligeable est A-5 (λ = 1e-15, λt = 10) et le puits A-6.
        REQUIRE_THAT(n[4], WithinRel(reference[4], 1e-12));
        REQUIRE_THAT(n[5], WithinRel(reference[5], 1e-13));
    }
}

TEST_CASE("cram: erreurs d'appel", "[cram][T1]") {
    const ChainFixture chain = make_linear_chain({1.0, 0.0});
    REQUIRE_THROWS_AS(solve_cram(chain.system, {1.0}, 1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(solve_cram(chain.system, {1.0, 0.0}, -1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(solve_cram(chain.system, {1.0, 0.0}, std::nan("")), std::invalid_argument);
}
