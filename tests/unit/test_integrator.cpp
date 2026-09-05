// T3 (ordre de convergence, contrôle rapide) et T5 (raideur) pour les intégrateurs.
//
// Le protocole complet des ordres de convergence (oracle multiprécision, ≥ 4 raffinements, normes
// L∞ et L2, figures) vit dans verification/ ; ici on fige un contrôle à deux raffinements sur le
// problème non raide du cahier des charges, λ = (1, 2, 3, 0), T = 4, référence = solution de
// Bateman en double (erreur ~1e-15, très inférieure aux erreurs de troncature mesurées, ≥ 1e-11).

#include <decaysolver/bateman.hpp>
#include <decaysolver/decay_system.hpp>
#include <decaysolver/integrator.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "chain_fixture.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using decaysolver::DecaySystem;
using decaysolver::integrate;
using decaysolver::NuclideLibrary;
using decaysolver::Scheme;
using decaysolver::solve_bateman;
namespace units = decaysolver::units;

namespace {

const NuclideLibrary& icrp_library() {
    static const NuclideLibrary library =
        NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    return library;
}

double max_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        worst = std::max(worst, std::abs(a[i] - b[i]));
    }
    return worst;
}

} // namespace

TEST_CASE("integrator: noms et ordres", "[integrator][T1]") {
    for (const Scheme scheme :
         {Scheme::euler_explicit, Scheme::euler_implicit, Scheme::rk4, Scheme::crank_nicolson}) {
        REQUIRE(decaysolver::scheme_from_string(decaysolver::to_string(scheme)) == scheme);
    }
    REQUIRE(decaysolver::theoretical_order(Scheme::rk4) == 4);
    REQUIRE(decaysolver::theoretical_order(Scheme::crank_nicolson) == 2);
    REQUIRE_THROWS_AS(decaysolver::scheme_from_string("euler"), std::invalid_argument);
}

TEST_CASE("integrator: ordre observe sur le probleme non raide", "[integrator][T3]") {
    const ChainFixture chain = make_linear_chain({1.0, 2.0, 3.0, 0.0});
    const std::vector<double> n0{1.0, 0.0, 0.0, 0.0};
    const double t_end = 4.0;
    const std::vector<double> reference = solve_bateman(chain.system, n0, t_end);

    const Scheme scheme = GENERATE(Scheme::euler_explicit, Scheme::euler_implicit,
                                   Scheme::crank_nicolson, Scheme::rk4);
    INFO(decaysolver::to_string(scheme));
    // h = 4/256 = 1/64 : h·λ_max = 0,047, régime asymptotique ; erreurs entre 1e-3 et 1e-11,
    // loin de l'arrondi.
    const double error_coarse =
        max_abs_error(integrate(chain.system, n0, t_end, 256, scheme), reference);
    const double error_fine =
        max_abs_error(integrate(chain.system, n0, t_end, 512, scheme), reference);
    const double observed_order = std::log2(error_coarse / error_fine);
    REQUIRE_THAT(observed_order,
                 WithinAbs(static_cast<double>(decaysolver::theoretical_order(scheme)), 0.1));
}

TEST_CASE("integrator: proprietes elementaires", "[integrator][T1]") {
    const ChainFixture chain = make_linear_chain({1.0, 2.0, 0.0});
    const std::vector<double> n0{1.0, 0.0, 0.0};
    SECTION("zéro pas rejeté, tailles et temps invalides rejetés") {
        REQUIRE_THROWS_AS(integrate(chain.system, n0, 1.0, 0, Scheme::rk4), std::invalid_argument);
        REQUIRE_THROWS_AS(integrate(chain.system, {1.0}, 1.0, 10, Scheme::rk4),
                          std::invalid_argument);
        REQUIRE_THROWS_AS(integrate(chain.system, n0, -1.0, 10, Scheme::rk4),
                          std::invalid_argument);
    }
    SECTION("t = 0 rend N₀ pour tout schéma") {
        for (const Scheme scheme : {Scheme::euler_explicit, Scheme::euler_implicit, Scheme::rk4,
                                    Scheme::crank_nicolson}) {
            const std::vector<double> n = integrate(chain.system, n0, 0.0, 3, scheme);
            REQUIRE_THAT(n[0], WithinRel(1.0, 0.0));
            REQUIRE_THAT(n[1], WithinAbs(0.0, 0.0));
        }
    }
    SECTION("conservation du nombre d'atomes à chaque pas (colonnes de A de somme nulle)") {
        for (const Scheme scheme : {Scheme::euler_explicit, Scheme::euler_implicit, Scheme::rk4,
                                    Scheme::crank_nicolson}) {
            const std::vector<double> n = integrate(chain.system, n0, 3.0, 50, scheme);
            REQUIRE_THAT(std::accumulate(n.begin(), n.end(), 0.0), WithinRel(1.0, 1e-13));
        }
    }
}

TEST_CASE("integrator: raideur, chaine du Ra-226 avec un pas d'un an", "[integrator][T5][data]") {
    // Rn-222 (3,8 j) comme graine : λ_Rn · h ≈ 66 pour h = 1 an. Loin, très loin de la limite de
    // stabilité des schémas explicites (h·λ ≤ 2).
    const DecaySystem system = DecaySystem::build(icrp_library(), {"Rn-222"});
    const std::size_t rn = system.index_of("Rn-222");
    std::vector<double> n0(system.size(), 0.0);
    n0[rn] = 1.0;
    const double h = units::seconds_per_year;

    SECTION("Euler implicite : positif, conservatif, L-stable") {
        const std::vector<double> n = integrate(system, n0, 10.0 * h, 10, Scheme::euler_implicit);
        for (std::size_t i = 0; i < n.size(); ++i) {
            INFO(system.names()[i]);
            REQUIRE(n[i] >= 0.0);
        }
        // Conservation au niveau des arrondis d'ICRP-107 sur les rapports d'embranchement
        // (Bi-210 : 1 + 1,3e-6), pas de l'arithmétique : l'exactitude à 1e-13 est testée sur
        // une chaîne synthétique dans « propriétés élémentaires ».
        REQUIRE_THAT(std::accumulate(n.begin(), n.end(), 0.0), WithinRel(1.0, 1e-5));
        // Le mode raide est amorti à chaque pas d'un facteur 1/(1 + hλ) ≈ 1/67.
        REQUIRE(n[rn] < 1e-15);
    }
    SECTION(
        "Crank–Nicolson : A-stable mais non L-stable, concentration négative dès le premier pas") {
        // R(z) = (1 + z/2)/(1 − z/2) → −1 quand z → −∞ : le mode raide est réfléchi, pas amorti.
        const std::vector<double> n = integrate(system, n0, h, 1, Scheme::crank_nicolson);
        REQUIRE(n[rn] < 0.0);
        const double z = -system.decay_constants_per_s()[rn] * h;
        REQUIRE_THAT(n[rn], WithinRel((1.0 + z / 2.0) / (1.0 - z / 2.0), 1e-12));
    }
    SECTION("Euler explicite : divergence (D5)") {
        // Facteur d'amplification 1 − hλ ≈ −65 par pas : |N| explose et change de signe.
        const std::vector<double> n = integrate(system, n0, 10.0 * h, 10, Scheme::euler_explicit);
        REQUIRE(std::abs(n[rn]) > 1e15);
    }
}

TEST_CASE("integrator: Euler implicite vs Bateman, mode inventaire type",
          "[integrator][T2][data]") {
    // Cs-137 (30 a) sur 30 ans avec un pas d'un mois : erreur relative de troncature ≈ h λ / 2
    // par pas, cumulée ≈ (hλ)·(λT)/2 ≈ 2e-3 · 0,69 / 2 ≈ 7e-4. On vérifie qu'on est dans cet ordre.
    const DecaySystem system = DecaySystem::build(icrp_library(), {"Cs-137"});
    const std::vector<double> n0{1.0, 0.0, 0.0};
    const double t_end = 30.0 * units::seconds_per_year;
    const std::vector<double> exact = solve_bateman(system, n0, t_end);
    const std::vector<double> approx = integrate(system, n0, t_end, 360, Scheme::euler_implicit);
    REQUIRE_THAT(approx[0], WithinRel(exact[0], 2e-3));
    REQUIRE_THAT(approx[2], WithinRel(exact[2], 2e-3));
    const std::vector<double> rk4 = integrate(system, n0, t_end, 360, Scheme::rk4);
    REQUIRE_THAT(rk4[0], WithinRel(exact[0], 1e-10));
}
