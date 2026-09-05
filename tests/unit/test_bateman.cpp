// T2 (solutions analytiques), T4 (invariants) et T5 (cas dégénérés) pour la solution de Bateman.
//
// Valeurs de référence du jeu D : verification/scripts/oracle_degenerate.py (mpmath, 50 chiffres).
// Tolérances : 1e-13 relatif sur les cas sains (une dizaine d'opérations flottantes, chacune à
// 1 ulp ≈ 2e-16) ; 1e-12 sur les cas quasi dégénérés corrigés, où la série de Taylor ajoute une
// troncature contrôlée par epsilon.

#include <decaysolver/bateman.hpp>
#include <decaysolver/decay_system.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "chain_fixture.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using decaysolver::DecaySystem;
using decaysolver::divided_difference_exp;
using decaysolver::divided_difference_exp_naive;
using decaysolver::NuclideLibrary;
using decaysolver::solve_bateman;
namespace units = decaysolver::units;

namespace {

const NuclideLibrary& icrp_library() {
    static const NuclideLibrary library =
        NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    return library;
}

std::vector<double> unit_seed(std::size_t size, std::size_t index) {
    std::vector<double> n0(size, 0.0);
    n0[index] = 1.0;
    return n0;
}

// Égalité bit à bit de deux doubles, sans comparaison flottante.
bool same_bits(double a, double b) {
    return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

} // namespace

TEST_CASE("bateman: differences divisees de l'exponentielle", "[bateman][T2]") {
    const double t = 4.0;
    SECTION("un nœud : exp(x t)") {
        REQUIRE_THAT(divided_difference_exp({-1.0}, t), WithinRel(std::exp(-t), 1e-15));
    }
    SECTION("deux nœuds bien séparés : formule fermée") {
        const double expected = (std::exp(-2.0 * t) - std::exp(-1.0 * t)) / (-2.0 + 1.0);
        REQUIRE_THAT(divided_difference_exp({-1.0, -2.0}, t), WithinRel(expected, 1e-14));
        REQUIRE_THAT(divided_difference_exp_naive({-1.0, -2.0}, t), WithinRel(expected, 1e-14));
    }
    SECTION("nœuds confondus : dérivées, tᵏ e^{xt} / k!") {
        REQUIRE_THAT(divided_difference_exp({-1.0, -1.0}, t), WithinRel(t * std::exp(-t), 1e-14));
        REQUIRE_THAT(divided_difference_exp({-1.0, -1.0, -1.0}, t),
                     WithinRel(t * t / 2.0 * std::exp(-t), 1e-14));
        REQUIRE_THROWS_AS(divided_difference_exp_naive({-1.0, -1.0}, t), std::domain_error);
    }
    SECTION("ordre des nœuds indifférent") {
        REQUIRE_THAT(divided_difference_exp({-3.0, -1.0, -2.0}, t),
                     WithinRel(divided_difference_exp({-1.0, -2.0, -3.0}, t), 1e-15));
    }
    SECTION("t = 0 : exp[x₀..x_k] vaut 1 pour k = 0 et exactement 0 sinon") {
        REQUIRE(same_bits(divided_difference_exp({-1.0}, 0.0), 1.0));
        REQUIRE(same_bits(divided_difference_exp({-1.0, -2.0}, 0.0), 0.0));
        REQUIRE(same_bits(divided_difference_exp({-1.0, -2.0, -3.0}, 0.0), 0.0));
    }
}

TEST_CASE("bateman: chaine lineaire vs formules fermees", "[bateman][T2]") {
    SECTION("un nucléide : N₀ e^{−λt}") {
        const ChainFixture chain = make_linear_chain({0.3, 0.0});
        const std::vector<double> n = solve_bateman(chain.system, {2.0, 0.0}, 5.0);
        REQUIRE_THAT(n[0], WithinRel(2.0 * std::exp(-0.3 * 5.0), 1e-15));
        REQUIRE_THAT(n[1], WithinRel(2.0 * (1.0 - std::exp(-0.3 * 5.0)), 1e-14));
    }
    SECTION("deux nucléides radioactifs : N₁ = N₀ λ₀/(λ₁−λ₀) (e^{−λ₀t} − e^{−λ₁t})") {
        const double l0 = 0.5;
        const double l1 = 2.0;
        const double t = 1.5;
        const ChainFixture chain = make_linear_chain({l0, l1});
        const std::vector<double> n = solve_bateman(chain.system, {1.0, 0.0, 0.0}, t);
        REQUIRE_THAT(n[1],
                     WithinRel(l0 / (l1 - l0) * (std::exp(-l0 * t) - std::exp(-l1 * t)), 1e-14));
        REQUIRE_THAT(n[0] + n[1] + n[2], WithinRel(1.0, 1e-15)); // conservation
    }
    SECTION("D4 — (1, 2, 3), t = 4, contrôle sain") {
        const ChainFixture chain = make_linear_chain({1.0, 2.0, 3.0});
        const std::vector<double> n = solve_bateman(chain.system, unit_seed(4, 0), 4.0);
        REQUIRE_THAT(n[1], WithinRel(0.017980176260831668455, 1e-14));
        REQUIRE_THAT(n[2], WithinRel(0.017650857845282484826, 1e-14));
    }
}

TEST_CASE("bateman: cas degeneres D1-D3 corriges", "[bateman][T5]") {
    const double t = 4.0;
    SECTION("D1 — (1, 1) exactement dégénéré : N₁ = λ t e^{−λt}") {
        const ChainFixture chain = make_linear_chain({1.0, 1.0});
        const std::vector<double> n = solve_bateman(chain.system, unit_seed(3, 0), t);
        REQUIRE_THAT(n[1], WithinRel(0.073262555554936721175, 1e-14));
    }
    SECTION("D2 — (3, 3−1e-7, 3+1e-7)") {
        const ChainFixture chain = make_linear_chain({3.0, 3.0 - 1e-7, 3.0 + 1e-7});
        const std::vector<double> n = solve_bateman(chain.system, unit_seed(4, 0), t);
        REQUIRE_THAT(n[1], WithinRel(0.00007373056298605013124, 1e-12));
        REQUIRE_THAT(n[2], WithinRel(0.00044238327469352735308, 1e-12));
    }
    SECTION("D3 — (3, 3−1e-11, 3+1e-11)") {
        const ChainFixture chain = make_linear_chain({3.0, 3.0 - 1e-11, 3.0 + 1e-11});
        const std::vector<double> n = solve_bateman(chain.system, unit_seed(4, 0), t);
        REQUIRE_THAT(n[1], WithinRel(0.000073730548241413128069, 1e-12));
        REQUIRE_THAT(n[2], WithinRel(0.00044238328943815649166, 1e-12));
    }
}

TEST_CASE("bateman: formule naive, perte de precision pour lambdas quasi degeneres",
          "[bateman][T5][verification][known-limitation]") {
    // Limitation documentée, pas masquée : la somme fermée de Bateman en double perd environ
    // log10(λ / |λ_l − λ_k|) chiffres. Ce test fige les ordres de grandeur du cahier des charges.
    const double t = 4.0;
    SECTION("D2 — écart 1e-7 : erreur relative ~2e-3") {
        const double prefactor = 3.0 * (3.0 - 1e-7);
        const double naive =
            prefactor * divided_difference_exp_naive({-3.0, -(3.0 - 1e-7), -(3.0 + 1e-7)}, t);
        const double error = std::abs(naive / 0.00044238327469352735308 - 1.0);
        REQUIRE(error > 1e-4);
        REQUIRE(error < 1e-2);
    }
    SECTION("D3 — écart 1e-11 : le résultat vaut exactement 0") {
        const double naive =
            divided_difference_exp_naive({-3.0, -(3.0 - 1e-11), -(3.0 + 1e-11)}, t);
        REQUIRE(same_bits(naive, 0.0));
    }
    SECTION("D4 — écarts d'ordre 1 : la formule naïve est exacte") {
        const double naive = 1.0 * 2.0 * divided_difference_exp_naive({-1.0, -2.0, -3.0}, t);
        REQUIRE_THAT(naive, WithinRel(0.017650857845282484826, 1e-14));
    }
}

TEST_CASE("bateman: chaine reelle Sr-90 / Y-90", "[bateman][T2][data]") {
    const DecaySystem system = DecaySystem::build(icrp_library(), {"Sr-90"});
    const std::vector<double>& lambda = system.decay_constants_per_s();
    const double l_sr = lambda[0];
    const double l_y = lambda[1];

    SECTION("30 jours : formule à deux corps") {
        const double t = 30.0 * units::seconds_per_day;
        const std::vector<double> n = solve_bateman(system, {1.0, 0.0, 0.0}, t);
        REQUIRE_THAT(n[0], WithinRel(std::exp(-l_sr * t), 1e-15));
        REQUIRE_THAT(
            n[1],
            WithinRel(l_sr / (l_y - l_sr) * (std::exp(-l_sr * t) - std::exp(-l_y * t)), 1e-13));
    }
    SECTION("1 an : équilibre séculaire, A_Y / A_Sr = λ_Y / (λ_Y − λ_Sr)") {
        const double t = units::seconds_per_year;
        const std::vector<double> n = solve_bateman(system, {1.0, 0.0, 0.0}, t);
        const double activity_ratio = (l_y * n[1]) / (l_sr * n[0]);
        REQUIRE_THAT(activity_ratio, WithinRel(l_y / (l_y - l_sr), 1e-12));
    }
}

TEST_CASE("bateman: invariants physiques", "[bateman][T4][data]") {
    SECTION("N(0) = N₀ bit à bit, filles exactement nulles") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Cs-137"});
        const std::vector<double> n0{0.123456789, 0.0, 0.0};
        const std::vector<double> n = solve_bateman(system, n0, 0.0);
        REQUIRE(same_bits(n[0], n0[0]));
        REQUIRE(same_bits(n[1], 0.0));
        REQUIRE(same_bits(n[2], 0.0));
    }
    SECTION("positivité et conservation, Sr-90 sur 10 000 ans") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Sr-90"});
        for (const double years : {0.01, 1.0, 100.0, 10000.0}) {
            const std::vector<double> n =
                solve_bateman(system, {1.0, 0.0, 0.0}, years * units::seconds_per_year);
            for (const double value : n) {
                REQUIRE(value >= 0.0);
            }
            REQUIRE_THAT(std::accumulate(n.begin(), n.end(), 0.0), WithinRel(1.0, 1e-14));
        }
    }
    SECTION("semi-groupe Φ(t₁+t₂) = Φ(t₂)∘Φ(t₁) sur la chaîne du Ra-226") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Ra-226"});
        const double t1 = 100.0 * units::seconds_per_year;
        const double t2 = 1000.0 * units::seconds_per_year;
        const std::vector<double> n0 = unit_seed(system.size(), 0);
        const std::vector<double> direct = solve_bateman(system, n0, t1 + t2);
        const std::vector<double> composed =
            solve_bateman(system, solve_bateman(system, n0, t1), t2);
        for (std::size_t i = 0; i < system.size(); ++i) {
            INFO(system.names()[i]);
            REQUIRE_THAT(composed[i], WithinRel(direct[i], 1e-10) || WithinAbs(direct[i], 1e-250));
        }
        // Conservation : pas de fission dans la chaîne Ra-226 → Pb-206, mais les rapports
        // d'embranchement ICRP-107 sont arrondis (Bi-210 : 1 + 1,3e-6) ; l'exactitude à 1e-14
        // est testée ci-dessus sur Sr-90, dont les rapports valent exactement 1.
        REQUIRE_THAT(std::accumulate(direct.begin(), direct.end(), 0.0), WithinRel(1.0, 1e-5));
    }
    SECTION("entrées invalides") {
        const DecaySystem system = DecaySystem::build(icrp_library(), {"Sr-90"});
        REQUIRE_THROWS_AS(solve_bateman(system, {1.0}, 1.0), std::invalid_argument);
        REQUIRE_THROWS_AS(solve_bateman(system, {1.0, 0.0, 0.0}, -1.0), std::invalid_argument);
    }
}
