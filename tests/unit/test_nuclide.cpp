// T1 — noms de nucléides, modes de décroissance, propriétés d'un nucléide.

#include <decaysolver/nuclide.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("nuclide: normalisation des noms", "[nuclide][T1]") {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"Cs-137", "Cs-137"},     {"Cs137", "Cs-137"},    {"cs137", "Cs-137"},
        {"CS-137", "Cs-137"},     {"137Cs", "Cs-137"},    {"137cs", "Cs-137"},
        {"Cs 137", "Cs-137"},     {"Ag-108m", "Ag-108m"}, {"Ag108m", "Ag-108m"},
        {"ag108M", "Ag-108m"},    {"U-235m", "U-235m"},   {"U235", "U-235"},
        {"235U", "U-235"},        {"H-3", "H-3"},         {"3H", "H-3"},
        {"Am-242m1", "Am-242m1"},
    };
    for (const auto& [input, expected] : cases) {
        INFO("entrée : '" << input << "'");
        REQUIRE(decaysolver::canonical_name(input) == expected);
    }
}

TEST_CASE("nuclide: graphies rejetees", "[nuclide][T1]") {
    for (const char* bad : {"", "Cs", "137", "Cs-", "Cesium-137", "Cs-1370", "Cs-0", "Cs-137x",
                            "Cs-137mm", "137mBa", "Cs-13-7", "-"}) {
        INFO("entrée : '" << bad << "'");
        REQUIRE_THROWS_AS(decaysolver::canonical_name(bad), std::invalid_argument);
    }
}

TEST_CASE("nuclide: modes de decroissance", "[nuclide][T1]") {
    using decaysolver::DecayMode;
    for (const DecayMode mode :
         {DecayMode::alpha, DecayMode::beta_minus, DecayMode::beta_plus_ec,
          DecayMode::isomeric_transition, DecayMode::spontaneous_fission, DecayMode::stable}) {
        REQUIRE(decaysolver::decay_mode_from_string(decaysolver::to_string(mode)) == mode);
    }
    REQUIRE(decaysolver::decay_mode_from_string("beta-") == DecayMode::beta_minus);
    REQUIRE_THROWS_AS(decaysolver::decay_mode_from_string("beta"), std::invalid_argument);
    REQUIRE_THROWS_AS(decaysolver::decay_mode_from_string("ALPHA"), std::invalid_argument);
}

TEST_CASE("nuclide: proprietes", "[nuclide][T1]") {
    using decaysolver::DecayMode;
    namespace units = decaysolver::units;

    SECTION("nucléide stable") {
        const decaysolver::Nuclide stable{"Ba-137", std::numeric_limits<double>::infinity(), {}};
        REQUIRE(stable.is_stable());
        REQUIRE_THAT(stable.decay_constant_per_s(), WithinAbs(0.0, 0.0));
        REQUIRE(stable.primary_mode() == DecayMode::stable);
    }
    SECTION("Cs-137 : deux voies β⁻, mode principal β⁻, λ = ln2 / T½") {
        const decaysolver::Nuclide cs137{"Cs-137",
                                         30.1671 * units::seconds_per_year,
                                         {{DecayMode::beta_minus, "Ba-137m", 0.94399},
                                          {DecayMode::beta_minus, "Ba-137", 0.056005}}};
        REQUIRE_FALSE(cs137.is_stable());
        REQUIRE(cs137.primary_mode() == DecayMode::beta_minus);
        REQUIRE_THAT(cs137.decay_constant_per_s(),
                     WithinRel(units::decay_constant_per_s(30.1671 * units::seconds_per_year),
                               std::numeric_limits<double>::epsilon()));
    }
    SECTION("Pu-241 : β⁻ à 99,998 % et α à 2,45e-5 → mode principal β⁻") {
        const decaysolver::Nuclide pu241{
            "Pu-241",
            14.35 * units::seconds_per_year,
            {{DecayMode::alpha, "U-237", 2.45e-5}, {DecayMode::beta_minus, "Am-241", 0.99998}}};
        REQUIRE(pu241.primary_mode() == DecayMode::beta_minus);
    }
}
