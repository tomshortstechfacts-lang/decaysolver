// T1 (lecture, durées) et T2 (vieillissement vs formules fermées) pour le mode inventaire.

#include <decaysolver/inventory.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using decaysolver::age_inventory;
using decaysolver::DataError;
using decaysolver::DaughterPolicy;
using decaysolver::Inventory;
using decaysolver::NuclideLibrary;
using decaysolver::read_inventory;
using decaysolver::ValueKind;
namespace units = decaysolver::units;

namespace {

const NuclideLibrary& icrp_library() {
    static const NuclideLibrary library =
        NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    return library;
}

Inventory inventory_from(const std::string& text, ValueKind kind = ValueKind::fraction) {
    std::istringstream in(text);
    return read_inventory(in, kind);
}

const decaysolver::AgedEntry& find(const decaysolver::AgedInventory& aged,
                                   const std::string& name) {
    for (const decaysolver::AgedEntry& entry : aged.entries) {
        if (entry.nuclide == name) {
            return entry;
        }
    }
    throw std::out_of_range(name);
}

} // namespace

TEST_CASE("inventory: lecture des fichiers d'entree", "[inventory][T1]") {
    SECTION("virgule décimale, notation scientifique, en-tête et commentaires") {
        const Inventory inv = inventory_from("# spectre de test\n"
                                             "Radionucleide;T0\n"
                                             "Be10;1,43E-07\n"
                                             "Cs-137 ; 1.10e-2\n"
                                             "\n"
                                             "pu241\t0.121\n",
                                             ValueKind::fraction);
        REQUIRE(inv.entries.size() == 3);
        REQUIRE(inv.entries[0].nuclide == "Be-10");
        REQUIRE_THAT(inv.entries[0].value, WithinRel(1.43e-7, 1e-15));
        REQUIRE(inv.entries[1].nuclide == "Cs-137");
        REQUIRE(inv.entries[2].nuclide == "Pu-241");
        REQUIRE(inv.kind == ValueKind::fraction);
    }
    SECTION("rejets : doublon, valeur négative, nom invalide, ligne sans séparateur, vide") {
        REQUIRE_THROWS_AS(inventory_from("Cs137;1\nCs-137;2\n"), DataError);
        REQUIRE_THROWS_AS(inventory_from("Cs137;-1\n"), DataError);
        REQUIRE_THROWS_AS(inventory_from("Cesium;1\n"), DataError);
        REQUIRE_THROWS_AS(inventory_from("Cs137 1\n"), DataError);
        REQUIRE_THROWS_AS(inventory_from("# rien\n"), DataError);
        REQUIRE_THROWS_AS(inventory_from("nuclide;valeur\nCs137;abc\n"), DataError);
    }
}

TEST_CASE("inventory: durees", "[inventory][T1]") {
    REQUIRE_THAT(decaysolver::parse_duration_s("+6a"),
                 WithinRel(6.0 * units::seconds_per_year, 1e-15));
    REQUIRE_THAT(decaysolver::parse_duration_s("30 j"),
                 WithinRel(30.0 * units::seconds_per_day, 1e-15));
    REQUIRE_THAT(decaysolver::parse_duration_s("1.5h"), WithinRel(5400.0, 1e-15));
    REQUIRE_THAT(decaysolver::parse_duration_s("90min"), WithinRel(5400.0, 1e-15));
    REQUIRE_THAT(decaysolver::parse_duration_s("0s"), WithinAbs(0.0, 0.0));
    REQUIRE_THROWS_AS(decaysolver::parse_duration_s("6 ans"), std::invalid_argument);
    REQUIRE_THROWS_AS(decaysolver::parse_duration_s("-1a"), std::invalid_argument);
    REQUIRE_THROWS_AS(decaysolver::parse_duration_s("a"), std::invalid_argument);
}

TEST_CASE("inventory: vieillissement d'un seul nucleide", "[inventory][T2][data]") {
    const Inventory inv = inventory_from("Co60;1000\n", ValueKind::activity_bq);
    const double t = 5.2713 * units::seconds_per_year; // une demi-vie ICRP-107
    const decaysolver::AgedInventory aged =
        age_inventory(icrp_library(), inv, t, DaughterPolicy::input_only);
    REQUIRE(aged.entries.size() == 1);
    REQUIRE_THAT(aged.entries[0].activity, WithinRel(500.0, 1e-12));
    REQUIRE_THAT(aged.entries[0].fraction, WithinRel(1.0, 1e-15));
    REQUIRE_THAT(aged.beta_gamma, WithinRel(500.0, 1e-12));
    REQUIRE_THAT(aged.alpha, WithinAbs(0.0, 0.0));
}

TEST_CASE("inventory: convention des filles", "[inventory][T2][data]") {
    const Inventory inv = inventory_from("Cs137;1\n", ValueKind::activity_bq);
    const double t = 10.0 * units::seconds_per_year;
    SECTION("input-only : seul Cs-137 est listé") {
        const decaysolver::AgedInventory aged =
            age_inventory(icrp_library(), inv, t, DaughterPolicy::input_only);
        REQUIRE(aged.entries.size() == 1);
        REQUIRE(aged.entries[0].nuclide == "Cs-137");
    }
    SECTION("all : Ba-137m apparaît à l'équilibre, 94,4 % de l'activité du Cs-137") {
        const decaysolver::AgedInventory aged =
            age_inventory(icrp_library(), inv, t, DaughterPolicy::all);
        REQUIRE(aged.entries.size() == 2); // Ba-137 stable exclu
        const double a_cs = find(aged, "Cs-137").activity;
        const double a_ba = find(aged, "Ba-137m").activity;
        // Équilibre séculaire : A_fille = b · A_mère · λ_f/(λ_f − λ_m), avec λ_f ≫ λ_m.
        const double l_cs = icrp_library().get("Cs-137").decay_constant_per_s();
        const double l_ba = icrp_library().get("Ba-137m").decay_constant_per_s();
        REQUIRE_THAT(a_ba, WithinRel(0.94399 * a_cs * l_ba / (l_ba - l_cs), 1e-10));
        REQUIRE_THAT(find(aged, "Cs-137").fraction + find(aged, "Ba-137m").fraction,
                     WithinRel(1.0, 1e-15));
    }
}

TEST_CASE("inventory: croissance de l'Am-241 a partir du Pu-241", "[inventory][T2][data]") {
    // Cas type des spectres de déchets : Pu-241 (β⁻, 14,35 a) alimente Am-241 (α, 432,2 a).
    const Inventory inv = inventory_from("Pu241;0.121\nAm241;3.93e-4\n", ValueKind::fraction);
    const double t = 6.0 * units::seconds_per_year;
    const decaysolver::AgedInventory aged =
        age_inventory(icrp_library(), inv, t, DaughterPolicy::input_only);

    const double l_pu = icrp_library().get("Pu-241").decay_constant_per_s();
    const double l_am = icrp_library().get("Am-241").decay_constant_per_s();
    const double b = 0.99998; // rapport β⁻ ICRP-107
    const double a_pu = 0.121 * std::exp(-l_pu * t);
    const double a_am =
        3.93e-4 * std::exp(-l_am * t) +
        b * 0.121 * l_am / (l_am - l_pu) * (std::exp(-l_pu * t) - std::exp(-l_am * t));
    REQUIRE_THAT(find(aged, "Pu-241").activity, WithinRel(a_pu, 1e-13));
    REQUIRE_THAT(find(aged, "Am-241").activity, WithinRel(a_am, 1e-12));
    REQUIRE_THAT(aged.alpha, WithinRel(a_am, 1e-12)); // Am-241 α ; Pu-241 compte en β-γ
    REQUIRE_THAT(aged.beta_gamma, WithinRel(a_pu, 1e-13));
    REQUIRE_THAT(aged.total, WithinRel(a_pu + a_am, 1e-13));
}

TEST_CASE("inventory: entrees rejetees", "[inventory][T1][data]") {
    REQUIRE_THROWS_AS(
        age_inventory(icrp_library(), inventory_from("Pb206;1\n"), 1.0, DaughterPolicy::input_only),
        DataError);
    REQUIRE_THROWS_AS(
        age_inventory(icrp_library(), inventory_from("Xx99;1\n"), 1.0, DaughterPolicy::input_only),
        std::out_of_range);
    REQUIRE_THROWS_AS(age_inventory(icrp_library(), inventory_from("Cs137;1\n"), -1.0,
                                    DaughterPolicy::input_only),
                      std::invalid_argument);
}

TEST_CASE("inventory: sortie CSV avec provenance", "[inventory][T1][data]") {
    const Inventory inv = inventory_from("Cs137;1\n", ValueKind::activity_bq);
    const decaysolver::AgedInventory aged =
        age_inventory(icrp_library(), inv, units::seconds_per_year, DaughterPolicy::input_only);
    std::ostringstream out;
    decaysolver::write_aged_inventory(out, aged, icrp_library(), "test");
    const std::string text = out.str();
    REQUIRE(text.find("# decaysolver: ") != std::string::npos);
    REQUIRE(text.find("# data: source: ICRP") != std::string::npos);
    REQUIRE(text.find("# daughters: input-only") != std::string::npos);
    REQUIRE(text.find("nuclide;activity;fraction;primary_mode\nCs-137;") != std::string::npos);
}
