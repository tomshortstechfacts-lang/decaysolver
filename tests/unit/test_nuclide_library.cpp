// T1 — chargement et validation de la bibliothèque de nucléides.
//
// Deux familles : un CSV minimal écrit dans le test (chaque invariant violé doit être rejeté), et
// le fichier réel data/nuclides_icrp107.csv (chargement complet, valeurs de contrôle, invariants
// vérifiés sur tout le jeu).

#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

using Catch::Matchers::WithinRel;
using decaysolver::DataError;
using decaysolver::DecayMode;
using decaysolver::NuclideLibrary;
namespace units = decaysolver::units;

namespace {

const std::string header =
    "nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction\n";

NuclideLibrary parse_text(const std::string& body) {
    std::istringstream in("# provenance de test\n" + header + body);
    return NuclideLibrary::parse(in);
}

} // namespace

TEST_CASE("library: chaine minimale A -> B -> C(stable)", "[library][T1]") {
    const NuclideLibrary lib = parse_text("Sr-90;28.79;y;beta-;Y-90;1\n"
                                          "Y-90;64.0;h;beta-;Zr-90;1\n"
                                          "Zr-90;stable;;stable;;\n");
    REQUIRE(lib.size() == 3);
    REQUIRE(lib.contains("sr90"));
    REQUIRE_FALSE(lib.contains("Cs-137"));
    REQUIRE(lib.provenance().size() == 1);
    REQUIRE(lib.names() == std::vector<std::string>{"Sr-90", "Y-90", "Zr-90"});

    const decaysolver::Nuclide& y90 = lib.get("90Y");
    REQUIRE_THAT(y90.half_life_s, WithinRel(64.0 * units::seconds_per_hour, 1e-15));
    REQUIRE(y90.branches.size() == 1);
    REQUIRE(y90.branches[0].daughter == "Zr-90");
    REQUIRE(lib.get("Zr-90").is_stable());
    REQUIRE_THROWS_AS(lib.get("Cs-137"), std::out_of_range);
}

TEST_CASE("library: invariants violes rejetes", "[library][T1]") {
    SECTION("somme des rapports d'embranchement différente de 1") {
        REQUIRE_THROWS_AS(parse_text("Cs-137;30.1671;y;beta-;Ba-137m;0.94399\n"
                                     "Ba-137m;2.552;m;IT;Ba-137;1\n"
                                     "Ba-137;stable;;stable;;\n"),
                          DataError);
    }
    SECTION("écart à 1 sous la tolérance accepté (arrondis ICRP-107)") {
        REQUIRE_NOTHROW(parse_text("Cs-137;30.1671;y;beta-;Ba-137m;0.94399\n"
                                   "Cs-137;30.1671;y;beta-;Ba-137;0.056005\n"
                                   "Ba-137m;2.552;m;IT;Ba-137;1\n"
                                   "Ba-137;stable;;stable;;\n"));
    }
    SECTION("fille absente de la bibliothèque") {
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;y;beta-;Y-90;1\n"), DataError);
    }
    SECTION("cycle de filiation") {
        REQUIRE_THROWS_AS(parse_text("A-1;1;s;beta-;B-1;1\n"
                                     "B-1;1;s;beta-;C-1;1\n"
                                     "C-1;1;s;beta-;A-1;1\n"),
                          DataError);
    }
    SECTION("demi-vie négative") {
        REQUIRE_THROWS_AS(parse_text("Sr-90;-28.79;y;beta-;Y-90;1\nY-90;stable;;stable;;\n"),
                          DataError);
    }
    SECTION("unité inconnue") {
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;an;beta-;Y-90;1\nY-90;stable;;stable;;\n"),
                          DataError);
    }
    SECTION("rapport d'embranchement nul ou > 1") {
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;y;beta-;Y-90;0\nY-90;stable;;stable;;\n"),
                          DataError);
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;y;beta-;Y-90;1.5\nY-90;stable;;stable;;\n"),
                          DataError);
    }
    SECTION("stable avec une voie, radioactif sans voie") {
        REQUIRE_THROWS_AS(parse_text("Zr-90;stable;;beta-;Y-90;1\nY-90;stable;;stable;;\n"),
                          DataError);
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;y;stable;;\n"), DataError);
    }
    SECTION("fission spontanée : fille vide obligatoire, produits non suivis") {
        REQUIRE_NOTHROW(parse_text("U-238;4.468;y;alpha;Th-234;1\n"
                                   "U-238;4.468;y;SF;;5.45e-07\n"
                                   "Th-234;stable;;stable;;\n"));
        REQUIRE_THROWS_AS(parse_text("U-238;4.468;y;SF;Th-234;1\nTh-234;stable;;stable;;\n"),
                          DataError);
    }
    SECTION("en-tête absent ou champs manquants") {
        std::istringstream no_header("Sr-90;28.79;y;beta-;Y-90;1\n");
        REQUIRE_THROWS_AS(NuclideLibrary::parse(no_header), DataError);
        REQUIRE_THROWS_AS(parse_text("Sr-90;28.79;y;beta-;Y-90\n"), DataError);
    }
    SECTION("nombre illisible") {
        REQUIRE_THROWS_AS(parse_text("Sr-90;28,79;y;beta-;Y-90;1\nY-90;stable;;stable;;\n"),
                          DataError);
    }
}

TEST_CASE("library: fichier ICRP-107 du depot", "[library][T1][data]") {
    const std::filesystem::path path =
        std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv";
    const NuclideLibrary lib = NuclideLibrary::load(path);

    SECTION("contenu attendu") {
        REQUIRE(lib.size() == 139);
        REQUIRE_FALSE(lib.provenance().empty());
        for (const char* name : {"Be-10", "C-14", "Cs-137", "Ba-137m", "Pu-241", "Am-241", "Cm-244",
                                 "Pb-206", "Zr-90", "Rn-217"}) {
            INFO(name);
            REQUIRE(lib.contains(name));
        }
    }
    SECTION("valeurs de contrôle (ICRP-107, converties avec 365,25 j/an)") {
        REQUIRE_THAT(lib.get("Cs-137").half_life_s,
                     WithinRel(30.1671 * units::seconds_per_year, 1e-15));
        REQUIRE_THAT(lib.get("Mn-54").half_life_s,
                     WithinRel(312.12 * units::seconds_per_day, 1e-15));
        REQUIRE_THAT(lib.get("Ba-137m").half_life_s,
                     WithinRel(2.552 * units::seconds_per_minute, 1e-15));
        REQUIRE(lib.get("Pu-241").primary_mode() == DecayMode::beta_minus);
        REQUIRE(lib.get("Am-241").primary_mode() == DecayMode::alpha);
        REQUIRE(lib.get("Ag-108m").primary_mode() == DecayMode::beta_plus_ec);
        REQUIRE(lib.get("Ba-137").is_stable());
        REQUIRE(lib.get("Pu-241").branches.size() == 2);
    }
    SECTION("invariants sur tout le jeu") {
        for (const std::string& name : lib.names()) {
            const decaysolver::Nuclide& nuclide = lib.get(name);
            INFO(name);
            if (nuclide.is_stable()) {
                REQUIRE(nuclide.branches.empty());
                continue;
            }
            double sum = 0.0;
            for (const decaysolver::DecayBranch& branch : nuclide.branches) {
                sum += branch.branching_fraction;
                if (branch.mode != DecayMode::spontaneous_fission) {
                    REQUIRE(lib.contains(branch.daughter));
                }
            }
            REQUIRE(std::abs(sum - 1.0) <= NuclideLibrary::branching_tolerance);
            REQUIRE(nuclide.decay_constant_per_s() > 0.0);
        }
    }
}
