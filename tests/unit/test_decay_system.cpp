// T1 — construction du système dN/dt = A·N depuis la bibliothèque : fermeture, ordre
// topologique, termes de production, second membre.

#include <decaysolver/decay_system.hpp>
#include <decaysolver/nuclide_library.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using decaysolver::DecaySystem;
using decaysolver::NuclideLibrary;

namespace {

const NuclideLibrary& icrp_library() {
    static const NuclideLibrary library =
        NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    return library;
}

} // namespace

TEST_CASE("system: chaine Sr-90 -> Y-90 -> Zr-90", "[system][T1][data]") {
    const DecaySystem system = DecaySystem::build(icrp_library(), {"sr90"});

    REQUIRE(system.size() == 3);
    REQUIRE(system.names() == std::vector<std::string>{"Sr-90", "Y-90", "Zr-90"});
    REQUIRE(system.index_of("Y-90") == 1);
    REQUIRE(system.contains("Zr-90"));
    REQUIRE_FALSE(system.contains("Cs-137"));
    REQUIRE_THROWS_AS(system.index_of("Cs-137"), std::out_of_range);

    const std::vector<double>& lambdas = system.decay_constants_per_s();
    REQUIRE(lambdas[0] > 0.0);
    REQUIRE(lambdas[1] > lambdas[0]); // Y-90 (64 h) décroît bien plus vite que Sr-90 (28,8 a)
    REQUIRE_THAT(lambdas[2], WithinAbs(0.0, 0.0));

    // Y-90 est produit par Sr-90 avec un rapport 1 ; Zr-90 par Y-90.
    REQUIRE(system.productions()[0].empty());
    REQUIRE(system.productions()[1].size() == 1);
    REQUIRE(system.productions()[1][0].parent == 0);
    REQUIRE_THAT(system.productions()[1][0].branching_fraction, WithinRel(1.0, 1e-15));
    REQUIRE(system.daughters()[0].size() == 1);
    REQUIRE(system.daughters()[0][0].parent == 1); // ici « parent » désigne l'indice de la fille
    REQUIRE(system.daughters()[2].empty());

    // Second membre à N = (1, 0, 0) : (−λ_Sr, +λ_Sr, 0).
    const std::vector<double> rate = system.rate({1.0, 0.0, 0.0});
    REQUIRE_THAT(rate[0], WithinRel(-lambdas[0], 1e-15));
    REQUIRE_THAT(rate[1], WithinRel(lambdas[0], 1e-15));
    REQUIRE_THAT(rate[2], WithinAbs(0.0, 0.0));
    REQUIRE_THROWS_AS(system.rate({1.0, 0.0}), std::invalid_argument);
}

TEST_CASE("system: ordre topologique sur un graphe branche", "[system][T1][data]") {
    // Cs-137 → Ba-137m → Ba-137 et Cs-137 → Ba-137 : deux chemins vers Ba-137, qui doit venir
    // après Ba-137m.
    const DecaySystem system = DecaySystem::build(icrp_library(), {"Cs-137"});
    REQUIRE(system.size() == 3);
    REQUIRE(system.index_of("Cs-137") < system.index_of("Ba-137m"));
    REQUIRE(system.index_of("Ba-137m") < system.index_of("Ba-137"));
    REQUIRE(system.productions()[system.index_of("Ba-137")].size() == 2);

    // Tout parent précède ses filles, pour tout le jeu ICRP-107.
    const DecaySystem all = DecaySystem::build(icrp_library(), icrp_library().names());
    REQUIRE(all.size() == icrp_library().size());
    for (std::size_t j = 0; j < all.size(); ++j) {
        for (const decaysolver::Production& daughter : all.daughters()[j]) {
            REQUIRE(daughter.parent > j);
        }
    }
}

TEST_CASE("system: bibliotheque complete chargee et ordonnee", "[system][T1][data]") {
    const NuclideLibrary full = NuclideLibrary::load(std::filesystem::path(DECAYSOLVER_DATA_DIR) /
                                                     "nuclides_icrp107_full.csv");
    REQUIRE(full.size() == 1512);
    REQUIRE(full.contains("Bi-212n"));
    const DecaySystem all = DecaySystem::build(full, full.names());
    REQUIRE(all.size() == 1512);
    for (std::size_t j = 0; j < all.size(); ++j) {
        for (const decaysolver::Production& daughter : all.daughters()[j]) {
            REQUIRE(daughter.parent > j);
        }
    }
}

TEST_CASE("system: fission spontanee exclue, graine inconnue rejetee", "[system][T1][data]") {
    const DecaySystem system = DecaySystem::build(icrp_library(), {"U-238"});
    REQUIRE(system.daughters()[system.index_of("U-238")].size() == 1); // Th-234 seulement
    REQUIRE(system.contains("Pb-206"));
    REQUIRE_THROWS_AS(DecaySystem::build(icrp_library(), {"Xx-99"}), std::out_of_range);
}
