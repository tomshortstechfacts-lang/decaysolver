// V1 — solution analytique vs oracle multiprécision (T2, CDC §4.3).
//
// Deux fichiers d'oracle, générés par des scripts versionnés (verification/scripts/) :
//  - oracle_lambda123.csv : chaîne (1, 2, 3, 0), T = 4 s ;
//  - oracle_ra226.csv     : chaîne réelle Ra-226 → Pb-206, toutes voies ICRP-107, à 30 j, 1 a, 100
//  a.
// L'oracle évalue la formule fermée à 50 chiffres ; decaysolver évalue des différences divisées
// en double. Deux implémentations, deux arithmétiques : l'accord n'est pas une tautologie.
//
// Tolérances : 1e-13 relatif sur la chaîne courte ; 1e-11 sur la chaîne du Ra-226 (15 nucléides,
// constantes étalées sur 14 ordres de grandeur, jusqu'à 15 nœuds par différence divisée) pour les
// populations au-dessus de 1e-30 ; en dessous, écart absolu 1e-40 (valeurs sans signification
// physique, la population totale valant 1).

#include <decaysolver/bateman.hpp>
#include <decaysolver/decay_system.hpp>
#include <decaysolver/nuclide_library.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "chain_fixture.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

struct OracleRow {
    std::string case_name;
    double t_s;
    std::string nuclide;
    double value;
};

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

} // namespace

TEST_CASE("oracle V1: chaine (1,2,3,0) a T = 4 s", "[bateman][T2][verification][V1]") {
    const std::vector<OracleRow> oracle =
        read_oracle(verification_dir / "V2_order_of_accuracy" / "oracle_lambda123.csv");
    const ChainFixture chain = make_linear_chain({1.0, 2.0, 3.0, 0.0});
    const std::vector<double> n =
        decaysolver::solve_bateman(chain.system, {1.0, 0.0, 0.0, 0.0}, 4.0);
    for (const OracleRow& row : oracle) {
        INFO(row.nuclide);
        REQUIRE_THAT(n[chain.system.index_of(row.nuclide)], WithinRel(row.value, 1e-13));
    }
}

TEST_CASE("oracle V1: chaine Ra-226 a 30 j, 1 a, 100 a", "[bateman][T2][verification][V1][data]") {
    const std::vector<OracleRow> oracle =
        read_oracle(verification_dir / "V1_analytic_bateman" / "oracle_ra226.csv");
    const decaysolver::NuclideLibrary library = decaysolver::NuclideLibrary::load(
        std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    const decaysolver::DecaySystem system = decaysolver::DecaySystem::build(library, {"Ra-226"});
    std::vector<double> n0(system.size(), 0.0);
    n0[system.index_of("Ra-226")] = 1.0;

    std::map<std::string, std::vector<double>> solutions;
    double worst_relative = 0.0;
    std::string worst_label;
    for (const OracleRow& row : oracle) {
        if (!solutions.contains(row.case_name)) {
            solutions[row.case_name] = decaysolver::solve_bateman(system, n0, row.t_s);
        }
        const double computed = solutions[row.case_name][system.index_of(row.nuclide)];
        INFO(row.case_name << " " << row.nuclide);
        if (row.value >= 1e-30) {
            REQUIRE_THAT(computed, WithinRel(row.value, 1e-11));
            const double relative = std::abs(computed / row.value - 1.0);
            if (relative > worst_relative) {
                worst_relative = relative;
                worst_label = row.case_name + " " + row.nuclide;
            }
        } else {
            REQUIRE_THAT(computed, WithinAbs(row.value, 1e-40));
        }
    }
    // Visible avec `decaysolver_unit_tests -s "[V1]"` ; repris dans le rapport de vérification.
    SUCCEED("écart relatif max vs oracle Ra-226 : " << worst_relative << " (" << worst_label
                                                    << ")");
}
