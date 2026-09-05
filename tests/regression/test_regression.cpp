// T6 — non-régression : sorties figées du mode inventaire, rejouées à chaque commit.
//
// Cas : les 35 nucléides de la liste standard de déclaration des déchets, tous à 1 Bq (valeurs
// synthétiques, aucun inventaire réel), vieillis de 10 ans, dans les deux conventions de filles.
// Référence : tests/regression/inventory_reference.csv, produite par decaysolver lui-même à la
// version indiquée dans son en-tête, puis figée. Tolérance 1e-12 relatif : un changement de
// données (demi-vie, voie) ou de numérique au-delà de ça doit être voulu, documenté dans le
// CHANGELOG (rubrique Numerics) et accompagné d'une régénération explicite de la référence
// (`decaysolver_regression --regenerate`).

#include <decaysolver/inventory.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <charconv>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using Catch::Matchers::WithinRel;
using decaysolver::DaughterPolicy;

namespace {

const char* const base_nuclides[] = {
    "Be-10",   "C-14",    "Cl-36",  "Ca-41",  "Mn-54",  "Fe-55",  "Co-60",  "Ni-59",  "Ni-63",
    "Zn-65",   "Se-79",   "Sr-90",  "Zr-93",  "Nb-94",  "Mo-93",  "Tc-99",  "Pd-107", "Ag-108m",
    "Ag-110m", "Sn-121m", "Sn-126", "Sb-125", "I-129",  "Cs-134", "Cs-135", "Cs-137", "Sm-151",
    "U-235",   "U-238",   "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Am-241", "Cm-244",
};

decaysolver::Inventory synthetic_inventory() {
    decaysolver::Inventory inventory{decaysolver::ValueKind::activity_bq, {}};
    for (const char* name : base_nuclides) {
        inventory.entries.push_back({name, 1.0});
    }
    return inventory;
}

using Reference = std::map<std::string, std::map<std::string, double>>; // policy -> nuclide -> A

Reference read_reference(const std::filesystem::path& path) {
    std::ifstream file(path);
    REQUIRE(file.good());
    Reference reference;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line.front() == '#' || line.rfind("policy;", 0) == 0) {
            continue;
        }
        std::stringstream stream(line);
        std::string policy;
        std::string nuclide;
        std::string value_text;
        std::getline(stream, policy, ';');
        std::getline(stream, nuclide, ';');
        std::getline(stream, value_text, ';');
        double value = 0.0;
        std::from_chars(value_text.data(), value_text.data() + value_text.size(), value);
        reference[policy][nuclide] = value;
    }
    REQUIRE(reference.size() == 2);
    return reference;
}

} // namespace

TEST_CASE("regression T6: inventaire synthetique de 35 nucleides a 10 ans",
          "[regression][T6][data]") {
    const decaysolver::NuclideLibrary library = decaysolver::NuclideLibrary::load(
        std::filesystem::path(DECAYSOLVER_DATA_DIR) / "nuclides_icrp107.csv");
    const Reference reference = read_reference(std::filesystem::path(DECAYSOLVER_REGRESSION_DIR) /
                                               "inventory_reference.csv");
    const double age_s = 10.0 * decaysolver::units::seconds_per_year;

    for (const DaughterPolicy policy : {DaughterPolicy::input_only, DaughterPolicy::all}) {
        const std::string policy_name(decaysolver::to_string(policy));
        const decaysolver::AgedInventory aged =
            decaysolver::age_inventory(library, synthetic_inventory(), age_s, policy);
        const std::map<std::string, double>& expected = reference.at(policy_name);
        INFO(policy_name);
        REQUIRE(aged.entries.size() == expected.size());
        for (const decaysolver::AgedEntry& entry : aged.entries) {
            INFO(entry.nuclide);
            REQUIRE(expected.contains(entry.nuclide));
            REQUIRE_THAT(entry.activity, WithinRel(expected.at(entry.nuclide), 1e-12));
        }
    }
}
