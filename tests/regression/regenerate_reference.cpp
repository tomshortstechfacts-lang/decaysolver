// Régénère tests/regression/inventory_reference.csv avec la version courante du code et des
// données. À lancer uniquement quand un changement de résultat est voulu ; le CHANGELOG doit
// alors porter une entrée Numerics. L'en-tête de provenance identifie la version qui a produit
// la référence.

#include <decaysolver/inventory.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/provenance.hpp>
#include <decaysolver/units.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

const char* const base_nuclides[] = {
    "Be-10",   "C-14",    "Cl-36",  "Ca-41",  "Mn-54",  "Fe-55",  "Co-60",  "Ni-59",  "Ni-63",
    "Zn-65",   "Se-79",   "Sr-90",  "Zr-93",  "Nb-94",  "Mo-93",  "Tc-99",  "Pd-107", "Ag-108m",
    "Ag-110m", "Sn-121m", "Sn-126", "Sb-125", "I-129",  "Cs-134", "Cs-135", "Cs-137", "Sm-151",
    "U-235",   "U-238",   "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Am-241", "Cm-244",
};

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage : decaysolver_regenerate_reference SORTIE.csv\n";
        return 2;
    }
    try {
        const decaysolver::NuclideLibrary library =
            decaysolver::NuclideLibrary::load(DECAYSOLVER_DEFAULT_LIBRARY);
        decaysolver::Inventory inventory{decaysolver::ValueKind::activity_bq, {}};
        for (const char* name : base_nuclides) {
            inventory.entries.push_back({name, 1.0});
        }
        std::ofstream out(argv[1]);
        if (!out) {
            std::cerr << "impossible d'écrire " << argv[1] << '\n';
            return 1;
        }
        out << decaysolver::provenance_header();
        for (const std::string& line : library.provenance()) {
            out << "# data:" << line.substr(1) << '\n';
        }
        out << "# cas : 35 nucléides de la liste standard, 1 Bq chacun, vieillis de 10 a\n";
        out << "policy;nuclide;activity_bq\n" << std::scientific << std::setprecision(17);
        const double age_s = 10.0 * decaysolver::units::seconds_per_year;
        for (const decaysolver::DaughterPolicy policy :
             {decaysolver::DaughterPolicy::input_only, decaysolver::DaughterPolicy::all}) {
            const decaysolver::AgedInventory aged =
                decaysolver::age_inventory(library, inventory, age_s, policy);
            for (const decaysolver::AgedEntry& entry : aged.entries) {
                out << decaysolver::to_string(policy) << ';' << entry.nuclide << ';'
                    << entry.activity << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "erreur : " << error.what() << '\n';
        return 1;
    }
}
