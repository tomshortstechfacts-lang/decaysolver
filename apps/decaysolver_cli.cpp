// Exécutable en ligne de commande. Séparé de la bibliothèque : il ne contient aucune numérique,
// seulement l'analyse des arguments, la lecture des fichiers et l'affichage.
//
//   decaysolver --version | --provenance
//   decaysolver age --input INV.csv --age 6a [--kind bq|fraction] [--daughters input-only|all]
//                   [--library data/nuclides_icrp107.csv] [--output OUT.csv]

#include <decaysolver/inventory.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/provenance.hpp>
#include <decaysolver/version.hpp>

#include <exception>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

namespace {

void print_usage(std::ostream& out) {
    out << "usage :\n"
           "  decaysolver --version | --provenance | --help\n"
           "  decaysolver age --input FICHIER --age DUREE [options]\n"
           "\n"
           "Vieillit un inventaire (lignes 'nuclide;valeur') de DUREE (ex. 6a, 30j, 12h, 90min, "
           "0s).\n"
           "options :\n"
           "  --kind bq|fraction         nature des valeurs d'entrée (défaut : bq)\n"
           "  --daughters input-only|all filles hors liste masquées (défaut) ou listées\n"
           "  --library FICHIER          bibliothèque de nucléides (défaut : "
        << DECAYSOLVER_DEFAULT_LIBRARY
        << ")\n"
           "  --output FICHIER           sortie CSV (défaut : sortie standard)\n";
}

struct AgeArguments {
    std::string input;
    std::string age;
    std::string kind = "bq";
    std::string daughters = "input-only";
    std::string library = DECAYSOLVER_DEFAULT_LIBRARY;
    std::string output;
};

AgeArguments parse_age_arguments(std::span<char*> args) {
    AgeArguments parsed;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string_view option{args[i]};
        if (i + 1 >= args.size()) {
            throw std::invalid_argument("option sans valeur : " + std::string(option));
        }
        const std::string value{args[i + 1]};
        ++i;
        if (option == "--input") {
            parsed.input = value;
        } else if (option == "--age") {
            parsed.age = value;
        } else if (option == "--kind") {
            parsed.kind = value;
        } else if (option == "--daughters") {
            parsed.daughters = value;
        } else if (option == "--library") {
            parsed.library = value;
        } else if (option == "--output") {
            parsed.output = value;
        } else {
            throw std::invalid_argument("option inconnue : " + std::string(option));
        }
    }
    if (parsed.input.empty() || parsed.age.empty()) {
        throw std::invalid_argument("--input et --age sont obligatoires");
    }
    return parsed;
}

int run_age(std::span<char*> args) {
    const AgeArguments arguments = parse_age_arguments(args);
    const decaysolver::ValueKind kind = decaysolver::value_kind_from_string(arguments.kind);
    const decaysolver::DaughterPolicy policy =
        decaysolver::daughter_policy_from_string(arguments.daughters);
    const double age_s = decaysolver::parse_duration_s(arguments.age);

    const decaysolver::NuclideLibrary library =
        decaysolver::NuclideLibrary::load(arguments.library);
    std::ifstream input_file(arguments.input);
    if (!input_file) {
        throw std::runtime_error("impossible d'ouvrir l'inventaire : " + arguments.input);
    }
    const decaysolver::Inventory inventory = decaysolver::read_inventory(input_file, kind);
    const decaysolver::AgedInventory aged =
        decaysolver::age_inventory(library, inventory, age_s, policy);

    const std::string description = arguments.input + " vieilli de " + arguments.age;
    if (arguments.output.empty()) {
        decaysolver::write_aged_inventory(std::cout, aged, library, description);
        return 0;
    }
    std::ofstream output_file(arguments.output);
    if (!output_file) {
        throw std::runtime_error("impossible d'écrire : " + arguments.output);
    }
    decaysolver::write_aged_inventory(output_file, aged, library, description);
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    const std::span<char*> args(argv, static_cast<std::size_t>(argc));
    if (args.size() < 2) {
        print_usage(std::cerr);
        return 2;
    }
    const std::string_view command{args[1]};
    try {
        if (command == "--version") {
            std::cout << decaysolver::version::string << '\n';
            return 0;
        }
        if (command == "--provenance") {
            std::cout << decaysolver::provenance_header();
            return 0;
        }
        if (command == "--help" || command == "-h") {
            print_usage(std::cout);
            return 0;
        }
        if (command == "age") {
            return run_age(args.subspan(2));
        }
        std::cerr << "commande inconnue : " << command << '\n';
        print_usage(std::cerr);
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "erreur : " << error.what() << '\n';
        return 1;
    }
}
