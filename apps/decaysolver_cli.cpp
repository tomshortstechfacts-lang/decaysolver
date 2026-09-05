// Exécutable en ligne de commande. Séparé de la bibliothèque : il ne contient aucune numérique,
// seulement l'analyse des arguments et l'affichage.
//
// Lot 0 : expose la version et l'en-tête de provenance. Les modes « chaîne » et « inventaire »
// arrivent avec le lot 1 (voir README, feuille de route).

#include <decaysolver/provenance.hpp>
#include <decaysolver/version.hpp>

#include <iostream>
#include <span>
#include <string_view>

namespace {

void print_usage(std::ostream& out) {
    out << "usage : decaysolver [--version | --provenance]\n";
}

} // namespace

int main(int argc, char** argv) {
    const std::span<char*> args(argv, static_cast<std::size_t>(argc));
    if (args.size() > 2) {
        std::cerr << "une seule option attendue\n";
        print_usage(std::cerr);
        return 2;
    }
    const std::string_view arg = args.size() == 2 ? std::string_view{args[1]} : std::string_view{};
    if (arg.empty() || arg == "--provenance") {
        std::cout << decaysolver::provenance_header();
        return 0;
    }
    if (arg == "--version") {
        std::cout << decaysolver::version::string << '\n';
        return 0;
    }
    if (arg == "--help" || arg == "-h") {
        print_usage(std::cout);
        return 0;
    }
    std::cerr << "option inconnue : " << arg << '\n';
    print_usage(std::cerr);
    return 2;
}
