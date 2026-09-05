#pragma once

// Assistant de test : chaîne linéaire A-1 → A-2 → … à constantes de décroissance choisies, passée
// par le vrai chemin (texte CSV → bibliothèque validée → système), pour que les tests numériques
// exercent le code de production et non un raccourci.

#include <decaysolver/decay_system.hpp>
#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <charconv>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

struct ChainFixture {
    // Le système référence les nucléides de la bibliothèque : elle doit lui survivre, d'où les
    // deux membres dans cet ordre.
    decaysolver::NuclideLibrary library;
    decaysolver::DecaySystem system;
};

// Représentation décimale la plus courte qui se relit exactement (std::to_chars), indépendante de
// la locale : la demi-vie écrite dans le CSV redonne λ à 1 ulp près.
inline std::string shortest(double value) {
    char buffer[32];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    return std::string(buffer, result.ptr);
}

// `lambdas_per_s` : constantes de A-1, A-2, … ; une constante nulle marque un nucléide stable et
// termine la chaîne. Si la dernière constante est non nulle, un puits stable est ajouté.
inline ChainFixture make_linear_chain(const std::vector<double>& lambdas_per_s) {
    std::ostringstream csv;
    csv << "nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction\n";
    std::size_t count = lambdas_per_s.size();
    for (std::size_t i = 0; i < lambdas_per_s.size(); ++i) {
        const std::string name = "A-" + std::to_string(i + 1);
        const std::string daughter = "A-" + std::to_string(i + 2);
        if (lambdas_per_s[i] <= 0.0) {
            csv << name << ";stable;;stable;;\n";
            count = i + 1;
            break;
        }
        const double half_life_s =
            decaysolver::units::half_life_s_from_decay_constant(lambdas_per_s[i]);
        csv << name << ';' << shortest(half_life_s) << ";s;beta-;" << daughter << ";1\n";
    }
    if (count == lambdas_per_s.size() && lambdas_per_s.back() > 0.0) {
        csv << "A-" << count + 1 << ";stable;;stable;;\n";
    }
    std::istringstream in(csv.str());
    ChainFixture fixture{decaysolver::NuclideLibrary::parse(in), decaysolver::DecaySystem{}};
    fixture.system = decaysolver::DecaySystem::build(fixture.library, {"A-1"});
    return fixture;
}
