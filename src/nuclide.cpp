#include <decaysolver/nuclide.hpp>
#include <decaysolver/units.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>

namespace decaysolver {

namespace {

struct ModeName {
    DecayMode mode;
    std::string_view text;
};

// Vocabulaire unique, partagé par la lecture et l'écriture.
constexpr ModeName mode_names[] = {
    {DecayMode::alpha, "alpha"},
    {DecayMode::beta_minus, "beta-"},
    {DecayMode::beta_plus_ec, "beta+/EC"},
    {DecayMode::isomeric_transition, "IT"},
    {DecayMode::spontaneous_fission, "SF"},
    {DecayMode::stable, "stable"},
};

bool is_letter(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) != 0;
}

bool is_digit(char c) {
    return std::isdigit(static_cast<unsigned char>(c)) != 0;
}

// Symbole chimique en casse canonique : première lettre majuscule, seconde minuscule.
std::string canonical_symbol(std::string_view symbol) {
    std::string out(symbol);
    out[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[0])));
    if (out.size() == 2) {
        out[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(out[1])));
    }
    return out;
}

} // namespace

DecayMode decay_mode_from_string(std::string_view text) {
    for (const ModeName& entry : mode_names) {
        if (entry.text == text) {
            return entry.mode;
        }
    }
    throw std::invalid_argument("mode de décroissance inconnu : '" + std::string(text) + "'");
}

std::string_view to_string(DecayMode mode) {
    for (const ModeName& entry : mode_names) {
        if (entry.mode == mode) {
            return entry.text;
        }
    }
    throw std::invalid_argument("mode de décroissance hors énumération");
}

bool Nuclide::is_stable() const {
    return std::isinf(half_life_s);
}

double Nuclide::decay_constant_per_s() const {
    return units::decay_constant_per_s(half_life_s);
}

DecayMode Nuclide::primary_mode() const {
    if (branches.empty()) {
        return DecayMode::stable;
    }
    const auto largest = std::max_element(branches.begin(), branches.end(),
                                          [](const DecayBranch& a, const DecayBranch& b) {
                                              return a.branching_fraction < b.branching_fraction;
                                          });
    return largest->mode;
}

std::string canonical_name(std::string_view raw) {
    const auto fail = [&raw] {
        return std::invalid_argument("nom de nucléide non reconnu : '" + std::string(raw) + "'");
    };
    // 1. On retire espaces et soulignés, et au plus un tiret : "Cs-137", "Cs 137" et "Cs137"
    //    deviennent identiques ; "Cs-13-7" est rejeté.
    std::string compact;
    int dashes = 0;
    for (const char c : raw) {
        if (c == '-') {
            ++dashes;
        } else if (c != ' ' && c != '_') {
            compact.push_back(c);
        }
    }
    if (compact.empty() || dashes > 1) {
        throw fail();
    }

    // 2. Deux formes : symbole puis masse (+ isomère), ou masse puis symbole.
    std::string_view symbol;
    std::string_view mass_text;
    std::string_view isomer;
    if (is_letter(compact.front())) {
        std::size_t i = 0;
        while (i < compact.size() && is_letter(compact[i])) {
            ++i;
        }
        std::size_t j = i;
        while (j < compact.size() && is_digit(compact[j])) {
            ++j;
        }
        symbol = std::string_view(compact).substr(0, i);
        mass_text = std::string_view(compact).substr(i, j - i);
        isomer = std::string_view(compact).substr(j);
    } else {
        std::size_t i = 0;
        while (i < compact.size() && is_digit(compact[i])) {
            ++i;
        }
        mass_text = std::string_view(compact).substr(0, i);
        symbol = std::string_view(compact).substr(i);
    }

    // 3. Contrôles : symbole de 1 ou 2 lettres, masse de 1 à 3 chiffres, isomère "m" ou "m<n>".
    if (symbol.empty() || symbol.size() > 2 ||
        !std::all_of(symbol.begin(), symbol.end(), is_letter)) {
        throw fail();
    }
    if (mass_text.empty() || mass_text.size() > 3) {
        throw fail();
    }
    int mass = 0;
    std::from_chars(mass_text.data(), mass_text.data() + mass_text.size(), mass);
    if (mass < 1 || mass > 300) {
        throw fail();
    }
    if (!isomer.empty()) {
        const bool well_formed =
            (isomer[0] == 'm' || isomer[0] == 'M') &&
            (isomer.size() == 1 || (isomer.size() == 2 && is_digit(isomer[1])));
        if (!well_formed) {
            throw fail();
        }
    }

    std::string out = canonical_symbol(symbol) + "-" + std::string(mass_text);
    if (!isomer.empty()) {
        out += 'm';
        out += isomer.substr(1);
    }
    return out;
}

} // namespace decaysolver
