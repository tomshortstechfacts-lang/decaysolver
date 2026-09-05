#include <decaysolver/bateman.hpp>
#include <decaysolver/decay_system.hpp>
#include <decaysolver/inventory.hpp>
#include <decaysolver/provenance.hpp>
#include <decaysolver/units.hpp>

#include <charconv>
#include <cmath>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace decaysolver {

namespace {

std::string_view trim(std::string_view text) {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t' || text.front() == '\r')) {
        text.remove_prefix(1);
    }
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t' || text.back() == '\r')) {
        text.remove_suffix(1);
    }
    return text;
}

// Lecture d'un nombre avec virgule ou point décimal, indépendante de la locale.
bool parse_number(std::string_view text, double& value) {
    std::string normalized(trim(text));
    for (char& c : normalized) {
        if (c == ',') {
            c = '.';
        }
    }
    const auto result =
        std::from_chars(normalized.data(), normalized.data() + normalized.size(), value);
    return result.ec == std::errc{} && result.ptr == normalized.data() + normalized.size();
}

// Coupe une ligne sur le premier séparateur reconnu : `;`, tabulation, ou `,` si la ligne n'a pas
// de `;` (la virgule peut aussi être décimale : "Cs137;1,43E-07" se coupe sur le `;`).
bool split_line(std::string_view line, std::string_view& left, std::string_view& right) {
    std::size_t cut = line.find(';');
    if (cut == std::string_view::npos) {
        cut = line.find('\t');
    }
    if (cut == std::string_view::npos) {
        cut = line.find(',');
    }
    if (cut == std::string_view::npos) {
        return false;
    }
    left = trim(line.substr(0, cut));
    right = trim(line.substr(cut + 1));
    return true;
}

} // namespace

DaughterPolicy daughter_policy_from_string(std::string_view text) {
    if (text == "input-only") {
        return DaughterPolicy::input_only;
    }
    if (text == "all") {
        return DaughterPolicy::all;
    }
    throw std::invalid_argument("convention des filles inconnue : '" + std::string(text) +
                                "' (attendu : input-only, all)");
}

std::string_view to_string(DaughterPolicy policy) {
    return policy == DaughterPolicy::input_only ? "input-only" : "all";
}

ValueKind value_kind_from_string(std::string_view text) {
    if (text == "bq") {
        return ValueKind::activity_bq;
    }
    if (text == "fraction") {
        return ValueKind::fraction;
    }
    throw std::invalid_argument("nature des valeurs inconnue : '" + std::string(text) +
                                "' (attendu : bq, fraction)");
}

std::string_view to_string(ValueKind kind) {
    return kind == ValueKind::activity_bq ? "bq" : "fraction";
}

Inventory read_inventory(std::istream& in, ValueKind kind) {
    Inventory inventory{kind, {}};
    std::set<std::string> seen;
    std::string raw;
    std::size_t line_number = 0;
    bool header_allowed = true;
    while (std::getline(in, raw)) {
        ++line_number;
        const std::string_view line = trim(raw);
        if (line.empty() || line.front() == '#') {
            continue;
        }
        std::string_view name_text;
        std::string_view value_text;
        if (!split_line(line, name_text, value_text)) {
            throw DataError("ligne " + std::to_string(line_number) +
                            " : attendu 'nuclide;valeur', lu '" + std::string(line) + "'");
        }
        double value = 0.0;
        if (!parse_number(value_text, value)) {
            if (header_allowed) { // première ligne non numérique : en-tête textuel toléré
                header_allowed = false;
                continue;
            }
            throw DataError("ligne " + std::to_string(line_number) + " : valeur illisible '" +
                            std::string(value_text) + "'");
        }
        header_allowed = false;
        if (std::isnan(value) || value < 0.0) {
            throw DataError("ligne " + std::to_string(line_number) + " : valeur négative ou NaN");
        }
        std::string name;
        try {
            name = canonical_name(name_text);
        } catch (const std::invalid_argument& error) {
            throw DataError("ligne " + std::to_string(line_number) + " : " + error.what());
        }
        if (!seen.insert(name).second) {
            throw DataError("ligne " + std::to_string(line_number) +
                            " : nucléide en double : " + name);
        }
        inventory.entries.push_back({name, value});
    }
    if (inventory.entries.empty()) {
        throw DataError("inventaire vide");
    }
    return inventory;
}

double parse_duration_s(std::string_view text) {
    std::string_view rest = trim(text);
    if (!rest.empty() && rest.front() == '+') {
        rest.remove_prefix(1);
    }
    double value = 0.0;
    const auto result = std::from_chars(rest.data(), rest.data() + rest.size(), value);
    if (result.ec != std::errc{} || result.ptr == rest.data()) {
        throw std::invalid_argument("durée illisible : '" + std::string(text) + "'");
    }
    const std::string_view unit = trim(std::string_view(
        result.ptr, static_cast<std::size_t>(rest.data() + rest.size() - result.ptr)));
    if (std::isnan(value) || value < 0.0) {
        throw std::invalid_argument("durée négative : '" + std::string(text) + "'");
    }
    if (unit == "s") {
        return units::to_seconds(value, units::TimeUnit::second);
    }
    if (unit == "min") {
        return units::to_seconds(value, units::TimeUnit::minute);
    }
    if (unit == "h") {
        return units::to_seconds(value, units::TimeUnit::hour);
    }
    if (unit == "j" || unit == "d") {
        return units::to_seconds(value, units::TimeUnit::day);
    }
    if (unit == "a" || unit == "y") {
        return units::to_seconds(value, units::TimeUnit::year);
    }
    throw std::invalid_argument("unité de durée inconnue : '" + std::string(unit) +
                                "' (attendu : s, min, h, j, a)");
}

AgedInventory age_inventory(const NuclideLibrary& library, const Inventory& inventory, double age_s,
                            DaughterPolicy policy) {
    if (std::isnan(age_s) || age_s < 0.0) {
        throw std::invalid_argument("âge négatif ou NaN");
    }
    std::vector<std::string> seeds;
    seeds.reserve(inventory.entries.size());
    for (const InventoryEntry& entry : inventory.entries) {
        const Nuclide& nuclide = library.get(entry.nuclide); // out_of_range si absent
        if (nuclide.is_stable()) {
            throw DataError(entry.nuclide + " : nucléide stable, activité indéfinie");
        }
        seeds.push_back(entry.nuclide);
    }
    const DecaySystem system = DecaySystem::build(library, seeds);
    const std::vector<double>& lambdas = system.decay_constants_per_s();

    // Activités → populations, évolution, populations → activités.
    std::vector<double> n0(system.size(), 0.0);
    for (const InventoryEntry& entry : inventory.entries) {
        const std::size_t i = system.index_of(entry.nuclide);
        n0[i] = entry.value / lambdas[i];
    }
    const std::vector<double> n = solve_bateman(system, n0, age_s);

    // Sélection et ordre de sortie : ordre d'entrée pour input-only, ordre topologique pour all.
    std::vector<std::size_t> selected;
    if (policy == DaughterPolicy::input_only) {
        for (const InventoryEntry& entry : inventory.entries) {
            selected.push_back(system.index_of(entry.nuclide));
        }
    } else {
        for (std::size_t i = 0; i < system.size(); ++i) {
            if (!system.nuclide(i).is_stable()) {
                selected.push_back(i);
            }
        }
    }

    AgedInventory aged{age_s, policy, inventory.kind, {}, 0.0, 0.0, 0.0};
    for (const std::size_t i : selected) {
        const double activity = lambdas[i] * n[i];
        const DecayMode mode = system.nuclide(i).primary_mode();
        aged.entries.push_back({system.names()[i], activity, 0.0, mode});
        aged.total += activity;
        if (mode == DecayMode::alpha) {
            aged.alpha += activity;
        } else {
            aged.beta_gamma += activity;
        }
    }
    for (AgedEntry& entry : aged.entries) {
        entry.fraction = aged.total > 0.0 ? entry.activity / aged.total : 0.0;
    }
    return aged;
}

void write_aged_inventory(std::ostream& out, const AgedInventory& aged,
                          const NuclideLibrary& library, std::string_view input_description) {
    out << provenance_header();
    for (const std::string& line : library.provenance()) {
        out << "# data:" << line.substr(1) << '\n';
    }
    out << "# input: " << input_description << '\n';
    out << "# input_kind: " << to_string(aged.kind) << '\n';
    out << "# age_s: " << std::setprecision(17) << aged.age_s
        << " (année julienne = 31 557 600 s)\n";
    out << "# daughters: " << to_string(aged.policy) << '\n';
    out << "# total_activity: " << std::setprecision(16) << aged.total << '\n';
    out << "# alpha_activity: " << aged.alpha << '\n';
    out << "# beta_gamma_activity: " << aged.beta_gamma << '\n';
    if (aged.alpha > 0.0) {
        out << "# beta_gamma_over_alpha: " << aged.beta_gamma / aged.alpha << '\n';
    }
    out << "nuclide;activity;fraction;primary_mode\n";
    // 17 chiffres significatifs : un double se relit exactement (aller-retour sans perte).
    out << std::scientific << std::setprecision(16);
    for (const AgedEntry& entry : aged.entries) {
        out << entry.nuclide << ';' << entry.activity << ';' << entry.fraction << ';'
            << to_string(entry.primary_mode) << '\n';
    }
}

} // namespace decaysolver
