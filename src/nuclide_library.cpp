#include <decaysolver/nuclide_library.hpp>
#include <decaysolver/units.hpp>

#include <charconv>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace decaysolver {

namespace {

constexpr std::string_view expected_header =
    "nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction";

std::string_view trim(std::string_view text) {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t' || text.front() == '\r')) {
        text.remove_prefix(1);
    }
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t' || text.back() == '\r')) {
        text.remove_suffix(1);
    }
    return text;
}

std::vector<std::string_view> split_fields(std::string_view line) {
    std::vector<std::string_view> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t end = line.find(';', start);
        if (end == std::string_view::npos) {
            fields.push_back(trim(line.substr(start)));
            return fields;
        }
        fields.push_back(trim(line.substr(start, end - start)));
        start = end + 1;
    }
}

// std::from_chars ne dépend pas de la locale : "30.1671" est lu de la même façon sur un poste
// configuré en français (où std::stod attendrait une virgule).
double parse_double(std::string_view text, std::size_t line_number, std::string_view what) {
    double value = 0.0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
        throw DataError("ligne " + std::to_string(line_number) + " : " + std::string(what) +
                        " illisible : '" + std::string(text) + "'");
    }
    return value;
}

units::TimeUnit parse_unit(std::string_view text, std::size_t line_number) {
    if (text == "s") {
        return units::TimeUnit::second;
    }
    if (text == "m") {
        return units::TimeUnit::minute;
    }
    if (text == "h") {
        return units::TimeUnit::hour;
    }
    if (text == "d") {
        return units::TimeUnit::day;
    }
    if (text == "y") {
        return units::TimeUnit::year;
    }
    throw DataError("ligne " + std::to_string(line_number) + " : unité de demi-vie inconnue : '" +
                    std::string(text) + "'");
}

} // namespace

NuclideLibrary NuclideLibrary::load(const std::filesystem::path& csv_path) {
    std::ifstream file(csv_path);
    if (!file) {
        throw DataError("impossible d'ouvrir le fichier de données : " + csv_path.string());
    }
    return parse(file);
}

NuclideLibrary NuclideLibrary::parse(std::istream& csv) {
    NuclideLibrary library;
    std::string raw_line;
    std::size_t line_number = 0;
    bool header_seen = false;

    while (std::getline(csv, raw_line)) {
        ++line_number;
        const std::string_view line = trim(raw_line);
        if (line.empty()) {
            continue;
        }
        if (line.front() == '#') {
            library.provenance_.emplace_back(line);
            continue;
        }
        if (!header_seen) {
            if (line != expected_header) {
                throw DataError("ligne " + std::to_string(line_number) +
                                " : en-tête inattendu, attendu '" + std::string(expected_header) +
                                "'");
            }
            header_seen = true;
            continue;
        }

        const std::vector<std::string_view> fields = split_fields(line);
        if (fields.size() != 6) {
            throw DataError("ligne " + std::to_string(line_number) + " : 6 champs attendus, " +
                            std::to_string(fields.size()) + " trouvés");
        }
        const std::string name = canonical_name(fields[0]);
        const DecayMode mode = decay_mode_from_string(fields[3]);

        // Première rencontre du nucléide : on crée l'entrée avec sa demi-vie.
        auto [it, inserted] = library.nuclides_.try_emplace(name);
        Nuclide& nuclide = it->second;
        if (inserted) {
            nuclide.name = name;
            if (fields[1] == "stable") {
                nuclide.half_life_s = std::numeric_limits<double>::infinity();
            } else {
                const double value = parse_double(fields[1], line_number, "demi-vie");
                nuclide.half_life_s = units::to_seconds(value, parse_unit(fields[2], line_number));
            }
        }

        if (mode == DecayMode::stable) {
            if (!nuclide.is_stable() || !nuclide.branches.empty()) {
                throw DataError("ligne " + std::to_string(line_number) + " : " + name +
                                " déclaré stable mais porte une demi-vie ou des voies");
            }
            continue;
        }
        if (nuclide.is_stable()) {
            throw DataError("ligne " + std::to_string(line_number) + " : " + name +
                            " stable ne peut pas avoir de voie de décroissance");
        }
        DecayBranch branch;
        branch.mode = mode;
        branch.daughter = fields[4].empty() ? std::string() : canonical_name(fields[4]);
        branch.branching_fraction = parse_double(fields[5], line_number, "rapport d'embranchement");
        nuclide.branches.push_back(branch);
    }

    if (!header_seen) {
        throw DataError("fichier de données vide ou sans ligne d'en-tête");
    }
    library.validate();
    return library;
}

void NuclideLibrary::validate() const {
    for (const auto& [name, nuclide] : nuclides_) {
        if (nuclide.is_stable()) {
            continue;
        }
        if (std::isnan(nuclide.half_life_s) || nuclide.half_life_s <= 0.0) {
            throw DataError(name + " : demi-vie non positive");
        }
        if (nuclide.branches.empty()) {
            throw DataError(name + " : nucléide radioactif sans voie de décroissance");
        }
        double sum = 0.0;
        for (const DecayBranch& branch : nuclide.branches) {
            if (!(branch.branching_fraction > 0.0) || branch.branching_fraction > 1.0) {
                throw DataError(name + " : rapport d'embranchement hors ]0 ; 1]");
            }
            sum += branch.branching_fraction;
            if (branch.mode == DecayMode::spontaneous_fission) {
                if (!branch.daughter.empty()) {
                    throw DataError(name + " : une voie de fission spontanée n'a pas de fille");
                }
            } else if (!nuclides_.contains(branch.daughter)) {
                throw DataError(name + " : fille absente de la bibliothèque : " + branch.daughter);
            }
        }
        if (std::abs(sum - 1.0) > branching_tolerance) {
            throw DataError(name + " : somme des rapports d'embranchement = " +
                            std::to_string(sum) + ", attendu 1");
        }
    }

    // Absence de cycle : parcours en profondeur depuis chaque nucléide. Un cycle de filiation
    // est physiquement impossible (la masse ou le numéro atomique change à chaque étape hors
    // transition isomérique, et une transition isomérique va toujours vers l'état fondamental),
    // mais un fichier de données peut en contenir un par erreur de saisie.
    std::set<std::string> finished;
    for (const auto& [start, unused] : nuclides_) {
        if (finished.contains(start)) {
            continue;
        }
        std::set<std::string> on_path;
        std::vector<std::pair<std::string, std::size_t>> stack{{start, 0}};
        on_path.insert(start);
        while (!stack.empty()) {
            auto& [current, next_branch] = stack.back();
            const std::vector<DecayBranch>& branches = nuclides_.at(current).branches;
            if (next_branch < branches.size()) {
                const std::string& daughter = branches[next_branch].daughter;
                ++next_branch;
                if (daughter.empty() || finished.contains(daughter)) {
                    continue;
                }
                if (on_path.contains(daughter)) {
                    throw DataError("cycle de filiation détecté passant par " + daughter);
                }
                on_path.insert(daughter);
                stack.emplace_back(daughter, 0);
            } else {
                on_path.erase(current);
                finished.insert(current);
                stack.pop_back();
            }
        }
    }
}

bool NuclideLibrary::contains(std::string_view name) const {
    return nuclides_.contains(canonical_name(name));
}

const Nuclide& NuclideLibrary::get(std::string_view name) const {
    const std::string key = canonical_name(name);
    const auto it = nuclides_.find(key);
    if (it == nuclides_.end()) {
        throw std::out_of_range("nucléide absent de la bibliothèque : " + key);
    }
    return it->second;
}

std::vector<std::string> NuclideLibrary::names() const {
    std::vector<std::string> out;
    out.reserve(nuclides_.size());
    for (const auto& [name, unused] : nuclides_) {
        out.push_back(name);
    }
    return out;
}

} // namespace decaysolver
