#include <decaysolver/decay_system.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace decaysolver {

DecaySystem DecaySystem::build(const NuclideLibrary& library,
                               const std::vector<std::string>& seeds) {
    // 1. Fermeture par filiation : parcours en profondeur depuis chaque nucléide de départ.
    std::set<std::string> visited;
    std::vector<std::string> stack;
    for (const std::string& seed : seeds) {
        stack.push_back(library.get(seed).name); // canonise et vérifie la présence
    }
    while (!stack.empty()) {
        const std::string current = stack.back();
        stack.pop_back();
        if (!visited.insert(current).second) {
            continue;
        }
        for (const DecayBranch& branch : library.get(current).branches) {
            if (!branch.daughter.empty()) {
                stack.push_back(branch.daughter);
            }
        }
    }

    // 2. Ordre topologique par parcours en profondeur : un nucléide est émis après toutes ses
    //    filles (ordre postfixe), puis la liste est inversée. La bibliothèque garantit déjà
    //    l'absence de cycle.
    std::vector<std::string> postorder;
    std::set<std::string> done;
    for (const std::string& start : visited) {
        if (done.contains(start)) {
            continue;
        }
        std::vector<std::pair<std::string, std::size_t>> path{{start, 0}};
        while (!path.empty()) {
            auto& [name, next] = path.back();
            const std::vector<DecayBranch>& branches = library.get(name).branches;
            if (next < branches.size()) {
                const std::string& daughter = branches[next].daughter;
                ++next;
                if (!daughter.empty() && !done.contains(daughter)) {
                    path.emplace_back(daughter, 0);
                }
            } else {
                if (done.insert(name).second) {
                    postorder.push_back(name);
                }
                path.pop_back();
            }
        }
    }
    std::reverse(postorder.begin(), postorder.end());

    // 3. Tables indexées.
    DecaySystem system;
    system.names_ = postorder;
    std::map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < system.names_.size(); ++i) {
        index[system.names_[i]] = i;
    }
    system.nuclides_.reserve(system.size());
    system.lambdas_.reserve(system.size());
    system.productions_.assign(system.size(), {});
    system.daughters_.assign(system.size(), {});
    for (std::size_t j = 0; j < system.size(); ++j) {
        const Nuclide& nuclide = library.get(system.names_[j]);
        system.nuclides_.push_back(&nuclide);
        system.lambdas_.push_back(nuclide.decay_constant_per_s());
        for (const DecayBranch& branch : nuclide.branches) {
            if (branch.daughter.empty()) {
                continue; // fission spontanée : les noyaux quittent le système
            }
            const std::size_t i = index.at(branch.daughter);
            system.productions_[i].push_back({j, branch.branching_fraction});
            system.daughters_[j].push_back({i, branch.branching_fraction});
        }
    }
    return system;
}

std::size_t DecaySystem::index_of(std::string_view name) const {
    const std::string key = canonical_name(name);
    const auto it = std::find(names_.begin(), names_.end(), key);
    if (it == names_.end()) {
        throw std::out_of_range("nucléide absent du système : " + key);
    }
    return static_cast<std::size_t>(it - names_.begin());
}

bool DecaySystem::contains(std::string_view name) const {
    return std::find(names_.begin(), names_.end(), canonical_name(name)) != names_.end();
}

std::vector<double> DecaySystem::rate(const std::vector<double>& populations) const {
    if (populations.size() != size()) {
        throw std::invalid_argument("rate : vecteur de taille " +
                                    std::to_string(populations.size()) + ", attendu " +
                                    std::to_string(size()));
    }
    std::vector<double> out(size());
    for (std::size_t i = 0; i < size(); ++i) {
        double value = -lambdas_[i] * populations[i];
        for (const Production& production : productions_[i]) {
            value += production.branching_fraction * lambdas_[production.parent] *
                     populations[production.parent];
        }
        out[i] = value;
    }
    return out;
}

} // namespace decaysolver
