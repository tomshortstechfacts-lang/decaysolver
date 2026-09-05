#pragma once

/// @file nuclide_library.hpp
/// Bibliothèque de nucléides chargée depuis un fichier CSV, validée au chargement.
///
/// Format du fichier (voir data/PROVENANCE.md) : lignes d'en-tête commençant par `#`, puis
/// `nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction`, une ligne par voie
/// de décroissance. Un nucléide stable a une seule ligne `X;stable;;stable;;`.
///
/// Invariants garantis après chargement (sinon DataError) :
///  - demi-vie strictement positive, ou nucléide stable ;
///  - rapports d'embranchement dans ]0 ; 1] et de somme 1 à `branching_tolerance` près ;
///  - chaque fille nommée existe dans la bibliothèque (aucune chaîne tronquée en silence) ;
///  - aucun cycle de filiation.

#include <decaysolver/nuclide.hpp>

#include <filesystem>
#include <istream>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace decaysolver {

/// Erreur de données : fichier mal formé ou invariant physique violé. Le message cite la ligne.
class DataError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NuclideLibrary {
public:
    /// Tolérance sur |Σ b − 1|. ICRP-107 arrondit chaque rapport indépendamment, et arrondit à
    /// 1.0 une voie principale que des voies mineures complètent : l'écart observé sur le jeu
    /// extrait atteint 6e-5 (Fr-223 : β⁻ 1.0 + α 6e-5). 5e-4 laisse une marge d'un facteur 8.
    /// Une voie oubliée de rapport inférieur à cette tolérance est indétectable par la somme :
    /// c'est le script d'extraction, qui rapporte tout écart > 1e-5 dans l'en-tête du fichier,
    /// et la comparaison à une seconde source qui la détectent (voir data/PROVENANCE.md).
    static constexpr double branching_tolerance = 5e-4;

    /// Charge et valide un fichier CSV.
    [[nodiscard]] static NuclideLibrary load(const std::filesystem::path& csv_path);
    /// Même chose depuis un flux déjà ouvert (utile pour les tests).
    [[nodiscard]] static NuclideLibrary parse(std::istream& csv);

    /// `name` peut être une graphie usuelle : elle est canonisée avant la recherche.
    [[nodiscard]] bool contains(std::string_view name) const;
    /// @throws std::out_of_range si le nucléide est absent.
    [[nodiscard]] const Nuclide& get(std::string_view name) const;

    [[nodiscard]] std::size_t size() const { return nuclides_.size(); }
    /// Noms canoniques, triés.
    [[nodiscard]] std::vector<std::string> names() const;
    /// Lignes `#` d'en-tête du fichier, telles quelles : elles portent la provenance des données
    /// et sont recopiées dans l'en-tête des sorties.
    [[nodiscard]] const std::vector<std::string>& provenance() const { return provenance_; }

private:
    std::map<std::string, Nuclide> nuclides_;
    std::vector<std::string> provenance_;

    void validate() const;
};

} // namespace decaysolver
