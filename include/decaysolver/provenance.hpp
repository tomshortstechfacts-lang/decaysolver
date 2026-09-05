#pragma once

/// @file provenance.hpp
/// En-tête de provenance placé en tête de chaque fichier de sortie du solveur.
///
/// Un résultat numérique sans sa provenance (version du code, compilateur, données) n'est pas
/// reproductible, donc pas vérifiable. Chaque ligne commence par `#` pour rester lisible par
/// tout lecteur CSV.

#include <string>

namespace decaysolver {

/// Horodatage UTC courant au format ISO 8601, par exemple `2026-09-05T13:04:00Z`.
[[nodiscard]] std::string utc_timestamp_iso8601();

/// Lignes `# clé: valeur` décrivant la version du code (avec SHA git), le compilateur, le type
/// de build et l'horodatage. Les informations sur les données nucléaires (source, version,
/// SHA-256) seront ajoutées avec la bibliothèque de nucléides.
[[nodiscard]] std::string provenance_header();

} // namespace decaysolver
