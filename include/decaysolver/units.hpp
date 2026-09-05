#pragma once

/// @file units.hpp
/// Unités de temps et conversion demi-vie <-> constante de décroissance.
///
/// Convention interne : toutes les durées sont en secondes, toutes les constantes de
/// décroissance en s⁻¹. Les noms portent l'unité (`half_life_s`, `lambda_per_s`) pour que
/// toute conversion implicite soit visible à la relecture.

namespace decaysolver::units {

inline constexpr double seconds_per_minute = 60.0;
inline constexpr double seconds_per_hour = 3'600.0;
inline constexpr double seconds_per_day = 86'400.0;

/// Année julienne : 365,25 j = 31 557 600 s (convention de l'UAI, employée par les
/// bibliothèques de données nucléaires évaluées, dont DDEP et ICRP-107).
/// L'écart avec l'année civile moyenne (365,2425 j) vaut 2e-5 en relatif, et avec 365 j
/// exactement 7e-4 : invisible à trois chiffres significatifs, mais à déclarer.
inline constexpr double seconds_per_year = 31'557'600.0;

enum class TimeUnit { second, minute, hour, day, year };

/// Convertit `value` exprimée dans `unit` en secondes.
[[nodiscard]] double to_seconds(double value, TimeUnit unit);

/// λ = ln 2 / T½.
/// @param half_life_s demi-vie en secondes, strictement positive, ou +∞ pour un nucléide
///        stable (retourne alors exactement 0).
/// @throws std::domain_error si la demi-vie est ≤ 0 ou NaN.
[[nodiscard]] double decay_constant_per_s(double half_life_s);

/// T½ = ln 2 / λ, inverse de decay_constant_per_s.
/// @param lambda_per_s constante de décroissance en s⁻¹, ≥ 0 et finie ; 0 donne +∞.
/// @throws std::domain_error si λ est < 0, NaN ou infinie.
[[nodiscard]] double half_life_s_from_decay_constant(double lambda_per_s);

} // namespace decaysolver::units
