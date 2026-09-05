#include <decaysolver/units.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace decaysolver::units {

double to_seconds(double value, TimeUnit unit) {
    switch (unit) {
    case TimeUnit::second:
        return value;
    case TimeUnit::minute:
        return value * seconds_per_minute;
    case TimeUnit::hour:
        return value * seconds_per_hour;
    case TimeUnit::day:
        return value * seconds_per_day;
    case TimeUnit::year:
        return value * seconds_per_year;
    }
    throw std::invalid_argument("unité de temps inconnue");
}

double decay_constant_per_s(double half_life_s) {
    if (std::isnan(half_life_s) || half_life_s <= 0.0) {
        throw std::domain_error("demi-vie invalide (attendu : > 0 s, ou +inf pour un nucléide "
                                "stable) : " +
                                std::to_string(half_life_s));
    }
    if (std::isinf(half_life_s)) {
        return 0.0; // nucléide stable
    }
    return std::numbers::ln2 / half_life_s;
}

double half_life_s_from_decay_constant(double lambda_per_s) {
    if (std::isnan(lambda_per_s) || std::isinf(lambda_per_s) || lambda_per_s < 0.0) {
        throw std::domain_error("constante de décroissance invalide (attendu : ≥ 0 et finie) : " +
                                std::to_string(lambda_per_s));
    }
    // Ici λ est fini et ≥ 0 ; le seul cas restant sans division est λ = 0 (stable).
    if (!(lambda_per_s > 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::numbers::ln2 / lambda_per_s;
}

} // namespace decaysolver::units
