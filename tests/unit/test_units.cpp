// T1 — tests unitaires des conversions d'unités (CDC §4.2).
//
// Tolérances : une division IEEE-754 est correctement arrondie, donc ln2 / T ne peut différer
// de la valeur attendue que par l'arrondi final, soit au plus 1 ulp (epsilon = 2^-52 en relatif).
// Un aller-retour enchaîne deux divisions : 2 ulp. Les constantes de conversion sont des entiers
// exactement représentables : comparaison exacte via static_assert sur des entiers.

#include <decaysolver/units.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace units = decaysolver::units;

namespace {
constexpr double one_ulp = std::numeric_limits<double>::epsilon();
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
} // namespace

static_assert(static_cast<long long>(units::seconds_per_year) == 31'557'600LL,
              "année julienne : 365,25 j × 86 400 s");
static_assert(static_cast<long long>(units::seconds_per_day) == 86'400LL);

TEST_CASE("units: demi-vie vers constante de decroissance", "[units][T1]") {
    SECTION("T½ = 1 s donne λ = ln 2") {
        REQUIRE_THAT(units::decay_constant_per_s(1.0), WithinRel(std::numbers::ln2, one_ulp));
    }
    SECTION("T½ = 30,05 a (Cs-137) donne λ ≈ 7,31e-10 s⁻¹") {
        // Valeur de contrôle calculée indépendamment : ln2 / (30.05 × 31 557 600).
        const double half_life_s = 30.05 * units::seconds_per_year;
        REQUIRE_THAT(units::decay_constant_per_s(half_life_s),
                     WithinRel(7.309e-10, 1e-3)); // 4 chiffres significatifs de la référence
    }
    SECTION("nucléide stable : T½ = +inf donne exactement λ = 0") {
        REQUIRE_THAT(units::decay_constant_per_s(infinity), WithinAbs(0.0, 0.0));
    }
    SECTION("entrées invalides rejetées") {
        REQUIRE_THROWS_AS(units::decay_constant_per_s(0.0), std::domain_error);
        REQUIRE_THROWS_AS(units::decay_constant_per_s(-1.0), std::domain_error);
        REQUIRE_THROWS_AS(units::decay_constant_per_s(quiet_nan), std::domain_error);
    }
}

TEST_CASE("units: constante de decroissance vers demi-vie", "[units][T1]") {
    SECTION("λ = 0 donne T½ = +inf") {
        const double half_life = units::half_life_s_from_decay_constant(0.0);
        REQUIRE(std::isinf(half_life));
        REQUIRE(half_life > 0.0);
    }
    SECTION("aller-retour T½ -> λ -> T½ à 2 ulp") {
        for (const double half_life_s :
             {1.0e-3, 1.0, 312.2 * units::seconds_per_day, 4.468e9 * units::seconds_per_year}) {
            const double lambda = units::decay_constant_per_s(half_life_s);
            REQUIRE_THAT(units::half_life_s_from_decay_constant(lambda),
                         WithinRel(half_life_s, 2 * one_ulp));
        }
    }
    SECTION("entrées invalides rejetées") {
        REQUIRE_THROWS_AS(units::half_life_s_from_decay_constant(-1.0), std::domain_error);
        REQUIRE_THROWS_AS(units::half_life_s_from_decay_constant(quiet_nan), std::domain_error);
        REQUIRE_THROWS_AS(units::half_life_s_from_decay_constant(infinity), std::domain_error);
    }
}

TEST_CASE("units: conversion vers les secondes", "[units][T1]") {
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::microsecond), WithinRel(1e-6, one_ulp));
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::millisecond), WithinRel(1e-3, one_ulp));
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::second), WithinRel(1.0, one_ulp));
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::minute), WithinRel(60.0, one_ulp));
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::hour), WithinRel(3'600.0, one_ulp));
    REQUIRE_THAT(units::to_seconds(2.0, units::TimeUnit::day), WithinRel(172'800.0, one_ulp));
    REQUIRE_THAT(units::to_seconds(1.0, units::TimeUnit::year), WithinRel(31'557'600.0, one_ulp));
}
