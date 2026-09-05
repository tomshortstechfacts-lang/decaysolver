"""Oracle multiprécision des équations de Bateman, indépendant du code C++.

Deux différences volontaires avec decaysolver, pour que la comparaison ait un sens :
  - la formule fermée de Bateman (somme alternée) est évaluée telle quelle, sans différences
    divisées ni série de Taylor : à 50 chiffres, l'annulation catastrophique coûte au plus une
    dizaine de chiffres, il en reste 40 ;
  - l'arithmétique est celle de mpmath, pas le `double`.

Les graphes branchés sont traités par superposition des chemins (même principe physique que le
code, implémentation distincte). Les constantes de décroissance égales sur un même chemin ne sont
pas admises ici (le cas dégénéré exact D1 a sa propre formule dans oracle_degenerate.py).

Convention d'année : julienne, 365,25 j, comme decaysolver.
"""

from __future__ import annotations

import csv
from pathlib import Path

from mpmath import mp, mpf, exp, log

SECONDS_PER_UNIT = {
    "s": mpf(1),
    "m": mpf(60),
    "h": mpf(3600),
    "d": mpf(86400),
    "y": mpf(86400) * mpf("365.25"),
}


class Nuclide:
    def __init__(self, name: str) -> None:
        self.name = name
        self.half_life_s: mpf | None = None  # None : stable
        self.branches: list[tuple[str, mpf]] = []  # (fille, rapport) ; fille "" pour la fission

    @property
    def decay_constant(self) -> mpf:
        if self.half_life_s is None:
            return mpf(0)
        return log(2) / self.half_life_s


def load_library(path: Path) -> dict[str, Nuclide]:
    """Lit data/nuclides_icrp107.csv avec les mêmes conventions que decaysolver."""
    library: dict[str, Nuclide] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if not row or row[0].startswith("#") or row[0] == "nuclide":
                continue
            name, value, unit, mode, daughter, fraction = row
            nuclide = library.setdefault(name, Nuclide(name))
            if value != "stable":
                nuclide.half_life_s = mpf(value) * SECONDS_PER_UNIT[unit]
            if mode != "stable":
                nuclide.branches.append((daughter, mpf(fraction)))
    return library


def bateman_path(lambdas: list[mpf], t: mpf) -> mpf:
    """Σ_m exp(−λ_m t) / Π_{l≠m}(λ_l − λ_m), constantes deux à deux distinctes."""
    total = mpf(0)
    for m, lam_m in enumerate(lambdas):
        denominator = mpf(1)
        for l, lam_l in enumerate(lambdas):
            if l != m:
                if lam_l == lam_m:
                    raise ValueError("constantes égales sur un chemin : oracle non applicable")
                denominator *= lam_l - lam_m
        total += exp(-lam_m * t) / denominator
    return total


def solve(library: dict[str, Nuclide], n0: dict[str, mpf], t: mpf) -> dict[str, mpf]:
    """Populations à t par superposition des chemins issus de chaque nucléide de n0."""
    result: dict[str, mpf] = {}

    def explore(path: list[str], branching: mpf, lambda_product: mpf, seed_n0: mpf) -> None:
        last = library[path[-1]]
        lambdas = [library[name].decay_constant for name in path]
        contribution = seed_n0 * branching * lambda_product * bateman_path(lambdas, t)
        result[last.name] = result.get(last.name, mpf(0)) + contribution
        for daughter, fraction in last.branches:
            if daughter:
                explore(path + [daughter], branching * fraction,
                        lambda_product * last.decay_constant, seed_n0)

    for seed, population in n0.items():
        if population != 0:
            explore([seed], mpf(1), mpf(1), population)
    return result


def write_oracle(path: Path, rows: list[tuple[str, mpf, str, mpf]], digits: int = 20) -> None:
    """CSV `case;t_s;nuclide;N` à `digits` chiffres significatifs, précision mpmath en en-tête."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# oracle mpmath, mp.dps = {mp.dps}, valeurs à {digits} chiffres significatifs\n")
        handle.write("# généré par verification/scripts/ (voir le script nommé dans le cas)\n")
        handle.write("case;t_s;nuclide;N\n")
        for case, t, nuclide, value in rows:
            handle.write(f"{case};{mp.nstr(t, 17)};{nuclide};{mp.nstr(value, digits)}\n")
