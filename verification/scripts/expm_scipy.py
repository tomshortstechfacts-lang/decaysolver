"""Troisième voie indépendante : exponentielle de matrice N(t) = expm(A t) N0 avec SciPy
(Padé + scaling-and-squaring, Al-Mohy & Higham 2009), en double.

On la compare à l'oracle mpmath sur les mêmes cas que le solveur (chaîne (1,2,3,0) et chaîne du
Ra-226 à 30 j, 1 a, 100 a). Comme le solveur, elle travaille en double ; contrairement à lui, elle
ne connaît ni la formule de Bateman ni les différences divisées : c'est une méthode générale
d'algèbre linéaire. L'accord des trois (oracle, solveur, expm) sur une chaîne raide à 14 ordres de
grandeur est la preuve croisée demandée par le cahier des charges.

Note : sur un système linéaire à coefficients constants, expm est exacte ; son erreur ne dépend
pas d'un pas de temps et ne fait donc pas l'objet d'une mesure d'ordre.

Usage : python verification/scripts/expm_scipy.py
Produit : verification/report/expm_comparison.md
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from oracle_common import load_library

ROOT = Path(__file__).resolve().parents[2]


def read_oracle(path: Path) -> dict[str, tuple[float, dict[str, float]]]:
    cases: dict[str, tuple[float, dict[str, float]]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if not row or row[0].startswith("#") or row[0] == "case":
                continue
            case, t_s, nuclide, value = row
            cases.setdefault(case, (float(t_s), {}))[1][nuclide] = float(value)
    return cases


def build_matrix(library, names: list[str]) -> np.ndarray:
    """A_ii = -λ_i, A_ij = b_ji λ_j (parents j avant filles i, ordre quelconque ici)."""
    index = {name: i for i, name in enumerate(names)}
    a = np.zeros((len(names), len(names)))
    for j, name in enumerate(names):
        nuclide = library[name]
        lam = float(nuclide.decay_constant)
        a[j, j] = -lam
        for daughter, fraction in nuclide.branches:
            if daughter:
                a[index[daughter], j] += float(fraction) * lam
    return a


def closure(library, seed: str) -> list[str]:
    seen, stack = [], [seed]
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.append(name)
        stack.extend(d for d, _ in library[name].branches if d)
    return seen


def compare(library, seed: str, cases, lines: list[str]) -> float:
    names = closure(library, seed)
    a = build_matrix(library, names)
    n0 = np.zeros(len(names))
    n0[names.index(seed)] = 1.0
    worst = 0.0
    for case, (t_s, oracle) in cases.items():
        n = expm(a * t_s) @ n0
        lines.append(f"\n### {case} (t = {t_s:.6g} s)\n")
        lines.append("| nucléide | oracle mpmath | expm SciPy | écart relatif |")
        lines.append("|---|---|---|---|")
        for name in names:
            ref = oracle.get(name, 0.0)
            val = n[names.index(name)]
            if ref >= 1e-30:
                rel = abs(val / ref - 1.0)
                worst = max(worst, rel)
                lines.append(f"| {name} | {ref:.12e} | {val:.12e} | {rel:.1e} |")
            else:
                lines.append(f"| {name} | {ref:.3e} | {val:.3e} | (< 1e-30, absolu {abs(val - ref):.1e}) |")
    return worst


def main() -> int:
    library = load_library(ROOT / "data" / "nuclides_icrp107.csv")
    lines = ["<!-- généré par verification/scripts/expm_scipy.py -->",
             "# Exponentielle de matrice (SciPy) vs oracle multiprécision", "",
             "Méthode : `scipy.linalg.expm` (Padé 13 + scaling-and-squaring), `double`. Même bibliothèque",
             "de données et même convention d'année que decaysolver (chargeur Python indépendant).", ""]

    # chaîne (1,2,3,0) : bibliothèque synthétique
    from mpmath import mp, mpf
    from oracle_common import Nuclide
    synthetic = {}
    for i, lam in enumerate([1.0, 2.0, 3.0, 0.0]):
        nuc = Nuclide(f"A-{i + 1}")
        if lam:
            nuc.half_life_s = mp.log(2) / mpf(lam)
            nuc.branches = [(f"A-{i + 2}", mpf(1))]
        synthetic[nuc.name] = nuc
    worst_123 = compare(synthetic, "A-1",
                        read_oracle(ROOT / "verification" / "V2_order_of_accuracy" / "oracle_lambda123.csv"), lines)
    worst_ra = compare(library, "Ra-226",
                       read_oracle(ROOT / "verification" / "V1_analytic_bateman" / "oracle_ra226.csv"), lines)

    lines.insert(5, f"**Écart relatif max** : chaîne (1,2,3,0) : {worst_123:.1e} ; chaîne du Ra-226 (15 nucléides, 3 échéances) : {worst_ra:.1e}.")
    out = ROOT / "verification" / "report" / "expm_comparison.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"écart max (1,2,3,0) : {worst_123:.2e} ; Ra-226 : {worst_ra:.2e} -> {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
