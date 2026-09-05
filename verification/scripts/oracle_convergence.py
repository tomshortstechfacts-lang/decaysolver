"""Oracle du problème de convergence (CDC §4.3) : chaîne λ = (1, 2, 3, 0) s⁻¹, N(0) = (1, 0, 0, 0),
T = 4 s. Produit verification/V2_order_of_accuracy/oracle_lambda123.csv, lu par l'exécutable
decaysolver_convergence et par le test V1.

Le problème est volontairement non raide (h·λ_max ≤ 0,75 dès 16 pas) : on mesure l'ordre nominal
des schémas, pas la réduction d'ordre des problèmes raides.

Usage : python verification/scripts/oracle_convergence.py
"""

from pathlib import Path

from mpmath import mp, mpf

from oracle_common import Nuclide, solve, write_oracle

mp.dps = 50
T = mpf(4)


def main() -> None:
    library = {}
    lambdas = [mpf(1), mpf(2), mpf(3), mpf(0)]
    for i, lam in enumerate(lambdas):
        nuclide = Nuclide(f"A-{i + 1}")
        if lam != 0:
            nuclide.half_life_s = mp.log(2) / lam
            nuclide.branches = [(f"A-{i + 2}", mpf(1))]
        library[nuclide.name] = nuclide
    populations = solve(library, {"A-1": mpf(1)}, T)
    rows = [("lambda123", T, name, populations[name]) for name in ["A-1", "A-2", "A-3", "A-4"]]
    out = Path(__file__).resolve().parents[1] / "V2_order_of_accuracy" / "oracle_lambda123.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_oracle(out, rows)
    total = sum(populations.values())
    print(f"{out.name} écrit ; conservation : Σ N = {mp.nstr(total, 25)}")
    for name in ["A-1", "A-2", "A-3", "A-4"]:
        print(f"  {name}: {mp.nstr(populations[name], 20)}")


if __name__ == "__main__":
    main()
