"""Oracle du jeu R (CDC §4.3) : chaîne réelle Ra-226 → … → Pb-206, toutes voies ICRP-107 incluses,
N(Ra-226)(0) = 1, filles nulles, à t = 30 j, 1 a et 100 a. Produit
verification/V1_analytic_bateman/oracle_ra226.csv, lu par le test V1.

Les constantes s'étalent sur 14 ordres de grandeur (Po-214 : 4,2e3 s⁻¹ ; Pb-206 : 0) : à 50
chiffres, la formule fermée de Bateman reste exacte, ce qui est précisément le contre-exemple du
cahier des charges (raideur ≠ annulation).

Usage : python verification/scripts/oracle_ra226.py
"""

from pathlib import Path

from mpmath import mp, mpf

from oracle_common import load_library, solve, write_oracle

mp.dps = 50
DAY = mpf(86400)
YEAR = mpf("365.25") * DAY
TIMES = {"30j": 30 * DAY, "1a": YEAR, "100a": 100 * YEAR}


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    library = load_library(root / "data" / "nuclides_icrp107.csv")
    rows = []
    for label, t in TIMES.items():
        populations = solve(library, {"Ra-226": mpf(1)}, t)
        total = sum(populations.values())
        print(f"t = {label:>4} : {len(populations)} nucléides, Σ N = {mp.nstr(total, 20)}")
        for name in sorted(populations, key=lambda n: -populations[n]):
            rows.append((f"ra226_{label}", t, name, populations[name]))
    out = root / "verification" / "V1_analytic_bateman" / "oracle_ra226.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_oracle(out, rows)
    print(f"{out.relative_to(root)} écrit ({len(rows)} lignes)")


if __name__ == "__main__":
    main()
