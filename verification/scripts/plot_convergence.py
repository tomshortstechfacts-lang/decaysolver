"""Figures log-log et tableaux Markdown des ordres de convergence, à partir du CSV produit par
l'exécutable decaysolver_convergence.

Usage : python verification/scripts/plot_convergence.py CONVERGENCE.csv
Produit : verification/V2_order_of_accuracy/figures/convergence_{Linf,L2}.png
          verification/report/convergence_tables.md
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
LABELS = {
    "euler-explicit": "Euler explicite (p = 1)",
    "euler-implicit": "Euler implicite (p = 1)",
    "crank-nicolson": "Crank–Nicolson (p = 2)",
    "rk4": "RK4 (p = 4)",
}


def read(path: Path) -> tuple[dict[str, list[dict]], list[str]]:
    data: dict[str, list[dict]] = {}
    comments: list[str] = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if not row:
                continue
            if row[0].startswith("#"):
                comments.append(row[0])
                continue
            if row[0] == "scheme":
                continue
            scheme, order, k, n_steps, h, e_inf, e_2, p_inf, p_2 = row
            data.setdefault(scheme, []).append({
                "order": int(order), "k": int(k), "n": int(n_steps), "h": float(h),
                "E_inf": float(e_inf), "E_2": float(e_2),
                "p_inf": float(p_inf) if p_inf else None, "p_2": float(p_2) if p_2 else None,
            })
    return data, comments


def figure(data: dict[str, list[dict]], key: str, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for scheme, rows in data.items():
        h = [r["h"] for r in rows]
        e = [r[key] for r in rows]
        line, = ax.loglog(h, e, marker="o", markersize=4, label=LABELS.get(scheme, scheme))
        # droite de pente théorique, calée sur le point le plus grossier
        p = rows[0]["order"]
        h_ref = [h[0], h[-1]]
        ax.loglog(h_ref, [e[0] * (x / h[0]) ** p for x in h_ref], linestyle="--", linewidth=0.8,
                  color=line.get_color())
    ax.axhline(1e-12, color="grey", linewidth=0.6, linestyle=":")
    ax.text(ax.get_xlim()[0] * 1.5, 1.6e-12, "plancher d'arrondi retenu (1e-12)", fontsize=7, color="grey")
    ax.set_xlabel("pas de temps h (s)")
    ax.set_ylabel(title)
    ax.set_title("Ordres de convergence, λ = (1, 2, 3, 0) s⁻¹, T = 4 s")
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def tables(data: dict[str, list[dict]], comments: list[str], out: Path) -> None:
    lines = ["<!-- généré par verification/scripts/plot_convergence.py, ne pas éditer -->", ""]
    for scheme, rows in data.items():
        lines.append(f"### {LABELS.get(scheme, scheme)}")
        lines.append("")
        lines.append("| k | pas | h (s) | E_L∞ | E_L2 | p_obs (L∞) | p_obs (L2) |")
        lines.append("|---|---|---|---|---|---|---|")
        for r in rows:
            p_inf = f"{r['p_inf']:.3f}" if r["p_inf"] is not None else "—"
            p_2 = f"{r['p_2']:.3f}" if r["p_2"] is not None else "—"
            lines.append(f"| {r['k']} | {r['n']} | {r['h']:.3e} | {r['E_inf']:.3e} | {r['E_2']:.3e} | {p_inf} | {p_2} |")
        lines.append("")
    lines.append("Synthèse de l'exécutable :")
    lines.append("")
    lines.extend(f"    {c}" for c in comments if "suite valide" in c or "resultat global" in c)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    data, comments = read(Path(sys.argv[1]))
    fig_dir = ROOT / "verification" / "V2_order_of_accuracy" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    figure(data, "E_inf", "erreur L∞ à T", fig_dir / "convergence_Linf.png")
    figure(data, "E_2", "erreur L2 à T", fig_dir / "convergence_L2.png")
    report_dir = ROOT / "verification" / "report"
    report_dir.mkdir(parents=True, exist_ok=True)
    tables(data, comments, report_dir / "convergence_tables.md")
    print("figures et tableaux écrits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
