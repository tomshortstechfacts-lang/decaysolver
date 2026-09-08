"""Contrôle indépendant de la méthode CRAM (lot 5), côté Python.

Deux choses, sans rien partager avec src/cram.cpp hormis les coefficients publiés :

1. **Scalaire** : r_k(x) = α₀ Π (1 + 2 Re[α_i/(x − θ_i)]) contre e^x calculé par mpmath à 30
   chiffres, pour x de −10⁻⁸ à −10¹², ordres 16 et 48. Tableau des erreurs absolue et relative,
   et figure log-log. C'est le contrôle le plus direct des coefficients : une faute de frappe sur
   un seul chiffre se voit à 10⁻⁵ ou pire.
2. **Matriciel** : N(t) sur la chaîne du Ra-226 (15 nucléides) à 30 j, 1 a, 100 a, par CRAM avec
   une résolution LU dense (`numpy.linalg.solve`, pas la substitution avant du C++), contre
   l'oracle mpmath et contre `scipy.linalg.expm` (Padé + scaling-and-squaring). Le point à
   montrer : expm perd quatre chiffres à 100 ans (‖At‖ ≈ 10¹³), CRAM n'en perd aucun.

Usage : python verification/scripts/cram_python.py
Produit : verification/report/cram_comparison.md,
          verification/V3_cram/figures/cram_scalar_error.png
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import numpy as np
from mpmath import mp, mpf
from scipy.linalg import expm

from expm_scipy import build_matrix, closure, read_oracle
from oracle_common import load_library

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
mp.dps = 30


def read_coefficients(path: Path) -> dict[str, tuple[float, np.ndarray, np.ndarray]]:
    """Relit les tables de src/cram_coefficients.cpp : la source de vérité est le fichier C++, pas
    une copie Python qui pourrait diverger."""
    text = path.read_text(encoding="utf-8")
    tables = {}
    for name in ("cram16", "cram48"):
        start = text.index(f"extern const CramCoefficients {name} = {{")
        end = text.index("};", start)
        body = text[start:end]
        alpha0 = float(re.search(r"=\s*\{\s*\n\s*([0-9.e+-]+),", body).group(1))
        pairs = re.findall(r"\{([+-][0-9.]+e[+-]\d+), ([+-][0-9.]+e[+-]\d+)\}", body)
        values = np.array([complex(float(a), float(b)) for a, b in pairs])
        half = len(values) // 2
        tables[name] = (alpha0, values[:half], values[half:])
    return tables


def cram_scalar(x: np.ndarray, alpha0: float, alpha: np.ndarray, theta: np.ndarray) -> np.ndarray:
    y = np.ones_like(x, dtype=float)
    for a, th in zip(alpha, theta):
        y = y + 2.0 * (a * y / (x - th)).real
    return alpha0 * y


def cram_matrix(a: np.ndarray, n0: np.ndarray, t: float, alpha0: float, alpha: np.ndarray,
                theta: np.ndarray) -> np.ndarray:
    identity = np.eye(len(n0))
    y = n0.astype(float)
    for a_i, th in zip(alpha, theta):
        y = y + 2.0 * (a_i * np.linalg.solve(a * t - th * identity, y)).real
    return alpha0 * y


def scalar_section(tables, lines: list[str], fig_path: Path) -> None:
    xs = -np.logspace(-8, 12, 4001)
    exact = np.array([float(mp.exp(mpf(x))) for x in xs])  # e^x exact ; 0.0 sous 1e-308, voulu
    lines.append("## 1. Scalaire : r_k(x) contre e^x sur l'axe négatif\n")
    lines.append("| x | e^x | erreur abs. k=16 | rel. k=16 | erreur abs. k=48 | rel. k=48 |")
    lines.append("|---|---|---|---|---|---|")
    samples = [-1e-6, -0.1, -1.0, -5.0, -10.0, -20.0, -30.0, -40.0, -60.0, -100.0, -1e3, -1e6, -1e12]
    for x in samples:
        e = float(mp.exp(mpf(x)))
        cells = [f"{x:.3g}", f"{e:.3e}"]
        for name in ("cram16", "cram48"):
            alpha0, alpha, theta = tables[name]
            r = float(cram_scalar(np.array([x]), alpha0, alpha, theta)[0])
            err = abs(r - e)
            rel = f"{err / e:.1e}" if e > 0 else "—"
            cells += [f"{err:.1e}", rel]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    window = 40  # enveloppe : maximum glissant, sinon les zéros exacts du double brouillent le tracé
    for name, style in (("cram16", "-"), ("cram48", "--")):
        alpha0, alpha, theta = tables[name]
        r = cram_scalar(xs, alpha0, alpha, theta)
        # Écart mesuré en multiprécision : |r (double) − e^x (mpmath)|, jamais exactement nul.
        err = np.array([float(abs(mpf(float(ri)) - mp.exp(mpf(x)))) for ri, x in zip(r, xs)])
        envelope = np.array([err[max(0, i - window):i + 1].max() for i in range(len(err))])
        axes[0].loglog(-xs, envelope, style, label=f"k = {name[4:]}")
        near = xs >= -150.0
        rel = envelope[near] / np.maximum(exact[near], 1e-300)
        axes[1].loglog(-xs[near], np.clip(rel, 1e-17, 1e3), style, label=f"k = {name[4:]}")
        lines.append(f"- k = {name[4:]} : erreur absolue max {err.max():.2e} (en x = {xs[err.argmax()]:.3g}), "
                     f"α₀ = {alpha0:.3e} (valeur de r_k en −∞, plancher de l'erreur absolue).")
    axes[0].set_xlabel("−x"); axes[0].set_ylabel("max glissant de |r_k(x) − e^x|")
    axes[0].set_title("erreur absolue, x de −10⁻⁸ à −10¹²"); axes[0].set_ylim(1e-50, 1e-12)
    axes[1].set_xlabel("−x"); axes[1].set_ylabel("|r_k(x) − e^x| / e^x")
    axes[1].set_title("erreur relative, x de −10⁻⁸ à −150"); axes[1].set_ylim(1e-17, 1e3)
    axes[1].axhline(1e-14, color="grey", lw=0.6, ls=":")
    for ax in axes:
        ax.grid(True, which="both", lw=0.3); ax.legend()
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=130)
    lines.append("")
    lines.append("![erreur scalaire](../V3_cram/figures/cram_scalar_error.png)")
    lines.append("")
    lines.append("Lecture : l'erreur *absolue* est bornée partout (c'est la définition de l'approximation de "
                 "Tchebychev sur ]−∞, 0]) ; l'erreur *relative* explose dès que e^x descend sous le plancher α₀. "
                 "Pour k = 16 cela arrive vers x = −30 (e^x ≈ 10⁻¹³) ; pour k = 48, vers x = −100 (10⁻⁴⁴). "
                 "Une population très petite devant la population initiale n'est donc résolue qu'en absolu.")
    lines.append("")


def matrix_section(tables, lines: list[str]) -> dict[str, float]:
    library = load_library(ROOT / "data" / "nuclides_icrp107.csv")
    names = closure(library, "Ra-226")
    a = build_matrix(library, names)
    n0 = np.zeros(len(names)); n0[names.index("Ra-226")] = 1.0
    cases = read_oracle(ROOT / "verification" / "V1_analytic_bateman" / "oracle_ra226.csv")
    lines.append("## 2. Matriciel : chaîne du Ra-226 (15 nucléides) contre l'oracle mpmath\n")
    lines.append("Écart maximal sur les 15 populations, par méthode. Relatif : populations ≥ 10⁻³⁰ ; "
                 "absolu : toutes (population totale = 1). CRAM ici : résolution LU dense NumPy, "
                 "indépendante de la substitution avant du C++.\n")
    lines.append("| cas | ‖At‖₁ | expm SciPy rel. | CRAM-16 rel. | CRAM-16 abs. | CRAM-48 rel. | CRAM-48 abs. |")
    lines.append("|---|---|---|---|---|---|---|")
    worst = {}
    for case, (t, oracle) in cases.items():
        ref = np.array([oracle.get(name, 0.0) for name in names])
        mask = ref >= 1e-30
        row = [case, f"{np.abs(a * t).sum(axis=0).max():.1e}"]
        n_expm = expm(a * t) @ n0
        row.append(f"{np.max(np.abs(n_expm[mask] / ref[mask] - 1)):.1e}")
        for name in ("cram16", "cram48"):
            alpha0, alpha, theta = tables[name]
            n = cram_matrix(a, n0, t, alpha0, alpha, theta)
            rel = np.max(np.abs(n[mask] / ref[mask] - 1)); absolute = np.max(np.abs(n - ref))
            worst[f"{name}_rel"] = max(worst.get(f"{name}_rel", 0.0), rel)
            worst[f"{name}_abs"] = max(worst.get(f"{name}_abs", 0.0), absolute)
            row += [f"{rel:.1e}", f"{absolute:.1e}"]
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    lines.append("Lecture : expm perd un chiffre par facteur 2 de ‖At‖ (mise à l'échelle puis élévations au "
                 "carré) : 10⁻⁷ à 30 j, 10⁻⁴ à 100 a. CRAM ne met rien à l'échelle : k/2 résolutions "
                 "linéaires, erreur indépendante de t. L'ordre 16 est limité *en relatif* par son plancher "
                 "absolu 2·10⁻¹⁶ sur les petites populations (Pb-206 à 30 j : 3·10⁻⁵) ; l'ordre 48 est à "
                 "l'arrondi près partout.")
    lines.append("")
    return worst


def main() -> int:
    tables = read_coefficients(ROOT / "src" / "cram_coefficients.cpp")
    assert len(tables["cram16"][1]) == 8 and len(tables["cram48"][1]) == 24
    lines = ["<!-- généré par verification/scripts/cram_python.py -->",
             "# CRAM : contrôle indépendant des coefficients et comparaison à expm", "",
             "Coefficients relus dans `src/cram_coefficients.cpp` (source : Pusa 2016, forme IPF). "
             "Référence scalaire : mpmath à 30 chiffres ; référence matricielle : oracle mpmath "
             "(`verification/V1_analytic_bateman/oracle_ra226.csv`).", ""]
    scalar_section(tables, lines, ROOT / "verification" / "V3_cram" / "figures" / "cram_scalar_error.png")
    worst = matrix_section(tables, lines)
    out = ROOT / "verification" / "report" / "cram_comparison.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"CRAM-16 : rel {worst['cram16_rel']:.1e}, abs {worst['cram16_abs']:.1e} ; "
          f"CRAM-48 : rel {worst['cram48_rel']:.1e}, abs {worst['cram48_abs']:.1e} -> {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
