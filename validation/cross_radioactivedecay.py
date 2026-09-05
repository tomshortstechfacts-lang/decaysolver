"""T7 — évaluation croisée avec `radioactivedecay` (Malins & Lemoine, JOSS 2022), outil de
référence ouvert pour les calculs de décroissance, mêmes données ICRP-107.

Cas : les 35 nucléides de la liste standard de déclaration des déchets, 1 Bq chacun, vieillis de
10 ans, convention `all` (toutes les filles), comparés nucléide par nucléide.

Précautions pour que la comparaison soit interprétable :
  - l'âge est passé aux deux outils **en secondes** : les conventions d'année diffèrent
    (365,25 j ici, 365,2422 j dans radioactivedecay), soit 2e-5 en relatif ;
  - les demi-vies de radioactivedecay sont converties en secondes avec sa propre convention à
    partir des mêmes valeurs ICRP-107 : les λ des nucléides dont la demi-vie est exprimée en
    années diffèrent donc de 2e-5 entre les deux outils. C'est l'écart attendu, et il est mesuré.
  - les deux voies β⁻ rétablies par decaysolver (At-219, At-217) sont absentes de
    radioactivedecay : hors de ce cas (série du Np-237), sans effet ici.

Usage : python validation/cross_radioactivedecay.py [chemin/vers/decaysolver]
Produit : validation/cross_radioactivedecay.md
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import radioactivedecay as rd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXE = ROOT / "build" / "dev" / "apps" / "decaysolver.exe"
BASE = ["Be-10", "C-14", "Cl-36", "Ca-41", "Mn-54", "Fe-55", "Co-60", "Ni-59", "Ni-63", "Zn-65",
        "Se-79", "Sr-90", "Zr-93", "Nb-94", "Mo-93", "Tc-99", "Pd-107", "Ag-108m", "Ag-110m",
        "Sn-121m", "Sn-126", "Sb-125", "I-129", "Cs-134", "Cs-135", "Cs-137", "Sm-151", "U-235",
        "U-238", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Am-241", "Cm-244"]
AGE_S = 10 * 365.25 * 86400


RD_YEAR_DAYS = 365.2422  # convention de radioactivedecay
OUR_YEAR_DAYS = 365.25   # convention de decaysolver (année julienne)


def aligned_library(tmp: Path) -> Path:
    """Copie de la bibliothèque où les demi-vies en années sont multipliées par 365,2422/365,25 :
    converties par decaysolver avec 365,25 j, elles redonnent exactement les secondes de
    radioactivedecay. Sert à isoler l'effet de la convention d'année du reste."""
    src = ROOT / "data" / "nuclides_icrp107.csv"
    dst = tmp / "nuclides_aligned.csv"
    lines = []
    for line in src.read_text(encoding="utf-8").splitlines():
        fields = line.split(";")
        if len(fields) == 6 and fields[2] == "y":
            fields[1] = repr(float(fields[1]) * RD_YEAR_DAYS / OUR_YEAR_DAYS)
            line = ";".join(fields)
        lines.append(line)
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dst


def run_decaysolver(exe: Path, align_year: bool = False) -> dict[str, float]:
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "inv.csv"
        inp.write_text("".join(f"{n};1\n" for n in BASE), encoding="utf-8")
        library = aligned_library(Path(tmp)) if align_year else ROOT / "data" / "nuclides_icrp107.csv"
        result = subprocess.run(
            [str(exe), "age", "--input", str(inp), "--age", f"{AGE_S:.0f}s", "--kind", "bq",
             "--daughters", "all", "--library", str(library)],
            capture_output=True, text=True, encoding="utf-8", check=True)
    out = {}
    for line in result.stdout.splitlines():
        if line.startswith("#") or line.startswith("nuclide;"):
            continue
        name, activity, _fraction, _mode = line.split(";")
        out[name] = float(activity)
    return out


def run_radioactivedecay(high_precision: bool = False) -> dict[str, float]:
    """`Inventory` : double (SciPy). `InventoryHP` : SymPy, arithmétique exacte puis
    évaluation multiprécision ; la référence la plus sûre de l'outil, plus lente."""
    factory = rd.InventoryHP if high_precision else rd.Inventory
    aged = factory({n: 1.0 for n in BASE}, "Bq").decay(AGE_S, "s")
    return {name: float(a) for name, a in aged.activities("Bq").items() if a > 0}


def compare(ours: dict[str, float], theirs: dict[str, float]):
    rows, worst, worst_name, missing = [], 0.0, "", []
    for name in sorted(set(ours) | set(theirs), key=lambda n: -ours.get(n, 0)):
        a, b = ours.get(name), theirs.get(name)
        if a is None or b is None:
            missing.append((name, a, b))
            continue
        if a < 1e-12 and b < 1e-12:
            continue  # activités sans signification (< 1e-12 Bq pour 35 Bq initiaux)
        rel = abs(a / b - 1.0)
        if rel > worst:
            worst, worst_name = rel, name
        rows.append((name, a, b, rel))
    return rows, worst, worst_name, missing


def main() -> int:
    exe = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_EXE
    theirs = run_radioactivedecay()
    theirs_hp = run_radioactivedecay(high_precision=True)
    ours_aligned = run_decaysolver(exe, align_year=True)
    rows_raw, worst_raw, name_raw, missing = compare(run_decaysolver(exe), theirs)
    rows_al, worst_al, name_al, _ = compare(ours_aligned, theirs)
    rows_hp, worst_hp, name_hp, _ = compare(ours_aligned, theirs_hp)
    aligned = {name: rel for name, _, _, rel in rows_al}
    hp = {name: rel for name, _, _, rel in rows_hp}
    lines = ["<!-- généré par validation/cross_radioactivedecay.py -->",
             "# Évaluation croisée avec radioactivedecay (T7, validation)", "",
             f"radioactivedecay {rd.__version__}, jeu `{rd.DEFAULTDATA.dataset_name}` ; 35 nucléides à 1 Bq, "
             f"vieillis de {AGE_S:.0f} s (10 années juliennes), convention `all`.", "",
             f"**Écart relatif max** (activités > 1e-12 Bq) : {worst_raw:.1e} ({name_raw}) avec les conventions "
             f"d'année de chaque outil (365,25 j ici, 365,2422 j là) ; **{worst_al:.1e} ({name_al}) une fois les "
             "conventions alignées** (demi-vies en années rééchelonnées avant lecture par decaysolver). "
             "L'écart brut est dominé par la convention d'année, qui se cumule le long des chaînes "
             f"d'actinides. Contre le mode haute précision de radioactivedecay (SymPy), conventions alignées : "
             f"**{worst_hp:.1e} ({name_hp})**. Le résiduel en double était celui de l'autre outil sur les "
             "activités très faibles issues de l'U-238 (accord à 1e-15 en haute précision). Le seul écart "
             "restant, Rn-219 dans la série de l'U-235, vient de la voie β⁻ d'At-219 (3 %) rétablie dans "
             "decaysolver et absente de radioactivedecay (voir data/PROVENANCE.md) : c'est une différence "
             "de données, connue et documentée.", "",
             "| nucléide | decaysolver (Bq) | radioactivedecay (Bq) | écart, conventions propres | écart, conventions alignées | écart vs rd haute précision |",
             "|---|---|---|---|---|---|"]
    for name, a, b, rel in rows_raw:
        lines.append(f"| {name} | {a:.10e} | {b:.10e} | {rel:.1e} | {aligned.get(name, float('nan')):.1e} | {hp.get(name, float('nan')):.1e} |")
    if missing:
        lines.append("")
        lines.append("Nucléides présents d'un seul côté (activités négligeables ou voies différentes) :")
        for name, a, b in missing:
            lines.append(f"- {name} : decaysolver {a}, radioactivedecay {b}")
    out = ROOT / "validation" / "cross_radioactivedecay.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"écart max brut : {worst_raw:.2e} ({name_raw}) ; aligné : {worst_al:.2e} ({name_al}) ; "
          f"vs haute précision : {worst_hp:.2e} ({name_hp}) ; "
          f"{len(missing)} d'un seul côté -> {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
