"""Extraction de la bibliothèque de nucléides depuis ICRP-107 (via le paquet radioactivedecay).

Produit data/nuclides_icrp107.csv : les 35 nucléides de la liste standard de déclaration des
déchets, plus la fermeture complète de leurs descendants (aucune chaîne tronquée), avec les
demi-vies dans leur valeur et leur unité d'origine. La conversion en secondes est faite par
decaysolver lui-même (année julienne 365,25 j), pas ici : la donnée committée est la donnée
source, pas une donnée dérivée.

Usage (depuis la racine du dépôt, environnement Python avec radioactivedecay) :
    python data/scripts/extract_icrp107.py          # 35 nucléides + descendants
    python data/scripts/extract_icrp107.py --full   # les 1512 nucléides du jeu

Le fichier produit est déterministe (tri par nom) : deux exécutions donnent le même SHA-256.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import sys
from pathlib import Path

import radioactivedecay as rd

# Liste standard de déclaration des déchets (35 radionucléides).
BASE_NUCLIDES = [
    "Be-10", "C-14", "Cl-36", "Ca-41", "Mn-54", "Fe-55", "Co-60", "Ni-59", "Ni-63", "Zn-65",
    "Se-79", "Sr-90", "Zr-93", "Nb-94", "Mo-93", "Tc-99", "Pd-107", "Ag-108m", "Ag-110m",
    "Sn-121m", "Sn-126", "Sb-125", "I-129", "Cs-134", "Cs-135", "Cs-137", "Sm-151", "U-235",
    "U-238", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Am-241", "Cm-244",
]

# Vocabulaire des modes dans le CSV. ICRP-107 regroupe β+ et capture électronique.
MODE_MAP = {
    "α": "alpha",
    "β-": "beta-",
    "β+": "beta+/EC",
    "EC": "beta+/EC",
    "β+ & EC": "beta+/EC",
    "IT": "IT",
    "SF": "SF",
    "unlisted": "unlisted",
}

# Tolérance du chargeur de decaysolver sur |Σ b − 1| ; au-delà, le fichier est refusé.
LOADER_TOLERANCE = 5e-4

HEADER = "nuclide;half_life_value;half_life_unit;mode;daughter;branching_fraction"

# Voies absentes du jeu radioactivedecay alors qu'ICRP-107 les contient : At-219 et At-217 n'y
# ont que leur voie α (0,97 et 0,99988). La voie β⁻ est rétablie avec, pour rapport, le
# complément à 1 de la voie α d'ICRP-107, et pour fille l'isobare Z+1 (Rn-219, Rn-217), ce que
# confirme ENSDF (API Livechart de l'AIEA, nds.iaea.org/relnsd, consultée le 2026-09-05 :
# At-219 β⁻ 6,4 %, At-217 β⁻ 0,007 % dans l'évaluation actuelle ; l'écart avec ICRP-107 illustre
# la dispersion des évaluations sur ces nucléides très mineurs).
MISSING_BRANCHES = {
    "At-219": ("β-", "Rn-219"),
    "At-217": ("β-", "Rn-217"),
}

# Seuil au-delà duquel un écart |Σ b − 1| est rapporté dans l'en-tête du fichier produit.
SUM_REPORT_THRESHOLD = 1e-5


def branches_of(data: rd.DecayData, name: str, full: bool = False) -> list[tuple[str, str, float]]:
    """Voies (mode, fille, rapport) d'un nucléide, voies manquantes rétablies.

    En mode --full, un déficit de somme au-delà de la tolérance du chargeur (nucléides exotiques
    dont le jeu redistribué omet une voie sans qu'on ait vérifié la fille sur une seconde source)
    devient une voie `unlisted` sans fille : le déficit reste visible dans le fichier et la somme
    des rapports vaut 1. Les nucléides concernés sont listés dans l'en-tête."""
    index = data.nuclide_dict[name]
    branches = list(zip(data.modes[index], data.progeny[index], (float(b) for b in data.bfs[index])))
    if name in MISSING_BRANCHES:
        mode, daughter = MISSING_BRANCHES[name]
        branches.append((mode, daughter, 1.0 - sum(b for _, _, b in branches)))
    deficit = 1.0 - sum(b for _, _, b in branches)
    if full and deficit > LOADER_TOLERANCE:
        branches.append(("unlisted", "", deficit))
    return branches


def closure(data: rd.DecayData, base: list[str]) -> list[str]:
    """Tous les nucléides atteignables depuis `base` par filiation (fission exclue)."""
    seen: set[str] = set()
    stack = list(base)
    while stack:
        name = stack.pop()
        if name in seen or name == "SF":
            continue
        seen.add(name)
        stack.extend(daughter for _, daughter, _ in branches_of(data, name))
    return sorted(seen)


def is_stable(data: rd.DecayData, name: str) -> bool:
    return data.hldata[data.nuclide_dict[name]][2] == "stable"


def rows_for(data: rd.DecayData, name: str) -> list[str]:
    index = data.nuclide_dict[name]
    value, unit, _readable = data.hldata[index]
    if is_stable(data, name):
        return [f"{name};stable;;stable;;"]
    lines = []
    for mode, daughter, fraction in branches_of(data, name, FULL_MODE[0]):
        daughter_out = "" if daughter == "SF" else daughter
        lines.append(f"{name};{float(value):.10g};{unit};{MODE_MAP[mode]};{daughter_out};{fraction:.10g}")
    return lines


FULL_MODE = [False]  # positionné par main() ; évite de propager un paramètre partout


def unlisted_nuclides(data: rd.DecayData, names: list[str]) -> list[tuple[str, float]]:
    """Nucléides recevant une voie `unlisted` en mode --full, avec le déficit."""
    out = []
    for name in names:
        if is_stable(data, name):
            continue
        deficit = 1.0 - sum(b for _, _, b in branches_of(data, name))
        if deficit > LOADER_TOLERANCE:
            out.append((name, deficit))
    return sorted(out, key=lambda item: -item[1])


def sum_deviations(data: rd.DecayData, names: list[str]) -> list[tuple[str, float]]:
    """Nucléides radioactifs dont |Σ b − 1| dépasse SUM_REPORT_THRESHOLD, par écart décroissant."""
    out = []
    for name in names:
        if is_stable(data, name):
            continue
        deviation = sum(b for _, _, b in branches_of(data, name)) - 1.0
        if abs(deviation) > SUM_REPORT_THRESHOLD:
            out.append((name, deviation))
    return sorted(out, key=lambda item: -abs(item[1]))


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    data = rd.DEFAULTDATA
    full = "--full" in sys.argv
    FULL_MODE[0] = full
    if full:
        # Bibliothèque complète : tous les nucléides du jeu ICRP-107 redistribué (1512), fermeture
        # incluse par construction. Fichier data/nuclides_icrp107_full.csv, à passer par --library.
        out_path = repo_root / "data" / "nuclides_icrp107_full.csv"
        names = sorted(str(n) for n in data.nuclides)
    else:
        out_path = repo_root / "data" / "nuclides_icrp107.csv"
        names = closure(data, BASE_NUCLIDES)

    header_lines = [
        "# decaysolver — bibliothèque de nucléides",
        "# source: ICRP Publication 107, Nuclear Decay Data for Dosimetric Calculations (2008),",
        "#         Copyright (c) 2008 A. Endo and K.F. Eckerman ; licence : data/LICENSE.ICRP-07",
        f"# extraction: paquet radioactivedecay {rd.__version__}, jeu '{data.dataset_name}'",
        f"# script: data/scripts/extract_icrp107.py ; date UTC: {dt.datetime.now(dt.timezone.utc):%Y-%m-%d}",
        (f"# contenu: bibliothèque complète, {len(names)} nucléides du jeu ICRP-107 redistribué"
         if full else
         f"# contenu: {len(BASE_NUCLIDES)} nucléides de base (liste standard de déclaration des déchets)"),
        ("#          (--full)" if full else
         f"#          + fermeture complète de leurs descendants = {len(names)} nucléides"),
        "# demi-vies: valeur et unité d'origine ICRP-107 (s, m, h, d, y) ; aucune conversion ici.",
        "#            decaysolver convertit avec 1 y = 365,25 j (année julienne) ;",
        "#            radioactivedecay utilise 365,2422 j : écart relatif 2e-5, à garder en tête",
        "#            lors de l'évaluation croisée.",
        "# incertitudes: ICRP-107 n'en fournit pas ; voir data/PROVENANCE.md.",
        "# fission spontanée (SF): voie conservée avec fille vide, produits de fission non suivis.",
        "# voies rétablies (absentes du jeu radioactivedecay, voir le script) : "
        + ", ".join(f"{n} {MISSING_BRANCHES[n][0]} vers {MISSING_BRANCHES[n][1]}" for n in MISSING_BRANCHES),
        "# écarts |somme des rapports - 1| > 1e-5, dus aux arrondis indépendants d'ICRP-107 :",
    ]
    for name, deviation in sum_deviations(data, names):
        if full and deviation < -LOADER_TOLERANCE:
            continue  # traité ci-dessous comme voie non répertoriée
        header_lines.append(f"#     {name}: {deviation:+.1e}")
    if full:
        header_lines.append("# voies non répertoriées (déficit de somme > 5e-4 dans le jeu redistribué, fille")
        header_lines.append("#   inconnue, noyaux sortant du système) : ")
        for name, deficit in unlisted_nuclides(data, names):
            header_lines.append(f"#     {name}: {deficit:.1e}")
    header_lines += [
        "# une ligne par voie de décroissance ; les stables ont une seule ligne 'stable'.",
        HEADER,
    ]
    body = [line for name in names for line in rows_for(data, name)]
    text = "\n".join(header_lines + body) + "\n"
    out_path.write_text(text, encoding="utf-8", newline="\n")

    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    out_path.with_suffix(".sha256").write_text(
        f"{digest}  {out_path.name}\n", encoding="utf-8", newline="\n")
    print(f"{out_path.relative_to(repo_root)} : {len(names)} nucléides, {len(body)} voies")
    print(f"SHA-256 : {digest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
