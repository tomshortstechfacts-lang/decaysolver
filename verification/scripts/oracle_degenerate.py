"""Oracle multiprécision pour les cas de vérification du jeu D (pathologie d'annulation).

Chaîne linéaire à trois membres, N(0) = (1, 0, 0), formule de Bateman évaluée avec mpmath à
50 chiffres : à cette précision, l'annulation catastrophique de la somme alternée ne coûte
qu'une dizaine de chiffres, il en reste largement assez pour servir de référence au `double`.
Le cas exactement dégénéré D1 utilise la forme confluente λ t e^{-λt}, pas la formule.

Usage : python verification/scripts/oracle_degenerate.py
Les valeurs imprimées sont celles reprises dans tests/unit/test_bateman.cpp.
"""

from mpmath import mp, mpf, exp

mp.dps = 50

CASES = {
    "D2": (mpf(3), mpf(3) - mpf("1e-7"), mpf(3) + mpf("1e-7")),
    "D3": (mpf(3), mpf(3) - mpf("1e-11"), mpf(3) + mpf("1e-11")),
    "D4": (mpf(1), mpf(2), mpf(3)),
}
T = mpf(4)


def bateman_last(lambdas, t):
    """N_k(t) pour N_0(0) = 1, chaîne linéaire, constantes toutes distinctes."""
    k = len(lambdas) - 1
    prefactor = mpf(1)
    for lam in lambdas[:-1]:
        prefactor *= lam
    total = mpf(0)
    for m, lam_m in enumerate(lambdas):
        denominator = mpf(1)
        for l, lam_l in enumerate(lambdas):
            if l != m:
                denominator *= lam_l - lam_m
        total += exp(-lam_m * t) / denominator
    return prefactor * total


def main() -> None:
    print(f"t = {T}, mp.dps = {mp.dps}")
    for name, lambdas in CASES.items():
        n2 = bateman_last(lambdas, T)
        n1 = bateman_last(lambdas[:2], T)
        print(f"{name}: lambda = {[mp.nstr(x, 15) for x in lambdas]}")
        print(f"     N1(t) = {mp.nstr(n1, 20)}")
        print(f"     N2(t) = {mp.nstr(n2, 20)}")
    # D1 : (1, 1), forme confluente N1 = lambda t exp(-lambda t)
    lam = mpf(1)
    print(f"D1: lambda = (1, 1), N1(t) = {mp.nstr(lam * T * exp(-lam * T), 20)}")


if __name__ == "__main__":
    main()
