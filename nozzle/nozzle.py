#!/usr/bin/env python3
"""Nozzle sizing math: scale a rocket nozzle to a test inlet and size
converging / converging-diverging nozzles. Read-only; prints results.
Standard library only. Formulas and sources are in README.md."""
import argparse
import math


# ---- isentropic flow (ideal gas, 1-D) ----

def area_ratio(M, g):
    """A/A* for Mach M."""
    return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1)))


def mach_from_area(ratio, g):
    """Supersonic Mach for a given A/A* (bisection; A/A* rises with M above 1)."""
    lo, hi = 1.0, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if area_ratio(mid, g) < ratio:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def p0_over_p(M, g):
    return (1 + (g - 1) / 2 * M * M) ** (g / (g - 1))


def mach_from_pressure(p0, p, g):
    return math.sqrt(((p0 / p) ** ((g - 1) / g) - 1) * 2 / (g - 1))


# ---- geometry ----

def converging_length(r_in, r_t, angle_deg, k_round=1.5):
    """Cone plus an upstream throat arc of radius k_round*r_t (Huzel & Huang)."""
    t = math.radians(angle_deg)
    R = k_round * r_t
    return ((r_in - r_t) + R * (1 / math.cos(t) - 1)) / math.tan(t)


def cone_length(r_a, r_b, angle_deg):
    return abs(r_b - r_a) / math.tan(math.radians(angle_deg))


def read_contour(path):
    """x, r columns (x = 0 at throat). Returns throat r, exit r, length."""
    pts = [tuple(map(float, l.split()[:2])) for l in open(path) if l.strip()]
    return pts[0][1], pts[-1][1], pts[-1][0] - pts[0][0]


# ---- report ----

def main():
    a = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    a.add_argument("--chamber", type=float, default=3.5, help="motor chamber dia, in")
    a.add_argument("--throat", type=float, default=1.5, help="motor throat dia, in")
    a.add_argument("--exit", type=float, default=3.112, help="motor exit dia, in")
    a.add_argument("--contour", help="x r z file; overrides --throat/--exit")
    a.add_argument("--inlet", type=float, default=0.622,
                   help="test inlet dia, in (1/2in Sch40 = 0.622, Sch80 = 0.546)")
    a.add_argument("--gamma", type=float, default=1.4, help="test gas gamma (air 1.4)")
    a.add_argument("--supply", type=float, default=100, help="test supply, psig")
    a.add_argument("--ambient", type=float, default=13.27, help="ambient, psia")
    a.add_argument("--cone", type=float, default=15, help="diverging cone half-angle, deg")
    o = a.parse_args()

    if o.contour:
        rt, re, L = read_contour(o.contour)
        o.throat, o.exit = 2 * rt, 2 * re
        print(f"contour {o.contour}: throat {o.throat:.4f} in, exit {o.exit:.4f} in, "
              f"length {L:.3f} in")

    cr = (o.chamber / o.throat) ** 2
    er = (o.exit / o.throat) ** 2
    s = o.inlet / o.chamber
    dt, de = o.throat * s, o.exit * s
    print(f"\nMOTOR   contraction Ac/At = {cr:.3f}   expansion Ae/At = {er:.3f}")
    print(f"\nSCALED to {o.inlet} in inlet (x{s:.4f})")
    print(f"  throat {dt:.4f} in   exit (same Ae/At) {de:.4f} in")

    g, p0 = o.gamma, o.supply + o.ambient
    M = mach_from_pressure(p0, o.ambient, g)
    er_ideal = area_ratio(M, g)
    de_ideal = dt * math.sqrt(er_ideal)
    pe = p0 / p0_over_p(mach_from_area(er, g), g)
    p0_design = o.ambient * p0_over_p(mach_from_area(er, g), g)
    print(f"\nTEST GAS  gamma {g}, supply {o.supply} psig ({p0:.2f} psia), ambient {o.ambient} psia")
    print(f"  ideal Ae/At {er_ideal:.3f}  -> exit {de_ideal:.4f} in  (exit Mach {M:.2f})")
    print(f"  motor Ae/At {er:.3f} on this gas: Pe {pe:.2f} psia, "
          f"needs {p0_design - o.ambient:.0f} psig to be ideally expanded")
    print("  supply psig   ideal Ae/At   exit dia in")
    for psig in (60, 90, 100, 120, 150):
        e = area_ratio(mach_from_pressure(psig + o.ambient, o.ambient, g), g)
        print(f"  {psig:>11}   {e:>11.3f}   {dt * math.sqrt(e):>11.4f}")

    print("\nCONVERGING SECTION (upstream throat radius 1.5 Rt)")
    print("  half-angle   motor in   test in   (sharp cone test in)")
    for ang in (30, 37.5, 45, 60):
        print(f"  {ang:>10}   {converging_length(o.chamber / 2, o.throat / 2, ang):>8.3f}"
              f"   {converging_length(o.inlet / 2, dt / 2, ang):>7.3f}"
              f"   {cone_length(o.inlet / 2, dt / 2, ang):>7.3f}")
    print(f"  throat radii, test: upstream 1.5 Rt = {1.5 * dt / 2:.4f} in, "
          f"downstream 0.382 Rt = {0.382 * dt / 2:.4f} in")

    lam = (1 + math.cos(math.radians(o.cone))) / 2
    print(f"\nDIVERGING CONE {o.cone} deg (divergence factor {lam:.4f})")
    print(f"  motor-ratio exit: {cone_length(dt / 2, de / 2, o.cone):.3f} in")
    print(f"  ideal-for-gas exit: {cone_length(dt / 2, de_ideal / 2, o.cone):.3f} in")


if __name__ == "__main__":
    main()
