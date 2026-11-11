"""
Injector face layout for the N2O / IPA unlike-doublet injector.

One ring of radially-oriented doublets: fuel hole outboard aimed inward,
ox hole inboard aimed outward, meeting at an impingement circle. The spray
fan from each element spreads tangentially, so the ring of fans forms an
annulus and keeps the spray off the ablative.

Each entry in LAYOUTS is drilled with N elements, sized so N_DESIGN open
elements match the HalfCatSim CdA, and plugged down in balanced patterns.

Outputs per layout, in output/<name>/: <name>.svg (face + element section),
<name>.csv (hole coordinates). Run with --png to also write <name>.png
(macOS, uses Quick Look).
"""
import csv
import math
import os
import subprocess
import sys
import types

# ---------------- INPUTS ----------------
CHAMBER_DIA_MM = 4.75 * 25.4   # chamber bore at the injector face
ANG_FU = 30.0  # deg off axis, fuel
ANG_OX = 17.2  # deg off axis, ox -> combined spray goes straight down the axis
IMPINGE_L = 14.0   # mm below face (~4 d_ox; short is good for flashing N2O)
R_IMPINGE = 30.0   # mm radius of impingement circle (~half the chamber radius)

# Plug thread envelope on the face (for spacing check). Set to the major
# diameter of whatever plug/insert thread you choose.
THREAD_OX = 5.0    # mm  e.g. M5 / #10
THREAD_FU = 3.0    # mm  e.g. M3 / #4
MIN_WALL = 2.0     # mm  minimum material between thread envelopes
FACE_PLATE_T = 8.0  # mm  injector face plate thickness (for hole L/D)

# HalfCatSim totals, for checking the chosen drills
CDA_OX, CDA_FU = 3.0198955e-05, 7.0045861e-06   # m^2
CD_OX, CD_FU = 0.40, 0.60

# name, elements drilled, elements open at design, ox drill, fuel drill,
# plug patterns {open count: [element indices to plug]}, ox on outer ring?
LAYOUTS = [
    dict(name="layout_12_elements", n=12, n_design=8,
         d_ox=3.454, d_fu=1.35,          # #29, 1.35 mm  (8-pair ideal 3.47 / 1.36)
         plugs={12: [], 10: [0, 6], 9: [0, 4, 8], 8: [0, 3, 6, 9], 6: [0, 2, 4, 6, 8, 10]},
         ox_outside=False),
    dict(name="layout_15_elements", n=15, n_design=12,
         d_ox=2.819, d_fu=1.10,          # #34, 1.10 mm  (12-pair ideal 2.83 / 1.11)
         plugs={15: [], 12: [0, 5, 10], 10: [0, 3, 6, 9, 12]},
         ox_outside=False),
    dict(name="layout_12_elements_OX_OUTSIDE", n=12, n_design=8,
         d_ox=3.454, d_fu=1.35,
         plugs={12: [], 10: [0, 6], 9: [0, 4, 8], 8: [0, 3, 6, 9], 6: [0, 2, 4, 6, 8, 10]},
         ox_outside=True),
    dict(name="layout_15_elements_OX_OUTSIDE", n=15, n_design=12,
         d_ox=2.819, d_fu=1.10,
         plugs={15: [], 12: [0, 5, 10], 10: [0, 3, 6, 9, 12]},
         ox_outside=True),
]
# ----------------------------------------

R_CH = CHAMBER_DIA_MM / 2
OUT_DIR = "output"


def chord(r, deg):
    return 2 * r * math.sin(math.radians(deg) / 2)


def area(d_mm):
    return math.pi * (d_mm / 1000) ** 2 / 4


def geometry(name):
    """All inputs plus one layout's derived geometry, for other scripts to import."""
    lay = next(l for l in LAYOUTS if l["name"] == name)
    sgn_ox = 1 if lay["ox_outside"] else -1
    g = {k: v for k, v in globals().items() if k.isupper() or k in ("area", "chord")}
    g.update(N_DRILLED=lay["n"], N_DESIGN=lay["n_design"], D_OX=lay["d_ox"], D_FU=lay["d_fu"],
             PLUG_SETS=lay["plugs"], OX_OUTSIDE=lay["ox_outside"], pitch=360.0 / lay["n"],
             r_ox=R_IMPINGE + sgn_ox * IMPINGE_L * math.tan(math.radians(ANG_OX)),
             r_fu=R_IMPINGE - sgn_ox * IMPINGE_L * math.tan(math.radians(ANG_FU)),
             OX_C="#d9480f", FU_C="#1971c2", PLUG_C="#868e96")
    return types.SimpleNamespace(**g)


def make_layout(name, n, n_design, d_ox, d_fu, plugs, ox_outside):
    out = os.path.join(OUT_DIR, name)
    os.makedirs(out, exist_ok=True)
    N_DRILLED, D_OX, D_FU, PLUG_SETS = n, d_ox, d_fu, plugs
    sgn_ox = 1 if ox_outside else -1   # which side of the impingement circle the ox hole sits
    r_ox = R_IMPINGE + sgn_ox * IMPINGE_L * math.tan(math.radians(ANG_OX))
    r_fu = R_IMPINGE - sgn_ox * IMPINGE_L * math.tan(math.radians(ANG_FU))
    r_in, r_out = min(r_ox, r_fu), max(r_ox, r_fu)
    thr_in, thr_out = (THREAD_FU, THREAD_OX) if ox_outside else (THREAD_OX, THREAD_FU)
    pitch = 360.0 / N_DRILLED

    # ---------- checks ----------
    checks = []
    checks.append(("ox ring radius", f"{r_ox:.2f} mm"))
    checks.append(("fuel ring radius", f"{r_fu:.2f} mm"))
    checks.append(("ox-fuel spacing in element", f"{r_out - r_in:.2f} mm"))
    wall_ox = chord(r_ox, pitch) - THREAD_OX
    wall_fu = chord(r_fu, pitch) - THREAD_FU
    wall_in = (r_out - r_in) - (THREAD_OX + THREAD_FU) / 2
    wall_out = R_CH - r_out - thr_out / 2
    checks.append(("wall between ox threads", f"{wall_ox:.2f} mm {'OK' if wall_ox >= MIN_WALL else 'TOO TIGHT'}"))
    checks.append(("wall between fuel threads", f"{wall_fu:.2f} mm {'OK' if wall_fu >= MIN_WALL else 'TOO TIGHT'}"))
    checks.append(("wall ox-to-fuel thread", f"{wall_in:.2f} mm {'OK' if wall_in >= MIN_WALL else 'TOO TIGHT'}"))
    checks.append(("outer thread to chamber wall", f"{wall_out:.2f} mm {'OK' if wall_out >= MIN_WALL else 'TOO TIGHT'}"))
    checks.append(("ox hole L/D", f"{FACE_PLATE_T / math.cos(math.radians(ANG_OX)) / D_OX:.1f}"))
    checks.append(("fuel hole L/D", f"{FACE_PLATE_T / math.cos(math.radians(ANG_FU)) / D_FU:.1f}"))

    print(f"===== {name} =====")
    print("Geometry checks")
    for k, v in checks:
        print(f"  {k:30s} {v}")

    print(f"\nOpen elements vs HalfCatSim design CdA ({n_design}-pair ideal = 100%)")
    print("  open | ox CdA % | fuel CdA % | O/F shift")
    for n in sorted(PLUG_SETS, reverse=True):
        ox = n * area(D_OX) * CD_OX / CDA_OX
        fu = n * area(D_FU) * CD_FU / CDA_FU
        print(f"  {n:4d} | {ox*100:7.1f}  | {fu*100:8.1f}   | {(ox/fu - 1)*100:+.1f}%")

    # ---------- coordinates ----------
    rows = []
    for i in range(N_DRILLED):
        th = math.radians(90 - i * pitch)  # element 0 at 12 o'clock, clockwise
        for prop, r, d, ang in (("OX", r_ox, D_OX, ANG_OX), ("FUEL", r_fu, D_FU, ANG_FU)):
            rows.append({
                "element": i, "propellant": prop,
                "x_mm": round(r * math.cos(th), 3), "y_mm": round(r * math.sin(th), 3),
                "radius_mm": round(r, 3), "clock_deg": round(i * pitch, 1),
                "drill_mm": d,
                "tilt_deg": ang,
                "tilt_direction": "outward (away from center)" if r == r_in else "inward (toward center)",
                f"plug_for_{n_design}": "PLUG" if i in PLUG_SETS[n_design] else "",
            })

    with open(os.path.join(out, f"{name}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    # ---------- SVG ----------
    S = 4.0   # px per mm for face view
    FW = 2 * (R_CH + 12) * S
    cx = cy = FW / 2
    OX_C, FU_C, PLUG_C = "#d9480f", "#1971c2", "#868e96"

    svg = []
    add = svg.append
    W = H = FW + 600   # square, so the Quick Look PNG needs no cropping
    add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="13">')
    add(f'<rect width="100%" height="100%" fill="white"/>')
    add(f'<text x="{cx:.0f}" y="24" text-anchor="middle" font-size="18" font-weight="bold">'
        f'{N_DRILLED} elements ({2*N_DRILLED} holes), looking upstream from chamber</text>')
    if ox_outside:
        add(f'<text x="{cx:.0f}" y="46" text-anchor="middle" font-size="15" font-weight="bold" fill="#c92a2a">'
            f'NOT RECOMMENDED: ox on the outer ring puts oxidizer-rich gas on the ablative</text>')

    def P(x, y):
        return cx + x * S, cy + 20 - y * S

    ccx, ccy = P(0, 0)
    add(f'<circle cx="{ccx}" cy="{ccy}" r="{R_CH*S}" fill="#f8f9fa" stroke="black" stroke-width="2"/>')
    add(f'<circle cx="{ccx}" cy="{ccy}" r="{R_IMPINGE*S}" fill="none" stroke="#2b8a3e" stroke-dasharray="6 4"/>')
    add(f'<text x="{ccx}" y="{ccy + 32}" text-anchor="middle" fill="#2b8a3e">impingement circle R{R_IMPINGE:g} (dashed)</text>')
    for r in (r_in, r_out):
        add(f'<circle cx="{ccx}" cy="{ccy}" r="{r*S}" fill="none" stroke="#ced4da" stroke-width="0.8"/>')
    add(f'<line x1="{ccx-8}" y1="{ccy}" x2="{ccx+8}" y2="{ccy}" stroke="black"/>')
    add(f'<line x1="{ccx}" y1="{ccy-8}" x2="{ccx}" y2="{ccy+8}" stroke="black"/>')

    for row in rows:
        x, y = P(row["x_mm"], row["y_mm"])
        plugged = bool(row[f"plug_for_{n_design}"])
        col = PLUG_C if plugged else (OX_C if row["propellant"] == "OX" else FU_C)
        thr = THREAD_OX if row["propellant"] == "OX" else THREAD_FU
        add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{thr/2*S:.1f}" fill="none" stroke="{col}" stroke-dasharray="3 2"/>')
        add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{row["drill_mm"]/2*S:.1f}" fill="{col}"/>')
        if plugged:
            k = thr / 2 * S * 0.7
            add(f'<line x1="{x-k:.1f}" y1="{y-k:.1f}" x2="{x+k:.1f}" y2="{y+k:.1f}" stroke="{PLUG_C}" stroke-width="2"/>')
            add(f'<line x1="{x-k:.1f}" y1="{y+k:.1f}" x2="{x+k:.1f}" y2="{y-k:.1f}" stroke="{PLUG_C}" stroke-width="2"/>')
        if row["radius_mm"] == round(r_out, 3):
            th = math.radians(90 - row["element"] * pitch)
            lx, ly = P((r_out + 9) * math.cos(th), (r_out + 9) * math.sin(th))
            add(f'<text x="{lx:.1f}" y="{ly+4:.1f}" text-anchor="middle" font-weight="bold">{row["element"]}</text>')
        # arrow showing jet direction (radial)
        th = math.radians(90 - row["element"] * pitch)
        sgn = 1 if row["radius_mm"] == round(r_in, 3) else -1
        r0 = row["radius_mm"]
        x2, y2 = P((r0 + sgn * 3.5) * math.cos(th), (r0 + sgn * 3.5) * math.sin(th))
        if not plugged:
            add(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{col}" stroke-width="1.5"/>')

    # legend + notes
    lx = FW + 10
    ly = 60
    def text(s, dy=20, **kw):
        nonlocal ly
        attrs = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
        add(f'<text x="{lx}" y="{ly}" {attrs}>{s}</text>')
        ly += dy

    text("LEGEND", font_weight="bold", font_size="15")
    add(f'<circle cx="{lx+8}" cy="{ly-5}" r="6" fill="{OX_C}"/>'); lx += 22
    text(f"N2O hole  {D_OX} mm, tilted {ANG_OX}° {'inward' if ox_outside else 'outward'}"); lx -= 22
    add(f'<circle cx="{lx+8}" cy="{ly-5}" r="4" fill="{FU_C}"/>'); lx += 22
    text(f"IPA hole  {D_FU} mm, tilted {ANG_FU}° {'outward' if ox_outside else 'inward'}"); lx -= 22
    add(f'<circle cx="{lx+8}" cy="{ly-5}" r="6" fill="{PLUG_C}"/>'); lx += 22
    text(f"plugged pair ({n_design}-open config: {', '.join(map(str, PLUG_SETS[n_design]))})"); lx -= 22
    add(f'<circle cx="{lx+8}" cy="{ly-5}" r="7" fill="none" stroke="black" stroke-dasharray="3 2"/>'); lx += 22
    text(f"thread envelope (ox {THREAD_OX:g}, fuel {THREAD_FU:g} mm)"); lx -= 22
    ly += 10
    text("PLUG PATTERNS (element numbers)", font_weight="bold", font_size="15")
    for n in sorted(PLUG_SETS, reverse=True):
        text(f"{n:>2} open: plug {', '.join(map(str, PLUG_SETS[n])) or 'none'}")
    ly += 10
    text("GEOMETRY", font_weight="bold", font_size="15")
    for k, v in checks:
        text(f"{k}: {v}", dy=18, font_size="12")

    # element section view
    ly += 20
    text("ELEMENT SECTION (radial plane)", font_weight="bold", font_size="15")
    s2 = 7.0
    ox0 = lx + 40
    fy = ly + 10
    ix = ox0 + (R_IMPINGE - r_in) * s2
    iy = fy + FACE_PLATE_T * s2 + IMPINGE_L * s2
    pt = FACE_PLATE_T * s2
    add(f'<rect x="{ox0-30}" y="{fy}" width="{(r_out-r_in)*s2+60}" height="{pt}" fill="#e9ecef" stroke="black"/>')
    add(f'<text x="{ox0-34}" y="{fy+pt/2+4}" text-anchor="end" font-size="11">face plate {FACE_PLATE_T:g} mm</text>')
    face_y = fy + pt
    (in_c, in_d, in_a), (out_c, out_d, out_a) = sorted(
        [(OX_C, D_OX, ANG_OX, r_ox), (FU_C, D_FU, ANG_FU, r_fu)], key=lambda t: t[3])[0][:3], \
        sorted([(OX_C, D_OX, ANG_OX, r_ox), (FU_C, D_FU, ANG_FU, r_fu)], key=lambda t: t[3])[1][:3]
    t_in, t_out = math.tan(math.radians(in_a)), math.tan(math.radians(out_a))
    add(f'<line x1="{ix - (iy-fy)*t_in:.1f}" y1="{fy}" x2="{ix:.1f}" y2="{iy:.1f}" stroke="{in_c}" stroke-width="{in_d*s2*0.5:.1f}" opacity="0.8"/>')
    add(f'<line x1="{ix + (iy-fy)*t_out:.1f}" y1="{fy}" x2="{ix:.1f}" y2="{iy:.1f}" stroke="{out_c}" stroke-width="{out_d*s2*0.5:.1f}" opacity="0.8"/>')
    add(f'<line x1="{ix:.1f}" y1="{iy:.1f}" x2="{ix:.1f}" y2="{iy+40:.1f}" stroke="#2b8a3e" stroke-width="2" stroke-dasharray="5 3"/>')
    add(f'<circle cx="{ix:.1f}" cy="{iy:.1f}" r="4" fill="#2b8a3e"/>')
    add(f'<text x="{ix+8:.1f}" y="{iy+4:.1f}" font-size="11" fill="#2b8a3e">impinge {IMPINGE_L:g} mm below face</text>')
    add(f'<text x="{ox0-30}" y="{iy+58:.1f}" font-size="11">← center</text>')
    add(f'<text x="{ox0+(r_out-r_in)*s2+30}" y="{iy+58:.1f}" font-size="11" text-anchor="end">wall →</text>')
    add('</svg>')

    with open(os.path.join(out, f"{name}.svg"), "w") as f:
        f.write("\n".join(svg))
    print(f"\nWrote {out}/{name}.svg and .csv\n")


def to_png(name, px=1800):
    """macOS-only: render the SVG to PNG with Quick Look."""
    out = os.path.join(OUT_DIR, name)
    subprocess.run(["qlmanage", "-t", "-s", str(px), "-o", out, os.path.join(out, f"{name}.svg")],
                   check=True, capture_output=True)
    os.replace(os.path.join(out, f"{name}.svg.png"), os.path.join(out, f"{name}.png"))


if __name__ == "__main__":
    for layout in LAYOUTS:
        make_layout(**layout)
        if "--png" in sys.argv:
            to_png(layout["name"])
