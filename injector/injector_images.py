"""
PNG renders of the injector, using the geometry from injector_layout.py.

Outputs (in output/renders/):
  injector_face.png     face view, 8-open config, with spray fans
  plug_patterns.png     which elements to plug for 6/8/9/10/12 open
  element_section.png   one element cut along its radial plane
  injector_3d.png       3D view of all jets meeting on the impingement circle
  chamber_3d.png        3D view of the whole thrust chamber and nozzle
  chamber_mixing_3d.png local O/F through the chamber (model in chamber_mixing.py)
  chamber_mixing_ox_outside_3d.png  same, with the ox and fuel holes swapped
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

import chamber_mixing
import injector_layout

L = injector_layout.geometry("layout_12_elements")

OUT = os.path.join("output", "renders")
OX_C, FU_C, PLUG_C, IMP_C = L.OX_C, L.FU_C, L.PLUG_C, "#2b8a3e"
DPI = 200
os.makedirs(OUT, exist_ok=True)


def elem_angle(i):
    return math.radians(90 - i * L.pitch)   # element 0 at 12 o'clock, clockwise


def draw_face(ax, plugged, labels=True, fans=True):
    ax.add_patch(Circle((0, 0), L.R_CH, fc="#f1f3f5", ec="black", lw=2))
    ax.add_patch(Circle((0, 0), L.R_IMPINGE, fill=False, ec=IMP_C, ls="--", lw=1))
    for r in (L.r_ox, L.r_fu):
        ax.add_patch(Circle((0, 0), r, fill=False, ec="#ced4da", lw=0.6))
    for i in range(L.N_DRILLED):
        th = elem_angle(i)
        c, s = math.cos(th), math.sin(th)
        off = i in plugged
        if fans and not off:
            # spray fan spreads tangentially from the impingement point
            px, py = L.R_IMPINGE * c, L.R_IMPINGE * s
            fan_w = L.R_IMPINGE * math.radians(L.pitch) * 1.1
            ax.add_patch(Ellipse((px, py), fan_w, 3.5, angle=math.degrees(th) + 90,
                                 fc=IMP_C, alpha=0.25, ec="none"))
            ax.plot(px, py, "o", color=IMP_C, ms=3)
        for r, d, thr, col, sgn in ((L.r_ox, L.D_OX, L.THREAD_OX, OX_C, 1),
                                    (L.r_fu, L.D_FU, L.THREAD_FU, FU_C, -1)):
            x, y = r * c, r * s
            colr = PLUG_C if off else col
            ax.add_patch(Circle((x, y), thr / 2, fill=False, ec=colr, ls=(0, (2, 1.5)), lw=0.8))
            ax.add_patch(Circle((x, y), d / 2, fc=colr, ec="none"))
            if off:
                k = thr / 2 * 0.75
                ax.plot([x - k, x + k], [y - k, y + k], color=PLUG_C, lw=1.5)
                ax.plot([x - k, x + k], [y + k, y - k], color=PLUG_C, lw=1.5)
            else:
                ax.annotate("", xy=(x + sgn * 4 * c, y + sgn * 4 * s), xytext=(x, y),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2, mutation_scale=8))
        if labels:
            ax.text((L.r_fu + 8) * c, (L.r_fu + 8) * s, str(i), ha="center", va="center",
                    fontweight="bold", fontsize=10)
    ax.plot([-3, 3], [0, 0], "k", lw=0.8)
    ax.plot([0, 0], [-3, 3], "k", lw=0.8)
    lim = L.R_CH + 4
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")


# ---------- 1. face view ----------
fig, ax = plt.subplots(figsize=(8, 8))
draw_face(ax, L.PLUG_SETS[8])
ax.set_title("Injector face, looking upstream from chamber\n"
             f"{L.N_DRILLED} drilled, {L.N_DESIGN} open (plug 0, 3, 6, 9)", fontsize=13)
handles = [
    plt.Line2D([], [], marker="o", ls="", color=OX_C, ms=9, label=f"N2O  Ø{L.D_OX} mm, {L.ANG_OX}° outward"),
    plt.Line2D([], [], marker="o", ls="", color=FU_C, ms=5, label=f"IPA  Ø{L.D_FU} mm, {L.ANG_FU}° inward"),
    plt.Line2D([], [], marker="X", ls="", color=PLUG_C, ms=8, label="plugged pair"),
    plt.Line2D([], [], ls="--", color=IMP_C, label=f"impingement circle R{L.R_IMPINGE:g}, {L.IMPINGE_L:g} mm below face"),
    Polygon([[0, 0]], fc=IMP_C, alpha=0.25, label="spray fan (tangential)"),
]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=9, frameon=False)
ax.text(0, -L.R_CH + 3, f"chamber Ø{L.CHAMBER_DIA_MM:.1f} mm", ha="center", fontsize=9, color="#495057")
fig.savefig(f"{OUT}/injector_face.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- 2. plug patterns ----------
counts = sorted(L.PLUG_SETS, reverse=True)
fig, axes = plt.subplots(1, len(counts), figsize=(3.2 * len(counts), 3.9))
for ax, n in zip(axes, counts):
    draw_face(ax, L.PLUG_SETS[n], labels=True, fans=False)
    ox = n * L.area(L.D_OX) * L.CD_OX / L.CDA_OX
    fu = n * L.area(L.D_FU) * L.CD_FU / L.CDA_FU
    plug = ", ".join(map(str, L.PLUG_SETS[n])) or "none"
    ax.set_title(f"{n} open\nplug: {plug}", fontsize=11, fontweight="bold")
    ax.text(0, -L.R_CH - 9, f"ox CdA {ox*100:.0f}%   fuel CdA {fu*100:.0f}%",
            ha="center", fontsize=9)
fig.suptitle("Plug patterns (CdA as % of HalfCatSim design point)", fontsize=13, y=1.02)
fig.savefig(f"{OUT}/plug_patterns.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- 3. element section ----------
fig, ax = plt.subplots(figsize=(8, 7))
T, Z = L.FACE_PLATE_T, L.IMPINGE_L
tox, tfu = math.tan(math.radians(L.ANG_OX)), math.tan(math.radians(L.ANG_FU))
x_lo, x_hi = L.r_ox - 10, L.R_CH
# face plate (z from +T at manifold side to 0 at chamber face), chamber below
ax.add_patch(Rectangle((x_lo, 0), x_hi - x_lo, T, fc="#dee2e6", ec="black", lw=1.2))
ax.add_patch(Rectangle((L.R_CH, -Z - 14), 3, Z + 14 + T, fc="#495057", ec="black"))
ax.text(L.R_CH + 1.5, -Z - 16, "chamber\nwall", ha="center", va="top", fontsize=9)
ax.text(48, T / 2, f"face plate {T:g} mm", ha="center", va="center", fontsize=9)


def bore(x_face, ang, d, col, sgn):
    # hole through the plate, then free jet down to the impingement point
    t = math.tan(math.radians(ang))
    hw = d / 2 / math.cos(math.radians(ang))
    x_top = x_face - sgn * T * t
    ax.add_patch(Polygon([(x_top - hw, T), (x_top + hw, T), (x_face + hw, 0), (x_face - hw, 0)],
                         fc=col, alpha=0.35, ec=col, lw=1))
    ax.plot([x_face, L.R_IMPINGE], [0, -Z], color=col, lw=max(1.5, d * 2.2), alpha=0.75,
            solid_capstyle="butt")


bore(L.r_ox, L.ANG_OX, L.D_OX, OX_C, 1)
bore(L.r_fu, L.ANG_FU, L.D_FU, FU_C, -1)
ax.plot([L.R_IMPINGE, L.R_IMPINGE], [-Z, -Z - 14], color=IMP_C, lw=2, ls="--")
ax.annotate("", xy=(L.R_IMPINGE, -Z - 14), xytext=(L.R_IMPINGE, -Z - 10),
            arrowprops=dict(arrowstyle="-|>", color=IMP_C, lw=2))
ax.plot(L.R_IMPINGE, -Z, "o", color=IMP_C, ms=7)
ax.text(L.R_IMPINGE + 1.2, -Z - 0.5, f"impingement\nR{L.R_IMPINGE:g}, {Z:g} mm below face",
        color=IMP_C, fontsize=9, va="top")
ax.text(L.r_ox, T + 1, f"N2O Ø{L.D_OX}\n{L.ANG_OX}° outward", color=OX_C, ha="center", fontsize=9)
ax.text(L.r_fu, T + 1, f"IPA Ø{L.D_FU}\n{L.ANG_FU}° inward", color=FU_C, ha="center", fontsize=9)
ax.axhline(0, color="black", lw=0.6)
ax.text(x_lo + 1, -1.5, "injector face", fontsize=8, color="#495057", va="top")
ax.text(x_lo, -Z - 16, "← chamber axis", fontsize=9, va="top")
ax.set_xlim(x_lo - 1, L.R_CH + 6)
ax.set_ylim(-Z - 22, T + 7)
ax.set_aspect("equal")
ax.set_xlabel("radius from chamber axis (mm)")
ax.set_ylabel("height above face (mm)")
ax.set_title("Element section (radial plane)", fontsize=13)
ax.grid(alpha=0.25)
fig.savefig(f"{OUT}/element_section.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- 4. 3D jets ----------
fig = plt.figure(figsize=(9, 8))
ax = fig.add_subplot(projection="3d")
u = [2 * math.pi * k / 120 for k in range(121)]
ax.plot([L.R_CH * math.cos(a) for a in u], [L.R_CH * math.sin(a) for a in u], 0, color="black", lw=1.5)
ax.plot([L.R_CH * math.cos(a) for a in u], [L.R_CH * math.sin(a) for a in u], -45, color="#adb5bd", lw=0.8)
for a in u[::10]:
    ax.plot([L.R_CH * math.cos(a)] * 2, [L.R_CH * math.sin(a)] * 2, [0, -45], color="#dee2e6", lw=0.6)
ax.plot([L.R_IMPINGE * math.cos(a) for a in u], [L.R_IMPINGE * math.sin(a) for a in u], -Z,
        color=IMP_C, ls="--", lw=1)
for i in range(L.N_DRILLED):
    th = elem_angle(i)
    c, s = math.cos(th), math.sin(th)
    off = i in L.PLUG_SETS[8]
    px, py = L.R_IMPINGE * c, L.R_IMPINGE * s
    for r, col, w in ((L.r_ox, OX_C, 3), (L.r_fu, FU_C, 1.6)):
        ax.scatter(r * c, r * s, 0, color=PLUG_C if off else col, s=w * 12, depthshade=False)
        if not off:
            ax.plot([r * c, px], [r * s, py], [0, -Z], color=col, lw=w)
    if not off:
        ax.plot([px, px], [py, py], [-Z, -Z - 22], color=IMP_C, lw=1.5, alpha=0.7)
        # tangential fan as a thin translucent triangle
        tx, ty = -s, c
        sp = 9
        fan = [(px, py, -Z), (px + tx * sp, py + ty * sp, -Z - 22), (px - tx * sp, py - ty * sp, -Z - 22)]
        ax.add_collection3d(Poly3DCollection([fan], color=IMP_C, alpha=0.15))
ax.set_xlim(-L.R_CH, L.R_CH)
ax.set_ylim(-L.R_CH, L.R_CH)
ax.set_zlim(-45, 5)
ax.set_box_aspect((1, 1, 50 / (2 * L.R_CH)))
ax.view_init(elev=22, azim=-60)
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_zlabel("z (mm)")
ax.set_title("Jets and spray fans, 8-open config\n(orange N2O, blue IPA, grey plugged)", fontsize=13)
fig.savefig(f"{OUT}/injector_3d.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- 5. whole chamber 3D ----------
def draw_chamber(ax, title, g=L, fans=True, wall_alpha=0.18):
    """Wall, injector jets and labels in the injector_3d.png style; returns contour points and throat depth."""
    pts, zt = chamber_mixing.chamber_contour(g.R_CH)
    zs, rs = [p[0] for p in pts], [p[1] for p in pts]
    z_exit = zs[-1]
    # back half of the wall as a translucent surface, so the inside stays visible
    az0 = math.radians(-60)
    back = [az0 + math.pi / 2 + math.pi * k / 40 for k in range(41)]
    quads = []
    for j in range(len(pts) - 1):
        for k in range(len(back) - 1):
            quads.append([(rs[j + jj] * math.cos(back[k + kk]), rs[j + jj] * math.sin(back[k + kk]), -zs[j + jj])
                          for jj, kk in ((0, 0), (0, 1), (1, 1), (1, 0))])
    ax.add_collection3d(Poly3DCollection(quads, color="#ced4da", alpha=wall_alpha, lw=0))
    # wireframe: profile lines all round, rings at the key stations
    for a in u[::10]:
        ax.plot([r * math.cos(a) for r in rs], [r * math.sin(a) for r in rs], [-z for z in zs],
                color="#868e96", lw=0.6)
    for z, col, lw in ((0, "black", 1.5), (chamber_mixing.CHAMBER_LEN, "#495057", 1), (zt, "#c92a2a", 1.5),
                       (z_exit, "black", 1.2)):
        r = rs[min(range(len(zs)), key=lambda i: abs(zs[i] - z))]
        ax.plot([r * math.cos(a) for a in u], [r * math.sin(a) for a in u], -z, color=col, lw=lw)
    ax.plot([0, 0], [0, 0], [0, -z_exit], color="#adb5bd", lw=0.8, ls="-.")
    ax.plot([g.R_IMPINGE * math.cos(a) for a in u], [g.R_IMPINGE * math.sin(a) for a in u], -Z,
            color=IMP_C, ls="--", lw=1)
    for i in range(g.N_DRILLED):
        th = elem_angle(i)
        c, s = math.cos(th), math.sin(th)
        off = i in g.PLUG_SETS[8]
        px, py = g.R_IMPINGE * c, g.R_IMPINGE * s
        for r, col, w in ((g.r_ox, OX_C, 2.5), (g.r_fu, FU_C, 1.3)):
            ax.scatter(r * c, r * s, 0, color=PLUG_C if off else col, s=w * 6, depthshade=False)
            if not off:
                ax.plot([r * c, px], [r * s, py], [0, -Z], color=col, lw=w)
        if fans and not off:
            tx, ty, sp = -s, c, 9
            fan = [(px, py, -Z), (px + tx * sp, py + ty * sp, -Z - 30), (px - tx * sp, py - ty * sp, -Z - 30)]
            ax.add_collection3d(Poly3DCollection([fan], color=IMP_C, alpha=0.2))
    lab = dict(fontsize=9, color="#212529")
    ax.text(g.R_CH + 4, 0, 0, f"injector face  Ø{2 * g.R_CH:.1f}", **lab)
    ax.text(g.R_CH + 4, 0, -chamber_mixing.CHAMBER_LEN, f"end of cylinder  {chamber_mixing.CHAMBER_LEN:.0f} mm", **lab)
    ax.text(chamber_mixing.THROAT_R + 6, 0, -zt, f"throat  Ø{2 * chamber_mixing.THROAT_R:.1f} @ {zt:.0f} mm",
            color="#c92a2a", fontsize=9)
    ax.text(rs[-1] + 4, 0, -z_exit, f"exit  Ø{2 * rs[-1]:.1f} @ {z_exit:.0f} mm", **lab)
    lim = g.R_CH + 5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-z_exit - 5, 5)
    ax.set_box_aspect((1, 1, (z_exit + 10) / (2 * lim)))
    ax.view_init(elev=18, azim=-60)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)", labelpad=14)
    ax.set_title(title, fontsize=13)
    return pts, zt


fig = plt.figure(figsize=(8, 11))
ax = fig.add_subplot(projection="3d")
draw_chamber(ax, "Thrust chamber, 8-open config\n(orange N2O, blue IPA, grey plugged, red ring = throat)")
fig.savefig(f"{OUT}/chamber_3d.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- 6. mixing through the chamber ----------
OF_CMAP = LinearSegmentedColormap.from_list("of", [FU_C, "#a5d8ff", "#f8f9fa", "#ffc078", OX_C])
OF_SPAN = math.log(4)   # colour scale runs from design O/F / 4 to x 4


def of_rgba(fo, fu, area_ref):
    """Colour by log(O/F / design); fade out where little vapor has arrived yet."""
    tot = fo + fu
    of = fo / np.maximum(fu, 1e-12)
    t = np.clip(0.5 + 0.5 * np.log(np.maximum(of, 1e-6) / chamber_mixing.OF_DESIGN) / OF_SPAN, 0, 1)
    rgba = OF_CMAP(t)
    rgba[..., 3] = 0.92 * np.clip(tot / (0.35 * area_ref), 0, 1)
    return rgba


def mixing_render(layout_name, filename, title):
    g = injector_layout.geometry(layout_name)
    M = chamber_mixing.simulate(layout_name, n_open=8)

    def field_at(iz, xi, yi):
        n = chamber_mixing.NXY
        ix = np.clip(((xi + 1) / 2 * n).astype(int), 0, n - 1)
        iy = np.clip(((yi + 1) / 2 * n).astype(int), 0, n - 1)
        ref = M["mdot"] / M["mask"].sum()   # flux per pixel if fully vaporized and evenly mixed
        return of_rgba(M["flux_ox"][iz][ix, iy], M["flux_fu"][iz][ix, iy], ref)

    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(projection="3d")
    pts, zt = draw_chamber(ax, title + "\ncolour = local vapor O/F, faded where little has vaporized yet", g=g,
                           fans=False, wall_alpha=0.1)
    zc, nz = M["zc"], chamber_mixing.NZ
    # cut plane through the axis, square to the viewer
    phi = math.radians(-60 + 90)
    zi = np.arange(0, nz, 1)
    s = np.linspace(-1, 1, 61)
    S_, ZI = np.meshgrid(s, zi)
    Rz = np.interp(zc[ZI], [p[0] for p in pts], [p[1] for p in pts])
    cols = np.array([[field_at(k, np.array([v * math.cos(phi)]), np.array([v * math.sin(phi)]))[0]
                      for v in (s[:-1] + s[1:]) / 2] for k in zi[:-1]])
    ax.plot_surface(S_ * Rz * math.cos(phi), S_ * Rz * math.sin(phi), -zc[ZI], facecolors=cols,
                    shade=False, rstride=1, cstride=1, lw=0, antialiased=False)
    # cross-section slices
    for z_mm in (25, 75, 152, zt):
        k = int(np.argmin(abs(zc - z_mm)))
        Rk = np.interp(z_mm, [p[0] for p in pts], [p[1] for p in pts])
        rr, tt = np.meshgrid(np.linspace(0, 1, 25), np.linspace(0, 2 * math.pi, 97))
        rm, tm = (rr[:-1, :-1] + rr[1:, 1:]) / 2, (tt[:-1, :-1] + tt[1:, 1:]) / 2
        cols = field_at(k, rm * np.cos(tm), rm * np.sin(tm))
        cols[..., 3] = np.maximum(cols[..., 3], 0.15)
        ax.plot_surface(rr * Rk * np.cos(tt), rr * Rk * np.sin(tt), np.full_like(rr, -z_mm), facecolors=cols,
                        shade=False, rstride=1, cstride=1, lw=0, antialiased=False)
        if z_mm < zt:
            ax.text(-Rk * 0.7, -Rk - 6, -z_mm, f"{z_mm:.0f} mm", fontsize=8, color="#495057")
    # a light dusting of droplets: an instantaneous snapshot of the steady spray
    rng = np.random.default_rng(0)
    sn = M["snap"][rng.choice(len(M["snap"]), 2500, replace=False)]
    for t, col in ((0, OX_C), (1, FU_C)):
        m = sn[:, 4] == t
        ax.scatter(sn[m, 0] * 1000, sn[m, 1] * 1000, -sn[m, 2] * 1000, s=1.5 + 4 * sn[m, 3], color=col,
                   alpha=0.35, lw=0, depthshade=False)
    sm = plt.cm.ScalarMappable(cmap=OF_CMAP, norm=plt.Normalize(-1, 1))
    cb = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.02)
    ticks = [0.25, 0.5, 1, 2, 4]
    cb.set_ticks([math.log(t) / OF_SPAN for t in ticks])
    cb.set_ticklabels([f"{chamber_mixing.OF_DESIGN * t:.1f}" for t in ticks])
    cb.set_label(f"local O/F  (design {chamber_mixing.OF_DESIGN:.2f}; blue fuel-rich, orange ox-rich)")
    note = "Rupe mixing efficiency  " + "   ".join(
        f"{lab} {chamber_mixing.mixing_efficiency(M, z) * 100:.0f}%"
        for lab, z in (("75 mm", 75), ("end of cyl", chamber_mixing.CHAMBER_LEN), ("throat", zt)))
    note += (f"\nSMD N2O {M['smd'][0]:.0f} / IPA {M['smd'][1]:.0f} µm (doublet correlation)   "
             f"wall film unburnt at throat: IPA {M['film_throat'][1] * 100:.1f}%")
    fig.text(0.5, 0.13, note + f"\nturbulent diffusivity {chamber_mixing.D_TURB:g} m²/s (assumed; it sets how fast the streaks blur)",
             ha="center", fontsize=10, color="#343a40")
    fig.savefig(f"{OUT}/{filename}", dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return M


mixing_render("layout_12_elements", "chamber_mixing_3d.png",
              "Mixing through the chamber, 8-open config, ox inside")
mixing_render("layout_12_elements_OX_OUTSIDE", "chamber_mixing_ox_outside_3d.png",
              "Mixing through the chamber, 8-open config, ox and fuel swapped (ox outside)")

print("Wrote " + ", ".join(sorted(f"{OUT}/{f}" for f in os.listdir(OUT))))
