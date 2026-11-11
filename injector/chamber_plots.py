"""
PNG graphs of atomization, vaporization and mixing through the thrust chamber,
from the reduced-order model in chamber_mixing.py.

Outputs (in output/plots/):
  atomization.png         drop-size distributions, SMD vs element count and vs jet velocity
  drop_trajectories.png   side view of the spray: where drops go and how big they are
  vaporization.png        % vaporized and wall film along the chamber, with the contour
  mixing.png              Rupe mixing and c* efficiency along the chamber, O/F cross-sections
  sensitivity.png         how the answers move with drop size and turbulent diffusivity
  layout_comparison.png   12 vs 15 elements vs ox outside, side by side
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

import chamber_mixing as cm
import injector_layout

OUT = os.path.join("output", "plots")
DPI = 200
OX_C, FU_C = "#eb6834", "#2a78d6"                       # N2O, IPA (same hues as the renders)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]              # layouts, fixed order
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
OF_CMAP = LinearSegmentedColormap.from_list("of", ["#184f95", "#6da7ec", "#f0efec", "#ec8a87", "#a8201f"])
DESIGN = ("layout_12_elements", 8)
LAYOUTS = [("12 elements, 8 open", "layout_12_elements", 8),
           ("15 elements, 12 open", "layout_15_elements", 12),
           ("12 elements, ox outside", "layout_12_elements_OX_OUTSIDE", 8)]
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2, "legend.frameon": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def save(fig, name, note=None):
    if note:   # below everything; bbox_inches="tight" grows the canvas to fit it
        fig.text(0.01, -0.3 / fig.get_figheight(), note, fontsize=7.5, color=INK2, va="top", wrap=True)
    fig.savefig(f"{OUT}/{name}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def stations(ax, s, label=True):
    """Dotted lines at the end of the cylinder and the throat."""
    for z, lab in ((cm.CHAMBER_LEN, "end of cyl"), (s["zt"], "throat")):
        ax.axvline(z, color=INK2, lw=0.8, ls=":")
        if label:
            ax.text(z, 1.0, " " + lab, transform=ax.get_xaxis_transform(), fontsize=7.5, color=INK2, va="top")


def along_z(s, fn):
    """fn(sim, z) at every axial station, NaN before anything has vaporized."""
    out = np.full(len(s["zc"]), np.nan)
    for k, z in enumerate(s["zc"]):
        if s["flux_ox"][k].sum() + s["flux_fu"][k].sum() > 1e-6 * s["mdot"]:
            out[k] = fn(s, z)
    return out


def rr_pdf(d_um, smd):
    """Rosin-Rammler q=3 volume-fraction density per um, about this SMD."""
    X, q = 1.354 * smd, 3.0
    return q / X * (d_um / X) ** (q - 1) * np.exp(-(d_um / X) ** q)


base = cm.simulate(*DESIGN)
zc = base["zc"]

# ---------------------------------------------------------------- atomization
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), gridspec_kw=dict(wspace=0.32))
ax = axs[0]
d = np.linspace(1, 900, 600)
for t, (name, col) in enumerate((("N2O", OX_C), ("IPA", FU_C))):
    smd = base["smd"][t]
    ax.plot(d, rr_pdf(d, smd) * 100, color=col, label=f"{name}  SMD {smd:.0f} µm")
    ax.axvline(smd, color=col, lw=0.8, ls="--")
ax.set(title="Drop sizes at injection (design)", xlabel="drop diameter (µm)",
       ylabel="share of liquid volume (% per µm)", xlim=(0, 900), ylim=(0, None))
ax.legend()

ax = axs[1]
ns = np.arange(4, 21)
g = injector_layout.geometry(DESIGN[0])
a_ox, a_fu = g.CDA_OX / g.CD_OX, g.CDA_FU / g.CD_FU          # total geometric area, m^2
dia = lambda a: np.sqrt(4 * a / np.pi) * 1000
smds = np.array([cm.doublet_smd(dia(a_ox / n), base["vo"], dia(a_fu / n), base["vf"]) for n in ns])
ax.plot(ns, smds[:, 0] * cm.N2O_FLASH_SMD * cm.SMD_SCALE, color=OX_C, label="N2O")
ax.plot(ns, smds[:, 1] * cm.SMD_SCALE, color=FU_C, label="IPA")
for j, (lab, lay, n_open) in enumerate(LAYOUTS[:2]):
    n = int(lay.split("_")[1])
    ax.axvline(n_open, color=INK2, lw=0.8, ls=":")
    ax.text(n_open, 0.30 - 0.08 * j, f" {n}-el layout\n ({n_open} open)", transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=INK2, va="top")
ax.set(title="SMD vs number of open elements", xlabel="open doublets (same total area and ΔP)",
       ylabel="SMD (µm)", ylim=(0, None))
ax.legend()

ax = axs[2]
v = np.linspace(10, 80, 200)
d_ox, d_fu = base["d_ox"], base["d_fu"]
ax.axvspan(100 * 0.3048, 150 * 0.3048, color=GRID, alpha=0.7, lw=0)
ax.text(125 * 0.3048, 0.03, "correlation\nfit range", transform=ax.get_xaxis_transform(), ha="center",
        va="bottom", fontsize=7.5, color=INK2)
ax.plot(v, [cm.doublet_smd(d_ox, vv, d_fu, base["vf"])[0] * cm.N2O_FLASH_SMD * cm.SMD_SCALE for vv in v],
        color=OX_C, label="N2O, vs N2O jet velocity")
ax.plot(v, [cm.doublet_smd(d_ox, base["vo"], d_fu, vv)[1] * cm.SMD_SCALE for vv in v],
        color=FU_C, label="IPA, vs IPA jet velocity")
ax.plot(base["vo"], base["smd"][0], "o", ms=8, color=OX_C, mec="white", mew=2)
ax.plot(base["vf"], base["smd"][1], "o", ms=8, color=FU_C, mec="white", mew=2)
ax.annotate("design", (base["vo"], base["smd"][0]), xytext=(8, 6), textcoords="offset points", fontsize=8, color=INK2)
ax.annotate("design", (base["vf"], base["smd"][1]), xytext=(8, 6), textcoords="offset points", fontsize=8, color=INK2)
ax.set(title=f"SMD vs jet velocity ({d_ox:.2f} / {d_fu:.2f} mm holes)", xlabel="jet velocity (m/s, liquid density)",
       ylabel="SMD (µm)", ylim=(0, None))
ax.legend(loc="upper right")
save(fig, "atomization.png",
     f"Unlike-doublet correlation (Rocketdyne, NTRS 19700005245), ±50 %; N2O × {cm.N2O_FLASH_SMD} flash credit "
     f"(guess); drop size × {cm.SMD_SCALE}. Velocity held at the design ΔP when the element count changes.")

# ---------------------------------------------------------------- drop trajectories
pts = np.array(base["pts"])
sn = base["snap"]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True, gridspec_kw=dict(wspace=0.06))
rng = np.random.default_rng(0)
dmax = np.percentile(sn[:, 5], 99.5) * 1e6
for ax, t, name, cmap in ((axs[0], 0, "N2O", "Oranges"), (axs[1], 1, "IPA", "Blues")):
    m = sn[:, 4] == t
    pick = rng.choice(np.nonzero(m)[0], min(12000, m.sum()), replace=False)
    r = np.hypot(sn[pick, 0], sn[pick, 1]) * 1000
    dia = sn[pick, 5] * 1e6
    ramp = LinearSegmentedColormap.from_list(name, plt.get_cmap(cmap)(np.linspace(0.3, 1, 64)))   # skip near-white
    sc = ax.scatter(sn[pick, 2] * 1000, r, c=dia, s=1 + 10 * (dia / dmax) ** 2, cmap=ramp,
                    vmin=0, vmax=dmax, lw=0, alpha=0.6, rasterized=True)
    ax.fill_between(pts[:, 0], pts[:, 1], pts[:, 1].max() + 8, color="#d3d2cd", lw=0)
    ax.plot(pts[:, 0], pts[:, 1], color=INK, lw=1.5)
    film = base["film_z"][t] * 100
    ax.set(title=f"{name} drops (half-section, r vs depth)", xlabel="depth below injector face (mm)",
           xlim=(0, base["zt"] + 40), ylim=(0, pts[:, 1].max() + 8))
    stations(ax, base)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, shrink=0.85)
    cb.set_label(f"{name} drop diameter (µm)")
    ax.text(0.98, 0.05, f"hits wall {base['film_wall'][t] * 100:.0f}%\nfilm left at throat {base['film_throat'][t] * 100:.1f}%",
            transform=ax.transAxes, ha="right", fontsize=8, color=INK,
            bbox=dict(fc="white", ec=GRID, boxstyle="round,pad=0.3"))
axs[0].set_ylabel("radius (mm)")
save(fig, "drop_trajectories.png",
     "Snapshot of the steady spray, design config. Grey band is the chamber wall. "
     "No secondary breakup in the model yet, so large drops are probably overstated.")

# ---------------------------------------------------------------- vaporization
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.5, 6), sharex=True,
                              gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08))
for t, name, col, vap in ((0, "N2O", OX_C, base["vap_ox"]), (1, "IPA", FU_C, base["vap_fu"])):
    ax.plot(zc, vap * 100, color=col, label=f"{name} vaporized")
    ax.plot(zc, base["film_z"][t] * 100, color=col, ls="--", lw=1.5, label=f"{name} in wall film")
    ax.annotate(f"{name} {vap[base['it']] * 100:.0f}% at throat", (base["zt"], vap[base["it"]] * 100),
                xytext=(-6, -14) if t == 0 else (6, -16), textcoords="offset points", fontsize=8, color=INK,
                ha="right" if t == 0 else "left")
ax.set(ylabel="% of that propellant", ylim=(0, 105),
       title="Vaporization along the chamber (design config)")
stations(ax, base)
ax.legend(loc="center right")
ax2.fill_between(pts[:, 0], pts[:, 1], -pts[:, 1], color="#ecebe7", lw=0)
ax2.plot(pts[:, 0], pts[:, 1], color=INK, lw=1.2)
ax2.plot(pts[:, 0], -pts[:, 1], color=INK, lw=1.2)
rmax = pts[:, 1].max() + 5
ax2.set(xlabel="depth below injector face (mm)", ylabel="r (mm)", xlim=(0, base["z_exit"]), ylim=(-rmax, rmax))
stations(ax2, base, label=False)
save(fig, "vaporization.png",
     f"N2O flashes {base['chi'] * 100:.0f}% to vapor at the face. Drops: d²-law with convective correction; "
     f"wall film: ST {cm.ST_FILM}, gas {cm.T_FILM_GAS:.0f} K next to the film.")

# ---------------------------------------------------------------- mixing
em = along_z(base, cm.mixing_efficiency) * 100
eta = along_z(base, cm.cstar_efficiency) * 100
fig = plt.figure(figsize=(12, 7))
gs = fig.add_gridspec(2, 4, height_ratios=[1.1, 1], hspace=0.45, wspace=0.25)
ax = fig.add_subplot(gs[0, :])
ax.plot(zc, em, color=SERIES[0], label="Rupe mixing efficiency (vapor)")
ax.plot(zc, eta, color=SERIES[1], label="c* efficiency (vaporization × mixing)")
for y in (em, eta):
    ax.annotate(f"{y[base['it']]:.0f}% at throat", (base["zt"], y[base["it"]]), xytext=(-6, -14),
                textcoords="offset points", fontsize=8, color=INK, ha="right")
ax.set(xlabel="depth below injector face (mm)", ylabel="%", ylim=(0, 100), xlim=(0, base["zt"] + 20),
       title="Mixing along the chamber (design config)")
stations(ax, base)
ax.legend(loc="lower right")
xi = base["xi"]
norm = TwoSlopeNorm(vcenter=0, vmin=-1.2, vmax=1.2)
for j, (z, lab) in enumerate(((25, "25 mm"), (75, "75 mm"), (cm.CHAMBER_LEN, "end of cylinder"), (base["zt"], "throat"))):
    k = int(np.argmin(abs(zc - z)))
    fo, fu = base["flux_ox"][k], base["flux_fu"][k]
    tot = fo + fu
    of = np.log10(np.maximum(fo, 1e-12) / np.maximum(fu, 1e-12) / cm.OF_DESIGN)
    of = np.where((base["mask"] > 0) & (tot > 0.02 * tot.max()), of, np.nan)
    ax = fig.add_subplot(gs[1, j])
    im = ax.imshow(of.T, origin="lower", extent=(-1, 1, -1, 1), cmap=OF_CMAP, norm=norm, interpolation="nearest")
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, color=INK, lw=1))
    ax.set(title=f"{lab}\nØ{base['R'][k] * 2:.0f} mm", xticks=[], yticks=[])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
cb = fig.colorbar(im, ax=fig.axes[1:], shrink=0.8, pad=0.02, ticks=[-1, -0.5, 0, 0.5, 1])
cb.set_ticklabels([f"{cm.OF_DESIGN * 10 ** t:.1f}" for t in (-1, -0.5, 0, 0.5, 1)])
cb.set_label(f"local vapor O/F (design {cm.OF_DESIGN:.2f})\nblue fuel-rich · red ox-rich")
save(fig, "mixing.png",
     f"Turbulent diffusivity {cm.D_TURB:g} m²/s (assumed — mixing depends heavily on it). Cross-sections blank where "
     "little has vaporized. c* from rough N2O/IPA tables; no heat loss or kinetics.")

# ---------------------------------------------------------------- sensitivity
def run(**kw):
    old = {k: getattr(cm, k) for k in kw}
    for k, v in kw.items():
        setattr(cm, k, v)
    try:
        s = cm.simulate(*DESIGN)
    finally:
        for k, v in old.items():
            setattr(cm, k, v)
    it = s["it"]
    return dict(vap_fu=s["vap_fu"][it] * 100, film_fu=s["film_throat"][1] * 100,
                eta=cm.cstar_efficiency(s, s["zt"]) * 100, em=cm.mixing_efficiency(s, cm.CHAMBER_LEN) * 100)


scales = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
dturbs = [1e-3, 2.5e-3, 5e-3, 1e-2, 2e-2]
by_smd = [run(SMD_SCALE=x) for x in scales]
by_dt = [run(D_TURB=x) for x in dturbs]
fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True, gridspec_kw=dict(wspace=0.08))
for ax, xs, res, xlabel, cur in ((axs[0], scales, by_smd, "drop size × correlation", cm.SMD_SCALE),
                                 (axs[1], dturbs, by_dt, "turbulent diffusivity (m²/s)", cm.D_TURB)):
    for (key, lab), col in zip((("eta", "c* efficiency at throat"), ("vap_fu", "IPA vaporized at throat"),
                                ("em", "Rupe mixing, end of cyl")), SERIES):
        y = [r[key] for r in res]
        ax.plot(xs, y, color=col, marker="o", ms=5, mec="white", mew=1.5, label=lab)
    ax.plot(xs, [r["film_fu"] for r in res], color=INK2, ls="--", lw=1.5, marker="o", ms=4, label="IPA film unburnt at throat")
    ax.axvline(cur, color=INK2, lw=0.8, ls=":")
    ax.text(cur, 1.0, " current", transform=ax.get_xaxis_transform(), fontsize=7.5, color=INK2, va="top")
    ax.set(xlabel=xlabel, ylim=(0, 105))
    ax.margins(x=0.02)
axs[1].set_xscale("log")
axs[0].set(ylabel="%", title="Sensitivity to drop size")
axs[1].set_title("Sensitivity to turbulent diffusivity")
axs[0].legend(loc="center left", fontsize=8)
save(fig, "sensitivity.png", "Design config; one input changed at a time, everything else at its default.")

# ---------------------------------------------------------------- layout comparison
rows = []
for lab, lay, n_open in LAYOUTS:
    s = base if (lay, n_open) == DESIGN else cm.simulate(lay, n_open=n_open)
    it = s["it"]
    rows.append(dict(lab=lab, eta=cm.cstar_efficiency(s, s["zt"]) * 100, vap_fu=s["vap_fu"][it] * 100,
                     em=cm.mixing_efficiency(s, cm.CHAMBER_LEN) * 100, film=s["film_throat"][1] * 100,
                     smd=s["smd"]))
fig, axs = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw=dict(width_ratios=[2.2, 1], wspace=0.2))
ax = axs[0]
metrics = (("eta", "c* eff\nat throat"), ("vap_fu", "IPA vaporized\nat throat"), ("em", "Rupe mixing\nend of cyl"),
           ("film", "IPA film unburnt\nat throat"))
w = 0.26
for i, (row, col) in enumerate(zip(rows, SERIES)):
    xs = np.arange(len(metrics)) + (i - 1) * (w + 0.02)
    vals = [row[k] for k, _ in metrics]
    bars = ax.bar(xs, vals, w, color=col, label=row["lab"])
    ax.bar_label(bars, fmt="%.0f", fontsize=7.5, color=INK, padding=2)
ax.set(xticks=np.arange(len(metrics)), xticklabels=[m for _, m in metrics], ylabel="%", ylim=(0, 110),
       title="Layouts compared (each at its design open count)")
ax.legend(loc="upper right", fontsize=8)
ax.grid(axis="x", visible=False)
ax = axs[1]
for i, (row, col) in enumerate(zip(rows, SERIES)):
    xs = np.arange(2) + (i - 1) * (w + 0.02)
    bars = ax.bar(xs, row["smd"], w, color=col)
    ax.bar_label(bars, fmt="%.0f", fontsize=7.5, color=INK, padding=2)
ax.set(xticks=[0, 1], xticklabels=["N2O", "IPA"], ylabel="SMD (µm)", title="Drop size (SMD)", ylim=(0, None))
ax.grid(axis="x", visible=False)
save(fig, "layout_comparison.png", "Same flows and ΔP for every layout. Ox outside uses the 12-element holes with ox and fuel swapped.")

print("Wrote " + ", ".join(f"{OUT}/{f}" for f in sorted(os.listdir(OUT))))
