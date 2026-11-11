"""
Thrust chamber contour and a steady-state spray/mixing model, for the 3D renders
in injector_images.py and the graphs in chamber_plots.py.

Reduced-order, not CFD:
  - each open doublet sprays droplet parcels from its impingement point:
    momentum-weighted direction, tangential fan, and each droplet partly keeps
    its own jet's heading (more so for the jet with the higher rho v^2 d, and more
    overall the further the Rupe ratio is from 1)
  - drop sizes from an unlike-doublet correlation (hole sizes and jet velocities),
    Rosin-Rammler about that, Schiller-Naumann drag, d^2-law vaporization with a
    convective correction, random turbulent kicks; no secondary breakup yet
  - N2O flashes on injection: equilibrium vapor fraction at Pc goes straight to gas
  - liquid that reaches the wall runs downstream as a film and evaporates into
    the wall cells at a mass-transfer rate set by the local gas mass flux;
    film still liquid at the throat doesn't burn
  - gas follows the wall contour (r/R = const, no recirculation); axial velocity
    from how much has vaporized, isentropic in the nozzle; each propellant's
    vapor diffuses sideways with a turbulent diffusivity
"""
import math
import os

import numpy as np

import injector_layout

# ---------------- INPUTS ----------------
# Thrust chamber ('Motor Simulation' C15, C17, C19, C28) and bell nozzle
CHAMBER_LEN = 6.0 * 25.4       # mm, cylindrical section only
THROAT_R = 1.5 * 25.4 / 2
EXIT_R = 3.0 * 25.4 / 2        # only used for the 15 deg cone fallback
CONV_HALF_ANGLE = 30.0         # not in the spreadsheet; assumed
DIV_HALF_ANGLE = 15.0
NOZZLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "nozzle_xyz.txt")  # bell, inches from throat

# Design-point flows, same as doublet_sizing.py ('Motor Simulation' G29/K29)
MDOT_OX = 3.2396 * 0.45359237  # kg/s
MDOT_FU = 1.0181 * 0.45359237
RHO_OX, RHO_FU = 726.59, 786.0
N2O_TEMP_K = (80 - 32) / 1.8 + 273.15
PC_PSI = 316.9                 # O3

# Combustion gas (rough N2O/IPA values -- replace with CEA output if you have it)
GAMMA, MW, T0, MU_GAS = 1.22, 26.0, 3000.0, 9e-5
# local mixture fraction f = fuel / (ox + fuel)  ->  T (K), c* (m/s)
F_TABLE = [0.0, 0.05, 0.0909, 0.1316, 0.1667, 0.2, 0.25, 0.333, 0.5, 1.0]
T_TABLE = [1100, 1900, 2700, 3150, 3200, 3100, 2850, 2300, 1500, 600]
CSTAR_TABLE = [900, 1250, 1450, 1550, 1590, 1590, 1560, 1450, 1200, 600]

# Atomization / vaporization / mixing -- the big unknowns
SMD_SCALE = 1.0                # x the doublet_smd() correlation, both propellants
N2O_FLASH_SMD = 0.5            # x correlation for N2O, credit for flash atomization (guess, 0.3-1)
K_FU, K_OX = 2.0e-6, 4.0e-6    # m^2/s burning-rate constants
D_TURB = 5e-3                  # m^2/s turbulent diffusivity
FAN_SIGMA_DEG = 30.0
UNMIX = 0.25
RUPE_PENALTY = 1.0             # extra unmixing per unit |ln(Rupe ratio)|; calibrate with a patternator

# Wall film evaporates at m'' = ST_FILM * G * ln(1 + B): G local gas mass flux,
# B = cp (T_FILM_GAS - T_boil) / h_fg (Spalding transfer number)
ST_FILM = 0.003                # mass-transfer Stanton number, turbulent wall layer
T_FILM_GAS = 2000.0            # K, gas next to the film (fuel-rich wall layer)
IPA_TB, IPA_HFG = 457.0, 440e3 # K, J/kg: IPA saturation at ~Pc (Tr 0.9, Watson-scaled h_fg)

PARCELS = 600                  # per open element per propellant
NZ, NXY = 150, 72              # vapor grid: axial bins x (x/R, y/R) pixels
NA = 72                        # azimuthal strips for the wall film
# ----------------------------------------

OF_DESIGN = MDOT_OX / MDOT_FU


def chamber_contour(r_ch):
    """(depth below face, radius) mm: cylinder, 30 deg cone with a 1.5 Rt throat arc, then the
    bell from NOZZLE_FILE (15 deg cone to EXIT_R if it's missing)."""
    a, rcu = math.radians(CONV_HALF_ANGLE), 1.5 * THROAT_R
    r_arc = THROAT_R + rcu * (1 - math.cos(a))
    zt = CHAMBER_LEN + (r_ch - r_arc) / math.tan(a) + rcu * math.sin(a)
    pts = [(0, r_ch), (CHAMBER_LEN, r_ch)]
    pts += [(zt - rcu * math.sin(a * (1 - k / 12)), THROAT_R + rcu * (1 - math.cos(a * (1 - k / 12))))
            for k in range(13)]
    if os.path.exists(NOZZLE_FILE):
        with open(NOZZLE_FILE) as f:
            bell = [tuple(float(v) * 25.4 for v in line.split()[:2]) for line in f if len(line.split()) >= 2]
        pts += [(zt + bz, br) for bz, br in bell if bz > 0]
    else:
        pts.append((zt + (EXIT_R - THROAT_R) / math.tan(math.radians(DIV_HALF_ANGLE)), EXIT_R))
    return pts, zt


def flash_fraction(pc):
    """Equilibrium N2O vapor fraction flashing from N2O_TEMP_K down to pc (Pa), ESDU 91022 fits.
    Also returns the saturation temperature (K) and latent heat (J/kg) at pc."""
    tc = 309.57
    psat = lambda T: 7251e3 * math.exp(tc / T * sum(
        b * (1 - T / tc) ** e for b, e in ((-6.71893, 1), (1.35966, 1.5), (-1.3779, 2.5), (-4.051, 5))))
    hl = lambda T: -200 + sum(c * (1 - T / tc) ** e for c, e in
                              ((116.043, 1 / 3), (-917.225, 2 / 3), (794.779, 1), (-589.587, 4 / 3)))
    hv = lambda T: -200 + sum(c * (1 - T / tc) ** e for c, e in
                              ((440.055, 1 / 3), (-459.701, 2 / 3), (434.081, 1), (-485.338, 4 / 3)))
    lo, hi = 182.0, N2O_TEMP_K
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (lo, m) if psat(m) > pc else (m, hi)
    return (hl(N2O_TEMP_K) - hl(lo)) / (hv(lo) - hl(lo)), lo, (hv(lo) - hl(lo)) * 1e3


def doublet_smd(d_ox, v_ox, d_fu, v_fu):
    """Sauter mean diameters (um) of the ox and fuel drops from one unlike doublet.

    Mass median from the unlike-doublet correlation quoted in 'Correlation of
    Spray Dropsize Distribution and Injector Variables' (Rocketdyne, NTRS
    19700005245, p. 50):  D_f = 1.0e5 d_f^0.27 d_o^0.023 / (v_f^0.74 v_o^0.33),
    um, d in inches, v in ft/s; the ox side swaps subscripts. Fitted to cold
    simulants at 1 atm, ~0.063 in holes, 100-150 ft/s, so treat as +-50 %.
    Rosin-Rammler q=3 has mass median = 1.198 SMD.
    """
    di_o, di_f, fv_o, fv_f = d_ox / 25.4, d_fu / 25.4, v_ox / 0.3048, v_fu / 0.3048
    mm_f = 1.0e5 * di_f ** 0.27 * di_o ** 0.023 / (fv_f ** 0.74 * fv_o ** 0.33)
    mm_o = 1.0e5 * di_o ** 0.27 * di_f ** 0.023 / (fv_o ** 0.74 * fv_f ** 0.33)
    return mm_o / 1.198, mm_f / 1.198


def area_mach(ar, sup):
    """Mach number for area ratios ar (array), subsonic or supersonic branch."""
    g = GAMMA
    f = lambda M: (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1)))
    lo = np.where(sup, 1.0, 1e-4)
    hi = np.where(sup, 8.0, 1.0)
    for _ in range(60):
        m = (lo + hi) / 2
        go_up = np.where(sup, f(m) < ar, f(m) > ar)
        lo, hi = np.where(go_up, m, lo), np.where(go_up, hi, m)
    return np.where(ar <= 1.0, 1.0, (lo + hi) / 2)


def simulate(layout_name="layout_12_elements", n_open=8, seed=1):
    g = injector_layout.geometry(layout_name)
    rng = np.random.default_rng(seed)
    pts, zt_mm = chamber_contour(g.R_CH)
    cz = np.array([p[0] for p in pts]) / 1000
    cr = np.array([p[1] for p in pts]) / 1000
    zt, z_exit = zt_mm / 1000, cz[-1]
    R = lambda z: np.interp(z, cz, cr)
    dRdz = lambda z: (R(z + 1e-4) - R(z - 1e-4)) / 2e-4

    # gas grid
    dz = z_exit / NZ
    zc = (np.arange(NZ) + 0.5) * dz
    Rg = R(zc)
    it = int(np.argmin(Rg))
    mach = area_mach((Rg / (THROAT_R / 1000)) ** 2, np.arange(NZ) > it)
    pc = PC_PSI * 6894.757
    rgas = 8314.46 / MW
    a0, rho0 = math.sqrt(GAMMA * rgas * T0), pc / (rgas * T0)
    tfac = 1 + (GAMMA - 1) / 2 * mach ** 2
    u_isen = mach * a0 / np.sqrt(tfac)
    rho_g = rho0 * tfac ** (-1 / (GAMMA - 1))
    chi, t_sat_ox, hfg_ox = flash_fraction(pc)
    cp = GAMMA * rgas / (GAMMA - 1)
    B = np.array([cp * (T_FILM_GAS - t_sat_ox) / hfg_ox, cp * (T_FILM_GAS - IPA_TB) / IPA_HFG])

    # injector
    plugged = set(g.PLUG_SETS[n_open])
    a_ox, a_fu = g.area(g.D_OX), g.area(g.D_FU)
    mo = MDOT_OX * g.CD_OX * a_ox / g.CDA_OX     # per hole, injector dP held at design
    mf = MDOT_FU * g.CD_FU * a_fu / g.CDA_FU
    vo, vf = mo / (RHO_OX * a_ox), mf / (RHO_FU * a_fu)
    rupe = (RHO_OX * vo ** 2 * g.D_OX) / (RHO_FU * vf ** 2 * g.D_FU)
    # Rupe: mixing is best at rho v^2 d ratio = 1, so total unmixing grows away from it
    unmix = UNMIX * (1 + RUPE_PENALTY * abs(math.log(rupe)))
    pen = (min(1.0, 2 * unmix * rupe / (1 + rupe)), min(1.0, 2 * unmix / (1 + rupe)))
    smd_ox, smd_fu = doublet_smd(g.D_OX, vo, g.D_FU, vf)
    smd_ox, smd_fu = smd_ox * N2O_FLASH_SMD * SMD_SCALE, smd_fu * SMD_SCALE
    mdot_ox, mdot_fu = 0.0, 0.0
    P0, V0, D0, M0, TY = [], [], [], [], []
    flash = []
    for i in range(g.N_DRILLED):
        if i in plugged:
            continue
        th = math.radians(90 - i * g.pitch)
        er = np.array([math.cos(th), math.sin(th), 0.0])
        et = np.array([-math.sin(th), math.cos(th), 0.0])
        P = np.array([g.R_IMPINGE * er[0], g.R_IMPINGE * er[1], g.IMPINGE_L]) / 1000
        own = []
        for r_h in (g.r_ox, g.r_fu):
            d = P - np.array([r_h * er[0], r_h * er[1], 0]) / 1000
            own.append(d / np.linalg.norm(d))
        Mv = mo * vo * own[0] + mf * vf * own[1]
        Rd = Mv / np.linalg.norm(Mv)
        vs = 0.8 * (mo * vo + mf * vf) / (mo + mf)
        flash.append((P, mo * chi))
        mdot_ox += mo
        mdot_fu += mf
        for t, (m_liq, smd) in enumerate(((mo * (1 - chi), smd_ox), (mf, smd_fu))):
            n = PARCELS
            sig = math.radians(FAN_SIGMA_DEG)
            psi = np.clip(rng.normal(0, sig, n), -2.3 * sig, 2.3 * sig)[:, None]
            dirs = (Rd * np.cos(psi) + et * np.sin(psi) + er * rng.normal(0, 0.06, (n, 1))
                    + pen[t] * (own[t] - Rd))
            dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
            P0.append(P + dirs * 0.003 * rng.random((n, 1)))
            V0.append(dirs * vs * (0.85 + 0.3 * rng.random((n, 1))))
            # Rosin-Rammler q=3 volume distribution; equal-mass parcels, so sample it directly
            d = 1.354 * smd * 1e-6 * np.cbrt(-np.log(1 - rng.random(n) * 0.9999))
            D0.append(np.clip(d, 4e-6, 3.5 * smd * 1e-6))
            M0.append(np.full(n, m_liq / n))
            TY.append(np.full(n, t))
    P0, V0 = np.vstack(P0), np.vstack(V0)
    D0, M0, TY = np.concatenate(D0), np.concatenate(M0), np.concatenate(TY)
    mdot = mdot_ox + mdot_fu

    def pix(x, y, z):
        Rz = R(z)
        ix = np.clip(((x / Rz + 1) / 2 * NXY).astype(int), 0, NXY - 1)
        iy = np.clip(((y / Rz + 1) / 2 * NXY).astype(int), 0, NXY - 1)
        return np.clip((z / dz).astype(int), 0, NZ - 1), ix, iy

    # Iterate: gas velocity depends on how much has vaporized, which depends on the gas velocity
    frac_vap = np.clip(zc / 0.1, 0.05, 1)
    for _ in range(3):
        u_gas = np.where(np.arange(NZ) <= it,
                         np.maximum(5, u_isen * np.minimum(1.5, frac_vap / max(frac_vap[it], 0.05))), u_isen)
        src = np.zeros((2, NZ, NXY, NXY))
        for P, m in flash:
            iz, ix, iy = pix(P[0:1], P[1:2], P[2:3])
            src[0, iz, ix, iy] += m
        film_dep = np.zeros((2, NZ, NA))    # kg/s of liquid landing on each wall strip
        liq_out = np.zeros(2)               # kg/s of drops leaving the nozzle (or still flying)
        p, v, d = P0.copy(), V0.copy(), D0.copy()
        alive = np.ones(len(d), bool)
        snap = []
        dt = 5e-6
        rho_l = np.where(TY == 0, RHO_OX, RHO_FU)
        k0 = np.where(TY == 0, K_OX, K_FU)
        for step in range(4000):
            idx = np.nonzero(alive)[0]
            if len(idx) == 0:
                break
            x, y, z = p[idx, 0], p[idx, 1], p[idx, 2]
            r = np.hypot(x, y)
            Rz = R(z)
            dd, ty = d[idx], TY[idx]
            frac = (dd / D0[idx]) ** 3
            out = z >= z_exit
            wall = (r >= Rz * 0.995) & ~out
            if wall.any():   # liquid on the wall joins the film on that strip
                iz = np.clip((z[wall] / dz).astype(int), 0, NZ - 1)
                ia = (np.arctan2(y[wall], x[wall]) % (2 * math.pi) / (2 * math.pi) * NA).astype(int) % NA
                np.add.at(film_dep, (ty[wall], iz, ia), M0[idx[wall]] * frac[wall])
            np.add.at(liq_out, ty[out], M0[idx[out]] * frac[out])
            alive[idx[out | wall]] = False
            keep = ~(out | wall)
            idx, x, y, z, r, Rz, dd, ty = idx[keep], x[keep], y[keep], z[keep], r[keep], Rz[keep], dd[keep], ty[keep]
            if step % 20 == 0:
                snap.append(np.column_stack([x, y, z, dd / D0[idx], ty, dd]))
            ug = np.interp(z, zc, u_gas)
            rho = np.interp(z, zc, rho_g)
            ur = ug * (r / Rz) * dRdz(z)
            rs = np.maximum(r, 1e-9)
            gv = np.column_stack([ur * x / rs, ur * y / rs, ug])
            vv = v[idx]
            inside = z < zt
            q = 0.1 * np.maximum(ug, 20) * math.sqrt(2 * dt / 2e-4) * inside
            vv += rng.normal(size=vv.shape) * q[:, None] * np.array([1, 1, 0.5])
            rel = np.linalg.norm(gv - vv, axis=1)
            re = rho * rel * dd / MU_GAS
            tau = rho_l[idx] * dd ** 2 / (18 * MU_GAS)
            e = np.exp(-dt * (1 + 0.15 * re ** 0.687) / tau)[:, None]
            vv = gv + (vv - gv) * e
            v[idx] = vv
            p[idx] += vv * dt
            d2 = dd ** 2 - k0[idx] * (1 + 0.27 * np.sqrt(re)) * dt
            nd = np.sqrt(np.maximum(d2, 0))
            gone = d2 <= 4e-12
            nd[gone] = 0
            dm = M0[idx] * (dd ** 3 - nd ** 3) / D0[idx] ** 3
            iz, ix, iy = pix(x, y, z)
            np.add.at(src, (ty, iz, ix, iy), dm)
            d[idx] = nd
            alive[idx[gone]] = False
        np.add.at(liq_out, TY[alive], M0[alive] * (d[alive] / D0[alive]) ** 3)

        # Wall film: march each strip downstream, evaporating into the wall cells at
        # m'' = ST_FILM G ln(1+B), with G the gas mass flux from the last iteration
        G = np.maximum(frac_vap, 0.02) * mdot / (math.pi * Rg ** 2)
        ds = dz * np.sqrt(1 + dRdz(zc) ** 2)
        tha = (np.arange(NA) + 0.5) / NA * 2 * math.pi
        film = np.zeros((2, NA))
        film_wall = film_dep.sum(axis=(1, 2))
        film_z = np.zeros((2, NZ))          # kg/s of liquid film flowing past each station
        for iz in range(NZ):
            film += film_dep[:, iz]
            ev = np.minimum(film, ST_FILM * G[iz] * np.log1p(B)[:, None] * 2 * math.pi * Rg[iz] / NA * ds[iz])
            film -= ev
            film_z[:, iz] = film.sum(axis=1)
            _, ix, iy = pix(0.999 * Rg[iz] * np.cos(tha), 0.999 * Rg[iz] * np.sin(tha), np.full(NA, zc[iz]))
            for t in (0, 1):
                np.add.at(src[t, iz], (ix, iy), ev[t])
            if iz == it:
                film_throat = film.sum(axis=1)
        cum = np.cumsum(src.sum(axis=(2, 3)), axis=1)    # vaporized kg/s up to each z, per propellant
        frac_vap = cum.sum(axis=0) / mdot

    # March each propellant's vapor downstream with sideways turbulent diffusion
    xs = (np.arange(NXY) + 0.5) / NXY * 2 - 1
    XI, YI = np.meshgrid(xs, xs, indexing="ij")
    mask = (XI ** 2 + YI ** 2 <= 1).astype(float)

    def blur(a, sig_px):
        if sig_px < 0.05:
            return a
        k = np.arange(-int(3 * sig_px + 1), int(3 * sig_px + 1) + 1)
        w = np.exp(-k ** 2 / (2 * sig_px ** 2))
        w /= w.sum()
        for ax in (0, 1):
            a = sum(wi * np.roll(a, ki, axis=ax) for wi, ki in zip(w, k))
        return a

    flux = np.zeros((2, NZ, NXY, NXY))
    cur = np.zeros((2, NXY, NXY))
    for iz in range(NZ):
        cur += src[:, iz]
        if zc[iz] < zt:
            sig = math.sqrt(2 * D_TURB * dz / max(u_gas[iz], 10)) / Rg[iz] * NXY / 2
            norm = np.maximum(blur(mask, sig), 1e-9)
            cur = np.array([blur(c * mask, sig) / norm * mask for c in cur])
            cur *= (src[:, :iz + 1].sum(axis=(1, 2, 3)) / np.maximum(cur.sum(axis=(1, 2)), 1e-12))[:, None, None]
        flux[:, iz] = cur

    snap = np.vstack(snap)
    return dict(pts=pts, zt=zt_mm, z_exit=z_exit * 1000, zc=zc * 1000, R=Rg * 1000, xi=xs,
                flux_ox=flux[0], flux_fu=flux[1], mask=mask, mdot=mdot,
                vap_ox=cum[0] / mdot_ox, vap_fu=cum[1] / mdot_fu, chi=chi, rupe=rupe,
                snap=snap, u_gas=u_gas, smd=(smd_ox, smd_fu),
                film_wall=film_wall / (mdot_ox, mdot_fu), film_throat=film_throat / (mdot_ox, mdot_fu),
                liq_out=liq_out / (mdot_ox, mdot_fu), film_z=film_z / np.array([[mdot_ox], [mdot_fu]]),
                film_dep=film_dep.sum(axis=2) / np.array([[mdot_ox], [mdot_fu]]),
                mdot_ox=mdot_ox, mdot_fu=mdot_fu, vo=vo, vf=vf, it=it, d_ox=g.D_OX, d_fu=g.D_FU)


def mixing_efficiency(sim, z_mm):
    """Rupe's Em at depth z_mm, from the local vapor mixture fractions (1 = perfectly mixed)."""
    k = int(np.argmin(abs(sim["zc"] - z_mm)))
    fo, fu = sim["flux_ox"][k].ravel(), sim["flux_fu"][k].ravel()
    m = fo + fu
    w = m / m.sum()
    r = fo / np.maximum(m, 1e-15)          # local oxidizer mass fraction
    R_ = fo.sum() / m.sum()
    lean, rich = r < R_, r >= R_
    return 1 - (w[lean] * (R_ - r[lean])).sum() / R_ - (w[rich] * (r[rich] - R_)).sum() / (1 - R_)


def cstar_efficiency(sim, z_mm):
    """c* efficiency if the gas at depth z_mm went through the throat as is: each cell burns at its
    own mixture fraction (F_TABLE/CSTAR_TABLE), unvaporized propellant gives nothing. Covers
    vaporization and mixing losses only (no heat loss or kinetics)."""
    k = int(np.argmin(abs(sim["zc"] - z_mm)))
    fo, fu = sim["flux_ox"][k].ravel(), sim["flux_fu"][k].ravel()
    m = fo + fu
    f = fu / np.maximum(m, 1e-15)
    f_design = 1 / (1 + OF_DESIGN)
    return (m * np.interp(f, F_TABLE, CSTAR_TABLE)).sum() / (sim["mdot"] * np.interp(f_design, F_TABLE, CSTAR_TABLE))


if __name__ == "__main__":
    s = simulate()
    izc = int(CHAMBER_LEN / 1000 / (s["z_exit"] / 1000 / NZ))
    it = int(np.argmin(s["R"]))
    print(f"N2O flashed at face {s['chi']*100:.0f}%   Rupe {s['rupe']:.2f}")
    print(f"SMD N2O {s['smd'][0]:.0f} um  IPA {s['smd'][1]:.0f} um")
    print(f"liquid hitting wall: N2O {s['film_wall'][0]*100:.1f}%  IPA {s['film_wall'][1]*100:.1f}%;"
          f"  film still liquid at throat: N2O {s['film_throat'][0]*100:.1f}%  IPA {s['film_throat'][1]*100:.1f}%")
    print(f"drops leaving nozzle: N2O {s['liq_out'][0]*100:.1f}%  IPA {s['liq_out'][1]*100:.1f}%")
    print(f"vaporized at end of cylinder: N2O {s['vap_ox'][izc]*100:.1f}%  IPA {s['vap_fu'][izc]*100:.1f}%")
    print(f"vaporized at throat: N2O {s['vap_ox'][it]*100:.1f}%  IPA {s['vap_fu'][it]*100:.1f}%")
    print(f"Rupe mixing efficiency at end of cylinder {mixing_efficiency(s, CHAMBER_LEN)*100:.0f}%")
    print(f"gas velocity end of cylinder {s['u_gas'][izc]:.0f} m/s, snapshot points {len(s['snap'])}")
