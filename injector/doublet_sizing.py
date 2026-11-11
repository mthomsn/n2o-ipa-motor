"""
Unlike-doublet injector sizing for N2O / IPA, built on HalfCatSim results.

HalfCatSim models both propellants with the single-phase orifice equation
    m_dot = CdA * sqrt(2 * rho * dP)
and uses a low Cd (~0.4) on the N2O side to stand in for two-phase flashing.
We keep the TOTAL CdA your sim says works and split it into N doublet pairs,
so the hole pattern flows the same as the motor you simulated.

Edit the INPUTS block to match the current version of your spreadsheet.
"""
import math

IN2_TO_M2 = 0.00064516
LBM_TO_KG = 0.45359237

# ---------------- INPUTS (from 'Motor Simulation' tab) ----------------
CDA_OX = 3.0198955e-05      # m^2   G28  Oxidizer CdA
CDA_FU = 7.0045861e-06      # m^2   K28  Fuel CdA
MDOT_OX = 3.2396 * LBM_TO_KG  # kg/s  G29  (avg, first 10 steps)
MDOT_FU = 1.0181 * LBM_TO_KG  # kg/s  K29
RHO_OX = 726.59             # kg/m^3 G10  liquid N2O @ 80 F
RHO_FU = 786.0              # kg/m^3 K11  IPA
CHAMBER_DIA_MM = 4.75 * 25.4  # C27

# Discharge coefficients you expect for the NEW small drilled holes.
# These are the biggest unknowns -- measure them with a cold flow if at all possible.
CD_OX = 0.40   # matches HalfCatSim's N2O assumption
CD_FU = 0.60   # sharp-edged drilled hole, liquid, L/D ~ 3-5

FUEL_ANGLE_DEG = 30.0  # fuel jet angle off the injector axis
MIN_WALL_MM = 2.0      # minimum material between hole edges
# ----------------------------------------------------------------------


def dia_mm(area_m2):
    return math.sqrt(4 * area_m2 / math.pi) * 1000


A_OX_TOTAL = CDA_OX / CD_OX   # geometric area needed
A_FU_TOTAL = CDA_FU / CD_FU

# Momentum per jet (liquid density = lowest velocity, conservative for N2O)
v_ox = MDOT_OX / (RHO_OX * A_OX_TOTAL)
v_fu = MDOT_FU / (RHO_FU * A_FU_TOTAL)
mom_ratio = (MDOT_OX * v_ox) / (MDOT_FU * v_fu)

# Ox angle that puts the combined spray straight down the chamber axis
ox_angle = math.degrees(math.asin(min(1.0, math.sin(math.radians(FUEL_ANGLE_DEG)) / mom_ratio)))

# Resultant spray angle if both jets are at the same angle
th = math.radians(FUEL_ANGLE_DEG)
beta_equal = math.degrees(math.atan((mom_ratio - 1) * math.sin(th) / ((mom_ratio + 1) * math.cos(th))))

print(f"Total geometric area   ox {A_OX_TOTAL*1e6:7.2f} mm^2   fuel {A_FU_TOTAL*1e6:6.2f} mm^2")
print(f"Jet velocity (liquid)  ox {v_ox:5.1f} m/s       fuel {v_fu:5.1f} m/s")
print(f"Momentum ratio ox/fuel {mom_ratio:.2f}")
print(f"  equal {FUEL_ANGLE_DEG:.0f}/{FUEL_ANGLE_DEG:.0f} deg angles -> spray tilts {beta_equal:.1f} deg toward fuel side")
print(f"  fuel at {FUEL_ANGLE_DEG:.0f} deg, ox at {ox_angle:.1f} deg -> spray goes straight down the axis")
print()
print(" pairs | ox hole        | fuel hole      | ox:fu dia | impinge dist (5-7 d_ox)")
print(" ------+----------------+----------------+-----------+------------------------")
for n in range(2, 21):
    d_o = dia_mm(A_OX_TOTAL / n)
    d_f = dia_mm(A_FU_TOTAL / n)
    if d_f < 0.5:
        break
    print(f"  {n:3d}  | {d_o:5.2f} mm {d_o/25.4:.4f}\" | {d_f:5.2f} mm {d_f/25.4:.4f}\" |"
          f"   {d_o/d_f:4.2f}    | {5*d_o:4.1f} - {7*d_o:4.1f} mm")
