# Nozzle math

`nozzle.py` sizes a converging (C) and a converging-diverging (CD) test nozzle by
scaling the HCR motor nozzle down to a ½" pipe. It uses only the Python standard
library, only reads input files, and prints its results.

```
python3 nozzle.py                              # defaults: 3.5/1.5/3.112 in motor, 1/2" Sch40, air, 100 psig
python3 nozzle.py --contour ../nozzle_xyz.txt  # take throat/exit from the contour file
python3 nozzle.py --inlet 0.546 --supply 90    # Sch80 pipe, 90 psig air
python3 nozzle.py -h                           # all options
```

## Inputs and where they come from

| Value | Default | Source |
|---|---|---|
| Chamber dia | 3.5 in | HalfCatSim motor spreadsheet (not included), Motor Simulation tab |
| Throat dia | 1.5 in | same sheet; matches `nozzle_xyz.txt` (r = 0.75 at x = 0) |
| Exit dia | 3.112 in | last point of `nozzle_xyz.txt` (r = 1.5562); the sheet itself says 3.0 |
| Test inlet | 0.622 in | ½" NPS Sch 40 pipe ID (Sch 80 = 0.546) |
| Ambient | 13.27 psia | sheet, at 2800 ft |
| γ | 1.4 | air |

Units are inches; `nozzle_xyz.txt` is assumed to be in inches because its throat
matches the sheet's 1.5 in throat.

## The math

All of this is standard 1-D isentropic ideal-gas nozzle theory. It's in any
compressible-flow or rocket-propulsion textbook; see Sources below.

**1. Area ratios** (area goes as diameter squared)

    contraction  Ac/At = (Dc/Dt)²       3.5/1.5   → 5.444
    expansion    Ae/At = (De/Dt)²       3.112/1.5 → 4.305

**2. Scaling to the pipe.** Keep both ratios and set the chamber diameter equal to
the pipe ID:

    s  = D_inlet / D_chamber            0.622/3.5 = 0.1777
    Dt = 1.5 × s = 0.267 in             De = 3.112 × s = 0.553 in

The C nozzle stops at the throat. The CD nozzle continues to De.

**3. Matching the expansion ratio to the gas.** A CD nozzle is "ideally
expanded" when its exit pressure equals ambient. From supply pressure P0 (absolute):

    exit Mach    M = sqrt( 2/(γ-1) · ((P0/Pa)^((γ-1)/γ) − 1) )
    area ratio   A/A* = (1/M) · [ (2/(γ+1)) · (1 + (γ-1)/2 · M²) ]^((γ+1)/(2(γ-1)))
    exit dia     De = Dt · sqrt(A/A*)

Going the other way (exit pressure for a given Ae/At), the script solves the
area-ratio equation for supersonic M by bisection, then

    Pe = P0 / (1 + (γ-1)/2 · M²)^(γ/(γ-1))

The motor's 4.305 ratio on 100 psig air gives Pe ≈ 3.0 psia, far below the
13.27 psia ambient. A rough rule (Summerfield criterion) is that flow separates
from the wall when Pe drops below about 0.35–0.4 × Pa (≈ 4.6–5.3 psia), so on
shop air that nozzle will separate. To use 4.305 properly on air you'd need
about 487 psig. For the C vs CD demo on shop air, size the exit to the ideal
ratio for your supply (≈ 1.77 → 0.355 in at 100 psig).

To get the γ ≈ 1.24 used earlier for the motor exhaust, I worked backwards from
the sheet's own results (Pc 316.9 psi, Pe 12.80 psi at Ae/At = 4) using the same
equations. That number is not used by the script.

**4. Converging section length.** You don't calculate a required length. You
pick a cone half-angle θ (30–45° is typical) and that sets the length. With a
rounded throat of upstream radius R = 1.5·Rt (Huzel & Huang):

    L = [ Rt(√(Ac/At) − 1) + R(sec θ − 1) ] / tan θ
      = [ (Rc − Rt) + R(sec θ − 1) ] / tan θ

For a sharp cone (R = 0) it's just (Rc − Rt)/tan θ.

Things to keep in mind:
- **Round the throat.** A sharp throat edge makes the flow pinch smaller than the
  hole (a vena contracta), which cuts the discharge coefficient from about 0.98
  to about 0.85. Typical radii are 1.5·Rt upstream and 0.382·Rt downstream (Rao
  bell convention).
- **Angle.** Steeper is shorter but loses a little. Shallower is smoother but
  longer. Past about 60° the nozzle behaves like an orifice.
- **Demo fairness.** Make the converging section identical on the C and CD
  nozzles (same angle, throat rounding and throat diameter) so the only
  difference is the diverging part. 45° is a standard countersink angle and easy
  to machine.
- **Real motor only.** The converging volume counts toward L*, and more wall
  means more area to heat.

**5. Diverging cone length and loss**

    L = (Re − Rt) / tan α
    divergence efficiency λ = (1 + cos α)/2     15° → 0.983

## Sources

- Anderson, *Modern Compressible Flow*, ch. 3 & 5 (isentropic relations, area–Mach relation).
- Sutton & Biblarz, *Rocket Propulsion Elements*, ch. 3 (nozzle theory, separation, divergence factor λ).
- Huzel & Huang, *Modern Engineering for Design of Liquid-Propellant Rocket Engines*,
  NASA SP-125 / AIAA, ch. 4 (converging length formula, throat radii, contraction ratios).
- Rao, "Exhaust Nozzle Contour for Optimum Thrust," *Jet Propulsion* 28 (1958) (bell nozzle, 0.382·Rt).
- NASA Glenn "Beginner's Guide to Rockets" pages on nozzles and isentropic flow (free online summary).
- HalfCatSim v1.3.11 (Half Cat Rocketry), the source of the motor inputs above.
