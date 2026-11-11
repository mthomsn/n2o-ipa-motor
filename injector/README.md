# N2O / IPA unlike-doublet injector

Sizing, layout, renders and a reduced-order spray/mixing model for an unlike-doublet
injector on an N2O / IPA motor. The design point comes from the HalfCatSim
'Motor Simulation' spreadsheet (cell references are noted next to each input).

The injector is one ring of radially oriented doublets: the fuel hole sits outboard
and aims inward, the ox hole sits inboard and aims outward, and the jets meet on an
impingement circle. Each element's spray fans out tangentially, so the ring of fans
forms an annulus that keeps the spray off the ablative. Elements are over-drilled
(for example 12 drilled, 8 open at design) and plugged down in balanced patterns.

## Scripts

| Script | What it does | Output |
| --- | --- | --- |
| `doublet_sizing.py` | Splits the HalfCatSim total CdA into N doublet pairs: hole sizes, jet velocities, momentum ratio, and the ox angle that sends the combined spray straight down the axis. | Prints to the terminal |
| `injector_layout.py` | Face layout for each entry in `LAYOUTS`: hole positions, drills, tilts, plug patterns and a spacing check. | `output/<layout>/<layout>.svg`, `.csv` (and `.png` with `--png`) |
| `chamber_mixing.py` | Chamber contour plus a steady-state spray, vaporization and mixing model. Run on its own for a summary of the design case. | Prints to the terminal |
| `injector_images.py` | Renders of the face, plug patterns, element section, 3D jets, 3D chamber and local O/F through the chamber. | `output/renders/*.png` |
| `chamber_plots.py` | Graphs from the mixing model: atomization, drop trajectories, vaporization, mixing, sensitivity, layout comparison. | `output/plots/*.png` |

## Running

Needs Python 3 with `numpy` and `matplotlib`.

```bash
pip install numpy matplotlib
```

Run from this folder, since output paths are relative:

```bash
python doublet_sizing.py
```

```bash
python injector_layout.py --png
```

```bash
python chamber_mixing.py
```

```bash
python injector_images.py
```

```bash
python chamber_plots.py
```

`--png` on `injector_layout.py` converts the SVGs with Quick Look, so it only works on macOS.

The bell nozzle contour is read from `../nozzle_xyz.txt` at the repo root (inches, measured from the
throat). If that file is missing, the model falls back to a 15° cone.

## Layouts

Defined in `LAYOUTS` in `injector_layout.py`:

| Name | Drilled | Open at design | Ox drill | Fuel drill | Plug options (open count) |
| --- | --- | --- | --- | --- | --- |
| `layout_12_elements` | 12 | 8 | #29 (3.454 mm) | 1.35 mm | 12, 10, 9, 8, 6 |
| `layout_15_elements` | 15 | 12 | #34 (2.819 mm) | 1.10 mm | 15, 12, 10 |
| `*_OX_OUTSIDE` variants | same | same | same | same | same, with the ox and fuel rings swapped |

Fuel jets are at 30° off axis and ox at 17.2°, which puts the combined spray on the
chamber axis. Impingement is 14 mm below the face on a 30 mm radius circle.

## Inputs to update

Each script has an `INPUTS` block at the top. Keep these in step with the spreadsheet:

- **CdA, mass flows, densities, chamber pressure**: copy from 'Motor Simulation'.
  `doublet_sizing.py`, `injector_layout.py` and `chamber_mixing.py` each hold their own copy.
- **Discharge coefficients** (`CD_OX = 0.40`, `CD_FU = 0.60`): the biggest unknowns for
  hole sizing. Measure them with a cold flow if you can.
- **Plug thread sizes** (`THREAD_OX`, `THREAD_FU`): set to the plugs you actually use,
  so the spacing check is right.
- **Combustion gas table** in `chamber_mixing.py`: rough N2O/IPA values. Replace with
  CEA output.
- **Atomization and mixing constants** (`SMD_SCALE`, `N2O_FLASH_SMD`, `K_FU`, `K_OX`,
  `D_TURB`, `RUPE_PENALTY`, `ST_FILM`): guesses. Calibrate against patternator or
  hot-fire data. `output/plots/sensitivity.png` shows how much the answers move with
  drop size and turbulent diffusivity.

## About the mixing model

This is a reduced-order model, not CFD. In outline:

- Each open doublet sprays droplet parcels from its impingement point. Direction is
  momentum-weighted, with a tangential fan and some memory of each jet's own heading
  (more when the Rupe ratio is far from 1).
- Drop sizes come from an unlike-doublet correlation with a Rosin-Rammler spread.
  Drag is Schiller-Naumann, and vaporization follows the d²-law with a convective
  correction. There is no secondary breakup yet.
- N2O flashes on injection: the equilibrium vapor fraction at Pc goes straight to gas.
- Liquid that reaches the wall becomes a film that runs downstream and evaporates.
  Film still liquid at the throat is counted as unburnt.
- Gas follows the wall contour with no recirculation. It is isentropic in the nozzle,
  and each propellant's vapor diffuses with a turbulent diffusivity.

Use it to compare layouts and spot trends, not for absolute c* efficiency numbers.
