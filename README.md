# Nyx — Light-Sail Dynamics Workbench

**Nyx** is a Python research workbench for exploring how idealised light-sail geometries respond to radiation pressure from a Gaussian laser beam. It combines a numerical physics model with a desktop interface for configuring runs, comparing configurations, inspecting time histories, and exporting results.

> **Research / educational software.** Nyx is an exploratory simulation, not a flight-readiness tool or a validated engineering package. Results depend on the model assumptions and material presets described below.

<p align="center">
  <em>Geometry · Radiation pressure · Rigid-body dynamics · Thermal estimates</em>
</p>

## At a glance

- **Three idealised geometries:** rectangular sheet, open hemispherical shell, and paraboloid.
- **Material presets:** Graphene, Aluminium, and Silicon Carbide.
- **Gaussian-beam configuration:** peak intensity, beam waist, wavelength, and lateral beam position.
- **Rigid-body response:** translation, orientation, force, and torque histories.
- **Thermal estimates:** local surface-temperature estimates and comparison with the configured material temperature limit.
- **Desktop workbench:** configure a run, inspect results and plots, and export data for further analysis.

## Screenshots

_Add a screenshot of the Nyx interface here when you're ready._

<!-- Example:
![Nyx interface](docs/images/nyx-interface.png)
-->

## Geometries and presets

| Geometry | Model | Main parameters |
|---|---|---|
| Rectangle | Flat rectangular sail | Width, height, thickness, resolution |
| Hemisphere | Open hemispherical shell — **not** a closed sphere | Radius, thickness, resolution |
| Paraboloid | Parabolic sail surface | Radius, shape coefficient `a`, thickness, resolution |

The supplied material presets are convenient model inputs, not certified material datasets. Check the values in `materials/` before using them in a study, and record any changes when reporting results.

## Getting started

### Requirements

- Python 3.10 or newer
- NumPy
- Matplotlib
- Tkinter (normally included with standard Python installations on Windows and many macOS/Linux installations)

From the project directory, install the Python dependencies:

```bash
python -m pip install numpy matplotlib
```

If your checkout includes a `requirements.txt`, you can instead use:

```bash
python -m pip install -r requirements.txt
```

### Launch the connected workbench

Run the launcher—not the interface module directly:

```bash
python NyxApp.py
```

`NyxApp.py` opens the interface in `Runner.py` and connects it to the simulation adapter in `nyx_runner.py`. Keep the physics module and its dependencies in the project structure expected by their imports. In the current setup, the physics module is named `SimMain.py`.

### First run

1. Select one geometry and one material.
2. Start with the default beam and initial-condition settings.
3. Use a short duration and a sensible timestep for the first test.
4. Click **Run selected configurations**.
5. Inspect the result row and plots, then export the data if needed.

Start with a low mesh resolution and short run to check a configuration before increasing computational cost. Smaller timesteps and finer surface discretisation can substantially increase runtime.

## Project layout

The exact folder layout can vary between checkouts, but the connected workbench uses these main components:

```text
Nyx/
├── NyxApp.py              # Application launcher; connects UI and runner
├── Runner.py              # Tkinter desktop interface
├── nyx_runner.py          # Adapts UI settings to the physics model
├── SimMain.py             # Sail geometry, force/torque and dynamics model
├── materials/             # Material definitions and presets
└── analysis/
    ├── experiment.py      # Experiment utilities
    └── plotting.py        # Analysis/plotting utilities
```

Use the actual package and import structure in your checkout if it differs from this illustrative layout. Avoid renaming package directories without updating imports.

## Model overview

Nyx models an idealised Gaussian laser beam with a fixed nominal propagation direction along **+Z**. The beam settings include peak intensity \(I_0\), waist \(w_0\), wavelength \(\lambda\), and lateral beam-centre offsets.

The sail model estimates radiation-pressure force and torque over the discretised surface, then advances translational and rotational state. The workbench records time histories so configurations can be inspected and compared. A local radiative-equilibrium calculation estimates surface temperature from absorbed illumination and the configured emissivity.

These are deliberately simplified models. Treat outputs as model estimates, not experimental measurements.

## Outputs

Depending on the selected run and export option, the workbench can report or save quantities such as:

- Sail geometry and material preset
- Mass
- Force and torque
- Position, velocity, orientation, and angular velocity histories
- Peak and average surface-temperature estimates
- Temperature margin relative to the preset's configured limit

The connected runner supports summary CSV, full JSON histories, and a CSV history for an individual result. Exported values should be checked against the selected configuration and model assumptions before being used in reports.

## Current scope and limitations

Nyx is intended for comparative, exploratory studies of idealised light-sail response. The current model does **not** provide:

- A full electromagnetic scattering or optical simulation
- Lateral heat conduction, transient heat capacity, or a multilayer thermal model
- Structural deformation, buckling, or finite-element stress analysis
- Gravity, orbital mechanics, relativistic dynamics, or a complete interstellar trajectory
- Attitude-control hardware or closed-loop guidance
- Engineering certification or proof of flight readiness

A temperature estimate below a configured limit does not, by itself, establish that a design is safe or flight-ready. Geometry, optical properties, emissivity, discretisation, and numerical integration assumptions all affect the result.

The interface includes parameter-sweep controls, but the connected runner may not yet execute every sweep workflow. Check the current implementation before relying on batch-sweep functionality.

## Reproducibility notes

For a meaningful comparison, record at least:

- Geometry and all shape parameters
- Material preset and any edited properties
- Beam intensity, waist, wavelength, and offset
- Initial position, velocity, orientation/tilt, and angular velocity
- Surface resolution, timestep, and total simulated duration
- Software version or commit hash

When comparing geometries, state the comparison basis explicitly—for example, equal projected aperture, equal mass, or equal surface area. These are different constraints and can lead to different conclusions.

## Troubleshooting

**The GUI opens, but Run reports that no callback is connected**

Launch with `python NyxApp.py`. Running `Runner.py` directly opens the interface without the connected simulation callback.

**Python cannot import `SimMain` or a material preset**

Check that the physics module and material package are present and that their names match the imports used in `nyx_runner.py`. Run the application from the project environment with its dependencies installed.

**A run fails with a traceback**

Read the complete traceback in the terminal. The last few lines usually identify the exception and the exact source line. Avoid changing unrelated physics code until the failing operation is identified.

**A run is slow**

Reduce surface resolution, shorten the duration, or increase the timestep cautiously. Confirm that the chosen timestep still resolves the dynamics you want to study.

## Contributing and feedback

Bug reports, reproducible test cases, and suggestions are welcome. When reporting a problem, include your Python version, operating system, the configuration used, and the complete traceback. For scientific issues, explain the expected physical behaviour and the assumptions behind that expectation.

## Citation

If you use Nyx in a report, presentation, or research project, cite the repository and the specific commit you used:

**Evyvaan Singh, _Nyx: Light-Sail Dynamics Workbench_, GitHub repository.**  
https://github.com/Evyvaan-Space-Lover/Nyx
### No AI was used for the code itself. AI was only used for the documentation. 
