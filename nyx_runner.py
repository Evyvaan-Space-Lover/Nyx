from pathlib import Path
import csv
import json
import math
import tkinter as tk
from tkinter import filedialog, messagebox

import traceback

import numpy as np
import SimMain as physics
from materials import Graphene, Aluminium, SiliconCarbide

MATERIALS = {
    "Graphene": Graphene,
    "Aluminium": Aluminium,
    "Silicon Carbide": SiliconCarbide,
}
GEOMETRIES = {
    "Rectangle": physics.RectangleLightSail,
    "Hemisphere": physics.SphereLightSail,
    "Paraboloid": physics.ParaboloidLightSail,
}


def _number(settings, key, cast=float):
    try:
        value = cast(settings[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid value for {key}: {settings.get(key)!r}") from exc
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{key} must be finite.")
    return value


def _step_count(settings, dt):
    if "steps" in settings:
        steps = _number(settings, "steps", int)
        if steps <= 0:
            raise ValueError("Number of steps must be a positive integer.")
    else:
        duration = _number(settings, "duration")
        if duration <= 0:
            raise ValueError("Duration must be positive.")
        steps = math.ceil(duration / dt)
    if steps > 200000:
        raise ValueError("This run exceeds 200,000 steps. Reduce the step count.")
    return steps


def _make_sail(geometry, material_name, s):
    material = MATERIALS[material_name]
    tx, ty = math.radians(_number(s, "theta_x")), math.radians(_number(s, "theta_y"))
    thickness = _number(s, "thickness")
    resolution = _number(s, "resolution", int)
    if thickness <= 0 or resolution < 2:
        raise ValueError("Thickness must be positive and resolution must be at least 2.")

    if geometry == "Rectangle":
        sail = physics.RectangleLightSail(
            _number(s, "width"), _number(s, "height"), resolution,
            0, 0, tx, ty, thickness=thickness, material=material)
    elif geometry == "Hemisphere":
        sail = physics.SphereLightSail(
            _number(s, "radius"), resolution, tx, ty,
            thickness=thickness, material=material)
    elif geometry == "Paraboloid":
        sail = physics.ParaboloidLightSail(
            _number(s, "radius"), _number(s, "a"), resolution, tx, ty,
            material=material, thickness=thickness)
    else:
        raise ValueError(f"Unsupported geometry: {geometry}")

    sail.positionOffset = np.array([_number(s, "pos_x"), _number(s, "pos_y"), _number(s, "pos_z")], dtype=float)
    sail.velocity = np.array([_number(s, "vel_x"), _number(s, "vel_y"), _number(s, "vel_z")], dtype=float)
    sail.angularVelocity = np.array([_number(s, "omega_x"), _number(s, "omega_y"), _number(s, "omega_z")], dtype=float)
    return sail


def _apply_beam(s):
    physics.I0 = _number(s, "I0")
    physics.w0 = _number(s, "w0")
    physics.lambdaLaser = _number(s, "lambdaLaser")
    if physics.I0 < 0 or physics.w0 <= 0 or physics.lambdaLaser <= 0:
        raise ValueError("I₀ must be non-negative; w₀ and wavelength must be positive.")
    physics.zR = math.pi * physics.w0**2 / physics.lambdaLaser
    physics.beamX = _number(s, "beam_x")
    physics.beamY = _number(s, "beam_y")


def _run_one(geometry, material_name, s, app, run_index, run_total):
    sail = _make_sail(geometry, material_name, s)
    dt = float(_number(s, "dt"))
    if dt <= 0:
        raise ValueError("Timestep must be positive.")
    steps = _step_count(s, dt)

    time = 0.0
    times = []
    for step in range(steps):
        sail.compute()
        if hasattr(app, "log_simulation_step"):
            app.log_simulation_step(
                step + 1,
                steps,
                time,
                sail.positionOffset,
                sail.angularVelocity,
                sail.Torque,
                sail.Force,
            )
        sail.update(dt, time)
        times.append(time)
        time += dt
        progress_interval = max(1, steps // 100)
        if step % progress_interval == 0 or step == steps - 1:
            completed_steps = (run_index - 1) * steps + step + 1
            total_steps = run_total * steps
            app.set_progress(
                completed_steps,
                total_steps,
                f"Run {run_index}/{run_total} — {geometry} / {material_name}; "
                f"iteration {step + 1}/{steps}",
            )
            root = getattr(app, "root", None)
            if root is not None:
                root.update_idletasks()

    # Recompute after the last state update so the displayed force/torque match final state.
    sail.compute()
    h = sail.history
    peak = max(h["maxTemperature"], default=0.0)
    margin = MATERIALS[material_name].maxTemperature - peak
    app.add_result_row(
        geometry,
        material_name,
        f"{sail.mass:.5g}",
        f"{np.linalg.norm(sail.Force):.5g}",
        f"{np.linalg.norm(sail.Torque):.5g}",
        f"{peak:.5g}",
        f"{margin:.5g}",
    )

    return {
        "geometry": geometry, "material": material_name, "mass": float(sail.mass),
        "final_force": np.asarray(sail.Force).tolist(), "final_torque": np.asarray(sail.Torque).tolist(),
        "history": {k: [np.asarray(v).tolist() if isinstance(v, np.ndarray) else float(v) if isinstance(v, (np.floating,)) else v for v in values]
                    for k, values in h.items()},
        "thermal_limit": float(MATERIALS[material_name].maxTemperature),
    }


def run_simulation(app, settings):
    try:
        if hasattr(app, "clear_step_log"):
            app.clear_step_log()
        if hasattr(app, "tabs") and hasattr(app, "steps_tab"):
            app.tabs.select(app.steps_tab)
        _apply_beam(settings)
        pairs = [(g, m) for g in settings["geometries"] for m in settings["materials"]]
        max_runs = _number(settings, "max_runs", int)
        if len(pairs) > max_runs:
            raise ValueError(f"Selected {len(pairs)} configurations, but Maximum run count is {max_runs}.")
        dt = _number(settings, "dt")
        if dt <= 0:
            raise ValueError("Timestep must be positive.")
        steps_per_run = _step_count(settings, dt)
        total_steps = len(pairs) * steps_per_run
        app.set_progress(0, total_steps, f"Starting Nyx — 0/{total_steps} iterations.")
        results = []
        for i, (geometry, material) in enumerate(pairs, start=1):
            results.append(_run_one(geometry, material, settings, app, i, len(pairs)))
        app.last_nyx_results = results
        app.refresh_plot()
        app.set_progress(
            total_steps,
            total_steps,
            f"Completed {len(pairs)} configuration(s), {total_steps} integration iterations.",
        )
        messagebox.showinfo(
            "Nyx complete",
            f"Finished {len(pairs)} simulation(s). Step-by-step output is on the "
            "Simulation steps tab; plots and exports are on Results & visualisation.\n\n"
            "To save results, choose an export button in Results & visualisation. "
            "CSV and JSON files are not saved automatically; you choose their location "
            "in the Save dialog.",
            parent=app.root,
        )
    except Exception as exc:
        traceback.print_exc()
        app.set_status("Run failed.")
        messagebox.showerror("Nyx run failed", f"{exc}\n\nCheck that SimMain.py and its dependencies are in the project folder.", parent=app.root)


def export_results(kind, settings, app):
    results = getattr(app, "last_nyx_results", [])
    if kind == "sweep_definition":
        path = filedialog.asksaveasfilename(parent=app.root, defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            Path(path).write_text(settings.get("sweep_text", ""), encoding="utf-8")
            messagebox.showinfo("Nyx export", f"Sweep definition saved to:\n{path}", parent=app.root)
        return
    if not results:
        messagebox.showwarning("Nyx export", "Run at least one simulation before exporting results.", parent=app.root)
        return
    if kind == "histories_json":
        path = filedialog.asksaveasfilename(parent=app.root, defaultextension=".json", filetypes=[("JSON", "*.json")])
        if path:
            Path(path).write_text(json.dumps(results, indent=2), encoding="utf-8")
            messagebox.showinfo("Nyx export", f"Full histories saved to:\n{path}", parent=app.root)
    elif kind == "summary_csv":
        path = filedialog.asksaveasfilename(parent=app.root, defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["geometry", "material", "mass", "final_force_N", "final_torque_Nm", "peak_temperature_K", "thermal_margin_K"])
                writer.writeheader()
                for r in results:
                    peak = max(r["history"]["maxTemperature"], default=0.0)
                    writer.writerow({"geometry":r["geometry"], "material":r["material"], "mass":r["mass"],
                        "final_force_N":np.linalg.norm(r["final_force"]), "final_torque_Nm":np.linalg.norm(r["final_torque"]),
                        "peak_temperature_K":peak, "thermal_margin_K":r["thermal_limit"]-peak})
            messagebox.showinfo("Nyx export", f"Summary CSV saved to:\n{path}", parent=app.root)
    elif kind == "selected_csv":
        path = filedialog.asksaveasfilename(parent=app.root, defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            r = results[-1]
            keys = ["time", "position", "velocity", "force", "torque", "orientation", "angularVelocity", "maxTemperature", "averageTemperature"]
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f); writer.writerow(keys)
                hist = r["history"]
                for i in range(len(hist["time"])):
                    writer.writerow([hist[k][i] for k in keys])
            messagebox.showinfo("Nyx export", f"Last run history CSV saved to:\n{path}", parent=app.root)
