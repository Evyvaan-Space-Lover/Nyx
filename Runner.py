import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

BG = "#0b1220"
PANEL = "#111c2c"
TEXT = "#e5edf7"
MUTED = "#9aaec5"
ACCENT = "#79b8ff"
GEOMETRIES = ("Rectangle", "Hemisphere", "Paraboloid")
MATERIALS = ("Graphene", "Aluminium", "Silicon Carbide")
PLOTS = (
    "Force magnitude", "Force components", "Torque magnitude", "Position",
    "Velocity", "Angular velocity", "Peak temperature", "Mean temperature",
    "Thermal margin",
)
MAX_STEP_LOG_LINES = 5000


class NyxInterface:
    def __init__(self, root, on_run=None, on_sweep=None, on_export=None):
        self.root = root
        self.on_run = on_run
        self.on_sweep = on_sweep
        self.on_export = on_export
        self.vars = {}
        self.geometry_vars = {}
        self.material_vars = {}
        self.root.title("NYX | Research Workbench")
        self.root.geometry("1420x900")
        self.root.minsize(1050, 680)
        self.root.configure(bg=BG)
        self._style()
        self._build()
        self.set_status("Ready — configure a run or parameter sweep.")

    def _style(self):
        s = ttk.Style()
        try:
            s.theme_use("clam")
        except tk.TclError:
            pass
        s.configure("TFrame", background=BG)
        s.configure("TLabel", background=BG, foreground=TEXT)
        s.configure("TCheckbutton", background=BG, foreground=TEXT)
        s.configure("TButton", padding=7)
        s.configure("TNotebook", background=BG, borderwidth=0)
        s.configure("TNotebook.Tab", padding=(13, 8))
        s.configure("TLabelframe", background=BG, foreground=ACCENT)
        s.configure("TLabelframe.Label", background=BG, foreground=ACCENT)
        s.configure("TEntry", fieldbackground=PANEL, foreground=TEXT)
        s.configure("TCombobox", fieldbackground=PANEL, foreground=TEXT)
        s.configure("Treeview", background=PANEL, fieldbackground=PANEL,
                    foreground=TEXT, rowheight=25)
        s.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))

    def _var(self, key, default):
        if key not in self.vars:
            self.vars[key] = tk.StringVar(value=str(default))
        return self.vars[key]

    def _entry(self, parent, label, key, default, row, col=0, width=14):
        ttk.Label(parent, text=label).grid(row=row, column=col * 2,
            sticky="w", padx=(7, 5), pady=4)
        e = ttk.Entry(parent, textvariable=self._var(key, default), width=width)
        e.grid(row=row, column=col * 2 + 1, sticky="ew", padx=(0, 9), pady=4)
        return e

    @staticmethod
    def _group(parent, title):
        return ttk.LabelFrame(parent, text=title, padding=9)

    def _build(self):
        header = ttk.Frame(self.root, padding=(16, 12))
        header.pack(fill="x")
        ttk.Label(header, text="NYX", font=("Segoe UI", 27, "bold"),
                  foreground=ACCENT).pack(side="left")
        ttk.Label(header, text="RESEARCH WORKBENCH  /  LIGHT-SAIL DYNAMICS",
                  foreground=MUTED).pack(side="left", padx=14, pady=(9, 0))
        self.status_var = tk.StringVar()
        ttk.Label(header, textvariable=self.status_var, foreground=MUTED).pack(
            side="right", pady=(9, 0))

        self.tabs = ttk.Notebook(self.root)
        self.tabs.pack(fill="both", expand=True, padx=12, pady=(0, 8))
        self.config_tab = ttk.Frame(self.tabs, padding=10)
        self.sweep_tab = ttk.Frame(self.tabs, padding=10)
        self.results_tab = ttk.Frame(self.tabs, padding=10)
        self.steps_tab = ttk.Frame(self.tabs, padding=10)
        self.tabs.add(self.config_tab, text="  Configuration  ")
        self.tabs.add(self.sweep_tab, text="  Parameter sweep  ")
        self.tabs.add(self.results_tab, text="  Results & visualisation  ")
        self.tabs.add(self.steps_tab, text="  Simulation steps  ")
        self._build_config()
        self._build_sweep()
        self._build_results()
        self._build_steps()

        footer = ttk.Frame(self.root, padding=(12, 4, 12, 12))
        footer.pack(fill="x")
        self.progress = ttk.Progressbar(footer, maximum=100, value=0)
        self.progress.pack(side="left", fill="x", expand=True, padx=(0, 12))
        self.progress_text = tk.StringVar(value="Idle")
        ttk.Label(footer, textvariable=self.progress_text, width=32).pack(side="right")

    def _build_config(self):
        left = ttk.Frame(self.config_tab)
        left.pack(side="left", fill="both", expand=True, padx=(0, 8))
        right = ttk.Frame(self.config_tab)
        right.pack(side="left", fill="both", expand=True, padx=(8, 0))

        box = self._group(left, "Geometries — select any combination")
        box.pack(fill="x", pady=(0, 9))
        for i, name in enumerate(GEOMETRIES):
            v = tk.BooleanVar(value=name == "Rectangle")
            self.geometry_vars[name] = v
            ttk.Checkbutton(box, text=name, variable=v).grid(
                row=0, column=i, sticky="w", padx=8, pady=3)
        ttk.Label(box, text="Hemisphere is treated as an open hemispherical shell.",
                  foreground=MUTED).grid(row=1, column=0, columnspan=3,
                                         sticky="w", padx=8, pady=4)

        box = self._group(left, "Materials — select any combination")
        box.pack(fill="x", pady=(0, 9))
        for i, name in enumerate(MATERIALS):
            v = tk.BooleanVar(value=name == "Graphene")
            self.material_vars[name] = v
            ttk.Checkbutton(box, text=name, variable=v).grid(
                row=0, column=i, sticky="w", padx=8, pady=3)

        box = self._group(left, "Sail geometry and initial conditions")
        box.pack(fill="x", pady=(0, 9))
        fields = [
            ("Rectangle width (m)", "width", "10"),
            ("Rectangle height (m)", "height", "10"),
            ("Curved radius (m)", "radius", "1"),
            ("Paraboloid a (1/m)", "a", "0.5"),
            ("Thickness (m)", "thickness", "1e-6"),
            ("Resolution", "resolution", "40"),
            ("Tilt X (degrees)", "theta_x", "0"),
            ("Tilt Y (degrees)", "theta_y", "0"),
            ("Position X (m)", "pos_x", "0"),
            ("Position Y (m)", "pos_y", "0"),
            ("Position Z (m)", "pos_z", "0"),
            ("Velocity X (m/s)", "vel_x", "0"),
            ("Velocity Y (m/s)", "vel_y", "0"),
            ("Velocity Z (m/s)", "vel_z", "0"),
            ("Angular velocity X (rad/s)", "omega_x", "0"),
            ("Angular velocity Y (rad/s)", "omega_y", "0"),
            ("Angular velocity Z (rad/s)", "omega_z", "0"),
        ]
        for i, (label, key, default) in enumerate(fields):
            self._entry(box, label, key, default, i // 2, i % 2)

        box = self._group(right, "Gaussian beam parameters")
        box.pack(fill="x", pady=(0, 9))
        for i, item in enumerate([
            ("Peak intensity I₀ (W/m²)", "I0", "1e9"),
            ("Waist w₀ (m)", "w0", "1"),
            ("Wavelength (m)", "lambdaLaser", "1e-6"),
            ("Beam centre X (m)", "beam_x", "0"),
            ("Beam centre Y (m)", "beam_y", "0"),
        ]):
            self._entry(box, *item, i // 2, i % 2)

        box = self._group(right, "Integration / batch limits")
        box.pack(fill="x", pady=(0, 9))
        self._entry(box, "Number of steps", "steps", "10", 0, 0)
        self._entry(box, "Timestep (s per step)", "dt", "0.01", 1, 0)
        self._entry(box, "Maximum run count", "max_runs", "500", 2, 0)
        ttk.Label(box, text="Simulated duration = number of steps × timestep.",
                  foreground=MUTED).grid(row=3, column=0, columnspan=2,
                                         sticky="w", padx=7, pady=5)

        row = ttk.Frame(right)
        row.pack(fill="x", pady=5)
        ttk.Button(row, text="Run selected configurations",
                   command=lambda: self._dispatch(self.on_run, "run")).pack(side="left")
        ttk.Button(row, text="Reset defaults", command=self.reset_defaults).pack(
            side="left", padx=8)

        box = self._group(right, "Workflow")
        box.pack(fill="both", expand=True)
        ttk.Label(box, justify="left", wraplength=490, foreground=MUTED,
            text=("This file is the interface only. It does not contain or import "
                  "the physics engine. Connect your existing simulator through the "
                  "callbacks to run experiments, update the results table and export data.")
        ).pack(anchor="nw", padx=6, pady=5)

    def _build_sweep(self):
        box = self._group(self.sweep_tab, "Cartesian parameter sweep")
        box.pack(fill="x", pady=(0, 10))
        ttk.Label(box, justify="left", foreground=MUTED, text=(
            "Enter one parameter per line: name = comma-separated values. Angles are degrees; other values use SI units.\n"
            "Supported: width, height, radius, a, thickness, resolution, theta_x, theta_y, beam_x, beam_y, I0, w0, lambdaLaser\n"
            "Example: 3 tilt values × 3 beam offsets = 9 combinations per geometry–material pair."
        )).pack(anchor="w", pady=(0, 8))
        self.sweep_text = tk.Text(box, height=12, bg=PANEL, fg=TEXT,
            insertbackground="white", relief="flat", font=("Consolas", 10),
            padx=8, pady=8)
        self.sweep_text.pack(fill="x")
        self.sweep_text.insert("1.0", "# Lines beginning with # are ignored\ntheta_x = 0, 15, 30\nbeam_x = 0, 0.25, 0.5\n")
        row = ttk.Frame(box)
        row.pack(fill="x", pady=8)
        ttk.Button(row, text="Run parameter sweep",
                   command=lambda: self._dispatch(self.on_sweep, "sweep")).pack(side="left")
        ttk.Button(row, text="Clear sweep",
                   command=lambda: self.sweep_text.delete("1.0", "end")).pack(side="left", padx=8)
        self.sweep_count_var = tk.StringVar(value="Sweep definition ready.")
        ttk.Label(box, textvariable=self.sweep_count_var, foreground=ACCENT).pack(anchor="w")

    def _build_results(self):
        export_bar = ttk.LabelFrame(
            self.results_tab, text="Export results", padding=(8, 6)
        )
        export_bar.pack(fill="x", pady=(0, 8))
        for title, kind in (
            ("Summary CSV", "summary_csv"),
            ("Full histories JSON", "histories_json"),
            ("Last run history CSV", "selected_csv"),
            ("Save sweep definition", "sweep_definition"),
        ):
            ttk.Button(
                export_bar,
                text=title,
                command=lambda export_kind=kind: self._export(export_kind),
            ).pack(side="left", padx=(0, 7))

        split = ttk.Panedwindow(self.results_tab, orient="horizontal")
        split.pack(fill="both", expand=True)
        left, right = ttk.Frame(split), ttk.Frame(split)
        split.add(left, weight=1)
        split.add(right, weight=3)
        cols = ("geometry", "material", "mass", "force", "torque", "peak", "margin")
        self.result_table = ttk.Treeview(left, columns=cols, show="headings", height=20)
        for col, title, width in zip(cols,
            ("Geometry", "Material", "Mass kg", "|F| N", "|τ| N·m", "Peak K", "Margin K"),
            (100, 110, 85, 90, 90, 80, 85)):
            self.result_table.heading(col, text=title)
            self.result_table.column(col, width=width, anchor="center")
        self.result_table.pack(fill="both", expand=True)
        row = ttk.Frame(left)
        row.pack(fill="x", pady=6)
        ttk.Button(row, text="Remove selected", command=self.remove_selected).pack(side="left")
        ttk.Button(row, text="Clear results", command=self.clear_results).pack(side="left", padx=6)

        toolbar = ttk.Frame(right)
        toolbar.pack(fill="x", pady=(0, 7))
        ttk.Label(toolbar, text="Plot").pack(side="left")
        self.plot_var = tk.StringVar(value=PLOTS[0])
        plot_selector = ttk.Combobox(toolbar, textvariable=self.plot_var, values=PLOTS,
                                     state="readonly", width=29)
        plot_selector.pack(side="left", padx=7)
        plot_selector.bind("<<ComboboxSelected>>", lambda _event: self.refresh_plot())
        ttk.Button(toolbar, text="Refresh plot", command=self.refresh_plot).pack(side="left")
        self.figure = Figure(figsize=(8, 5.5), dpi=100, tight_layout=True,
                             facecolor=BG)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_facecolor(PANEL)
        self.axes.set_title("Simulation results will appear here")
        self.canvas = FigureCanvasTkAgg(self.figure, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.details = tk.Text(right, height=5, bg=PANEL, fg=TEXT,
            insertbackground="white", relief="flat", font=("Consolas", 9), padx=8, pady=6)
        self.details.pack(fill="x", pady=(7, 0))
        self.details.insert("1.0", "Select a run to show its parameters and summary.")
        self.details.configure(state="disabled")

    def _build_steps(self):
        toolbar = ttk.Frame(self.steps_tab)
        toolbar.pack(fill="x", pady=(0, 7))
        ttk.Label(
            toolbar,
            text="Live step-by-step state (position, angular velocity, torque, force, time)",
            foreground=MUTED,
        ).pack(side="left")
        ttk.Button(toolbar, text="Clear log", command=self.clear_step_log).pack(side="right")

        log_frame = ttk.Frame(self.steps_tab)
        log_frame.pack(fill="both", expand=True)
        self.step_log = tk.Text(
            log_frame,
            bg=PANEL,
            fg=TEXT,
            insertbackground=TEXT,
            relief="flat",
            font=("Consolas", 10),
            wrap="none",
            state="disabled",
        )
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.step_log.yview)
        self.step_log.configure(yscrollcommand=scrollbar.set)
        self.step_log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self._step_log_line_count = 0
        self._step_log_truncated = False
        self.log_step_message("Run a simulation to see per-step position, angular velocity, torque, and force.")

    def collect_settings(self):
        """Collect current UI inputs for the simulator callback."""
        result = {key: var.get().strip() for key, var in self.vars.items()}
        result["geometries"] = [n for n, v in self.geometry_vars.items() if v.get()]
        result["materials"] = [n for n, v in self.material_vars.items() if v.get()]
        result["sweep_text"] = self.sweep_text.get("1.0", "end").strip()
        result["plot"] = self.plot_var.get()
        return result

    def _dispatch(self, callback, action):
        settings = self.collect_settings()
        if not settings["geometries"]:
            messagebox.showwarning("Configuration", "Select at least one geometry.", parent=self.root)
            return
        if not settings["materials"]:
            messagebox.showwarning("Configuration", "Select at least one material.", parent=self.root)
            return
        if callback is None:
            messagebox.showinfo("Interface only",
                f"No {action} callback is connected. The simulator is intentionally separate.",
                parent=self.root)
            return
        callback(settings)

    def _export(self, kind):
        if self.on_export is None:
            messagebox.showinfo("Interface only", f"Export slot: {kind}. Connect an export callback.", parent=self.root)
        else:
            self.on_export(kind, self.collect_settings())

    def set_status(self, text):
        self.status_var.set(text)

    def set_progress(self, current, total, message=None):
        self.progress.configure(maximum=max(int(total), 1), value=int(current))
        self.progress_text.set(f"{current:,} / {total:,} iterations")
        if message:
            self.set_status(message)

    def clear_step_log(self):
        self.step_log.configure(state="normal")
        self.step_log.delete("1.0", "end")
        self.step_log.configure(state="disabled")
        self._step_log_line_count = 0
        self._step_log_truncated = False

    def log_step_message(self, message):
        self.step_log.configure(state="normal")
        self.step_log.insert("end", f"{message}\n")
        self._step_log_line_count += 1
        if self._step_log_line_count > MAX_STEP_LOG_LINES:
            self.step_log.delete("1.0", "1001.0")
            self._step_log_line_count -= 1000
            if not self._step_log_truncated:
                self.step_log.insert(
                    "1.0",
                    "[Earlier step output trimmed; showing the most recent 5,000 lines.]\n",
                )
                self._step_log_line_count += 1
                self._step_log_truncated = True
        self.step_log.see("end")
        self.step_log.configure(state="disabled")

    def log_simulation_step(self, step, total_steps, time, position, omega, torque, force):
        self.log_step_message(
            f"Step {step}/{total_steps}\n"
            f"Position: {np.asarray(position)}\n"
            f"Omega: {np.asarray(omega)}\n"
            f"Torque: {np.asarray(torque)}\n"
            f"Force: {np.asarray(force)}\n"
            f"(Time: {time})\n{'=' * 45}"
        )

    def add_result_row(self, geometry, material, mass="—", force="—", torque="—", peak="—", margin="—"):
        """Backend hook: append a completed simulation to the visible results table."""
        return self.result_table.insert("", "end", values=(geometry, material, mass, force, torque, peak, margin))

    def show_plot(self, x, y, label="Run", xlabel="Time (s)", ylabel="Value"):
        """Backend hook: display one series from your existing simulator."""
        self.axes.clear()
        self.axes.plot(x, y, label=label)
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_title(self.plot_var.get())
        self._style_plot_axes()
        self.axes.legend(facecolor=BG, edgecolor=MUTED, labelcolor=TEXT)
        self.canvas.draw_idle()

    def refresh_plot(self):
        self.axes.clear()
        results = getattr(self, "last_nyx_results", [])
        plot_name = self.plot_var.get()
        history_keys = {
            "Force magnitude": ("force", "Force (N)", "magnitude"),
            "Force components": ("force", "Force (N)", "components"),
            "Torque magnitude": ("torque", "Torque (N·m)", "magnitude"),
            "Position": ("position", "Position (m)", "components"),
            "Velocity": ("velocity", "Velocity (m/s)", "components"),
            "Angular velocity": ("angularVelocity", "Angular velocity (rad/s)", "components"),
            "Peak temperature": ("maxTemperature", "Temperature (K)", "scalar"),
            "Mean temperature": ("averageTemperature", "Temperature (K)", "scalar"),
            "Thermal margin": ("maxTemperature", "Margin (K)", "margin"),
        }
        selected = history_keys.get(plot_name)
        if not results or selected is None:
            self.axes.set_title(
                f"{plot_name} — run a simulation to display results"
                if not results else f"{plot_name} — no plot available"
            )
            self.axes.set_xlabel("Time (s)")
            self.axes.set_ylabel("Value")
            self._style_plot_axes()
            self.canvas.draw_idle()
            return

        key, ylabel, mode = selected
        for result in results:
            history = result["history"]
            times = np.asarray(history["time"], dtype=float)
            values = np.asarray(history[key], dtype=float)
            run_label = f"{result['geometry']} / {result['material']}"
            if mode == "magnitude":
                self.axes.plot(times, np.linalg.norm(values, axis=1), label=run_label)
            elif mode == "margin":
                self.axes.plot(
                    times,
                    result["thermal_limit"] - values,
                    label=run_label,
                )
            elif mode == "components":
                for axis_index, axis_name in enumerate(("X", "Y", "Z")):
                    self.axes.plot(
                        times,
                        values[:, axis_index],
                        label=f"{run_label} {axis_name}",
                    )
            else:
                self.axes.plot(times, values, label=run_label)

        self.axes.set_title(plot_name)
        self.axes.set_xlabel("Time (s)")
        self.axes.set_ylabel(ylabel)
        self._style_plot_axes()
        self.axes.legend(facecolor=BG, edgecolor=MUTED, labelcolor=TEXT)
        self.canvas.draw_idle()

    def _style_plot_axes(self):
        self.axes.set_facecolor(PANEL)
        self.axes.tick_params(colors=TEXT)
        self.axes.xaxis.label.set_color(TEXT)
        self.axes.yaxis.label.set_color(TEXT)
        self.axes.title.set_color(TEXT)
        for spine in self.axes.spines.values():
            spine.set_color(MUTED)
        self.axes.grid(True, color=MUTED, alpha=0.25)

    def remove_selected(self):
        for item in self.result_table.selection():
            self.result_table.delete(item)

    def clear_results(self):
        for item in self.result_table.get_children():
            self.result_table.delete(item)
        self.refresh_plot()
        self.set_status("Results cleared.")

    def reset_defaults(self):
        defaults = {"width":"10", "height":"10", "radius":"1", "a":"0.5",
            "thickness":"1e-6", "resolution":"40", "theta_x":"0", "theta_y":"0",
            "pos_x":"0", "pos_y":"0", "pos_z":"0", "vel_x":"0", "vel_y":"0",
            "vel_z":"0", "omega_x":"0", "omega_y":"0", "omega_z":"0",
            "I0":"1e9", "w0":"1", "lambdaLaser":"1e-6", "beam_x":"0", "beam_y":"0",
            "steps":"10", "dt":"0.01", "max_runs":"500"}
        for key, value in defaults.items():
            self._var(key, value).set(value)
        for name, var in self.geometry_vars.items():
            var.set(name == "Rectangle")
        for name, var in self.material_vars.items():
            var.set(name == "Graphene")
        self.set_status("Defaults restored.")


if __name__ == "__main__":
    root = tk.Tk()
    app = NyxInterface(root)
    root.mainloop()
