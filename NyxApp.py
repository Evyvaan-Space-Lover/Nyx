import tkinter as tk 
from tkinter import messagebox
from Runner import NyxInterface
from nyx_runner import run_parameter_sweep, run_simulation, export_results

def main():
    root = tk.Tk()
    app = NyxInterface(
        root,
        on_run=lambda settings: run_simulation(app, settings),
        on_sweep=lambda settings: run_parameter_sweep(app, settings),
        on_export=lambda kind, settings: export_results(kind, settings, app),
    )
    root.mainloop()

if __name__=="__main__":
    main()