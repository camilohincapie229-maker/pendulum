# The draining pendulum

[Leer en español](README.md)

An interactive simulator of a physical pendulum made from a cylindrical bottle that loses liquid through an outlet at its base. The simulation combines Torricelli's law with the angular momentum balance for a variable-mass system. In the default case, the small-angle period rises from 1.5973 s to 1.6989 s and drops to 1.5981 s after draining.

## Run the simulator

Requires Python 3.8 or newer and a browser with JavaScript. The launcher uses the Python standard library.

```sh
git clone https://github.com/camilohincapie229-maker/pendulum.git
cd pendulum
python3 main.py
```

After checking the model, the launcher opens the simulator in your browser. If it does not open automatically, use the URL printed in the terminal. The starting port is 8000. Press `Ctrl+C` in the terminal to stop the server.

| Command | Action |
| --- | --- |
| `python3 main.py --solo-verificar` | Run the numerical checks and print the default-case summary. |
| `python3 main.py --sin-navegador` | Start the server and open the URL manually. |
| `python3 main.py --puerto 8080` | Look for an available port starting at 8080. |
| `python3 main.py --sin-verificar` | Start the server immediately. |

You can also open [web/index.html](web/index.html) directly in your browser. Keep all three files in `web/` together. The simulator interface and terminal output are in Spanish.

## Explore the model

1. Adjust the bottle geometry, fill level, string length, initial angle, and outlet.
2. Select **Iniciar** to start and watch the mass, center of mass, period, amplitude, and energy change.
3. Repeat with a **0 mm** outlet diameter to compare against the constant-mass case. Select **Valores por defecto** to restore the defaults.

The checks cover the period without drainage, derivatives of the inertial properties, the adiabatic approximation for amplitude, a comparison between two outlet sizes, and the point-mass limit.

## Repository layout and documentation

| Path | Contents |
| --- | --- |
| [main.py](main.py) | Launcher, local HTTP server, and numerical checks. |
| [fisica.py](fisica.py) | Python reference model and fourth-order Runge–Kutta integrator. |
| [web/](web/) | User interface, animation, plots, and JavaScript implementation. |
| [Lab report](documento/informe_laboratorio.pdf) · [LaTeX](documento/informe_laboratorio.tex) | Objectives, methodology, implementation, results, and discussion in Spanish. |
| [Theoretical background](documento/pendulo_vaciado.pdf) · [LaTeX](documento/pendulo_vaciado.tex) | Derivations and analysis of the physical model in Spanish. |

To build the documents from the repository root, use a LaTeX distribution with `latexmk`, pdfLaTeX, Spanish Babel, TikZ/PGFPlots, and the packages declared in the source files.

```sh
latexmk -pdf -outdir=documento/build documento/informe_laboratorio.tex
latexmk -pdf -outdir=documento/build documento/pendulo_vaciado.tex
```
