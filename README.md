# El péndulo que se vacía

[Read this in English](README.en.md)

Simulador interactivo de un péndulo físico formado por una botella cilíndrica que pierde líquido por un orificio en el fondo. La simulación combina el vaciado según la ley de Torricelli con el balance de momento angular de un sistema de masa variable. En el caso base, el período de pequeña amplitud sube de 1,5973 s a 1,6989 s y baja a 1,5981 s tras el vaciado.

## Ejecutar el simulador

Requiere Python 3.8 o superior y un navegador con JavaScript. El lanzador usa la biblioteca estándar de Python.

```sh
git clone https://github.com/camilohincapie229-maker/pendulum.git
cd pendulum
python3 main.py
```

Tras verificar el modelo, el lanzador abre el simulador en el navegador. Si no se abre automáticamente, usa la URL que aparece en la terminal. El puerto inicial es 8000. Pulsa `Ctrl+C` en la terminal para detener el servidor.

| Comando | Acción |
| --- | --- |
| `python3 main.py --solo-verificar` | Ejecutar las comprobaciones numéricas y mostrar el resumen del caso base. |
| `python3 main.py --sin-navegador` | Iniciar el servidor y abrir la URL manualmente. |
| `python3 main.py --puerto 8080` | Buscar un puerto disponible a partir del 8080. |
| `python3 main.py --sin-verificar` | Iniciar el servidor directamente. |

También puedes abrir [web/index.html](web/index.html) directamente en el navegador. Mantén los tres archivos de `web/` juntos.

## Explorar el modelo

1. Ajusta la geometría de la botella, el llenado, la longitud de la cuerda, el ángulo inicial y el orificio.
2. Pulsa **Iniciar** y observa cómo cambian la masa, el centro de masa, el período, la amplitud y la energía.
3. Repite con un diámetro de orificio de **0 mm** para comparar con el caso de masa constante. Usa **Valores por defecto** para recuperar el caso base.

Las comprobaciones cubren el período sin vaciado, las derivadas de las propiedades inerciales, la aproximación adiabática de la amplitud, la comparación entre dos orificios y el límite de masa puntual.

## Estructura y documentación

| Ruta | Contenido |
| --- | --- |
| [main.py](main.py) | Lanzador, servidor HTTP local y comprobaciones numéricas. |
| [fisica.py](fisica.py) | Modelo de referencia en Python e integrador Runge–Kutta de cuarto orden. |
| [web/](web/) | Interfaz, animación, gráficas e implementación del modelo en JavaScript. |
| [Informe de laboratorio](documento/informe_laboratorio.pdf) · [LaTeX](documento/informe_laboratorio.tex) | Objetivo, metodología, implementación, resultados y discusión. |
| [Desarrollo teórico](documento/pendulo_vaciado.pdf) · [LaTeX](documento/pendulo_vaciado.tex) | Derivaciones y análisis del modelo físico. |

Para compilar los documentos desde la raíz, instala una distribución de LaTeX con `latexmk`, pdfLaTeX, Babel para español, TikZ/PGFPlots y los paquetes declarados en los archivos fuente.

```sh
latexmk -pdf -outdir=documento/build documento/informe_laboratorio.tex
latexmk -pdf -outdir=documento/build documento/pendulo_vaciado.tex
```
