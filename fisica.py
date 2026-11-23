# -*- coding: utf-8 -*-
"""
fisica.py -- Modelo de referencia del pendulo fisico que se vacia.

Implementacion en Python puro (solo biblioteca estandar) de las mismas
ecuaciones que ejecuta el simulador web en web/simulador.js.  Sirve para:

  * validar numericamente la implementacion en JavaScript,
  * reproducir las figuras del documento LaTeX,
  * usarse como modulo desde otros scripts.

ECUACIONES (ver documento/pendulo_vaciado.pdf)
---------------------------------------------
Geometria:      A = pi R^2,   s = Lc + H,   db = Lc + H/2
Vaciado:        h(t) = h0 (1 - t/tv)^2,   tv = A/(Cd a) sqrt(2 h0/g)
Momento estatico:
                S(h) = mb db + rho A h (s - h/2)
Momento de inercia respecto del pivote:
                I(h) = Ib + mb db^2 + rho A h [R^2/4 + h^2/12 + (s - h/2)^2]
Frenado por evacuacion de momento angular:
                Gamma = rho A |dh/dt| [ h(2s - h) - R^2/4 ]
Ecuacion de movimiento:
                I(t) theta'' + Gamma(t) theta' + g S(t) sin(theta) = 0
"""

import math

G_STD = 9.81      # m/s^2
RHO_AGUA = 1000.0  # kg/m^3


class Parametros:
    """Parametros geometricos y fisicos del sistema."""

    def __init__(self, H=0.250, R=0.040, mb=0.040, Lc=0.500,
                 h0=0.250, da=0.005, Cd=0.60, theta0_deg=12.0,
                 rho=RHO_AGUA, g=G_STD):
        self.H = H                      # altura interior del cilindro [m]
        self.R = R                      # radio interior [m]
        self.mb = mb                    # masa de la botella vacia [kg]
        self.Lc = Lc                    # longitud de la cuerda [m]
        self.h0 = min(h0, H)            # altura inicial de liquido [m]
        self.da = da                    # diametro del orificio [m]
        self.Cd = Cd                    # coeficiente de descarga [-]
        self.theta0 = math.radians(theta0_deg)
        self.rho = rho
        self.g = g

    # ---------------- magnitudes derivadas (constantes) ----------------
    @property
    def A(self):
        """Area de la seccion interior del cilindro [m^2]."""
        return math.pi * self.R ** 2

    @property
    def a(self):
        """Area del orificio de salida [m^2]."""
        return math.pi * self.da ** 2 / 4.0

    @property
    def s(self):
        """Distancia pivote -> fondo (orificio) [m]."""
        return self.Lc + self.H

    @property
    def db(self):
        """Distancia pivote -> CM de la botella vacia [m]."""
        return self.Lc + self.H / 2.0

    @property
    def Ib(self):
        """Inercia propia de la botella (cascaron cilindrico) [kg m^2]."""
        return self.mb * (self.R ** 2 / 2.0 + self.H ** 2 / 12.0)

    @property
    def Ibot(self):
        """Inercia de la botella vacia respecto del pivote [kg m^2]."""
        return self.Ib + self.mb * self.db ** 2

    @property
    def tv(self):
        """Tiempo de vaciado completo [s]; inf si no hay orificio."""
        if self.a <= 0.0 or self.h0 <= 0.0 or self.Cd <= 0.0:
            return float('inf')
        return (self.A / (self.Cd * self.a)) * math.sqrt(2.0 * self.h0 / self.g)


# ----------------------------------------------------------------------
#  Funciones del modelo
# ----------------------------------------------------------------------
def altura(t, p):
    """Altura de la columna de liquido h(t) [m] -- ley de Torricelli."""
    tv = p.tv
    if tv == float('inf'):
        return p.h0
    u = 1.0 - t / tv
    return p.h0 * u * u if u > 0.0 else 0.0


def altura_punto(t, p):
    """Derivada dh/dt [m/s] (negativa o cero)."""
    tv = p.tv
    if tv == float('inf'):
        return 0.0
    u = 1.0 - t / tv
    return -2.0 * p.h0 * u / tv if u > 0.0 else 0.0


def momento_estatico(h, p):
    """S(h) = M*l : momento de primer orden respecto del pivote [kg m]."""
    return p.mb * p.db + p.rho * p.A * h * (p.s - h / 2.0)


def inercia(h, p):
    """I_O(h): momento de inercia respecto del pivote [kg m^2]."""
    return (p.Ibot + p.rho * p.A * h *
            (p.R ** 2 / 4.0 + h * h / 12.0 + (p.s - h / 2.0) ** 2))


def masa_total(h, p):
    """M(h) = mb + masa de liquido [kg]."""
    return p.mb + p.rho * p.A * h


def gamma_frenado(h, hp, p):
    """Gamma = dI/dt + (dm_out/dt) s^2  [kg m^2/s].  Positivo = frena."""
    return p.rho * p.A * abs(hp) * (h * (2.0 * p.s - h) - p.R ** 2 / 4.0)


def centro_masa(h, p):
    """Distancia pivote -> centro de masa del conjunto [m]."""
    return momento_estatico(h, p) / masa_total(h, p)


def longitud_equivalente(h, p):
    """L_eq = I_O/S : longitud del pendulo simple isocrono [m]."""
    return inercia(h, p) / momento_estatico(h, p)


def periodo_lineal(h, p):
    """Periodo de oscilaciones infinitesimales [s]."""
    return 2.0 * math.pi * math.sqrt(longitud_equivalente(h, p) / p.g)


def periodo_amplitud(h, p, amp):
    """Periodo con correccion de amplitud finita (serie hasta amp^4)."""
    return periodo_lineal(h, p) * (1.0 + amp ** 2 / 16.0 +
                                   11.0 * amp ** 4 / 3072.0)


def energias(t, th, om, p):
    """Devuelve (cinetica, potencial, mecanica) [J].

    U se mide respecto de la posicion de equilibrio instantanea (theta=0),
    de modo que U = g S (1 - cos theta).  No se incluye el termino
    (1/2) M l'^2 del arrastre longitudinal del CM: es funcion solo del
    tiempo y no interviene en la ecuacion de Euler-Lagrange.
    """
    h = altura(t, p)
    K = 0.5 * inercia(h, p) * om * om
    U = p.g * momento_estatico(h, p) * (1.0 - math.cos(th))
    return K, U, K + U


def aceleracion(t, th, om, p):
    """theta'' a partir de la ecuacion de movimiento corregida."""
    h = altura(t, p)
    hp = altura_punto(t, p)
    I = inercia(h, p)
    S = momento_estatico(h, p)
    Gm = gamma_frenado(h, hp, p)
    return -(Gm * om + p.g * S * math.sin(th)) / I


# ----------------------------------------------------------------------
#  Integrador Runge-Kutta 4 con deteccion de periodos y amplitudes
# ----------------------------------------------------------------------
def integrar(p, t_final=None, dt=1.0e-3, muestreo=0.02):
    """Integra el sistema y devuelve un diccionario con las series.

    Claves: t, theta, omega, h, mw, cm, Leq, K, U, E,
            t_periodo, periodo, t_amp, amplitud
    """
    if t_final is None:
        t_final = 1.05 * p.tv if p.tv != float('inf') else 20.0 * periodo_lineal(p.h0, p)

    t, th, om = 0.0, p.theta0, 0.0
    n = int(round(t_final / dt))
    cada = max(1, int(round(muestreo / dt)))

    out = {k: [] for k in ('t', 'theta', 'omega', 'h', 'mw', 'cm', 'Leq',
                           'K', 'U', 'E', 't_periodo', 'periodo',
                           't_amp', 'amplitud')}
    cruce_previo = None
    amp_actual = abs(th)

    def registrar():
        h = altura(t, p)
        K, U, E = energias(t, th, om, p)
        out['t'].append(t);          out['theta'].append(th)
        out['omega'].append(om);     out['h'].append(h)
        out['mw'].append(p.rho * p.A * h)
        out['cm'].append(centro_masa(h, p))
        out['Leq'].append(longitud_equivalente(h, p))
        out['K'].append(K); out['U'].append(U); out['E'].append(E)

    registrar()
    for i in range(n):
        th_ant = th
        # --- RK4 ---
        k1th, k1om = om, aceleracion(t, th, om, p)
        k2th, k2om = om + 0.5 * dt * k1om, aceleracion(t + 0.5 * dt, th + 0.5 * dt * k1th, om + 0.5 * dt * k1om, p)
        k3th, k3om = om + 0.5 * dt * k2om, aceleracion(t + 0.5 * dt, th + 0.5 * dt * k2th, om + 0.5 * dt * k2om, p)
        k4th, k4om = om + dt * k3om,       aceleracion(t + dt,       th + dt * k3th,       om + dt * k3om, p)
        th += dt * (k1th + 2 * k2th + 2 * k3th + k4th) / 6.0
        om += dt * (k1om + 2 * k2om + 2 * k3om + k4om) / 6.0
        t = (i + 1) * dt

        amp_actual = max(amp_actual, abs(th))
        # cruce por cero ascendente -> medida de periodo y amplitud
        if th_ant < 0.0 <= th:
            frac = -th_ant / (th - th_ant)
            t_cruce = t - dt + frac * dt
            if cruce_previo is not None:
                out['periodo'].append(t_cruce - cruce_previo)
                out['t_periodo'].append(0.5 * (t_cruce + cruce_previo))
                out['amplitud'].append(amp_actual)
                out['t_amp'].append(0.5 * (t_cruce + cruce_previo))
            cruce_previo = t_cruce
            amp_actual = 0.0

        if (i + 1) % cada == 0:
            registrar()
    return out


# ----------------------------------------------------------------------
#  Prediccion analitica WKB (regimen adiabatico)
# ----------------------------------------------------------------------
def integral_frenado(p, n=20000):
    """K = int_0^{tv} (Gamma/I) dt, integral puramente geometrica en h.

    K = rho A int_0^{h0} [h(2s-h) - R^2/4] / I(h)  dh
    Es independiente del area del orificio, de Cd y de g.
    """
    if p.h0 <= 0.0:
        return 0.0
    acc, dh = 0.0, p.h0 / n
    for i in range(n):                      # regla del trapecio
        h1, h2 = i * dh, (i + 1) * dh
        f1 = (h1 * (2 * p.s - h1) - p.R ** 2 / 4.0) / inercia(h1, p)
        f2 = (h2 * (2 * p.s - h2) - p.R ** 2 / 4.0) / inercia(h2, p)
        acc += 0.5 * (f1 + f2) * dh
    return p.rho * p.A * acc


def amplitud_wkb(p):
    """Factor de amplitud tras un vaciado completo (prediccion adiabatica).

    Theta_f/Theta_0 = sqrt(w0/wf) * exp(-K/2)
                    = (S0 I0 / Sf If)^(1/4) * exp(-(s^2/2) int dm_out/I)
    """
    w0 = math.sqrt(p.g * momento_estatico(p.h0, p) / inercia(p.h0, p))
    wf = math.sqrt(p.g * momento_estatico(0.0, p) / inercia(0.0, p))
    return math.sqrt(w0 / wf) * math.exp(-0.5 * integral_frenado(p))


def instante_periodo_maximo(p):
    """Busca h* tal que (s-h*)^2 + R^2/4 = (s-h*) L_eq(h*)  [Sec. 9.2].

    Devuelve (h_estrella, t_estrella, periodo_maximo) o None si no existe
    extremo interior en 0 < h < h0.
    """
    def f(h):
        return ((p.s - h) ** 2 + p.R ** 2 / 4.0) / (p.s - h) - longitud_equivalente(h, p)

    lo, hi = 1e-9, p.h0
    if f(lo) * f(hi) > 0.0:
        return None
    for _ in range(200):                     # biseccion
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    h_est = 0.5 * (lo + hi)
    tv = p.tv
    t_est = tv * (1.0 - math.sqrt(h_est / p.h0)) if tv != float('inf') else float('inf')
    return h_est, t_est, periodo_lineal(h_est, p)
