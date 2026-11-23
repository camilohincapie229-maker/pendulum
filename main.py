#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py -- Lanzador del simulador «El péndulo que se vacía».

Hace tres cosas:
  1. Verifica numéricamente el modelo con el módulo de referencia fisica.py.
  2. Sirve la carpeta web/ en un servidor HTTP local.
  3. Imprime el enlace y abre el navegador.

Uso
---
    python3 main.py                  # verifica, sirve y abre el navegador
    python3 main.py --puerto 8080    # puerto concreto
    python3 main.py --sin-navegador  # no abrir el navegador automáticamente
    python3 main.py --solo-verificar # solo las pruebas, sin servidor

Requisitos: Python 3.8 o superior. No hace falta instalar nada.
"""

import argparse
import functools
import http.server
import math
import os
import socket
import socketserver
import sys
import threading
import webbrowser

RAIZ = os.path.dirname(os.path.abspath(__file__))
ARCHIVOS_WEB = ('index.html', 'estilo.css', 'simulador.js')

sys.path.insert(0, RAIZ)
import fisica  # noqa: E402


def encontrar_carpeta_web():
    """Localiza la carpeta que contiene index.html, estilo.css y simulador.js.

    Comprueba primero RAIZ/web (estructura recomendada). Si no está, y por
    error los tres archivos quedaron sueltos directamente junto a main.py
    (típico al descargarlos uno por uno desde un chat, sin conservar
    subcarpetas), los sirve igualmente desde RAIZ para no bloquear al
    usuario. Devuelve la ruta encontrada, o None si no aparecen en ningún
    sitio razonable.
    """
    candidatas = [os.path.join(RAIZ, 'web'), RAIZ]
    for carpeta in candidatas:
        if all(os.path.isfile(os.path.join(carpeta, f)) for f in ARCHIVOS_WEB):
            return carpeta
    return None


# ----------------------------------------------------------------------
#  Presentación en consola
# ----------------------------------------------------------------------
class C:
    """Códigos ANSI de color (se desactivan si la salida no es un terminal)."""
    activo = sys.stdout.isatty() and os.name != 'nt'
    AZUL = '\033[94m' if activo else ''
    VERDE = '\033[92m' if activo else ''
    ROJO = '\033[91m' if activo else ''
    AMAR = '\033[93m' if activo else ''
    GRIS = '\033[90m' if activo else ''
    NEG = '\033[1m' if activo else ''
    FIN = '\033[0m' if activo else ''


def titulo(texto):
    print('\n' + C.NEG + texto + C.FIN)
    print(C.GRIS + '-' * 66 + C.FIN)


# ----------------------------------------------------------------------
#  Verificación del modelo
# ----------------------------------------------------------------------
def verificar():
    """Ejecuta las pruebas de consistencia física. Devuelve True si pasan."""
    titulo('Verificación del modelo (fisica.py)')
    pruebas = []

    # -- 1. Sin vaciado: período de amplitud finita --------------------
    p = fisica.Parametros(da=0.0, theta0_deg=12.0)
    r = fisica.integrar(p, t_final=30.0, dt=5e-4, muestreo=1.0)
    medido = sum(r['periodo']) / len(r['periodo'])
    teorico = fisica.periodo_amplitud(p.h0, p, p.theta0)
    err = abs(medido / teorico - 1.0)
    pruebas.append((err < 2e-4,
                    'Sin vaciado: T = T0(1 + Θ²/16 + 11Θ⁴/3072)',
                    'medido %.6f s · teórico %.6f s · error %.1e' % (medido, teorico, err)))

    # -- 2. Identidades dS/dh y dI/dh ----------------------------------
    p = fisica.Parametros()
    h, e = 0.12, 1e-6
    dS = (fisica.momento_estatico(h + e, p) - fisica.momento_estatico(h - e, p)) / (2 * e)
    dI = (fisica.inercia(h + e, p) - fisica.inercia(h - e, p)) / (2 * e)
    dS_t = p.rho * p.A * (p.s - h)
    dI_t = p.rho * p.A * ((p.s - h) ** 2 + p.R ** 2 / 4)
    err = max(abs(dS / dS_t - 1), abs(dI / dI_t - 1))
    pruebas.append((err < 1e-6,
                    'Identidades dS/dh = ρA(s−h) y dI/dh = ρA[(s−h)²+R²/4]',
                    'error %.1e' % err))

    # -- 3. Ley adiabática frente a integración numérica ---------------
    p = fisica.Parametros(da=0.008, theta0_deg=12.0)
    r = fisica.integrar(p, t_final=0.995 * p.tv, dt=1e-3, muestreo=1.0)
    num = r['amplitud'][-1] / p.theta0
    wkb = fisica.amplitud_wkb(p)
    err = abs(num / wkb - 1.0)
    pruebas.append((err < 0.01,
                    'Ley adiabática Θ ∝ (S·I)^(−1/4)·exp(−½∫ṁ s²/I)',
                    'numérico %.5f · WKB %.5f · error %.2f %%' % (num, wkb, err * 100)))

    # -- 4. Universalidad respecto del área del orificio ---------------
    factores = []
    for da in (0.006, 0.012):
        p = fisica.Parametros(da=da, theta0_deg=12.0)
        r = fisica.integrar(p, t_final=0.995 * p.tv, dt=1e-3, muestreo=1.0)
        factores.append(r['amplitud'][-1] / p.theta0)
    err = abs(factores[0] / factores[1] - 1.0)
    pruebas.append((err < 0.02,
                    'El factor de amplitud no depende del área del orificio',
                    'd=6 mm → %.5f · d=12 mm → %.5f · dif %.2f %%'
                    % (factores[0], factores[1], err * 100)))

    # -- 5. Límite de masa puntual -------------------------------------
    p = fisica.Parametros(H=2e-3, R=1e-3, h0=2e-3, mb=1e-4, Lc=0.8,
                          da=2e-4, theta0_deg=12.0)
    r = fisica.integrar(p, t_final=0.995 * p.tv, dt=5e-4, muestreo=1.0)
    dT = abs(r['periodo'][-1] / r['periodo'][0] - 1.0)
    dA = abs(r['amplitud'][-1] / r['amplitud'][0] - 1.0)
    pruebas.append((dT < 5e-3 and dA < 5e-3,
                    'Masa puntual: período y amplitud constantes',
                    'ΔT %.3f %% · ΔΘ %.3f %%' % (dT * 100, dA * 100)))

    todo_ok = True
    for ok, texto, detalle in pruebas:
        marca = (C.VERDE + '  OK  ' + C.FIN) if ok else (C.ROJO + ' FALLA' + C.FIN)
        print('%s  %s' % (marca, texto))
        print('        ' + C.GRIS + detalle + C.FIN)
        todo_ok = todo_ok and ok
    return todo_ok


def resumen_caso_base():
    """Imprime las magnitudes características del caso por defecto."""
    titulo('Caso por defecto (botella de 1.26 L, orificio de 5 mm)')
    p = fisica.Parametros()
    ext = fisica.instante_periodo_maximo(p)
    filas = [
        ('Tiempo de vaciado  t_v', '%.2f s' % p.tv),
        ('Período inicial    T(0)', '%.4f s' % fisica.periodo_lineal(p.h0, p)),
        ('Período final      T(t_v)', '%.4f s' % fisica.periodo_lineal(0.0, p)),
    ]
    if ext:
        h_e, t_e, T_e = ext
        filas += [
            ('Período máximo', '%.4f s  (+%.2f %%)'
             % (T_e, 100 * (T_e / fisica.periodo_lineal(p.h0, p) - 1))),
            ('  ocurre en', 't = %.1f s = %.2f t_v,  h* = %.1f mm' % (t_e, t_e / p.tv, h_e * 1000)),
            ('  se cumple', 's − h* = %.4f m = L_eq' % (p.s - h_e)),
        ]
    filas += [
        ('Integral de frenado  K', '%.5f' % fisica.integral_frenado(p)),
        ('Factor de amplitud', '×%.4f  (independiente del orificio)' % fisica.amplitud_wkb(p)),
    ]
    for etiqueta, valor in filas:
        print('  %-28s %s%s%s' % (etiqueta, C.AZUL, valor, C.FIN))


# ----------------------------------------------------------------------
#  Servidor
# ----------------------------------------------------------------------
class Manejador(http.server.SimpleHTTPRequestHandler):
    """Sirve web/ sin caché y en silencio."""

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, must-revalidate')
        super().end_headers()

    def log_message(self, formato, *args):
        pass                      # no ensuciar la consola


class Servidor(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


def puerto_libre(preferido):
    """Devuelve el primer puerto disponible a partir del preferido."""
    for p in range(preferido, preferido + 40):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', p)) != 0:
                return p
    raise RuntimeError('No se encontró ningún puerto libre entre %d y %d'
                       % (preferido, preferido + 40))


def servir(puerto, abrir_navegador=True):
    carpeta_web = encontrar_carpeta_web()
    if carpeta_web is None:
        print(C.ROJO + 'ERROR: no encuentro index.html, estilo.css y simulador.js.' + C.FIN)
        print('Busqué en:')
        print('  - ' + os.path.join(RAIZ, 'web'))
        print('  - ' + RAIZ)
        print()
        print('Contenido de ' + RAIZ + ':')
        try:
            for nombre in sorted(os.listdir(RAIZ)):
                marca = '/' if os.path.isdir(os.path.join(RAIZ, nombre)) else ''
                print('    ' + nombre + marca)
        except OSError:
            pass
        print()
        print('Es probable que al descargar los archivos uno por uno se haya')
        print('perdido la carpeta web/. Crea una carpeta "web" junto a main.py')
        print('y mueve dentro index.html, estilo.css y simulador.js.')
        print('La estructura completa está descrita en README.md.')
        return 1
    if carpeta_web == RAIZ:
        print(C.AMAR + 'Aviso: index.html, estilo.css y simulador.js están sueltos junto a')
        print('main.py en vez de dentro de una carpeta web/. Los sirvo igualmente,')
        print('pero conviene reorganizarlos como indica README.md.' + C.FIN)

    puerto = puerto_libre(puerto)
    manejador = functools.partial(Manejador, directory=carpeta_web)
    url = 'http://localhost:%d/' % puerto

    with Servidor(('127.0.0.1', puerto), manejador) as httpd:
        titulo('Servidor en marcha')
        print('  Abre este enlace en tu navegador:\n')
        print('      ' + C.NEG + C.VERDE + url + C.FIN + '\n')
        print(C.GRIS + '  Pulsa Ctrl+C para detener el servidor.' + C.FIN)
        if abrir_navegador:
            threading.Timer(0.8, lambda: webbrowser.open(url)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print('\n' + C.AMAR + 'Servidor detenido.' + C.FIN)
    return 0


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description='Simulador interactivo del péndulo que se vacía.')
    ap.add_argument('--puerto', type=int, default=8000,
                    help='puerto de escucha (por defecto 8000)')
    ap.add_argument('--sin-navegador', action='store_true',
                    help='no abrir el navegador automáticamente')
    ap.add_argument('--solo-verificar', action='store_true',
                    help='ejecutar solo las pruebas del modelo')
    ap.add_argument('--sin-verificar', action='store_true',
                    help='saltar las pruebas y arrancar el servidor directamente')
    args = ap.parse_args()

    print(C.NEG + C.AZUL + """
   ____  _____ _   _ ____  _   _ _     ___
  |  _ \\| ____| \\ | |  _ \\| | | | |   / _ \\    El péndulo que se vacía
  | |_) |  _| |  \\| | | | | | | | |  | | | |   péndulo físico de masa variable
  |  __/| |___| |\\  | |_| | |_| | |__| |_| |   vaciado tipo Torricelli
  |_|   |_____|_| \\_|____/ \\___/|_____\\___/
""" + C.FIN)

    if not args.sin_verificar:
        if not verificar():
            print('\n' + C.ROJO + 'Alguna prueba ha fallado: revisa fisica.py.' + C.FIN)
            return 2
        resumen_caso_base()

    if args.solo_verificar:
        print()
        return 0
    return servir(args.puerto, not args.sin_navegador)


if __name__ == '__main__':
    sys.exit(main())
