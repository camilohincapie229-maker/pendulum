/* =====================================================================
   El péndulo que se vacía — motor de simulación
   ---------------------------------------------------------------------
   Física idéntica a fisica.py (implementación de referencia en Python).

        I_O(t) θ¨ + Γ(t) θ˙ + g S(t) sen θ = 0

        h(t)  = h0 (1 − t/tv)²                 (Torricelli, tanque cilíndrico)
        tv    = A/(Cd·a) · √(2h0/g)
        S(h)  = mb·db + ρAh(s − h/2)
        I_O(h)= Ib + mb·db² + ρAh[R²/4 + h²/12 + (s − h/2)²]
        Γ     = ρA|h˙|[h(2s − h) − R²/4]       (evacuación de mom. angular)

   Sin dependencias externas: todo se dibuja sobre <canvas>.
   ===================================================================== */
'use strict';

const G = 9.81;

/* =====================================================================
   1. FÍSICA
   ===================================================================== */

/** Calcula las constantes geométricas derivadas de los parámetros. */
function derivar(p) {
  const A  = Math.PI * p.R * p.R;          // sección interior
  const a  = Math.PI * p.da * p.da / 4;    // área del orificio
  const s  = p.Lc + p.H;                   // pivote → fondo
  const db = p.Lc + p.H / 2;               // pivote → CM botella vacía
  const Ib = p.mb * (p.R * p.R / 2 + p.H * p.H / 12);   // cascarón cilíndrico
  const Ibot = Ib + p.mb * db * db;        // botella vacía respecto del pivote
  const tv = (a > 0 && p.h0 > 0 && p.Cd > 0)
           ? (A / (p.Cd * a)) * Math.sqrt(2 * p.h0 / G)
           : Infinity;
  return { A, a, s, db, Ib, Ibot, tv };
}

/** Altura de la columna de líquido h(t). */
function altura(t, p, d) {
  if (!isFinite(d.tv)) return p.h0;
  const u = 1 - t / d.tv;
  return u > 0 ? p.h0 * u * u : 0;
}

/** Derivada dh/dt (≤ 0). */
function alturaPunto(t, p, d) {
  if (!isFinite(d.tv)) return 0;
  const u = 1 - t / d.tv;
  return u > 0 ? -2 * p.h0 * u / d.tv : 0;
}

/** Momento estático S = M·ℓ respecto del pivote. */
function momEstatico(h, p, d) {
  return p.mb * d.db + p.rho * d.A * h * (d.s - h / 2);
}

/** Momento de inercia respecto del pivote. */
function inercia(h, p, d) {
  return d.Ibot + p.rho * d.A * h *
         (p.R * p.R / 4 + h * h / 12 + Math.pow(d.s - h / 2, 2));
}

/** Masa total instantánea. */
function masaTotal(h, p, d) { return p.mb + p.rho * d.A * h; }

/** Coeficiente de frenado Γ = İ_O + ṁ_out·s². */
function gammaFrenado(h, hp, p, d) {
  return p.rho * d.A * Math.abs(hp) * (h * (2 * d.s - h) - p.R * p.R / 4);
}

/** Aceleración angular a partir de la ecuación de movimiento corregida. */
function aceleracion(t, th, om, p, d) {
  const h  = altura(t, p, d);
  const hp = alturaPunto(t, p, d);
  const I  = inercia(h, p, d);
  const S  = momEstatico(h, p, d);
  const Gm = gammaFrenado(h, hp, p, d);
  return -(Gm * om + G * S * Math.sin(th)) / I;
}

/** Un paso de Runge–Kutta de cuarto orden. */
function pasoRK4(t, th, om, dt, p, d) {
  const k1t = om,                     k1o = aceleracion(t, th, om, p, d);
  const k2t = om + 0.5 * dt * k1o,    k2o = aceleracion(t + 0.5 * dt, th + 0.5 * dt * k1t, om + 0.5 * dt * k1o, p, d);
  const k3t = om + 0.5 * dt * k2o,    k3o = aceleracion(t + 0.5 * dt, th + 0.5 * dt * k2t, om + 0.5 * dt * k2o, p, d);
  const k4t = om + dt * k3o,          k4o = aceleracion(t + dt, th + dt * k3t, om + dt * k3o, p, d);
  return [th + dt * (k1t + 2 * k2t + 2 * k3t + k4t) / 6,
          om + dt * (k1o + 2 * k2o + 2 * k3o + k4o) / 6];
}

/** Longitud equivalente L_eq = I_O/S (centro de oscilación). */
function longEquiv(h, p, d) { return inercia(h, p, d) / momEstatico(h, p, d); }

/** Período lineal instantáneo. */
function periodoLineal(h, p, d) { return 2 * Math.PI * Math.sqrt(longEquiv(h, p, d) / G); }

/** Período con corrección de amplitud finita (serie hasta Θ⁴). */
function periodoAmplitud(h, p, d, amp) {
  return periodoLineal(h, p, d) * (1 + amp * amp / 16 + 11 * Math.pow(amp, 4) / 3072);
}

/** Integral de frenado K = ∫Γ/I dt, puramente geométrica en h. */
function integralFrenado(p, d, n) {
  n = n || 4000;
  if (p.h0 <= 0) return 0;
  const dh = p.h0 / n;
  let acc = 0;
  for (let i = 0; i < n; i++) {
    const h1 = i * dh, h2 = (i + 1) * dh;
    const f1 = (h1 * (2 * d.s - h1) - p.R * p.R / 4) / inercia(h1, p, d);
    const f2 = (h2 * (2 * d.s - h2) - p.R * p.R / 4) / inercia(h2, p, d);
    acc += 0.5 * (f1 + f2) * dh;
  }
  return p.rho * d.A * acc;
}

/** Predicción adiabática (WKB) del factor de amplitud tras vaciado completo. */
function factorAmplitudWKB(p, d) {
  const w0 = Math.sqrt(G * momEstatico(p.h0, p, d) / inercia(p.h0, p, d));
  const wf = Math.sqrt(G * momEstatico(0, p, d) / inercia(0, p, d));
  return Math.sqrt(w0 / wf) * Math.exp(-0.5 * integralFrenado(p, d));
}

/* =====================================================================
   2. GRÁFICAS (canvas puro, sin librerías)
   ===================================================================== */
class Grafica {
  constructor(canvas, opciones) {
    this.cv = canvas;
    this.ctx = canvas.getContext('2d');
    this.op = Object.assign({
      ejeX: 't (s)', ejeY: '', series: [], yMin: null, yMax: null, xMin: 0
    }, opciones || {});
    this.series = this.op.series;
  }

  ajustarTamano() {
    const dpr = window.devicePixelRatio || 1;
    const w = this.cv.clientWidth, h = this.cv.clientHeight;
    if (this.cv.width !== Math.round(w * dpr) || this.cv.height !== Math.round(h * dpr)) {
      this.cv.width = Math.round(w * dpr);
      this.cv.height = Math.round(h * dpr);
    }
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.w = w; this.h = h;
  }

  /** Escalón "bonito" para las marcas de los ejes. */
  static paso(rango, objetivo) {
    if (!(rango > 0)) return 1;
    const bruto = rango / objetivo;
    const exp = Math.floor(Math.log10(bruto));
    const base = Math.pow(10, exp);
    const n = bruto / base;
    const mult = n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10;
    return mult * base;
  }

  static formatear(v, paso) {
    const dec = Math.max(0, Math.min(6, -Math.floor(Math.log10(paso)) + (paso < 1 ? 0 : 0)));
    if (Math.abs(v) >= 1e5 || (v !== 0 && Math.abs(v) < 1e-3)) return v.toExponential(1);
    return v.toFixed(Math.max(0, dec));
  }

  dibujar() {
    this.ajustarTamano();
    const ctx = this.ctx, W = this.w, H = this.h;
    ctx.clearRect(0, 0, W, H);

    const m = { l: 56, r: 12, t: 26, b: 30 };
    const px = W - m.l - m.r, py = H - m.t - m.b;
    if (px <= 10 || py <= 10) return;

    // --- rangos ---
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    let hayDatos = false;
    for (const s of this.series) {
      for (let i = 0; i < s.datos.length; i++) {
        const d = s.datos[i];
        if (!isFinite(d[0]) || !isFinite(d[1])) continue;
        hayDatos = true;
        if (d[0] < xmin) xmin = d[0];
        if (d[0] > xmax) xmax = d[0];
        if (d[1] < ymin) ymin = d[1];
        if (d[1] > ymax) ymax = d[1];
      }
    }
    if (!hayDatos) {
      ctx.fillStyle = '#5c6879';
      ctx.font = '13px system-ui, sans-serif';
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText('Pulsa «Iniciar» para generar datos', W / 2, H / 2);
      return;
    }
    if (this.op.xMin !== null) xmin = Math.min(xmin, this.op.xMin);
    if (xmax - xmin < 1e-9) xmax = xmin + 1;
    if (this.op.yMin !== null) ymin = Math.min(ymin, this.op.yMin);
    if (this.op.yMax !== null) ymax = Math.max(ymax, this.op.yMax);
    let margenY = (ymax - ymin) * 0.10;
    if (margenY < 1e-12) margenY = Math.max(Math.abs(ymax) * 0.1, 1e-6);
    ymin -= margenY; ymax += margenY;

    const X = v => m.l + (v - xmin) / (xmax - xmin) * px;
    const Y = v => m.t + py - (v - ymin) / (ymax - ymin) * py;

    // --- rejilla ---
    const pasoX = Grafica.paso(xmax - xmin, 6), pasoY = Grafica.paso(ymax - ymin, 5);
    ctx.font = '10.5px system-ui, sans-serif';
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(255,255,255,.07)';
    ctx.fillStyle = '#7a8697';

    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (let v = Math.ceil(ymin / pasoY) * pasoY; v <= ymax; v += pasoY) {
      const y = Y(v);
      ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + px, y); ctx.stroke();
      ctx.fillText(Grafica.formatear(v, pasoY), m.l - 7, y);
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let v = Math.ceil(xmin / pasoX) * pasoX; v <= xmax; v += pasoX) {
      const x = X(v);
      ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + py); ctx.stroke();
      ctx.fillText(Grafica.formatear(v, pasoX), x, m.t + py + 6);
    }

    // marco
    ctx.strokeStyle = 'rgba(255,255,255,.16)';
    ctx.strokeRect(m.l + .5, m.t + .5, px, py);

    // etiquetas de eje
    ctx.fillStyle = '#8b97a8'; ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
    ctx.fillText(this.op.ejeX, m.l + px, H - 2);
    if (this.op.ejeY) {
      ctx.save(); ctx.translate(11, m.t + py / 2); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(this.op.ejeY, 0, 0); ctx.restore();
    }

    // --- series ---
    ctx.save();
    ctx.beginPath(); ctx.rect(m.l, m.t, px, py); ctx.clip();
    for (const s of this.series) {
      if (!s.datos.length) continue;
      ctx.strokeStyle = s.color; ctx.fillStyle = s.color;
      ctx.lineWidth = s.grosor || 1.7;
      ctx.setLineDash(s.trazo || []);
      if (s.puntos) {
        for (const d of s.datos) {
          ctx.beginPath();
          ctx.arc(X(d[0]), Y(d[1]), s.radio || 2.1, 0, 2 * Math.PI);
          ctx.fill();
        }
      } else {
        ctx.beginPath();
        let iniciado = false;
        // diezmado: no dibujar más de ~2 puntos por píxel
        const salto = Math.max(1, Math.floor(s.datos.length / (px * 2)));
        for (let i = 0; i < s.datos.length; i += salto) {
          const d = s.datos[i];
          if (!isFinite(d[1])) { iniciado = false; continue; }
          const x = X(d[0]), y = Y(d[1]);
          if (!iniciado) { ctx.moveTo(x, y); iniciado = true; } else ctx.lineTo(x, y);
        }
        const ult = s.datos[s.datos.length - 1];
        if (iniciado && isFinite(ult[1])) ctx.lineTo(X(ult[0]), Y(ult[1]));
        ctx.stroke();
      }
      ctx.setLineDash([]);
    }
    ctx.restore();

    // --- leyenda ---
    ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    let cx = m.l + 2;
    for (const s of this.series) {
      if (!s.nombre) continue;
      const anch = ctx.measureText(s.nombre).width;
      ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = 2.2;
      ctx.setLineDash(s.trazo || []);
      ctx.beginPath(); ctx.moveTo(cx, m.t - 13); ctx.lineTo(cx + 15, m.t - 13); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#aab6c6';
      ctx.fillText(s.nombre, cx + 20, m.t - 12.5);
      cx += anch + 34;
    }
  }
}

/* =====================================================================
   3. ESTADO DE LA APLICACIÓN
   ===================================================================== */
const POR_DEFECTO = {
  H: 0.250, D: 0.080, llenado: 100, rho: 1000,
  Lc: 0.500, mb: 0.040, th0: 12, da: 0.005, Cd: 0.60, vel: 6
};

const DT = 0.001;                 // paso de integración [s]
const MUESTREO = 0.04;            // intervalo de muestreo de las series [s]
const MAX_PUNTOS = 3000;          // tope por serie (se diezma al superarlo)

const app = {
  p: null, d: null,
  t: 0, th: 0, om: 0,
  corriendo: false,
  tMuestra: 0,
  ampActual: 0, ampMedida: 0, cruceAnt: null, periodoMedido: NaN,
  rastro: [],
  gotas: [],
  series: {
    periodoMed: [], periodoTeo: [],
    cm: [], leq: [], sup: [],
    mw: [], M: [],
    K: [], U: [], E: []
  }
};

/* ------------------------- lectura de controles ------------------------- */
const $ = id => document.getElementById(id);

function leerParametros() {
  const H  = +$('H').value;
  const R  = +$('D').value / 2;
  const h0 = H * (+$('llenado').value) / 100;
  return {
    H, R, h0,
    rho: +$('rho').value,
    Lc:  +$('Lc').value,
    mb:  +$('mb').value,
    da:  +$('da').value,
    Cd:  +$('Cd').value,
    theta0: (+$('th0').value) * Math.PI / 180
  };
}

function fmt(v, dec) { return (isFinite(v) ? v.toFixed(dec === undefined ? 3 : dec) : '—'); }

function actualizarEtiquetas() {
  const p = app.p, d = app.d;
  $('H_val').textContent   = (p.H * 100).toFixed(1) + ' cm';
  $('D_val').textContent   = (p.R * 200).toFixed(1) + ' cm';
  $('llenado_val').textContent = $('llenado').value + ' %';
  $('rho_val').textContent = p.rho + ' kg/m³';
  $('Lc_val').textContent  = (p.Lc * 100).toFixed(0) + ' cm';
  $('mb_val').textContent  = (p.mb * 1000).toFixed(0) + ' g';
  $('th0_val').textContent = $('th0').value + '°';
  $('da_val').textContent  = (p.da * 1000).toFixed(1) + ' mm';
  $('Cd_val').textContent  = p.Cd.toFixed(2);
  $('vel_val').textContent = (+$('vel').value) + '×';

  $('vol_val').textContent = (d.A * p.H * 1000).toFixed(2) + ' L';
  $('mw_ini').textContent  = (p.rho * d.A * p.h0).toFixed(3) + ' kg  ('
                           + (p.h0 * 100).toFixed(1) + ' cm)';
  $('s_val').textContent   = (d.s * 100).toFixed(1) + ' cm';
  $('a_val').textContent   = (d.a * 1e6).toFixed(1) + ' mm²';
  $('tv_val').textContent  = isFinite(d.tv) ? d.tv.toFixed(1) + ' s' : '∞ (sin orificio)';
}

/* ------------------------- reinicio ------------------------- */
function reiniciar() {
  app.p = leerParametros();
  app.d = derivar(app.p);
  app.t = 0;
  app.th = app.p.theta0;
  app.om = 0;
  app.tMuestra = 0;
  app.ampActual = Math.abs(app.th);
  app.ampMedida = Math.abs(app.th);
  app.cruceAnt = null;
  app.periodoMedido = NaN;
  app.rastro.length = 0;
  app.gotas.length = 0;
  for (const k in app.series) app.series[k].length = 0;
  actualizarEtiquetas();
  comprobarAvisos();
  muestrear();
  dibujarTodo();
}

function comprobarAvisos() {
  const p = app.p, d = app.d, av = $('aviso');
  const msgs = [];
  const T0 = periodoLineal(p.h0, p, d);
  if (isFinite(d.tv) && d.tv < 3 * T0) {
    msgs.push('Vaciado rápido: t<sub>v</sub> = ' + d.tv.toFixed(2) + ' s ≈ ' +
      (d.tv / T0).toFixed(1) + ' períodos. La aproximación adiabática deja de valer y el ' +
      'resultado depende de la fase inicial; la integración numérica sigue siendo válida.');
  }
  if (+$('th0').value > 35) {
    msgs.push('Amplitud grande: la corrección no lineal del período supera el 2.5 % y la ' +
      'hipótesis de superficie libre plana se degrada.');
  }
  if (p.da > 0 && p.da > 2 * p.R) {
    msgs.push('El orificio es mayor que la sección del envase: la ley de Torricelli no aplica.');
  }
  if (msgs.length) { av.innerHTML = msgs.join('<br>'); av.classList.remove('oculto'); }
  else av.classList.add('oculto');
}

/* ------------------------- muestreo de series ------------------------- */
function diezmar(arr) {
  if (arr.length <= MAX_PUNTOS) return;
  for (let i = 0, j = 0; j < arr.length; i++, j += 2) arr[i] = arr[j];
  arr.length = Math.ceil(arr.length / 2);
}

function muestrear() {
  const p = app.p, d = app.d, t = app.t;
  const h = altura(t, p, d);
  const I = inercia(h, p, d), S = momEstatico(h, p, d);
  const K = 0.5 * I * app.om * app.om;
  const U = G * S * (1 - Math.cos(app.th));
  const s = app.series;

  s.cm.push([t, momEstatico(h, p, d) / masaTotal(h, p, d)]);
  s.leq.push([t, I / S]);
  s.sup.push([t, d.s - h]);
  s.mw.push([t, p.rho * d.A * h]);
  s.M.push([t, masaTotal(h, p, d)]);
  s.K.push([t, K]); s.U.push([t, U]); s.E.push([t, K + U]);
  s.periodoTeo.push([t, periodoAmplitud(h, p, d, app.ampMedida)]);

  for (const k of ['cm', 'leq', 'sup', 'mw', 'M', 'K', 'U', 'E', 'periodoTeo']) diezmar(s[k]);
}

/* ------------------------- avance temporal ------------------------- */
function avanzar(dtReal) {
  const p = app.p, d = app.d;
  const factor = +$('vel').value;
  let objetivo = dtReal * factor;
  const maxPaso = 4500 * DT;
  if (objetivo > maxPaso) objetivo = maxPaso;      // evita bloqueos
  const n = Math.max(1, Math.round(objetivo / DT));

  for (let i = 0; i < n; i++) {
    const thAnt = app.th;
    const r = pasoRK4(app.t, app.th, app.om, DT, p, d);
    app.th = r[0]; app.om = r[1];
    app.t += DT;

    if (Math.abs(app.th) > app.ampActual) app.ampActual = Math.abs(app.th);

    // cruce por cero ascendente → período y amplitud del ciclo
    if (thAnt < 0 && app.th >= 0) {
      const frac = -thAnt / (app.th - thAnt);
      const tc = app.t - DT + frac * DT;
      if (app.cruceAnt !== null) {
        app.periodoMedido = tc - app.cruceAnt;
        app.ampMedida = app.ampActual;
        app.series.periodoMed.push([0.5 * (tc + app.cruceAnt), app.periodoMedido]);
        diezmar(app.series.periodoMed);
      }
      app.cruceAnt = tc;
      app.ampActual = 0;
    }

    if (app.t - app.tMuestra >= MUESTREO) { app.tMuestra = app.t; muestrear(); }
  }

  // gotas del chorro (puramente decorativas)
  const h = altura(app.t, p, d);
  if (h > 0 && d.a > 0) {
    const u = p.Cd * Math.sqrt(2 * G * h);
    if (app.gotas.length < 60 && Math.random() < 0.85) {
      app.gotas.push({ r: d.s, ang: app.th, v: u, vida: 1 });
    }
  }
  const dtG = Math.min(dtReal * factor, 0.05);
  for (const g of app.gotas) { g.r += g.v * dtG; g.v += G * dtG; g.vida -= dtG * 0.9; }
  app.gotas = app.gotas.filter(g => g.vida > 0 && g.r < d.s * 2.2);
}

/* =====================================================================
   4. DIBUJO DE LA ESCENA
   ===================================================================== */
function dibujarEscena() {
  const cv = $('escena'), ctx = cv.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) {
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const p = app.p, d = app.d;
  const h = altura(app.t, p, d);

  // escala: que quepa el sistema completo con margen
  const alcance = d.s * 1.16;
  const esc = Math.min((H - 56) / alcance, (W * 0.44) / (d.s * Math.sin(Math.max(Math.abs(app.th), 0.30)) + p.R * 1.6));
  const ox = W / 2, oy = 34;                  // pivote en pantalla

  const aX = r => ox + r * Math.sin(app.th) * esc;   // punto del eje a distancia r
  const aY = r => oy + r * Math.cos(app.th) * esc;

  // ---- techo ----
  ctx.fillStyle = '#222c3a';
  ctx.fillRect(ox - 130, oy - 16, 260, 13);
  ctx.strokeStyle = '#3a4759'; ctx.lineWidth = 1;
  for (let x = -130; x < 130; x += 10) {
    ctx.beginPath(); ctx.moveTo(ox + x, oy - 3); ctx.lineTo(ox + x + 9, oy - 16); ctx.stroke();
  }
  ctx.strokeStyle = '#4a5a70'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(ox - 130, oy - 3); ctx.lineTo(ox + 130, oy - 3); ctx.stroke();

  // ---- vertical de referencia y arco del ángulo ----
  ctx.strokeStyle = 'rgba(255,255,255,.16)'; ctx.lineWidth = 1; ctx.setLineDash([4, 5]);
  ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox, oy + d.s * esc * 0.62); ctx.stroke();
  ctx.setLineDash([]);
  const rArc = Math.min(p.Lc * esc * 0.62, 66);
  if (Math.abs(app.th) > 0.005) {
    ctx.strokeStyle = 'rgba(77,163,255,.75)'; ctx.lineWidth = 1.7;
    ctx.beginPath();
    ctx.arc(ox, oy, rArc, Math.PI / 2, Math.PI / 2 - app.th, app.th > 0);
    ctx.stroke();
    ctx.fillStyle = '#4da3ff'; ctx.font = '12px system-ui, sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const am = -app.th / 2;
    ctx.fillText((app.th * 180 / Math.PI).toFixed(1) + '°',
                 ox + (rArc + 17) * Math.sin(am), oy + (rArc + 17) * Math.cos(am));
  }

  // ---- rastro del fondo de la botella ----
  if ($('chk_rastro').checked) {
    app.rastro.push([aX(d.s), aY(d.s)]);
    if (app.rastro.length > 130) app.rastro.shift();
    ctx.lineWidth = 1.5;
    for (let i = 1; i < app.rastro.length; i++) {
      ctx.strokeStyle = 'rgba(77,163,255,' + (i / app.rastro.length * 0.34).toFixed(3) + ')';
      ctx.beginPath();
      ctx.moveTo(app.rastro[i - 1][0], app.rastro[i - 1][1]);
      ctx.lineTo(app.rastro[i][0], app.rastro[i][1]);
      ctx.stroke();
    }
  } else app.rastro.length = 0;

  // ---- gotas del chorro ----
  for (const g of app.gotas) {
    const gx = ox + g.r * Math.sin(g.ang) * esc, gy = oy + g.r * Math.cos(g.ang) * esc;
    ctx.fillStyle = 'rgba(93,170,225,' + Math.max(0, g.vida * 0.75).toFixed(3) + ')';
    ctx.beginPath(); ctx.ellipse(gx, gy, 1.9, 3.4, -app.th, 0, 2 * Math.PI); ctx.fill();
  }

  // ---- cuerda ----
  ctx.strokeStyle = '#cfd8e4'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(aX(p.Lc), aY(p.Lc)); ctx.stroke();

  // ---- botella (sistema rotado) ----
  ctx.save();
  ctx.translate(ox, oy);
  ctx.rotate(-app.th);           // eje +y del canvas apunta hacia el fondo
  const R = p.R * esc, y1 = p.Lc * esc, y2 = d.s * esc;
  const yAgua = y2 - h * esc;

  // agua
  if (h > 1e-6) {
    const grad = ctx.createLinearGradient(-R, yAgua, R, y2);
    grad.addColorStop(0, 'rgba(61,143,214,.88)');
    grad.addColorStop(1, 'rgba(41,110,175,.88)');
    ctx.fillStyle = grad;
    ctx.fillRect(-R + 1.5, yAgua, 2 * R - 3, y2 - yAgua);
    ctx.strokeStyle = '#7fc2f0'; ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.moveTo(-R + 1.5, yAgua); ctx.lineTo(R - 1.5, yAgua); ctx.stroke();
  }

  // paredes
  ctx.strokeStyle = '#e4eaf2'; ctx.lineWidth = 2.1; ctx.lineJoin = 'round';
  ctx.beginPath();
  ctx.moveTo(-R, y1); ctx.lineTo(R, y1);            // tapa
  ctx.stroke();
  const rOri = Math.min(p.da / 2 * esc, R * 0.8);
  ctx.beginPath();
  ctx.moveTo(-R, y1); ctx.lineTo(-R, y2); ctx.lineTo(-rOri, y2);
  ctx.moveTo(R, y1);  ctx.lineTo(R, y2);  ctx.lineTo(rOri, y2);
  ctx.stroke();

  // marcadores de CM y centro de oscilación
  if ($('chk_centros').checked) {
    const yCM = (momEstatico(h, p, d) / masaTotal(h, p, d)) * esc;
    const yLeq = longEquiv(h, p, d) * esc;
    ctx.strokeStyle = '#f0a93b'; ctx.lineWidth = 1.4; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(-R - 16, yLeq); ctx.lineTo(R + 16, yLeq); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#e5646b';
    ctx.beginPath(); ctx.arc(0, yCM, 4.2, 0, 2 * Math.PI); ctx.fill();
    ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(0, yCM, 4.2, 0, 2 * Math.PI); ctx.stroke();
  }
  ctx.restore();

  // ---- leyenda de la escena ----
  if ($('chk_centros').checked) {
    ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    const bx = 12, by = H - 34;
    ctx.fillStyle = '#e5646b'; ctx.beginPath(); ctx.arc(bx + 5, by, 4, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = '#9aa6b6'; ctx.fillText('centro de masa ℓ(t)', bx + 16, by);
    ctx.strokeStyle = '#f0a93b'; ctx.lineWidth = 1.4; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(bx, by + 16); ctx.lineTo(bx + 11, by + 16); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#9aa6b6';
    ctx.fillText('centro de oscilación L' + '\u2091' + 'q', bx + 16, by + 16);
  }

  // ---- barra de vaciado ----
  const frac = p.h0 > 0 ? h / p.h0 : 0;
  const bw = 96, bh = 7, bx2 = W - bw - 14, by2 = H - 20;
  ctx.fillStyle = 'rgba(255,255,255,.10)'; ctx.fillRect(bx2, by2, bw, bh);
  ctx.fillStyle = '#3d8fd6'; ctx.fillRect(bx2, by2, bw * frac, bh);
  ctx.fillStyle = '#8b97a8'; ctx.font = '10.5px system-ui, sans-serif';
  ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
  ctx.fillText('líquido restante ' + (frac * 100).toFixed(0) + ' %', bx2 + bw, by2 - 3);
}

/* =====================================================================
   5. GRÁFICAS Y LECTURAS
   ===================================================================== */
let gPeriodo, gCM, gMasa, gEnergia;

function crearGraficas() {
  const s = app.series;
  gPeriodo = new Grafica($('g_periodo'), {
    ejeY: 'T (s)',
    series: [
      { nombre: 'T teórico', color: '#4da3ff', datos: s.periodoTeo },
      { nombre: 'T medido',  color: '#f0a93b', datos: s.periodoMed, puntos: true, radio: 2.3 }
    ]
  });
  gCM = new Grafica($('g_cm'), {
    ejeY: 'distancia al pivote (m)',
    series: [
      { nombre: 'ℓ (centro de masa)', color: '#e5646b', datos: s.cm },
      { nombre: 'L_eq (centro de oscilación)', color: '#f0a93b', datos: s.leq },
      { nombre: 's − h (superficie libre)', color: '#4da3ff', datos: s.sup, trazo: [5, 4] }
    ]
  });
  gMasa = new Grafica($('g_masa'), {
    ejeY: 'masa (kg)', yMin: 0,
    series: [
      { nombre: 'm_w (líquido)', color: '#3d8fd6', datos: s.mw },
      { nombre: 'M (total)', color: '#b07fe0', datos: s.M, trazo: [5, 4] }
    ]
  });
  gEnergia = new Grafica($('g_energia'), {
    ejeY: 'energía (J)', yMin: 0,
    series: [
      { nombre: 'E = K + U', color: '#4fc98a', datos: s.E, grosor: 2.1 },
      { nombre: 'K cinética', color: '#f0a93b', datos: s.K, grosor: 1.3 },
      { nombre: 'U potencial', color: '#4da3ff', datos: s.U, grosor: 1.3 }
    ]
  });
}

function actualizarLecturas() {
  const p = app.p, d = app.d, t = app.t;
  const h = altura(t, p, d), hp = alturaPunto(t, p, d);
  const I = inercia(h, p, d), S = momEstatico(h, p, d), M = masaTotal(h, p, d);
  const K = 0.5 * I * app.om * app.om;
  const U = G * S * (1 - Math.cos(app.th));

  // ε = |ω̇|/ω² por diferencias centradas sobre la frecuencia instantánea
  const dtE = 1e-3;
  const w = (tt) => { const hh = altura(tt, p, d); return Math.sqrt(G * momEstatico(hh, p, d) / inercia(hh, p, d)); };
  const w0 = w(t);
  const eps = Math.abs((w(t + dtE) - w(Math.max(0, t - dtE))) / (dtE * (t > dtE ? 2 : 1))) / (w0 * w0);

  $('reloj_t').textContent = 't = ' + t.toFixed(2) + ' s';
  $('l_theta').textContent = (app.th * 180 / Math.PI).toFixed(2) + '°';
  $('l_amp').textContent   = (app.ampMedida * 180 / Math.PI).toFixed(2) + '°';
  $('l_T').textContent     = isFinite(app.periodoMedido) ? app.periodoMedido.toFixed(4) + ' s'
                                                         : periodoLineal(h, p, d).toFixed(4) + ' s';
  $('l_h').textContent     = (h * 100).toFixed(2) + ' cm';
  $('l_mw').textContent    = (p.rho * d.A * h).toFixed(4) + ' kg';
  $('l_M').textContent     = M.toFixed(4) + ' kg';
  $('l_cm').textContent    = (S / M).toFixed(4) + ' m';
  $('l_Leq').textContent   = (I / S).toFixed(4) + ' m';
  $('l_I').textContent     = I.toFixed(5) + ' kg·m²';
  $('l_gam').textContent   = (gammaFrenado(h, hp, p, d) / I).toExponential(2) + ' s⁻¹';
  $('l_E').textContent     = (K + U).toFixed(5) + ' J';
  $('l_eps').textContent   = eps.toExponential(2);
}

function dibujarGraficas() {
  gPeriodo.dibujar(); gCM.dibujar(); gMasa.dibujar(); gEnergia.dibujar();
}

function dibujarTodo() { dibujarEscena(); actualizarLecturas(); dibujarGraficas(); }

/* =====================================================================
   6. BUCLE PRINCIPAL
   ===================================================================== */
let ultimo = 0, contador = 0;

function bucle(ahora) {
  requestAnimationFrame(bucle);
  if (!ultimo) ultimo = ahora;
  const dtReal = Math.min((ahora - ultimo) / 1000, 0.05);
  ultimo = ahora;

  if (app.corriendo) avanzar(dtReal);
  dibujarEscena();
  if (++contador % 4 === 0) { actualizarLecturas(); dibujarGraficas(); }
}

/* =====================================================================
   7. VERIFICACIÓN DEL MODELO (se ejecuta al cargar)
   ===================================================================== */
function integrarSimple(p, d, tFin, dt) {
  let t = 0, th = p.theta0, om = 0;
  let cruceAnt = null, amp = Math.abs(th), periodos = [], amps = [];
  const n = Math.round(tFin / dt);
  for (let i = 0; i < n; i++) {
    const thAnt = th;
    const r = pasoRK4(t, th, om, dt, p, d);
    th = r[0]; om = r[1]; t += dt;
    if (Math.abs(th) > amp) amp = Math.abs(th);
    if (thAnt < 0 && th >= 0) {
      const tc = t - dt + (-thAnt / (th - thAnt)) * dt;
      if (cruceAnt !== null) { periodos.push(tc - cruceAnt); amps.push(amp); }
      cruceAnt = tc; amp = 0;
    }
  }
  return { periodos, amps, th, om, t };
}

function verificar() {
  const res = [];

  // --- Prueba 1: sin vaciado, período de amplitud finita ---
  {
    const p = { H: .25, R: .04, h0: .25, rho: 1000, Lc: .5, mb: .04, da: 0, Cd: .6,
                theta0: 12 * Math.PI / 180 };
    const d = derivar(p);
    const r = integrarSimple(p, d, 30, 5e-4);
    const medido = r.periodos.reduce((a, b) => a + b, 0) / r.periodos.length;
    const teorico = periodoAmplitud(p.h0, p, d, p.theta0);
    const err = Math.abs(medido / teorico - 1);
    res.push({ ok: err < 2e-4,
      texto: 'Sin vaciado, el período reproduce T₀(1 + Θ²/16 + 11Θ⁴/3072)',
      detalle: 'error ' + (err * 100).toExponential(1) + ' %' });
  }

  // --- Prueba 2: identidades dS/dh y dI/dh ---
  {
    const p = { H: .25, R: .04, h0: .25, rho: 1000, Lc: .5, mb: .04, da: .005, Cd: .6, theta0: 0 };
    const d = derivar(p);
    const h = 0.12, e = 1e-6;
    const dS = (momEstatico(h + e, p, d) - momEstatico(h - e, p, d)) / (2 * e);
    const dI = (inercia(h + e, p, d) - inercia(h - e, p, d)) / (2 * e);
    const dSt = p.rho * d.A * (d.s - h);
    const dIt = p.rho * d.A * (Math.pow(d.s - h, 2) + p.R * p.R / 4);
    const err = Math.max(Math.abs(dS / dSt - 1), Math.abs(dI / dIt - 1));
    res.push({ ok: err < 1e-6,
      texto: 'Identidades dS/dh = ρA(s−h) y dI/dh = ρA[(s−h)²+R²/4]',
      detalle: 'error ' + err.toExponential(1) });
  }

  // --- Prueba 3: la predicción WKB reproduce la amplitud final ---
  {
    const p = { H: .25, R: .04, h0: .25, rho: 1000, Lc: .5, mb: .04, da: .008, Cd: .6,
                theta0: 12 * Math.PI / 180 };
    const d = derivar(p);
    const r = integrarSimple(p, d, d.tv * 0.995, 1e-3);
    const medido = r.amps[r.amps.length - 1] / p.theta0;
    const wkb = factorAmplitudWKB(p, d);
    const err = Math.abs(medido / wkb - 1);
    res.push({ ok: err < 0.01,
      texto: 'La ley adiabática Θ ∝ (S·I)^(−1/4)·e^(−½∫ṁ s²/I) reproduce la amplitud final',
      detalle: 'numérico ' + medido.toFixed(4) + ' · WKB ' + wkb.toFixed(4) });
  }

  // --- Prueba 4: universalidad respecto del área del orificio ---
  {
    const base = { H: .25, R: .04, h0: .25, rho: 1000, Lc: .5, mb: .04, Cd: .6,
                   theta0: 12 * Math.PI / 180 };
    const f = [];
    for (const da of [0.006, 0.012]) {
      const p = Object.assign({}, base, { da });
      const d = derivar(p);
      const r = integrarSimple(p, d, d.tv * 0.995, 1e-3);
      f.push(r.amps[r.amps.length - 1] / p.theta0);
    }
    const err = Math.abs(f[0] / f[1] - 1);
    res.push({ ok: err < 0.02,
      texto: 'El factor total de amplitud no depende del área del orificio',
      detalle: 'd=6 mm → ' + f[0].toFixed(4) + ' · d=12 mm → ' + f[1].toFixed(4) });
  }

  // --- Prueba 5: masa puntual ⇒ período y amplitud constantes ---
  {
    // R y H muy pequeños frente a Lc: toda la masa concentrada en el fondo
    const p = { H: 2e-3, R: 1e-3, h0: 2e-3, rho: 1000, Lc: .8, mb: 1e-4, da: 2e-4, Cd: .6,
                theta0: 12 * Math.PI / 180 };
    const d = derivar(p);
    const r = integrarSimple(p, d, d.tv * 0.995, 5e-4);
    const dT = Math.abs(r.periodos[r.periodos.length - 1] / r.periodos[0] - 1);
    const dA = Math.abs(r.amps[r.amps.length - 1] / r.amps[0] - 1);
    res.push({ ok: dT < 5e-3 && dA < 5e-3,
      texto: 'Límite de masa puntual: período y amplitud permanecen constantes',
      detalle: 'ΔT ' + (dT * 100).toFixed(2) + ' % · ΔΘ ' + (dA * 100).toFixed(2) + ' %' });
  }

  const ul = $('lista_validacion');
  ul.innerHTML = '';
  for (const r of res) {
    const li = document.createElement('li');
    li.className = r.ok ? 'ok' : 'fallo';
    li.innerHTML = '<span class="marca">' + (r.ok ? '✓' : '✗') + '</span>' +
                   '<span>' + r.texto + '</span>' +
                   '<span class="detalle">' + r.detalle + '</span>';
    ul.appendChild(li);
  }
}

/* =====================================================================
   8. ARRANQUE
   ===================================================================== */
function inicializar() {
  app.p = leerParametros();
  app.d = derivar(app.p);
  crearGraficas();
  reiniciar();

  const ctrls = ['H', 'D', 'llenado', 'rho', 'Lc', 'mb', 'th0', 'da', 'Cd'];
  for (const id of ctrls) {
    $(id).addEventListener('input', () => {
      const estaba = app.corriendo;
      app.corriendo = false;
      reiniciar();
      if (estaba) { app.corriendo = true; $('btn_play').textContent = '❚❚ Pausar'; }
    });
  }
  $('vel').addEventListener('input', actualizarEtiquetas);

  $('btn_play').addEventListener('click', () => {
    app.corriendo = !app.corriendo;
    $('btn_play').textContent = app.corriendo ? '❚❚ Pausar' : '▶ Iniciar';
    $('btn_play').classList.toggle('primario', !app.corriendo);
  });
  $('btn_reset').addEventListener('click', () => {
    app.corriendo = false;
    $('btn_play').textContent = '▶ Iniciar';
    $('btn_play').classList.add('primario');
    reiniciar();
  });
  $('btn_defecto').addEventListener('click', () => {
    for (const k in POR_DEFECTO) if ($(k)) $(k).value = POR_DEFECTO[k];
    app.corriendo = false;
    $('btn_play').textContent = '▶ Iniciar';
    $('btn_play').classList.add('primario');
    reiniciar();
  });
  window.addEventListener('resize', dibujarTodo);

  setTimeout(verificar, 60);     // no bloquear el primer pintado
  requestAnimationFrame(bucle);
}

document.addEventListener('DOMContentLoaded', inicializar);
