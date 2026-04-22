const els = {
  canvas: document.getElementById('simCanvas'),
  playPauseBtn: document.getElementById('playPauseBtn'),
  resetBtn: document.getElementById('resetBtn'),
  successPresetBtn: document.getElementById('successPresetBtn'),
  failPresetBtn: document.getElementById('failPresetBtn'),
  B: document.getElementById('B'),
  sigma: document.getElementById('sigma'),
  L: document.getElementById('L'),
  N: document.getElementById('N'),
  dt: document.getElementById('dt'),
  substeps: document.getElementById('substeps'),
  viewMode: document.getElementById('viewMode'),
  BValue: document.getElementById('BValue'),
  sigmaValue: document.getElementById('sigmaValue'),
  LValue: document.getElementById('LValue'),
  NValue: document.getElementById('NValue'),
  dtValue: document.getElementById('dtValue'),
  substepsValue: document.getElementById('substepsValue'),
  timeReadout: document.getElementById('timeReadout'),
  energyReadout: document.getElementById('energyReadout'),
  dxReadout: document.getElementById('dxReadout'),
  pairReadout: document.getElementById('pairReadout'),
  messageReadout: document.getElementById('messageReadout')
};

const TWO_PI = 2 * Math.PI;
const PHI_PADDING = 1.0;

const ctx = els.canvas.getContext('2d');
let running = false;
let sim = null;
let rafId = null;

function params() {
  return {
    B: Number(els.B.value),
    sigma: Number(els.sigma.value),
    L: Number(els.L.value),
    N: Number(els.N.value),
    dt: Number(els.dt.value),
    substeps: Number(els.substeps.value),
    viewMode: els.viewMode.value
  };
}

function syncLabels() {
  const p = params();
  els.BValue.textContent = p.B.toFixed(1);
  els.sigmaValue.textContent = p.sigma.toFixed(1);
  els.LValue.textContent = p.L.toFixed(0);
  els.NValue.textContent = p.N.toFixed(0);
  els.dtValue.textContent = p.dt.toFixed(3);
  els.substepsValue.textContent = p.substeps.toFixed(0);
}

function initializePhiSectorWindow(phi) {
  let minPhi = Infinity;
  let maxPhi = -Infinity;
  for (let i = 0; i < phi.length; i++) {
    minPhi = Math.min(minPhi, phi[i]);
    maxPhi = Math.max(maxPhi, phi[i]);
  }
  const smin = Math.floor(minPhi / TWO_PI);
  const smax = Math.floor(maxPhi / TWO_PI);
  return { smin, smax };
}

function updatePhiSectorWindow() {
  let minPhi = Infinity;
  let maxPhi = -Infinity;
  for (let i = 0; i < sim.phi.length; i++) {
    minPhi = Math.min(minPhi, sim.phi[i]);
    maxPhi = Math.max(maxPhi, sim.phi[i]);
  }
  const currentYMin = TWO_PI * sim.phiSectorWindow.smin - PHI_PADDING;
  const currentYMax = TWO_PI * (sim.phiSectorWindow.smax + 1) + PHI_PADDING;

  if (minPhi < currentYMin || maxPhi > currentYMax) {
    const newSMin = Math.floor(minPhi / TWO_PI);
    const newSMax = Math.floor(maxPhi / TWO_PI);
    sim.phiSectorWindow.smin = Math.min(sim.phiSectorWindow.smin, newSMin);
    sim.phiSectorWindow.smax = Math.max(sim.phiSectorWindow.smax, newSMax);
  }
}

function getPhiAxisBounds() {
  updatePhiSectorWindow();
  return {
    ymin: TWO_PI * sim.phiSectorWindow.smin - PHI_PADDING,
    ymax: TWO_PI * (sim.phiSectorWindow.smax + 1) + PHI_PADDING
  };
}

function getPhiSectorLines() {
  const lines = [];
  for (let n = sim.phiSectorWindow.smin; n <= sim.phiSectorWindow.smax + 1; n++) {
    lines.push(n * TWO_PI);
  }
  return lines;
}

function initialize() {
  const p = params();
  const dx = p.L / p.N;
  const x = new Float64Array(p.N);
  const phi = new Float64Array(p.N);
  const pi = new Float64Array(p.N);
  for (let i = 0; i < p.N; i++) {
    x[i] = -p.L / 2 + i * dx;
    pi[i] = p.B * Math.exp(-((x[i] / p.sigma) ** 2));
  }
  sim = {
    ...p,
    dx,
    x,
    phi,
    pi,
    t: 0,
    phiSectorWindow: initializePhiSectorWindow(phi)
  };
  updateReadout('initialized');
  draw();
}

function laplacian(phi, dx) {
  const n = phi.length;
  const out = new Float64Array(n);
  const inv = 1 / (dx * dx);
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const im = i === 0 ? n - 1 : i - 1;
    out[i] = (phi[ip] - 2 * phi[i] + phi[im]) * inv;
  }
  return out;
}

function step() {
  const { phi, pi, dx, dt } = sim;
  const n = phi.length;
  const a = laplacian(phi, dx);
  const piHalf = new Float64Array(n);
  const phiNew = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    piHalf[i] = pi[i] + 0.5 * dt * (a[i] - Math.sin(phi[i]));
    phiNew[i] = phi[i] + dt * piHalf[i];
  }
  const aNew = laplacian(phiNew, dx);
  const piNew = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    piNew[i] = piHalf[i] + 0.5 * dt * (aNew[i] - Math.sin(phiNew[i]));
  }
  sim.phi = phiNew;
  sim.pi = piNew;
  sim.t += dt;
}

function energyDensity() {
  const { phi, pi, dx } = sim;
  const n = phi.length;
  const ed = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const grad = (phi[ip] - phi[i]) / dx;
    ed[i] = 0.5 * pi[i] * pi[i] + 0.5 * grad * grad + (1 - Math.cos(phi[i]));
  }
  return ed;
}

function totalEnergy(ed) {
  let sum = 0;
  for (let i = 0; i < ed.length; i++) sum += ed[i];
  return sum * sim.dx;
}

function estimatePair() {
  const { phi, dx, L } = sim;
  const n = phi.length;
  const grad = new Float64Array(n);
  let max1 = -1, idx1 = -1, max2 = -1, idx2 = -1;
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const im = i === 0 ? n - 1 : i - 1;
    grad[i] = Math.abs((phi[ip] - phi[im]) / (2 * dx));
  }
  for (let i = 0; i < n; i++) {
    const v = grad[i];
    if (v > max1) {
      max2 = max1; idx2 = idx1;
      max1 = v; idx1 = i;
    } else if (v > max2 && Math.abs(i - idx1) > Math.floor(4 / dx)) {
      max2 = v; idx2 = i;
    }
  }
  if (idx1 < 0 || idx2 < 0) return { label: 'none' };
  const i = Math.min(idx1, idx2), j = Math.max(idx1, idx2);
  const mid = [], out = [];
  for (let k = i; k <= j; k++) mid.push(phi[k]);
  for (let k = 0; k < i; k++) out.push(phi[k]);
  for (let k = j + 1; k < n; k++) out.push(phi[k]);
  if (mid.length < 5 || out.length < 5) return { label: 'none' };
  const median = (arr) => {
    const a = [...arr].sort((u, v) => u - v);
    return a[Math.floor(a.length / 2)];
  };
  const jump = Math.abs(median(mid) - median(out));
  let sep = (j - i) * dx;
  sep = Math.min(sep, L - sep);
  const looksLikePair = jump > 4.5 && sep > 5 && sep < 0.48 * L;
  return { label: looksLikePair ? `pair-like (Δphi=${jump.toFixed(2)}, sep=${sep.toFixed(2)})` : `no clear pair (Δphi=${jump.toFixed(2)})` };
}

function updateReadout(message) {
  const ed = energyDensity();
  els.timeReadout.textContent = sim.t.toFixed(2);
  els.energyReadout.textContent = totalEnergy(ed).toFixed(4);
  els.dxReadout.textContent = sim.dx.toFixed(4);
  els.pairReadout.textContent = estimatePair().label;
  els.messageReadout.textContent = message;
}

function drawAxes(x0, y0, w, h, ymin, ymax, title, dashed = []) {
  ctx.strokeStyle = '#cbd5e1';
  ctx.lineWidth = 1;
  ctx.strokeRect(x0, y0, w, h);
  ctx.fillStyle = '#0f172a';
  ctx.font = '14px sans-serif';
  ctx.fillText(title, x0, y0 - 10);
  ctx.fillStyle = '#475569';
  ctx.font = '12px sans-serif';
  ctx.fillText(ymax.toFixed(1), 8, y0 + 12);
  ctx.fillText(ymin.toFixed(1), 8, y0 + h);
  ctx.fillText((-sim.L / 2).toFixed(0), x0 - 4, y0 + h + 18);
  ctx.fillText((sim.L / 2).toFixed(0), x0 + w - 24, y0 + h + 18);
  dashed.forEach((yVal) => {
    if (yVal < ymin || yVal > ymax) return;
    const yy = y0 + h - ((yVal - ymin) / (ymax - ymin)) * h;
    ctx.setLineDash([6, 6]);
    ctx.strokeStyle = '#94a3b8';
    ctx.beginPath();
    ctx.moveTo(x0, yy);
    ctx.lineTo(x0 + w, yy);
    ctx.stroke();
    ctx.setLineDash([]);
  });
}

function drawCurve(arr, x0, y0, w, h, ymin, ymax, color) {
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < arr.length; i++) {
    const xx = x0 + (i / (arr.length - 1)) * w;
    const yy = y0 + h - ((arr[i] - ymin) / (ymax - ymin)) * h;
    if (i === 0) ctx.moveTo(xx, yy);
    else ctx.lineTo(xx, yy);
  }
  ctx.stroke();
}

function draw() {
  ctx.clearRect(0, 0, els.canvas.width, els.canvas.height);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, els.canvas.width, els.canvas.height);
  const pad = 44;
  const gap = 28;
  const mode = sim.viewMode;
  const panels = mode === 'both' ? 2 : 1;
  const panelH = (els.canvas.height - 2 * pad - (panels - 1) * gap) / panels;
  const plotW = els.canvas.width - 2 * pad;
  const ed = energyDensity();
  if (mode === 'phi' || mode === 'both') {
    const phiBounds = getPhiAxisBounds();
    drawAxes(pad, pad, plotW, panelH, phiBounds.ymin, phiBounds.ymax, 'Field profile phi(x,t)', getPhiSectorLines());
    drawCurve(sim.phi, pad, pad, plotW, panelH, phiBounds.ymin, phiBounds.ymax, '#2563eb');
  }
  if (mode === 'energy' || mode === 'both') {
    const y0 = mode === 'both' ? pad + panelH + gap : pad;
    let edMax = 0.1;
    for (let i = 0; i < ed.length; i++) edMax = Math.max(edMax, ed[i]);
    drawAxes(pad, y0, plotW, panelH, 0, edMax * 1.08, 'Energy density', []);
    drawCurve(ed, pad, y0, plotW, panelH, 0, edMax * 1.08, '#dc2626');
  }
}

function tick() {
  if (running) {
    for (let k = 0; k < sim.substeps; k++) step();
    updateReadout('running');
    draw();
  }
  rafId = requestAnimationFrame(tick);
}

function applyPreset(kind) {
  if (kind === 'success') {
    els.B.value = 4.5;
    els.sigma.value = 1.5;
    els.L.value = 120;
    els.N.value = 400;
    els.dt.value = 0.02;
    els.substeps.value = 4;
  } else {
    els.B.value = 2.0;
    els.sigma.value = 2.0;
    els.L.value = 120;
    els.N.value = 400;
    els.dt.value = 0.02;
    els.substeps.value = 4;
  }
  syncLabels();
  initialize();
}

els.playPauseBtn.addEventListener('click', () => {
  running = !running;
  els.playPauseBtn.textContent = running ? 'Pause' : 'Play';
  updateReadout(running ? 'running' : 'paused');
});

els.resetBtn.addEventListener('click', () => {
  running = false;
  els.playPauseBtn.textContent = 'Play';
  initialize();
});

els.successPresetBtn.addEventListener('click', () => applyPreset('success'));
els.failPresetBtn.addEventListener('click', () => applyPreset('fail'));

[els.B, els.sigma, els.L, els.N, els.dt, els.substeps, els.viewMode].forEach(el => {
  el.addEventListener('input', () => {
    syncLabels();
    running = false;
    els.playPauseBtn.textContent = 'Play';
    initialize();
  });
});

syncLabels();
initialize();
tick();
