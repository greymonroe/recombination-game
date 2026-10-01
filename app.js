'use strict';

// ---------- genome ----------
// Positions are in centimorgans; one crossover on a bivalent is worth 50 cM.
// Each marker is a SNP: [site name, position, blue-parent base, red-parent base].
const CHRS = [
  { name: 'Chromosome 1', len: 120, cen: 62, markers: [['A', 6, 'G', 'T'], ['B', 16, 'C', 'T'], ['C', 44, 'A', 'G'], ['D', 112, 'T', 'C']] },
  { name: 'Chromosome 2', len: 80, cen: 48, markers: [['E', 8, 'G', 'A'], ['F', 30, 'A', 'C'], ['G', 72, 'C', 'G']] },
];
const MK = [];
CHRS.forEach((ch, c) => ch.markers.forEach(([n, pos, b0, b1]) => MK.push({ n, c, pos, base: [b0, b1] })));
const TOTAL = CHRS.reduce((s, ch) => s + ch.len, 0);
const COL = ['#2a78d6', '#e34948'], HET = '#b9b6ad';
const MAX_ROWS = 800;

const $ = (s) => document.querySelector(s);
const rnd = (n) => Math.floor(Math.random() * n);
const pct = (x, len) => (100 * x / len).toFixed(2) + '%';

// ---------- haplotypes ----------
// A chromatid is a list of segments [{a: allele 0|1, e: end position}].
const alleleAt = (segs, x) => (segs.find((s) => x < s.e) || segs[segs.length - 1]).a;
const pure = (a) => CHRS.map((ch) => [{ a, e: ch.len }]);

function poisson(mean) {
  let l = Math.exp(-mean), k = 0, p = 1;
  do { k++; p *= Math.random(); } while (p > l);
  return k - 1;
}

// Fast path for bulk simulation: one chromatid drawn directly.
function simGamete() {
  return CHRS.map((ch) => {
    const n = poisson(ch.len / 100);
    const bps = Array.from({ length: n }, () => Math.random() * ch.len).sort((p, q) => p - q);
    let a = rnd(2);
    const segs = [];
    for (const b of bps) { segs.push({ a, e: b }); a = 1 - a; }
    segs.push({ a, e: ch.len });
    return segs;
  });
}

// Full meiosis for the animation: crossovers on the four-chromatid bivalent,
// then one chromatid of each chromosome goes into the gamete. Same model as
// simGamete, with the intermediate steps kept for drawing.
function simMeiosis() {
  const chr = CHRS.map((ch) => {
    const xos = Array.from({ length: poisson(ch.len / 50) }, () => ({ x: 1 + Math.random() * (ch.len - 2), i: rnd(2), j: 2 + rnd(2) }));
    return { xos, pick: rnd(4), ...resolveBivalent(ch, xos) };
  });
  return { chr, g: chr.map((b) => b.rows[b.pick]) };
}

const makeKid = (g1, g2) => ({ g: [g1, g2], geno: MK.map((m) => alleleAt(g1[m.c], m.pos) + alleleAt(g2[m.c], m.pos)) });

// Resolve crossovers on a bivalent into four chromatids. Strands 0,1 start as
// allele 0 (sisters) and 2,3 as allele 1. Each chromatid keeps its centromere
// and switches partner strand at every crossover it takes part in.
function resolveBivalent(ch, xos) {
  const al = (s) => (s < 2 ? 0 : 1);
  const sorted = xos.map((o, k) => ({ ...o, k })).sort((p, q) => p.x - q.x);
  const right = sorted.filter((o) => o.x >= ch.cen);
  const left = sorted.filter((o) => o.x < ch.cen).reverse();
  const links = xos.map(() => []);
  const rows = [];
  for (let c = 0; c < 4; c++) {
    const pieces = [];
    let cur = c, to = ch.cen;
    for (const o of left) {
      if (cur !== o.i && cur !== o.j) continue;
      pieces.unshift({ a: al(cur), e: to });
      to = o.x; cur = cur === o.i ? o.j : o.i; links[o.k].push(c);
    }
    pieces.unshift({ a: al(cur), e: to });
    cur = c;
    for (const o of right) {
      if (cur !== o.i && cur !== o.j) continue;
      pieces.push({ a: al(cur), e: o.x });
      cur = cur === o.i ? o.j : o.i; links[o.k].push(c);
    }
    pieces.push({ a: al(cur), e: ch.len });
    const segs = [];
    for (const p of pieces) {
      const last = segs[segs.length - 1];
      if (last && last.a === p.a) last.e = p.e; else segs.push({ ...p });
    }
    rows.push(segs);
  }
  return { rows, links };
}

// ---------- recombination frequency from F2 genotypes ----------
function pairCounts(pop, i, j) {
  const t = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const k of pop) t[k.geno[i]][k.geno[j]]++;
  return t;
}
// Maximum likelihood by EM: double heterozygotes are either two parental or
// two recombinant gametes, so their expected recombinant count depends on r.
function estRF(t) {
  const n = t.flat().reduce((s, v) => s + v, 0);
  if (!n) return null;
  const n1 = t[0][1] + t[1][0] + t[1][2] + t[2][1], n2 = t[0][2] + t[2][0], dh = t[1][1];
  let r = 0.25;
  for (let it = 0; it < 500; it++) {
    const r2 = (n1 + 2 * n2 + dh * 2 * r * r / ((1 - r) * (1 - r) + r * r)) / (2 * n);
    if (Math.abs(r2 - r) < 1e-9) { r = r2; break; }
    r = r2;
  }
  return Math.min(r, 0.5);
}
const haldane = (d) => 0.5 * (1 - Math.exp(-2 * d / 100));
const expectedRF = (i, j) => (MK[i].c === MK[j].c ? haldane(Math.abs(MK[i].pos - MK[j].pos)) : 0.5);
function rfMatrix(pop) {
  const m = MK.map(() => MK.map(() => null));
  for (let i = 0; i < MK.length; i++) for (let j = i + 1; j < MK.length; j++) m[i][j] = m[j][i] = estRF(pairCounts(pop, i, j));
  return m;
}

// ---------- chromosome drawing (HTML) ----------
function gradient(segs, len) {
  let from = 0;
  const stops = segs.map((s) => { const g = `${COL[s.a]} ${pct(from, len)} ${pct(s.e, len)}`; from = s.e; return g; });
  return `linear-gradient(90deg, ${stops.join(', ')})`;
}
const mkHead = (c) => `<div class="mk-head">${CHRS[c].markers.map(([n, p]) => `<span style="left:${pct(p, CHRS[c].len)}">${n}</span>`).join('')}</div>`;

const chrLab = (c, short) => `<div class="lab">${short ? 'Chr' : 'Chromosome'} ${c + 1}</div>`;
// calls: base calls stacked under each marker, one per chromosome drawn
// (top base = top chromosome), colored by which parent the base came from.
function karyo(haps, { head = false, short = false, calls = false } = {}) {
  return `<div class="karyo">${CHRS.map((ch, c) => `<div class="chrom" style="flex:${ch.len}">
    ${chrLab(c, short)}${head ? mkHead(c) : ''}
    <div class="bars">${haps.map((h) => `<div class="bar" style="background:${gradient(h[c], ch.len)}"><i class="cen" style="left:${pct(ch.cen, ch.len)}"></i></div>`).join('')}
    ${ch.markers.map(([, p]) => `<i class="tick" style="left:${pct(p, ch.len)}"></i>`).join('')}</div>
    ${calls ? `<div class="calls">${MK.filter((m) => m.c === c).map((m) => `<span class="callcol" style="left:${pct(m.pos, ch.len)}">${haps.map((h) => { const a = alleleAt(h[c], m.pos); return `<b class="a${a}">${m.base[a]}</b>`; }).join('')}</span>`).join('')}</div>` : ''}</div>`).join('')}</div>`;
}
const genoName = (m, g) => (g === 0 ? `${m.base[0]}/${m.base[0]}` : g === 1 ? `${m.base[0]}/${m.base[1]}` : `${m.base[1]}/${m.base[1]}`);
function chips(haps) {
  return `<div class="chips">${MK.map((m, i) => {
    const sep = i && MK[i - 1].c !== m.c ? '<span class="sep"></span>' : '';
    const calls = haps.map((h) => { const a = alleleAt(h[m.c], m.pos); return `<b class="a${a}">${m.base[a]}</b>`; }).join('<i>/</i>');
    return `${sep}<span class="chip"><span class="site">${m.n}</span><span class="call">${calls}</span></span>`;
  }).join('')}</div>`;
}

$('#cross').innerHTML = `<div class="cross">
  <div><div class="who">Blue parent</div>${karyo([pure(0), pure(0)], { short: true })}${chips([pure(0)])}</div><div class="x">×</div>
  <div><div class="who">Red parent</div>${karyo([pure(1), pure(1)], { short: true })}${chips([pure(1)])}</div>
  <div class="arrow">↓</div>
  <div class="f1"><div class="who">F1 (selfed to make the offspring)</div>${karyo([pure(0), pure(1)], { head: true, calls: true })}</div></div>`;

// ---------- heat map + detail ----------
function mix(t) {
  const a = [241, 240, 250], b = [74, 58, 167];
  return `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * t)).join(',')})`;
}
function renderHeat(el, detailEl, pop, st) {
  const rf = rfMatrix(pop), S = 44, L = 24, T = 38, GAPC = 8;
  const off = (i) => i * S + (MK[i].c ? GAPC : 0);
  const size = L + off(MK.length - 1) + S;
  let s = `<svg class="heat" viewBox="0 0 ${size} ${T + size - L}" role="img" aria-label="Recombination frequency between each pair of markers">`;
  CHRS.forEach((ch, c) => {
    const idx = MK.map((m, i) => (m.c === c ? i : -1)).filter((i) => i >= 0);
    const x1 = L + off(idx[0]) + 3, x2 = L + off(idx[idx.length - 1]) + S - 3;
    s += `<text x="${(x1 + x2) / 2}" y="11" text-anchor="middle" font-size="10" fill="#8a8880">CHR ${c + 1}</text><line x1="${x1}" x2="${x2}" y1="16" y2="16" stroke="#c9c6bc"/>`;
  });
  MK.forEach((m, i) => {
    s += `<text x="${L + off(i) + S / 2}" y="${T - 7}" text-anchor="middle" font-size="13" font-weight="700" fill="#16161a">${m.n}</text>`;
    s += `<text x="${L / 2}" y="${T + off(i) + S / 2 + 5}" text-anchor="middle" font-size="13" font-weight="700" fill="#16161a">${m.n}</text>`;
  });
  for (let i = 0; i < MK.length; i++) for (let j = 0; j < MK.length; j++) {
    const x = L + off(j) + 1, y = T + off(i) + 1, r = rf[i][j];
    if (i === j) { s += `<rect x="${x}" y="${y}" width="${S - 2}" height="${S - 2}" rx="5" fill="#eeede8"/>`; continue; }
    const sel = st.sel && ((st.sel[0] === i && st.sel[1] === j) || (st.sel[0] === j && st.sel[1] === i));
    s += `<g class="cell" data-i="${i}" data-j="${j}"><rect x="${x}" y="${y}" width="${S - 2}" height="${S - 2}" rx="5" fill="${r == null ? '#f6f5f1' : mix(r / 0.5)}" ${sel ? 'stroke="#16161a" stroke-width="2.5"' : ''}/>
      <text x="${x + S / 2 - 1}" y="${y + S / 2 + 3}" text-anchor="middle" font-size="12" font-weight="600" fill="${r == null ? '#8a8880' : r > 0.27 ? '#fff' : '#16161a'}">${r == null ? '–' : Math.round(r * 100) + '%'}</text></g>`;
  }
  el.innerHTML = s + '</svg><div class="scale">0% <i></i> 50%</div>';
  el.onclick = (e) => {
    const g = e.target.closest('.cell');
    if (!g) return;
    st.sel = [+g.dataset.i, +g.dataset.j].sort((p, q) => p - q);
    renderHeat(el, detailEl, pop, st);
  };
  if (!st.sel) { detailEl.innerHTML = ''; return; }
  const [i, j] = st.sel, a = MK[i], b = MK[j], t = pairCounts(pop, i, j), r = rf[i][j];
  const exp = a.c === b.c
    ? `<p class="note">These markers are ${Math.abs(a.pos - b.pos)} cM apart on chromosome ${a.c + 1}. With unlimited offspring the estimate would settle near ${(100 * expectedRF(i, j)).toFixed(1)}%.</p>`
    : '<p class="note">These markers are on different chromosomes, so they assort independently. With unlimited offspring the estimate would settle at 50%.</p>';
  detailEl.innerHTML = `<div class="detail"><div>Markers <b>${a.n}</b> and <b>${b.n}</b> · ${pop.length} offspring</div>
    <div class="big">${r == null ? '–' : (100 * r).toFixed(1) + '%'}</div>${exp}
    <table><tr><th></th>${[0, 1, 2].map((g) => `<th>${b.n}: ${genoName(b, g)}</th>`).join('')}</tr>
    ${[0, 1, 2].map((g) => `<tr><th>${a.n}: ${genoName(a, g)}</th>${t[g].map((v) => `<td>${v}</td>`).join('')}</tr>`).join('')}</table>
    <p class="note">Number of offspring with each combination of genotypes at the two markers.</p></div>`;
}

// ---------- population ----------
const sim = { pop: [], view: 'hap', heat: {} };
const GAP = 14;

function drawPop() {
  const scroll = $('#sim-scroll'), cv = $('#sim-canvas'), n = sim.pop.length;
  $('#sim-empty').hidden = n > 0; cv.hidden = n === 0;
  $('#sim-n').textContent = n ? `n = ${n}` : '';
  const shown = sim.pop.slice(-MAX_ROWS);
  $('#sim-shown').textContent = n > MAX_ROWS ? `Showing the most recent ${MAX_ROWS} offspring. All ${n} are used for the estimates below.` : n ? 'Each row is one offspring.' : '';
  if (!n) return;
  const W = scroll.clientWidth, rh = shown.length <= 10 ? 24 : shown.length <= 40 ? 14 : shown.length <= 150 ? 8 : 5;
  const dpr = Math.min(window.devicePixelRatio || 1, 2), gapRow = rh >= 14 ? 5 : rh >= 8 ? 2 : 1;
  cv.width = Math.round(W * dpr); cv.height = shown.length * rh * dpr; cv.style.height = shown.length * rh + 'px';
  const ctx = cv.getContext('2d');
  ctx.scale(dpr, dpr);
  const cw = CHRS.map((ch) => (W - GAP) * ch.len / TOTAL), x0 = [0, cw[0] + GAP];
  const inner = rh >= 14 ? 1 : 0, strip = (rh - gapRow - inner) / 2;
  shown.forEach((k, row) => {
    const y = row * rh;
    CHRS.forEach((ch, c) => {
      const sc = cw[c] / ch.len;
      if (sim.view === 'hap') {
        k.g.forEach((g, h) => {
          let from = 0;
          for (const s of g[c]) { ctx.fillStyle = COL[s.a]; ctx.fillRect(x0[c] + from * sc, y + h * (strip + inner), (s.e - from) * sc, strip); from = s.e; }
        });
      } else {
        const cuts = [...new Set([...k.g[0][c], ...k.g[1][c]].map((s) => s.e))].sort((p, q) => p - q);
        let from = 0;
        for (const e of cuts) {
          const z = alleleAt(k.g[0][c], from) + alleleAt(k.g[1][c], from);
          ctx.fillStyle = z === 1 ? HET : COL[z / 2];
          ctx.fillRect(x0[c] + from * sc, y, (e - from) * sc, rh - gapRow); from = e;
        }
      }
    });
  });
  ctx.fillStyle = 'rgba(255,255,255,.85)';
  MK.forEach((m) => ctx.fillRect(x0[m.c] + m.pos * cw[m.c] / CHRS[m.c].len - 0.5, 0, 1, shown.length * rh));
  scroll.scrollTop = scroll.scrollHeight;
}

function renderScatter() {
  const rf = rfMatrix(sim.pop), W = 360, H = 250, l = 38, r = 250, u0 = 272, u1 = 348, top = 14, bot = 204;
  const X = (d) => l + (r - l) * d / 110, Y = (v) => bot - (bot - top) * v / 0.55;
  let s = `<svg class="scatter" viewBox="0 0 ${W} ${H}" role="img" aria-label="Recombination frequency against distance between markers">`;
  for (const v of [0, 0.1, 0.2, 0.3, 0.4, 0.5]) {
    s += `<line x1="${l}" x2="${u1}" y1="${Y(v)}" y2="${Y(v)}" stroke="${v === 0.5 ? '#b9b6ad' : '#e9e7e0'}" ${v === 0.5 ? 'stroke-dasharray="3 3"' : ''}/><text x="${l - 6}" y="${Y(v) + 3.5}" text-anchor="end" font-size="10">${v * 100}%</text>`;
  }
  for (const d of [0, 25, 50, 75, 100]) s += `<text x="${X(d)}" y="${bot + 14}" text-anchor="middle" font-size="10">${d}</text>`;
  s += `<text x="${(l + r) / 2}" y="${bot + 30}" text-anchor="middle" font-size="10.5">Distance between markers (cM)</text>
    <text x="${(u0 + u1) / 2}" y="${bot + 14}" text-anchor="middle" font-size="10">different</text><text x="${(u0 + u1) / 2}" y="${bot + 26}" text-anchor="middle" font-size="10">chromosomes</text>
    <rect x="${u0}" y="${top}" width="${u1 - u0}" height="${bot - top}" fill="#f1f0fa" opacity=".6" rx="4"/>`;
  let path = '';
  for (let d = 0; d <= 110; d += 2) path += `${d ? 'L' : 'M'}${X(d).toFixed(1)} ${Y(haldane(d)).toFixed(1)}`;
  s += `<path d="${path}" fill="none" stroke="#8a8880" stroke-width="1.5" stroke-dasharray="4 3"/><text x="${X(82)}" y="${Y(haldane(82)) + 16}" font-size="10">expected</text>`;
  let u = 0;
  // alternate label side in order of distance so neighbours don't collide
  const dists = [];
  for (let i = 0; i < MK.length; i++) for (let j = i + 1; j < MK.length; j++) if (MK[i].c === MK[j].c) dists.push(Math.abs(MK[i].pos - MK[j].pos));
  dists.sort((p, q) => p - q);
  for (let i = 0; i < MK.length; i++) for (let j = i + 1; j < MK.length; j++) {
    const v = rf[i][j];
    if (v == null) continue;
    const same = MK[i].c === MK[j].c, d = Math.abs(MK[i].pos - MK[j].pos);
    const x = same ? X(d) : u0 + 9 + (u1 - u0 - 18) * (u++ / 11), y = Y(v), name = `${MK[i].n}–${MK[j].n}`;
    const tip = `${name}: ${(100 * v).toFixed(1)}%` + (same ? ` (${d} cM apart)` : ' (different chromosomes)');
    s += `<g data-tip="${tip}"><circle cx="${x}" cy="${y}" r="13" fill="transparent"/><circle cx="${x}" cy="${y}" r="5" fill="#4a3aa7" stroke="#fcfcfb" stroke-width="2"/>`;
    if (same) s += `<text x="${x}" y="${dists.indexOf(d) % 2 ? y + 17 : y - 9}" text-anchor="middle" font-size="9.5" font-weight="600">${name}</text>`;
    s += '</g>';
  }
  $('#sim-scatter').innerHTML = s + '</svg>' + (sim.pop.length ? '' : '<p class="empty">Simulate offspring to see the points.</p>');
}

function renderSim() {
  drawPop();
  renderHeat($('#sim-heat'), $('#sim-detail'), sim.pop, sim.heat);
  renderScatter();
}
function renderLegend() {
  $('#sim-legend').innerHTML = sim.view === 'hap'
    ? '<span><i class="sw p1"></i> from blue parent</span><span><i class="sw p2"></i> from red parent</span>'
    : '<span><i class="sw p1"></i> homozygous blue</span><span><i class="sw het"></i> heterozygous</span><span><i class="sw p2"></i> homozygous red</span>';
}

// ---------- animated single offspring ----------
// Phases: 0 paired homologs, 1 crossovers, 2 one chromatid chosen per
// chromosome, 3 gametes fly together, 4 offspring with base calls.
const M_Y = [0, 11, 28, 39], M_H = 9;
const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
let anim = null;

function meiosisCol(m, w) {
  return `<div class="mcol" data-w="${w}"><h3>${w ? 'Pollen' : 'Egg'}</h3>${CHRS.map((ch, c) => `<div class="mbiv" data-c="${c}" style="width:${pct(ch.len, CHRS[0].len)}">${chrLab(c)}${mkHead(c)}
    <div class="strands">${[0, 1, 2, 3].map((r) => `<div class="strand" data-r="${r}" style="top:${M_Y[r]}px;background:${COL[r < 2 ? 0 : 1]}"></div>`).join('')}
    ${ch.markers.map(([, p]) => `<i class="mline" style="left:${pct(p, ch.len)}"></i>`).join('')}
    ${m.chr[c].xos.map((o, k) => {
      const [r1, r2] = m.chr[c].links[k].slice().sort((p, q) => p - q), y1 = M_Y[r1] + M_H / 2, y2 = M_Y[r2] + M_H / 2;
      return `<div class="xo" style="left:${pct(o.x, ch.len)};top:${y1}px;height:${y2 - y1}px"><b>✕</b></div>`;
    }).join('')}</div></div>`).join('')}</div>`;
}

function animate() {
  const st = $('#stage'), ms = [simMeiosis(), simMeiosis()], kid = makeKid(ms[0].g, ms[1].g), num = sim.pop.length + 1;
  const nxo = ms.reduce((s, m) => s + m.chr.reduce((t, b) => t + b.xos.length, 0), 0);
  const caps = [
    'Meiosis in the F1: each chromosome has been copied into two sister chromatids, and the blue and red homologs pair up.',
    nxo ? `Crossovers (✕) swap segments between non-sister chromatids: ${nxo} in these two meioses.` : 'No crossovers happened in these two meioses.',
    'Each gamete receives one chromatid of each chromosome.',
    'The egg and the pollen come together.',
    `Offspring ${num}: the base calls under each marker are read from the two chromosomes it inherited (egg on top, pollen below).`,
  ];
  st.innerHTML = `<div class="meio">${meiosisCol(ms[0], 0)}${meiosisCol(ms[1], 1)}</div>
    <p class="hint" id="cap"></p>
    <div class="kid"><div class="id">Offspring ${num}</div>${karyo(kid.g, { head: true, calls: true })}</div>`;
  const timers = [];
  let finished = false;
  const strand = (w, c, r) => st.querySelector(`.mcol[data-w="${w}"] .mbiv[data-c="${c}"] .strand[data-r="${r}"]`);
  const setPhase = (p) => {
    st.dataset.phase = p;
    $('#cap').textContent = caps[p];
    if (p >= 1) ms.forEach((m, w) => CHRS.forEach((ch, c) => m.chr[c].rows.forEach((segs, r) => { strand(w, c, r).style.background = gradient(segs, ch.len); })));
    if (p >= 2) ms.forEach((m, w) => CHRS.forEach((ch, c) => strand(w, c, m.chr[c].pick).classList.add('picked')));
  };
  const finish = () => {
    if (finished) return;
    finished = true; timers.forEach(clearTimeout); anim = null;
    st.querySelectorAll('.fly').forEach((f) => f.remove());
    setPhase(4);
    sim.pop.push(kid); renderSim();
  };
  const fly = () => {
    setPhase(3);
    const base = st.getBoundingClientRect();
    ms.forEach((m, w) => CHRS.forEach((ch, c) => {
      const a = strand(w, c, m.chr[c].pick).getBoundingClientRect();
      const b = st.querySelectorAll('.kid .chrom')[c].querySelectorAll('.bar')[w].getBoundingClientRect();
      const f = document.createElement('div');
      f.className = 'fly';
      f.style.cssText = `left:${a.left - base.left}px;top:${a.top - base.top}px;width:${a.width}px;height:${a.height}px;background:${gradient(m.chr[c].rows[m.chr[c].pick], ch.len)}`;
      st.appendChild(f);
      f.getBoundingClientRect(); // commit the start position before transitioning
      f.style.transform = `translate(${b.left - a.left}px, ${b.top - a.top}px) scale(${b.width / a.width}, ${b.height / a.height})`;
    }));
  };
  anim = { finish };
  setPhase(0);
  if (reduceMotion) { finish(); return; }
  timers.push(setTimeout(() => setPhase(1), 900), setTimeout(() => setPhase(2), 2300), setTimeout(fly, 3500), setTimeout(finish, 4600));
}

$('#sim-head').innerHTML = CHRS.map((ch, c) => `<div class="chrom" style="flex:${ch.len}">${chrLab(c)}${mkHead(c)}</div>`).join('');
// A click during an animation completes the offspring in progress first.
$('#sim-one').onclick = () => { if (anim) anim.finish(); animate(); };
$('#sim-many').onclick = () => {
  if (anim) anim.finish();
  for (let k = +$('#sim-x').value; k > 0; k--) sim.pop.push(makeKid(simGamete(), simGamete()));
  renderSim();
};
$('#sim-x').onchange = (e) => { $('#sim-x-label').textContent = e.target.value; };
$('#sim-reset').onclick = () => {
  if (anim) anim.finish();
  sim.pop = []; sim.heat.sel = null;
  $('#stage').innerHTML = '<p class="empty">Simulate one offspring to watch the F1 make an egg and a pollen grain.</p>';
  delete $('#stage').dataset.phase;
  renderSim();
};
document.querySelectorAll('.seg button').forEach((b) => {
  b.onclick = () => {
    sim.view = b.dataset.view;
    document.querySelectorAll('.seg button').forEach((o) => o.classList.toggle('on', o === b));
    renderLegend(); drawPop();
  };
});

// ---------- tooltip, resize ----------
const tip = $('#tip');
let tipTimer;
function showTip(e) {
  const g = e.target.closest && e.target.closest('[data-tip]');
  if (!g) { tip.hidden = true; return; }
  tip.textContent = g.dataset.tip; tip.hidden = false;
  const w = tip.offsetWidth;
  tip.style.left = Math.max(8, Math.min(window.innerWidth - w - 8, e.clientX - w / 2)) + 'px';
  tip.style.top = Math.max(8, e.clientY - tip.offsetHeight - 14) + 'px';
  clearTimeout(tipTimer);
  if (e.type === 'click') tipTimer = setTimeout(() => { tip.hidden = true; }, 3000);
}
document.addEventListener('pointermove', (e) => { if (e.pointerType === 'mouse') showTip(e); });
document.addEventListener('click', showTip);
window.addEventListener('scroll', () => { tip.hidden = true; }, { passive: true });
let rz;
window.addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(drawPop, 150); });

renderLegend(); renderSim();
