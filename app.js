'use strict';

// ---------- genome ----------
// Positions are in centimorgans; one crossover on a bivalent is worth 50 cM.
const CHRS = [
  { name: 'Chromosome 1', len: 120, cen: 62, markers: [['A', 6], ['B', 16], ['C', 44], ['D', 112]] },
  { name: 'Chromosome 2', len: 80, cen: 48, markers: [['E', 8], ['F', 30], ['G', 72]] },
];
const MK = [];
CHRS.forEach((ch, c) => ch.markers.forEach(([n, pos]) => MK.push({ n, c, pos })));
const TOTAL = CHRS.reduce((s, ch) => s + ch.len, 0);
const COL = ['#2a78d6', '#eb6834'], HET = '#b9b6ad';
const MAX_XO = 4, MAX_ROWS = 800;

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

function karyo(haps, { head = false, lab = false } = {}) {
  return `<div class="karyo">${CHRS.map((ch, c) => `<div class="chrom" style="flex:${ch.len}">
    ${lab ? `<div class="lab">Chr ${c + 1}</div>` : ''}${head ? mkHead(c) : ''}
    <div class="bars">${haps.map((h) => `<div class="bar" style="background:${gradient(h[c], ch.len)}"><i class="cen" style="left:${pct(ch.cen, ch.len)}"></i></div>`).join('')}
    ${ch.markers.map(([, p]) => `<i class="tick" style="left:${pct(p, ch.len)}"></i>`).join('')}</div></div>`).join('')}</div>`;
}
const genoName = (m, g) => (g === 0 ? m.n + m.n : g === 1 ? m.n + m.n.toLowerCase() : m.n.toLowerCase().repeat(2));
function chips(haps) {
  return `<div class="chips">${MK.map((m, i) => {
    const sep = i && MK[i - 1].c !== m.c ? '<span class="sep"></span>' : '';
    if (haps.length === 1) {
      const a = alleleAt(haps[0][m.c], m.pos);
      return `${sep}<span class="chip a${a}">${a ? m.n.toLowerCase() : m.n}</span>`;
    }
    const g = alleleAt(haps[0][m.c], m.pos) + alleleAt(haps[1][m.c], m.pos);
    return `${sep}<span class="chip ${g === 1 ? 'h' : 'a' + g / 2}">${genoName(m, g)}</span>`;
  }).join('')}</div>`;
}

$('#cross').innerHTML = `<div class="cross">
  <div><div class="who">Blue parent</div>${karyo([pure(0), pure(0)])}</div><div class="x">×</div>
  <div><div class="who">Orange parent</div>${karyo([pure(1), pure(1)])}</div>
  <div class="arrow">↓</div>
  <div class="f1"><div class="who">F1 (selfed to make the offspring)</div>${karyo([pure(0), pure(1)], { head: true })}${chips([pure(0), pure(1)])}</div></div>`;

// ---------- heat map + detail ----------
function mix(t) {
  const a = [241, 240, 250], b = [74, 58, 167];
  return `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * t)).join(',')})`;
}
function renderHeat(el, detailEl, pop, st, showExpected) {
  const rf = rfMatrix(pop), S = 44, L = 24, T = 38, GAP = 8;
  const off = (i) => i * S + (MK[i].c ? GAP : 0);
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
    renderHeat(el, detailEl, pop, st, showExpected);
  };
  if (!st.sel) { detailEl.innerHTML = ''; return; }
  const [i, j] = st.sel, a = MK[i], b = MK[j], t = pairCounts(pop, i, j), r = rf[i][j];
  const same = a.c === b.c;
  let exp = '';
  if (showExpected) {
    exp = same
      ? `<p class="note">These markers are ${Math.abs(a.pos - b.pos)} cM apart on chromosome ${a.c + 1}. With unlimited offspring the estimate would settle near ${(100 * expectedRF(i, j)).toFixed(1)}%.</p>`
      : '<p class="note">These markers are on different chromosomes, so they assort independently. With unlimited offspring the estimate would settle at 50%.</p>';
  }
  detailEl.innerHTML = `<div class="detail"><div>Markers <b>${a.n}</b> and <b>${b.n}</b> · ${pop.length} offspring</div>
    <div class="big">${r == null ? '–' : (100 * r).toFixed(1) + '%'}</div>${exp}
    <table><tr><th></th>${[0, 1, 2].map((g) => `<th>${genoName(b, g)}</th>`).join('')}</tr>
    ${[0, 1, 2].map((g) => `<tr><th>${genoName(a, g)}</th>${t[g].map((v) => `<td>${v}</td>`).join('')}</tr>`).join('')}</table>
    <p class="note">Number of offspring with each combination of genotypes at the two markers.</p></div>`;
}

// ---------- Activity 1: by hand ----------
const KEY = 'recomb-hand-v1';
const hand = { pop: [], step: 0, phase: 'place', xos: [[], []], products: null, egg: null, pollen: null, heat: {} };
try { hand.pop = (JSON.parse(localStorage.getItem(KEY)) || []).map(([g1, g2]) => makeKid(g1, g2)); } catch (e) { hand.pop = []; }
const saveHand = () => { try { localStorage.setItem(KEY, JSON.stringify(hand.pop.map((k) => k.g))); } catch (e) { /* storage unavailable */ } };
const ROW_Y = [0, 17, 51, 68], ROW_H = 14;

function bivalentHTML(c) {
  const ch = CHRS[c], { rows, links } = resolveBivalent(ch, hand.xos[c]);
  return `<div class="biv-wrap"><div class="chrom"><div class="lab">${ch.name}</div>${mkHead(c)}</div>
    <div class="bivalent" data-c="${c}"><div class="strands">
    ${rows.map((segs, r) => `<div class="strand" style="top:${ROW_Y[r]}px;background:${gradient(segs, ch.len)}"><i class="cenp" style="left:${pct(ch.cen, ch.len)};top:2px"></i></div>`).join('')}
    ${ch.markers.map(([, p]) => `<i class="mline" style="left:${pct(p, ch.len)}"></i>`).join('')}
    ${hand.xos[c].map((o, k) => {
      const [r1, r2] = links[k].sort((p, q) => p - q), y1 = ROW_Y[r1] + ROW_H / 2, y2 = ROW_Y[r2] + ROW_H / 2;
      return `<div class="xo" style="left:${pct(o.x, ch.len)};top:${y1}px;height:${y2 - y1}px"><b data-k="${k}" role="button" aria-label="Remove crossover">✕</b></div>`;
    }).join('')}</div></div></div>`;
}

function gameteCard(g, pickLabel) {
  return `${pickLabel ? `<span class="pick">${pickLabel}</span>` : ''}${karyo([g])}${chips([g])}`;
}

function renderStage() {
  const which = hand.step === 0 ? 'egg' : 'pollen';
  document.querySelectorAll('#steps li').forEach((li) => {
    const s = +li.dataset.s;
    li.className = s === hand.step ? 'on' : s < hand.step ? 'done' : '';
  });
  const st = $('#stage');
  if (hand.step === 2) {
    st.innerHTML = `<h3>Fertilization: egg + pollen = offspring ${hand.pop.length}</h3>
      <div class="pair"><div><h3>Egg</h3><div class="gam static">${gameteCard(hand.egg)}</div></div><div><h3>Pollen</h3><div class="gam static">${gameteCard(hand.pollen)}</div></div></div>
      <h3>Offspring genotype at each marker</h3>
      <div class="kid new">${karyo(hand.pop[hand.pop.length - 1].g, { head: true, lab: true })}${chips(hand.pop[hand.pop.length - 1].g)}</div>
      <p class="hint">The offspring got one chromosome of each pair from the egg and one from the pollen. It has been added to your population below.</p>
      <div class="controls"><button class="btn primary" id="again">Make another offspring</button></div>`;
    $('#again').onclick = () => { hand.step = 0; hand.phase = 'place'; hand.xos = [[], []]; renderStage(); st.scrollIntoView({ behavior: 'smooth', block: 'start' }); };
    return;
  }
  if (hand.phase === 'place') {
    const n = hand.xos[0].length + hand.xos[1].length;
    st.innerHTML = `<h3>Make the ${which}: meiosis in the F1</h3>
      <p class="note">Each chromosome has been copied into two sister chromatids, and the blue and orange homologs have paired up. <b>Tap a chromosome to place a crossover</b> between two non-sister chromatids. Tap ✕ to remove one.</p>
      ${bivalentHTML(0)}${bivalentHTML(1)}
      <p class="hint">${n ? `${n} crossover${n > 1 ? 's' : ''} placed. Look at which markers now sit on a recombinant chromatid.` : 'No crossovers yet. Chromosome pairs usually have at least one.'}</p>
      <div class="controls"><button class="btn primary" id="divide">Divide the cell →</button><button class="btn ghost" id="randxo">Random crossovers</button><button class="btn ghost" id="clearxo" ${n ? '' : 'disabled'}>Clear</button></div>
      <button class="btn ghost small" id="skip">Skip: let the computer make the ${which}</button>`;
    st.querySelectorAll('.bivalent').forEach((el) => {
      el.onclick = (e) => {
        const c = +el.dataset.c, list = hand.xos[c], b = e.target.closest('b[data-k]');
        if (b) list.splice(+b.dataset.k, 1);
        else if (list.length < MAX_XO) {
          const rect = el.querySelector('.strands').getBoundingClientRect();
          const x = Math.min(CHRS[c].len - 1, Math.max(1, (e.clientX - rect.left) / rect.width * CHRS[c].len));
          list.push({ x, i: rnd(2), j: 2 + rnd(2) });
        }
        renderStage();
      };
    });
    $('#randxo').onclick = () => {
      hand.xos = CHRS.map((ch) => Array.from({ length: Math.min(MAX_XO, poisson(ch.len / 50)) }, () => ({ x: 1 + Math.random() * (ch.len - 2), i: rnd(2), j: 2 + rnd(2) })));
      renderStage();
    };
    $('#clearxo').onclick = () => { hand.xos = [[], []]; renderStage(); };
    $('#skip').onclick = () => choose(simGamete());
    $('#divide').onclick = () => {
      // Meiosis I: each homolog pair orients at random (independent assortment).
      // Meiosis II: sister chromatids separate.
      const perChr = CHRS.map((ch, c) => {
        const { rows } = resolveBivalent(ch, hand.xos[c]);
        const top = rnd(2), halves = [[rows[0], rows[1]], [rows[2], rows[3]]].map((h) => (rnd(2) ? h : [h[1], h[0]]));
        return top ? [halves[1], halves[0]] : halves;
      });
      hand.products = [0, 1].map((cell) => [0, 1].map((m) => perChr.map((p) => p[cell][m])));
      hand.phase = 'pick'; renderStage();
    };
    return;
  }
  st.innerHTML = `<h3>Four gametes: pick one to be the ${which}</h3>
    <p class="note"><b>Meiosis I</b> pulled the homologs into two cells. <b>Meiosis II</b> split the sister chromatids, giving four gametes with one copy of each chromosome.</p>
    ${hand.products.map((cell, k) => `<div class="gcell"><div class="lab">From meiosis I cell ${k + 1}</div>${cell.map((g, m) => `<button class="gam" data-k="${k}" data-m="${m}">${gameteCard(g, 'Pick')}</button>`).join('')}</div>`).join('')}
    <div class="controls"><button class="btn ghost" id="back">← Change crossovers</button></div>`;
  st.querySelectorAll('.gam').forEach((b) => { b.onclick = () => choose(hand.products[+b.dataset.k][+b.dataset.m]); });
  $('#back').onclick = () => { hand.phase = 'place'; renderStage(); };
}

function choose(g) {
  if (hand.step === 0) { hand.egg = g; hand.step = 1; } else {
    hand.pollen = g; hand.step = 2;
    hand.pop.push(makeKid(hand.egg, hand.pollen)); saveHand(); renderHandPop();
  }
  hand.phase = 'place'; hand.xos = [[], []];
  renderStage();
  $('#steps').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function renderHandPop() {
  const n = hand.pop.length;
  $('#hand-n').textContent = n ? `n = ${n}` : '';
  $('#hand-pop').innerHTML = n
    ? hand.pop.map((k, i) => `<div class="kid"><div class="id">Offspring ${i + 1}</div>${karyo(k.g)}${chips(k.g)}</div>`).reverse().join('')
    : '<p class="empty">No offspring yet. Make an egg and a pollen grain above.</p>';
  renderHeat($('#hand-heat'), $('#hand-detail'), hand.pop, hand.heat, false);
}
$('#hand-reset').onclick = () => {
  if (hand.pop.length && !confirm('Remove all of your offspring and start over?')) return;
  hand.pop = []; hand.step = 0; hand.phase = 'place'; hand.xos = [[], []]; hand.heat.sel = null;
  saveHand(); renderStage(); renderHandPop();
};

// ---------- Activity 2: simulation ----------
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
  renderHeat($('#sim-heat'), $('#sim-detail'), sim.pop, sim.heat, true);
  renderScatter();
}
function addSim(n) { for (let k = 0; k < n; k++) sim.pop.push(makeKid(simGamete(), simGamete())); renderSim(); }
function renderLegend() {
  $('#sim-legend').innerHTML = sim.view === 'hap'
    ? '<span><i class="sw p1"></i> from blue parent</span><span><i class="sw p2"></i> from orange parent</span>'
    : '<span><i class="sw p1"></i> homozygous blue</span><span><i class="sw het"></i> heterozygous</span><span><i class="sw p2"></i> homozygous orange</span>';
}
$('#sim-head').innerHTML = CHRS.map((ch, c) => `<div style="flex:${ch.len};min-width:0">${mkHead(c)}</div>`).join('');
$('#sim-one').onclick = () => addSim(1);
$('#sim-many').onclick = () => addSim(+$('#sim-x').value);
$('#sim-x').onchange = (e) => { $('#sim-x-label').textContent = e.target.value; };
$('#sim-reset').onclick = () => { sim.pop = []; sim.heat.sel = null; renderSim(); };
document.querySelectorAll('.seg button').forEach((b) => {
  b.onclick = () => {
    sim.view = b.dataset.view;
    document.querySelectorAll('.seg button').forEach((o) => o.classList.toggle('on', o === b));
    renderLegend(); drawPop();
  };
});

// ---------- tabs, tooltip, resize ----------
document.querySelectorAll('.tab').forEach((t) => {
  t.onclick = () => {
    document.querySelectorAll('.tab').forEach((o) => { o.classList.toggle('on', o === t); o.setAttribute('aria-selected', o === t); });
    $('#tab-hand').hidden = t.dataset.tab !== 'hand';
    $('#tab-sim').hidden = t.dataset.tab !== 'sim';
    if (t.dataset.tab === 'sim') renderSim();
  };
});
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
window.addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => { if (!$('#tab-sim').hidden) drawPop(); }, 150); });

renderStage(); renderHandPop(); renderLegend();
