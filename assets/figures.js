/* figures.js: the home page figure, the three interactive figures, the rotating contact link, the bio notes, the writings list, the photo carousel, the back link, and the footer clock.
   All data here is synthetic and generated in the browser. */
(function () {
  "use strict";

  var $ = function (s, el) { return (el || document).querySelector(s); };
  var $$ = function (s, el) { return Array.prototype.slice.call((el || document).querySelectorAll(s)); };
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var SVGNS = "http://www.w3.org/2000/svg";

  // utils
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(R) { var u = 1 - R(), v = R(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function svg(tag, attrs, parent) {
    var el = document.createElementNS(SVGNS, tag);
    for (var k in attrs) el.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(el);
    return el;
  }
  function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }
  function signed(x, d) { var s = Math.abs(x).toFixed(d); return (x > 0 ? "+" : x < 0 ? "−" : "") + s; }
  function trim(x) { return String(+x.toFixed(1)); }
  function put(fig, sel, html) { var el = $(sel, fig); if (el) el.innerHTML = html; }

  // Bridson's Poisson-disc sampling
  function poisson(W, H, r, R) {
    var k = 24, cell = r / Math.SQRT2, gw = Math.ceil(W / cell), gh = Math.ceil(H / cell);
    var grid = new Int32Array(gw * gh).fill(-1), pts = [], active = [];
    function add(x, y) { pts.push([x, y]); active.push(pts.length - 1); grid[Math.floor(y / cell) * gw + Math.floor(x / cell)] = pts.length - 1; }
    add(R() * W, R() * H);
    while (active.length) {
      var ai = Math.floor(R() * active.length), p = pts[active[ai]], found = false;
      for (var t = 0; t < k; t++) {
        var a = R() * Math.PI * 2, d = r * (1 + R());
        var x = p[0] + Math.cos(a) * d, y = p[1] + Math.sin(a) * d;
        if (x < 0 || y < 0 || x >= W || y >= H) continue;
        var gx = Math.floor(x / cell), gy = Math.floor(y / cell), ok = true;
        for (var yy = Math.max(0, gy - 2); yy <= Math.min(gh - 1, gy + 2) && ok; yy++) {
          for (var xx = Math.max(0, gx - 2); xx <= Math.min(gw - 1, gx + 2); xx++) {
            var j = grid[yy * gw + xx];
            if (j >= 0) { var q = pts[j], dx = q[0] - x, dy = q[1] - y; if (dx * dx + dy * dy < r * r) { ok = false; break; } }
          }
        }
        if (ok) { add(x, y); found = true; break; }
      }
      if (!found) active.splice(ai, 1);
    }
    return pts;
  }

  // tissue model

  var TYPES = [
    { id: "pod", name: "Podocyte", v: "--c-pod" },
    { id: "mes", name: "Mesangial", v: "--c-mes" },
    { id: "endo", name: "Endothelial", v: "--c-endo" },
    { id: "pec", name: "Parietal epithelial", v: "--c-pec" },
    { id: "pt", name: "Proximal tubule", v: "--c-pt" },
    { id: "dt", name: "Distal tubule", v: "--c-dt" },
    { id: "fib", name: "Fibroblast", v: "--c-fib" },
    { id: "imm", name: "Immune", v: "--c-imm" }
  ];
  var TI = {};
  TYPES.forEach(function (t, i) { TI[t.id] = i; });

  function buildTissue(seed, W, H, r, nGloms) {
    var R = mulberry32(seed), gloms = [], tubs = [], tries = 0;
    nGloms = nGloms || 5;
    while (gloms.length < nGloms && tries++ < 6000) {
      var gr = 54 + R() * 18;
      var gx = gr + 30 + R() * (W - 2 * gr - 60), gy = gr + 24 + R() * (H - 2 * gr - 48);
      if (gloms.every(function (g) { return Math.hypot(g.x - gx, g.y - gy) > g.r + gr + 80; }))
        gloms.push({ x: gx, y: gy, r: gr, thr: 0.08 + R() * 0.8, dir: R() * Math.PI * 2 });
    }
    tries = 0;
    var nT = Math.round(40 * (W * H) / 620000);
    while (tubs.length < nT && tries++ < 12000) {
      var rx = 20 + R() * 22, ry = rx * (0.55 + R() * 0.5), tx = R() * W, ty = R() * H, rot = R() * Math.PI, rr = Math.max(rx, ry);
      if (gloms.every(function (g) { return Math.hypot(g.x - tx, g.y - ty) > g.r + 30 + rr; }) &&
          tubs.every(function (t) { return Math.hypot(t.x - tx, t.y - ty) > t.rr + rr + 8; }))
        tubs.push({ x: tx, y: ty, rx: rx, ry: ry, cos: Math.cos(-rot), sin: Math.sin(-rot), rr: rr, kind: R() < 0.62 ? "pt" : "dt" });
    }
    var hot = [];
    for (var h = 0; h < Math.max(2, Math.round(4 * W * H / 620000)); h++) hot.push({ x: R() * W, y: R() * H, r: 60 + R() * 70 });

    var pts = poisson(W, H, r, R), cells = [];
    for (var i = 0; i < pts.length; i++) {
      var x = pts[i][0], y = pts[i][1], type = null, g = -1, ang = 0, inTub = false;
      for (var gi = 0; gi < gloms.length; gi++) {
        var G = gloms[gi], d = Math.hypot(x - G.x, y - G.y);
        if (d < G.r + 17) {
          ang = Math.atan2(y - G.y, x - G.x);
          var rel = d / G.r;
          if (d > G.r + 8) type = "pec";
          else if (d > G.r) type = null;                 // Bowman's space
          else type = rel > 0.8 ? "pod" : rel < 0.3 ? "mes" : (R() < 0.55 ? "endo" : "mes");
          g = gi;
          break;
        }
      }
      if (g < 0) {
        for (var ti = 0; ti < tubs.length; ti++) {
          var T = tubs[ti], dx = x - T.x, dy = y - T.y;
          var u = (dx * T.cos - dy * T.sin) / T.rx, v = (dx * T.sin + dy * T.cos) / T.ry, q = Math.hypot(u, v);
          if (q < 1) { inTub = true; type = q < 0.5 ? null : T.kind; break; }
        }
        if (!inTub) {
          var inHot = hot.some(function (hh) { return Math.hypot(x - hh.x, y - hh.y) < hh.r * (0.5 + R() * 0.6); });
          type = inHot && R() < 0.5 ? "imm" : (R() < 0.14 ? "endo" : "fib");
        }
      }
      if (type) cells.push({ x: x, y: y, t: TI[type], g: type === "pec" ? -1 : g, ang: ang, j: R() });
    }
    return { W: W, H: H, gloms: gloms, cells: cells };
  }

  function readColors() {
    var c = TYPES.map(function (t) { return cssVar(t.v); });
    c.scl = cssVar("--c-scl");
    return c;
  }

  function paint(ctx, tissue, scale, dot, colorOf, alphaOf) {
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    var cells = tissue.cells;
    for (var i = 0; i < cells.length; i++) {
      var c = cells[i], a = alphaOf ? alphaOf(c) : 1;
      if (a < 0.02) continue;
      ctx.globalAlpha = a;
      ctx.fillStyle = colorOf(c);
      ctx.beginPath();
      ctx.arc(c.x * scale, c.y * scale, dot * (c.state ? c.state.r : 1), 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
  }

  // hero (home)

  function Hero(fig) {
    var plate = $("[data-plate]", fig);
    var canvas = document.createElement("canvas");
    canvas.setAttribute("role", "img");
    canvas.setAttribute("aria-label", "Simulated kidney tissue section drawn as grey dots. Cell types appear in color under the cursor.");
    plate.appendChild(canvas);
    var ctx = canvas.getContext("2d");
    var narrow = window.innerWidth < 600;
    var tissue = narrow ? buildTissue(4021, 1000, 560, 11, 4) : buildTissue(4021, 1000, 300, 8.5, 3);
    var grey, color, dpr = 1, scale = 1, cssW = 0;
    var lens = { x: 0, y: 0, r: 0, target: 0 }, raf = 0, hideTimer;

    function layers() {
      var colors = readColors(), ink = cssVar("--ink");
      grey = document.createElement("canvas"); color = document.createElement("canvas");
      grey.width = color.width = canvas.width; grey.height = color.height = canvas.height;
      var dot = Math.max(1.3, 2.9 * scale) * dpr;
      // monochrome layer: glomeruli darkest, tubules mid, interstitium lightest
      paint(grey.getContext("2d"), tissue, scale * dpr, dot, function () { return ink; }, function (c) {
        return c.g >= 0 || c.t === 3 ? 0.42 : (c.t === 4 || c.t === 5) ? 0.27 : 0.14;
      });
      paint(color.getContext("2d"), tissue, scale * dpr, dot, function (c) { return colors[c.t]; });
    }
    function size() {
      var w = plate.clientWidth;
      if (!w || w === cssW) return;
      cssW = w; dpr = Math.min(2, window.devicePixelRatio || 1); scale = w / tissue.W;
      var h = Math.round(tissue.H * scale);
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); canvas.style.height = h + "px";
      layers(); draw();
    }
    function draw() {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(grey, 0, 0);
      if (lens.r > 0.5) {
        ctx.save();
        ctx.beginPath(); ctx.arc(lens.x * dpr, lens.y * dpr, lens.r * dpr, 0, Math.PI * 2); ctx.clip();
        ctx.drawImage(color, 0, 0);
        ctx.restore();
        ctx.strokeStyle = cssVar("--faint"); ctx.lineWidth = dpr;
        ctx.beginPath(); ctx.arc(lens.x * dpr, lens.y * dpr, lens.r * dpr, 0, Math.PI * 2); ctx.stroke();
      }
    }
    function step() {
      raf = 0;
      var d = lens.target - lens.r;
      lens.r = reduceMotion || Math.abs(d) < 0.3 ? lens.target : lens.r + d * 0.2;
      draw();
      if (lens.r !== lens.target) raf = requestAnimationFrame(step);
    }
    function kick() { if (!raf) raf = requestAnimationFrame(step); }
    function at(e) { var r = canvas.getBoundingClientRect(); lens.x = e.clientX - r.left; lens.y = e.clientY - r.top; }

    canvas.addEventListener("pointermove", function (e) { at(e); lens.target = 70; kick(); if (e.pointerType !== "mouse") { clearTimeout(hideTimer); hideTimer = setTimeout(function () { lens.target = 0; kick(); }, 1400); } });
    canvas.addEventListener("pointerdown", function (e) { at(e); lens.target = 70; kick(); clearTimeout(hideTimer); hideTimer = setTimeout(function () { lens.target = 0; kick(); }, 1400); });
    canvas.addEventListener("pointerleave", function (e) { if (e.pointerType === "mouse") { lens.target = 0; kick(); } });
    new ResizeObserver(size).observe(plate);
    size();
  }

  // tissue figure

  function Tissue(fig) {
    var plate = $("[data-plate]", fig), readout = $("[data-readout]", fig);
    var narrow = window.innerWidth < 560;
    var tissue = buildTissue(20261108, 1000, 620, narrow ? 11 : 9);
    var canvas = document.createElement("canvas");
    canvas.setAttribute("role", "img");
    canvas.setAttribute("aria-label", "Simulated kidney tissue section with glomeruli, tubules, and other cells colored by type.");
    plate.appendChild(canvas);
    var ctx = canvas.getContext("2d"), colors = readColors();
    var injury = 0, focus = -1, pinned = -1, hover = null, raf = 0, dpr = 1, scale = 1, cssW = 0;
    var scanStart = reduceMotion ? null : performance.now() + 200;
    var BUCKET = 24, bucket = {};
    tissue.cells.forEach(function (c, i) { var k = Math.floor(c.x / BUCKET) + "," + Math.floor(c.y / BUCKET); (bucket[k] = bucket[k] || []).push(i); });

    function angDist(a, b) { var d = Math.abs(a - b) % (2 * Math.PI); return d > Math.PI ? 2 * Math.PI - d : d; }

    // Glomerulus i is affected once injury passes its threshold (focal);
    // inside it, a sector around a random direction grows with injury (segmental).
    function targets() {
      var halves = tissue.gloms.map(function (G) { var e = clamp((injury - G.thr) * 2.2, 0, 1); return e > 0 ? 0.25 + e * 0.75 * Math.PI * 0.9 : 0; });
      var left = 0, affected = halves.filter(function (h) { return h > 0; }).length;
      var f = pinned >= 0 ? pinned : focus;
      tissue.cells.forEach(function (c) {
        var lost = false, scl = false;
        if (c.g >= 0 && halves[c.g] > 0 && angDist(c.ang, tissue.gloms[c.g].dir) < halves[c.g] + (c.j - 0.5) * 0.35) {
          if (c.t === 0) lost = true; else scl = true;
        }
        if (c.t === 0 && !lost) left++;
        var a = lost ? 0 : scl ? 0.55 : 1;
        if (f >= 0 && c.t !== f) a *= 0.1;
        c.target = { a: a, sc: scl, r: lost ? 0.2 : 1 };
        if (!c.state) c.state = { a: a, sc: scl, r: c.target.r };
      });
      put(fig, "[data-stat=pods]", left);
      put(fig, "[data-stat=gloms]", affected);
    }

    function step() {
      raf = 0;
      var moving = false, k = reduceMotion ? 1 : 0.16;
      tissue.cells.forEach(function (c) {
        var s = c.state, t = c.target, da = t.a - s.a, dr = t.r - s.r;
        if (Math.abs(da) > 0.004 || Math.abs(dr) > 0.004) { s.a += da * k; s.r += dr * k; moving = true; } else { s.a = t.a; s.r = t.r; }
        s.sc = t.sc;
      });
      var scan = Infinity;
      if (scanStart !== null) {
        var p = clamp((performance.now() - scanStart) / 1400, 0, 1);
        p = 1 - Math.pow(1 - p, 3);
        scan = p * tissue.W;
        if (p < 1) moving = true; else scanStart = null;
      }
      var dot = Math.max(1.4, 3.1 * scale) * dpr;
      paint(ctx, tissue, scale * dpr, dot,
        function (c) { return c.state.sc ? colors.scl : colors[c.t]; },
        function (c) { return c.x > scan ? 0 : c.state.a; });
      if (scan < tissue.W) {
        ctx.fillStyle = cssVar("--accent"); ctx.globalAlpha = 0.6;
        ctx.fillRect(scan * scale * dpr, 0, Math.max(1, dpr), canvas.height); ctx.globalAlpha = 1;
      }
      if (hover) {
        ctx.strokeStyle = cssVar("--ink"); ctx.lineWidth = dpr;
        ctx.beginPath(); ctx.arc(hover.x * scale * dpr, hover.y * scale * dpr, dot + 3.5 * dpr, 0, Math.PI * 2); ctx.stroke();
      }
      if (moving) raf = requestAnimationFrame(step);
    }
    function kick() { if (!raf) raf = requestAnimationFrame(step); }
    function update() { targets(); kick(); }

    function size() {
      var w = plate.clientWidth;
      if (!w || w === cssW) return;
      cssW = w; dpr = Math.min(2, window.devicePixelRatio || 1); scale = w / tissue.W;
      var h = Math.round(tissue.H * scale);
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); canvas.style.height = h + "px";
      kick();
    }

    function nearest(x, y) {
      var bx = Math.floor(x / BUCKET), by = Math.floor(y / BUCKET), best = null, bd = 400;
      for (var yy = by - 1; yy <= by + 1; yy++) for (var xx = bx - 1; xx <= bx + 1; xx++) {
        (bucket[xx + "," + yy] || []).forEach(function (i) {
          var c = tissue.cells[i], d = (c.x - x) * (c.x - x) + (c.y - y) * (c.y - y);
          if (d < bd && c.state && c.state.a > 0.3) { bd = d; best = c; }
        });
      }
      return best;
    }
    canvas.addEventListener("pointermove", function (e) {
      var r = canvas.getBoundingClientRect(), c = nearest((e.clientX - r.left) / scale, (e.clientY - r.top) / scale);
      hover = c;
      if (c) {
        var where = c.g >= 0 ? "glomerulus " + (c.g + 1) : c.t === 3 ? "Bowman's capsule" : (c.t === 4 || c.t === 5) ? "tubule" : "interstitium";
        readout.innerHTML = "<b>" + (c.state.sc ? "Scarred tissue" : TYPES[c.t].name) + "</b>, " + where;
      } else readout.textContent = "";
      kick();
    });
    canvas.addEventListener("pointerleave", function () { hover = null; readout.textContent = ""; kick(); });

    var counts = TYPES.map(function () { return 0; });
    tissue.cells.forEach(function (c) { counts[c.t]++; });
    var ul = $("[data-legend]", fig);
    var btns = TYPES.map(function (t, i) {
      var li = document.createElement("li"), b = document.createElement("button");
      b.type = "button"; b.setAttribute("aria-pressed", "false");
      b.innerHTML = '<i style="--sw:var(' + t.v + ')"></i><span>' + t.name + "</span><span>" + counts[i].toLocaleString() + "</span>";
      b.addEventListener("pointerenter", function (e) { if (e.pointerType === "mouse") { focus = i; update(); } });
      b.addEventListener("pointerleave", function () { focus = -1; update(); });
      b.addEventListener("click", function () {
        pinned = pinned === i ? -1 : i;
        btns.forEach(function (bb, j) { bb.setAttribute("aria-pressed", String(j === pinned)); });
        update();
      });
      li.appendChild(b); ul.appendChild(li);
      return b;
    });

    var input = $("[data-injury]", fig), out = $("[data-injury-out]", fig);
    input.addEventListener("input", function () { injury = input.value / 100; out.textContent = input.value + "%"; update(); });
    new ResizeObserver(size).observe(plate);
    size();
    update();
  }

  // spectra figure

  var LIB = [
    [57.0699, 8], [69.0335, 12], [81.0335, 18], [95.0491, 22], [107.0491, 35], [121.0648, 64],
    [133.0648, 28], [147.0804, 100], [159.0804, 41], [175.0754, 57], [189.0910, 19], [203.1067, 33],
    [217.0859, 14], [245.0808, 26], [263.0914, 48], [281.1020, 11]
  ];

  function Spectra(fig) {
    var R = mulberry32(147), QUERY = [];
    LIB.forEach(function (p, i) {
      if (i === 2 || i === 12) return;                       // two fragments not observed
      QUERY.push({ mz0: p[0], noise: gauss(R) * 2.2, I: clamp(p[1] * (0.72 + R() * 0.56), 3, 100), lib: i });
    });
    [[88.04, 9], [139.98, 6], [196.6, 12], [233.2, 7]].forEach(function (p) { QUERY.push({ mz0: p[0], noise: 0, I: p[1], lib: -1 }); });
    var mx = Math.max.apply(null, QUERY.map(function (q) { return q.I; }));
    QUERY.forEach(function (q) { q.I = q.I / mx * 100; });

    var state = { drift: 0, tol: 10, rt: 2, model: false, recal: false, last: null };
    var W = 680, H = 300, MID = 150, GAP = 10, X0 = 44, X1 = 312, PAD = 14, EH = 72, EMAX = 40;
    function x(mz) { return PAD + (mz - X0) / (X1 - X0) * (W - 2 * PAD); }
    function ex(ppm) { return PAD + (ppm + EMAX) / (2 * EMAX) * (W - 2 * PAD); }

    var plateEl = $("[data-plate]", fig), errEl = $("[data-err]", fig);
    var plot = svg("svg", { viewBox: "0 0 " + W + " " + H, role: "img", "aria-label": "Mirror plot of a query spectrum above a library spectrum." }, plateEl);
    svg("line", { x1: PAD, x2: W - PAD, y1: MID - GAP, y2: MID - GAP, class: "axis" }, plot);
    svg("line", { x1: PAD, x2: W - PAD, y1: MID + GAP, y2: MID + GAP, class: "axis" }, plot);
    for (var t = 50; t <= 300; t += 50) { var tx = svg("text", { x: x(t), y: MID + 4, "text-anchor": "middle" }, plot); tx.textContent = t; }
    var gq = svg("g", {}, plot), gl = svg("g", {}, plot), qLines = [], lLines = [];
    var tip = document.createElement("div"); tip.className = "spec-tip"; plateEl.appendChild(tip);

    LIB.forEach(function (p, i) {
      var h = p[1] / 100 * (MID - GAP - 20);
      lLines.push(svg("line", { x1: x(p[0]), x2: x(p[0]), y1: MID + GAP, y2: MID + GAP + h, class: "pk pk-l" }, gl));
      var hit = svg("line", { x1: x(p[0]), x2: x(p[0]), y1: MID + GAP, y2: MID + GAP + h, class: "pk-hit" }, gl);
      hit.addEventListener("pointerenter", function () { showTip("l", i, x(p[0])); });
      hit.addEventListener("pointerleave", hideTip);
    });
    QUERY.forEach(function (q, i) {
      var h = q.I / 100 * (MID - GAP - 20);
      qLines.push(svg("line", { x1: x(q.mz0), x2: x(q.mz0), y1: MID - GAP, y2: MID - GAP - h, class: "pk pk-q" + (q.lib < 0 ? " pk-n" : "") }, gq));
      var hit = svg("line", { x1: x(q.mz0), x2: x(q.mz0), y1: MID - GAP, y2: MID - GAP - h, class: "pk-hit" }, gq);
      hit.addEventListener("pointerenter", function () { showTip("q", i, x(q.mz0)); });
      hit.addEventListener("pointerleave", hideTip);
    });

    var err = svg("svg", { viewBox: "0 0 " + W + " " + EH, role: "img", "aria-label": "Mass error of each true fragment in parts per million." }, errEl);
    var band = svg("rect", { y: 6, height: EH - 28, class: "err-band" }, err);
    var curve = svg("path", { class: "err-curve", d: "" }, err);
    svg("line", { x1: ex(0), x2: ex(0), y1: 4, y2: EH - 20, class: "axis" }, err);
    for (var pp = -40; pp <= 40; pp += 10) svg("text", { x: ex(pp), y: EH - 5, "text-anchor": "middle" }, err).textContent = pp > 0 ? "+" + pp : pp;
    var dots = [];
    QUERY.forEach(function (q, i) { if (q.lib >= 0) dots.push({ q: i, el: svg("circle", { cx: ex(0), cy: 13 + (i % 5) * 8, r: 3, class: "err-dot" }, err) }); });

    function match(tol, offset, model) {
      var pairs = [];
      QUERY.forEach(function (q, qi) {
        LIB.forEach(function (l, li) {
          var ppm = (q.mz - l[0]) / l[0] * 1e6 - offset;
          if (Math.abs(ppm) <= tol) pairs.push({ qi: qi, li: li, ppm: ppm, w: Math.sqrt(q.I * l[1]) });
        });
      });
      pairs.sort(function (a, b) { return b.w - a.w; });
      var usedQ = {}, usedL = {}, ms = [], sigma = tol / 3, num = 0, nq = 0, nl = 0;
      pairs.forEach(function (p) { if (!usedQ[p.qi] && !usedL[p.li]) { usedQ[p.qi] = usedL[p.li] = 1; ms.push(p); } });
      QUERY.forEach(function (q) { nq += q.I; });
      LIB.forEach(function (l) { nl += l[1]; });
      ms.forEach(function (m) { m.g = model ? Math.exp(-0.5 * Math.pow(m.ppm / sigma, 2)) : 1; num += m.w * m.g; });
      return { matches: ms, cos: num / Math.sqrt(nq * nl), usedQ: usedQ, usedL: usedL };
    }
    function median(a) { if (!a.length) return 0; var s = a.slice().sort(function (p, q) { return p - q; }), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }

    function showTip(side, i, px) {
      var r = plot.getBoundingClientRect(), text;
      if (side === "q") {
        var q = QUERY[i], m = state.last.matches.filter(function (mm) { return mm.qi === i; })[0];
        text = "m/z " + q.mz.toFixed(4) + (q.lib < 0 ? ", noise" : m ? ", error " + signed(m.ppm, 1) + " ppm" + (state.model ? ", weight " + m.g.toFixed(2) : "") : ", no match");
      } else {
        text = "Library m/z " + LIB[i][0].toFixed(4) + (state.last.usedL[i] ? ", matched" : ", no match");
      }
      tip.textContent = text;
      tip.style.left = clamp(px / W * r.width, 110, r.width - 110) + "px";
      tip.style.top = (side === "q" ? 30 : r.height - 26) + "px";
      tip.classList.add("on");
    }
    function hideTip() { tip.classList.remove("on"); }

    function render() {
      QUERY.forEach(function (q) { q.mz = q.mz0 * (1 + (state.drift + q.noise) / 1e6); });
      var offset = state.recal ? median(match(60, 0, false).matches.map(function (m) { return m.ppm; })) : 0;
      var res = state.last = match(state.tol, offset, state.model);
      qLines.forEach(function (ln, i) { ln.classList.toggle("m", !!res.usedQ[i]); });
      lLines.forEach(function (ln, i) { ln.classList.toggle("m", !!res.usedL[i]); });

      var ink = cssVar("--ink"), faint = cssVar("--faint");
      dots.forEach(function (d) {
        var q = QUERY[d.q], l = LIB[q.lib], ppm = (q.mz - l[0]) / l[0] * 1e6 - offset;
        d.el.setAttribute("cx", ex(clamp(ppm, -EMAX - 2, EMAX + 2)));
        d.el.style.fill = res.usedQ[d.q] ? ink : faint;
      });
      band.setAttribute("x", ex(-state.tol));
      band.setAttribute("width", ex(state.tol) - ex(-state.tol));
      var sigma = state.tol / 3, dpath = "";
      for (var p = -EMAX; p <= EMAX; p += 0.5) dpath += (p === -EMAX ? "M" : "L") + ex(p).toFixed(1) + " " + (6 + (EH - 28) * (1 - Math.exp(-0.5 * Math.pow(p / sigma, 2)))).toFixed(1);
      curve.setAttribute("d", dpath);
      curve.style.opacity = state.model ? 0.85 : 0;

      var rtFit = Math.exp(-0.5 * Math.pow(state.rt / 5, 2));
      put(fig, "[data-stat=matched]", res.matches.length);
      put(fig, "[data-stat=cos]", res.cos.toFixed(3));
      put(fig, "[data-stat=score]", (res.cos * rtFit).toFixed(3));
      put(fig, "[data-offset]", state.recal ? "Drift estimate <b>" + signed(offset, 1) + " ppm</b>" : "");
    }

    function bind(key, fmt) {
      var input = $("[data-" + key + "]", fig), out = $("[data-" + key + "-out]", fig);
      input.addEventListener("input", function () { state[key] = parseFloat(input.value); out.textContent = fmt(state[key]); render(); });
    }
    bind("drift", function (v) { return signed(v, v % 1 ? 1 : 0) + " ppm"; });
    bind("tol", function (v) { return trim(v) + " ppm"; });
    bind("rt", function (v) { return trim(v) + " s"; });
    ["model", "recal"].forEach(function (k) {
      var c = $("[data-" + k + "]", fig);
      c.checked = false;
      c.addEventListener("change", function () { state[k] = c.checked; render(); });
    });
    render();
  }

  // exclusivity figure

  function Exclusivity(fig) {
    var state = { dev: 10, pte: true, orphan: false, ped: false, gain: false };
    var A = -15, B = 25;
    function pos(y) { return (clamp(y, A, B) - A) / (B - A); }

    function compute(s) {
      var remaining = 20 - s.dev;
      var credit = 0.5 * (s.dev - 1) + 1;                      // half of testing + 1-year review
      var ext = s.pte ? Math.max(0, Math.min(5, credit, 14 - remaining)) : 0;
      var ped = s.ped ? 0.5 : 0, gain = s.gain ? 5 : 0;
      var patentEnd = remaining + ext;
      var ends = { patent: patentEnd + ped, nce: 5 + gain + ped, orphan: s.orphan ? 7 + gain + ped : 0 };
      var bind = "Patent", entry = ends.patent;
      if (ends.nce > entry) { entry = ends.nce; bind = "NCE"; }
      if (s.orphan && ends.orphan > entry) { entry = ends.orphan; bind = "Orphan"; }
      return { remaining: remaining, ext: ext, patentEnd: patentEnd, ped: ped, gain: gain, entry: entry, bind: bind };
    }
    function seg(name, a, b) {
      var e = $('[data-seg="' + name + '"]', fig);
      e.style.setProperty("--a", pos(a));
      e.style.setProperty("--b", pos(Math.max(a, b)));
      e.style.opacity = b > a ? 1 : 0;
    }

    var axis = $("[data-axis]", fig), lines = $("[data-lines]", fig);
    for (var y = -15; y <= 25; y += 5) {
      var s = document.createElement("span");
      s.style.left = pos(y) * 100 + "%";
      s.textContent = y > 0 ? "+" + y : y < 0 ? "−" + Math.abs(y) : "0";
      axis.appendChild(s);
      if (y !== 0) { var i = document.createElement("i"); i.style.left = pos(y) * 100 + "%"; lines.appendChild(i); }
    }
    var zero = document.createElement("span"); zero.className = "xmark zero"; zero.setAttribute("data-label", "Approval"); zero.style.left = pos(0) * 100 + "%"; lines.appendChild(zero);
    var entry = document.createElement("span"); entry.className = "xmark entry"; entry.setAttribute("data-label", "Generic entry"); lines.appendChild(entry);

    function render() {
      var c = compute(state);
      seg("dev", -state.dev, 0);
      seg("patent", -state.dev, c.remaining);
      seg("pte", c.remaining, c.patentEnd);
      seg("patent-ped", c.patentEnd, c.patentEnd + c.ped);
      seg("nce", 0, 5);
      seg("nce-gain", 5, 5 + c.gain);
      seg("nce-ped", 5 + c.gain, 5 + c.gain + c.ped);
      seg("orphan", 0, state.orphan ? 7 : 0);
      seg("orphan-gain", 7, state.orphan ? 7 + c.gain : 7);
      seg("orphan-ped", 7 + c.gain, state.orphan ? 7 + c.gain + c.ped : 7 + c.gain);
      $('[data-row="orphan"]', fig).setAttribute("data-off", String(!state.orphan));
      entry.style.left = pos(c.entry) * 100 + "%";
      entry.classList.toggle("flip", pos(c.entry) > 0.6);
      put(fig, "[data-stat=entry]", c.entry.toFixed(1));
      put(fig, "[data-stat=bind]", c.bind);
    }

    var input = $("[data-dev]", fig), out = $("[data-dev-out]", fig);
    input.addEventListener("input", function () { state.dev = parseFloat(input.value); out.textContent = trim(state.dev) + " years"; render(); });
    $$("[data-opt]", fig).forEach(function (c) {
      var k = c.getAttribute("data-opt");
      c.checked = state[k];
      c.addEventListener("change", function () { state[k] = c.checked; render(); });
    });

    // Hovering a law in the list highlights the rows it created.
    var related = { orphan: ["orphan"], "patent nce": ["patent", "nce"], ped: ["patent", "nce", "orphan"], gain: ["nce", "orphan"] };
    $$("[data-hist] li[data-link]").forEach(function (li) {
      var rows = related[li.getAttribute("data-link")] || [];
      li.addEventListener("pointerenter", function () {
        fig.setAttribute("data-focus", "");
        $$(".xrow", fig).forEach(function (r) { r.classList.toggle("is-focus", rows.indexOf(r.getAttribute("data-row")) >= 0); });
      });
      li.addEventListener("pointerleave", function () { fig.removeAttribute("data-focus"); $$(".xrow", fig).forEach(function (r) { r.classList.remove("is-focus"); }); });
    });
    render();
  }

  // boot

  var MAP = { hero: Hero, tissue: Tissue, spectra: Spectra, exclusivity: Exclusivity };
  $$("[data-fig]").forEach(function (fig) {
    var fn = MAP[fig.getAttribute("data-fig")];
    if (fn) { try { fn(fig); } catch (e) { if (window.console) console.error(e); } }
  });

  // Opens or closes a block by animating its height, so whatever is below it moves with it.
  function slide(el, state) {
    clearTimeout(state.settle);
    var from = el.hidden ? 0 : el.getBoundingClientRect().height;
    el.hidden = false;
    el.style.height = "auto";
    var to = state.open ? el.getBoundingClientRect().height : 0;
    el.style.height = from + "px";
    void el.offsetWidth;
    el.style.height = to + "px";
    el.classList.toggle("on", state.open);
    state.settle = setTimeout(function () {
      el.style.height = "";
      el.hidden = !state.open;
    }, 260);
  }

  // Notes that open as small sub-paragraphs under a word in the bio.
  $$("[data-reveal]").forEach(function (bio) {
    var items = $$(".reveal", bio).map(function (btn) {
      return { btn: btn, note: document.getElementById(btn.getAttribute("aria-controls")), open: false, inside: false, closer: 0, settle: 0 };
    });

    function set(item, next) {
      clearTimeout(item.closer);
      if (next === item.open) return;
      item.open = next;
      slide(item.note, item);
      item.btn.setAttribute("aria-expanded", String(next));
    }
    function show(item) {
      items.forEach(function (other) { if (other !== item) set(other, false); });
      updateAge();
      set(item, true);
    }

    items.forEach(function (item) {
      function enter(e) { if (e.pointerType === "mouse") { item.inside = true; show(item); } }
      function leave(e) { if (e.pointerType === "mouse") { item.inside = false; item.closer = setTimeout(function () { set(item, false); }, 0); } }
      [item.btn, item.note].forEach(function (el) { el.addEventListener("pointerenter", enter); el.addEventListener("pointerleave", leave); });
      item.btn.addEventListener("click", function () { if (item.inside) return; if (item.open) set(item, false); else show(item); });
    });
    bio.addEventListener("keydown", function (e) {
      if (e.key !== "Escape") return;
      items.forEach(function (item) { if (item.open) { set(item, false); item.btn.focus(); } });
    });
  });

  // Whole years, then the days left over, since the camera's date (Providence time).
  function updateAge() {
    var el = $("[data-age]");
    if (!el) return;
    var parts = {};
    new Intl.DateTimeFormat("en-US", { timeZone: "America/New_York", year: "numeric", month: "numeric", day: "numeric" })
      .formatToParts(new Date()).forEach(function (p) { parts[p.type] = +p.value; });
    var today = Date.UTC(parts.year, parts.month - 1, parts.day);
    var years = parts.year - 2015;
    if (today < Date.UTC(parts.year, 1, 15)) years--;
    var days = Math.round((today - Date.UTC(2015 + years, 1, 15)) / 864e5);
    el.textContent = years + (years === 1 ? " year" : " years") + " and " + days + (days === 1 ? " day" : " days");
  }
  updateAge();

  // Writings: hovering an entry opens it for as long as the cursor stays; clicking keeps it open.
  $$("[data-days]").forEach(function (list) {
    var still = null, closers = [], days = $$(".day", list);
    days.forEach(function (day, i) {
      // Numbered from the oldest piece, which sits at the bottom of the list.
      var count = days.length - i, number = document.createElement("span");
      number.className = "day-num";
      number.textContent = (count < 10 ? "0" : "") + count + " /";
      $(".day-head", day).prepend(number, " ");
      var head = $(".day-head", day), body = $(".day-body", day);
      var state = { open: false, settle: 0 }, hovering = false, pinned = false;

      function apply(e) {
        var next = pinned || hovering;
        if (next === state.open) return;
        state.open = next;
        // Remember the cursor, so a row that slides under a still cursor doesn't open.
        if (e) still = { x: e.clientX, y: e.clientY };
        slide(body, state);
        head.setAttribute("aria-expanded", String(next));
      }
      day.addEventListener("pointerenter", function (e) {
        if (e.pointerType !== "mouse" || (still && still.x === e.clientX && still.y === e.clientY)) return;
        hovering = true; apply(e);
      });
      day.addEventListener("pointerleave", function (e) {
        if (e.pointerType !== "mouse") return;
        hovering = false; apply(e);
      });
      head.addEventListener("click", function () { pinned = !pinned; apply(); });
      body.addEventListener("click", function () {
        if (!pinned && String(window.getSelection()) === "") { pinned = true; apply(); }
      });
      closers.push(function (e) {
        if (e.type === "click" && day.contains(e.target)) return;
        pinned = false;
        if (e.type === "keydown") hovering = false;
        apply();
      });
    });
    document.addEventListener("click", function (e) { closers.forEach(function (close) { close(e); }); });
    document.addEventListener("keydown", function (e) {
      if (e.key !== "Escape") return;
      // An open piece takes this press, so the next one is what leaves the page.
      if ($(".day-head[aria-expanded=\"true\"]", list)) e.preventDefault();
      closers.forEach(function (close) { close(e); });
    });
  });

  // Photo carousel. Offsets are in photos from the center. Photos are one size until hovered.
  function Reel(reel) {
    var shots = $$(".shot", reel), imgs = shots.map(function (shot) { return $("img", shot); }), N = shots.length;
    var GAP = 8, SPEED = reduceMotion ? 0 : 0.12, GROW = 0.25, GRID_GROW = 0.05, REST = 0.8, DWELL = 320, ROWS = [2, 2, 2, 2];
    var SLIDE = reduceMotion ? 1 : 780, STAGGER = reduceMotion ? 0 : 56, TOTAL = SLIDE + STAGGER * (N - 1);
    var at = 4, v = SPEED, W = 0, last = 0, mouse = null, target = null, pending = -1, pendingSince = 0;
    var grid = false, gridHover = -1, change = null, frozen = false, slots = [], gridHeight = 0;
    var pose = shots.map(function () { return { x: 0, y: 0, s: REST }; });
    var lift = shots.map(function () { return 0; });

    function wrap(d) { return ((d + N / 2) % N + N) % N - N / 2; }
    function clamp(x) { return Math.max(0, Math.min(1, x)); }
    function glide(x) { x = clamp(x); return x * x * x * (x * (6 * x - 15) + 10); }
    function strip() { return !grid && !change && !frozen; }

    function layoutGrid() {
      var C = reel.clientWidth, top = 0, n = 0;
      slots = [];
      ROWS.forEach(function (count) {
        var w = (C - (count - 1) * GAP) / count, h = w * 0.75;
        for (var c = 0; c < count; c++, n++) {
          slots.push({ x: c * (w + GAP) + w / 2 - C / 2, y: top + h / 2 - W * 0.375, s: REST * w / W });
        }
        top += h + GAP;
      });
      gridHeight = top - GAP;
    }
    function stripPose(i, high) {
      var d = wrap(i - at), side = d < 0 ? -1 : 1;
      return { x: d * (W + GAP) + side * (high - lift[i]) * GROW * W / 2, y: 0, s: REST * (1 + GROW * lift[i]) };
    }
    // Every axis speeds up from rest and slows back to rest. To the grid x leads and y
    // follows, and back to the strip y leads and x follows, so photos travel in a curve.
    function moving(i, goal) {
      var m = clamp((change.t - change.order[i] * STAGGER) / SLIDE), from = change.from[i];
      var lead = glide(m * 1.25), trail = glide((m - 0.2) * 1.25);
      var ax = change.toGrid ? lead : trail, ay = change.toGrid ? trail : lead;
      return { x: from.x + (goal.x - from.x) * ax, y: from.y + (goal.y - from.y) * ay, s: from.s + (goal.s - from.s) * glide(m) };
    }
    function render() {
      var high = Math.max.apply(null, lift);
      shots.forEach(function (shot, i) {
        var goal = grid ? { x: slots[i].x, y: slots[i].y, s: slots[i].s * (1 + GRID_GROW * lift[i]) } : stripPose(i, high), p = change ? moving(i, goal) : goal;
        pose[i] = p;
        shot.style.transform = "translate3d(" + p.x.toFixed(2) + "px," + p.y.toFixed(2) + "px,0) scale(" + p.s.toFixed(4) + ")";
      });
      reel.classList.toggle("open", grid && !change);
      var strip = W * 0.75;
      if (change) reel.style.height = (change.height + ((grid ? gridHeight : strip) - change.height) * glide(change.t / TOTAL)).toFixed(2) + "px";
      else reel.style.height = grid ? gridHeight.toFixed(2) + "px" : "";
    }
    function beginChange() {
      var order = shots.map(function (_, i) { return i; });
      if (!grid) {
        // Photos that end up farthest from the middle go first, so the visible ones land last.
        var far = order.slice().sort(function (a, b) { return Math.abs(wrap(b - at)) - Math.abs(wrap(a - at)); });
        far.forEach(function (i, rank) { order[i] = rank; });
      }
      change = { t: 0, toGrid: grid, order: order, height: reel.getBoundingClientRect().height, from: pose.map(function (p) { return { x: p.x, y: p.y, s: p.s }; }) };
    }
    function shotAt(x, y) {
      var el = document.elementFromPoint(x, y), shot = el && el.closest && el.closest(".shot");
      return shot && reel.contains(shot) ? shots.indexOf(shot) : -1;
    }
    function pick(i) { target = i; new Image().src = imgs[i].getAttribute("data-full"); }
    // Resting on another photo brings it to the middle once the last slide has finished.
    function consider(now) {
      var hit = Math.abs(wrap(target - at)) < 0.02 ? shotAt(mouse.x, mouse.y) : -1;
      if (hit < 0 || hit === target) { pending = -1; return; }
      if (hit !== pending) { pending = hit; pendingSince = now; return; }
      if (now - pendingSince >= DWELL) { pending = -1; pick(hit); }
    }
    function frame(now) {
      var ms = Math.min(now - last || 0, 50), dt = ms / 1000;
      last = now;
      if (change) { change.t += ms; if (change.t >= TOTAL) change = null; }
      if (strip()) {
        if (mouse && target === null) { var hit = shotAt(mouse.x, mouse.y); if (hit >= 0) pick(hit); }
        else if (mouse) consider(now);
        var hold = target !== null || !!reel.querySelector(":focus-visible");
        v += ((hold ? 0 : SPEED) - v) * (1 - Math.exp(-dt / 0.1));
      } else v = 0;
      gridHover = grid && !change && !frozen && mouse ? shotAt(mouse.x, mouse.y) : -1;
      var ease = reduceMotion ? 1 : 1 - Math.exp(-dt / 0.14);
      if (target !== null && strip()) at += wrap(target - at) * ease;
      else if (strip()) at += v * dt;
      at = (at % N + N) % N;
      lift = lift.map(function (l, i) { return l + ((i === gridHover || (i === target && strip()) ? 1 : 0) - l) * ease; });
      render();
      requestAnimationFrame(frame);
    }

    function measure() { W = shots[0].offsetWidth * REST; layoutGrid(); render(); }
    reel.classList.add("live");
    new ResizeObserver(measure).observe(reel.parentNode);
    measure();
    imgs.forEach(function (img) {
      function show() { img.classList.add("in"); }
      if (img.complete) show(); else { img.addEventListener("load", show); img.addEventListener("error", show); }
    });

    reel.addEventListener("pointermove", function (e) { if (e.pointerType === "mouse") mouse = { x: e.clientX, y: e.clientY }; });
    reel.addEventListener("pointerleave", function (e) { if (e.pointerType === "mouse") { mouse = target = null; pending = -1; } });
    reel.addEventListener("pointerdown", function (e) {
      if (e.pointerType === "mouse" || !strip()) return;
      var hit = shotAt(e.clientX, e.clientY);
      if (hit >= 0) pick(hit);
    });
    ["pointerup", "pointercancel"].forEach(function (ev) { window.addEventListener(ev, function (e) { if (e.pointerType !== "mouse" && !frozen) target = null; }); });
    reel.addEventListener("focusin", function (e) { var i = imgs.indexOf(e.target); if (i >= 0 && e.target.matches(":focus-visible") && strip()) pick(i); });
    reel.addEventListener("focusout", function () { if (!mouse) target = null; });

    // Full-size view: the photo grows out of its place in the page and returns to it.
    var box = document.createElement("dialog"), big = document.createElement("img"), current = -1, closing = false, byKey = false;
    box.className = "lightbox";
    box.appendChild(big);
    document.body.appendChild(box);

    function flight(from, to) {
      return "translate(" + (from.left + from.width / 2 - to.left - to.width / 2) + "px," + (from.top + from.height / 2 - to.top - to.height / 2) + "px) scale(" + from.width / to.width + ")";
    }
    function open(i, keyboard) {
      if (box.open) return;
      byKey = keyboard;
      var img = imgs[i], full = img.getAttribute("data-full");
      current = i; frozen = true; mouse = target = null; pending = -1;
      big.src = img.currentSrc || img.src;
      big.alt = img.alt;
      document.documentElement.style.overflow = "hidden";
      box.showModal();
      var start = shots[i].getBoundingClientRect(), end = big.getBoundingClientRect();
      shots[i].style.visibility = "hidden";
      big.animate([{ transform: flight(start, end) }, { transform: "none" }], { duration: reduceMotion ? 1 : 360, easing: "cubic-bezier(.2,.7,.2,1)" });
      void box.offsetWidth;
      box.classList.add("on");
      var hi = new Image();
      hi.onload = function () { if (current === i) big.src = full; };
      hi.src = full;
    }
    function close() {
      if (!box.open || closing) return;
      closing = true;
      var shot = shots[current], end = big.getBoundingClientRect(), back = shot.getBoundingClientRect();
      box.classList.remove("on");
      big.animate([{ transform: "none" }, { transform: flight(back, end) }], { duration: reduceMotion ? 1 : 280, easing: "cubic-bezier(.4,0,.2,1)", fill: "forwards" }).onfinish = function (e) {
        shot.style.visibility = "";
        box.close();
        // Closing hands focus back to the photo. Only keyboard users need it there;
        // for everyone else it would leave a focus ring and hold the carousel still.
        if (!byKey) imgs[current].blur();
        e.target.cancel();
        big.removeAttribute("src");
        document.documentElement.style.overflow = "";
        closing = false; frozen = false; current = -1;
      };
    }
    box.addEventListener("click", function (e) { if (e.target !== big) close(); });
    box.addEventListener("cancel", function (e) { e.preventDefault(); close(); });
    reel.addEventListener("click", function (e) { var i = imgs.indexOf(e.target); if (i >= 0) open(i, false); });
    reel.addEventListener("keydown", function (e) {
      var i = imgs.indexOf(e.target);
      if (i >= 0 && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); open(i, true); }
    });

    var toggle = $("[data-reel-switch]");
    if (toggle) {
      toggle.parentNode.hidden = false;
      toggle.addEventListener("click", function () {
        grid = !grid;
        mouse = target = null; pending = -1;
        beginChange();
        toggle.textContent = grid ? "carousel" : "grid";
      });
    }
    requestAnimationFrame(frame);
    if (reduceMotion) setInterval(function () { if (strip() && target === null && !reel.querySelector(":focus-visible")) at = (Math.round(at) + 1) % N; }, 5000);
  }
  $$("[data-reel]").forEach(Reel);

  // Resume pop-up.
  $$("[data-resume]").forEach(function (link) {
    var box = document.createElement("dialog");
    box.className = "lightbox popup";
    box.setAttribute("aria-label", "R\u00e9sum\u00e9");
    box.innerHTML = '<div class="resume-sheet" tabindex="0"><img src="assets/resume.png" width="1760" height="2278" alt="R\u00e9sum\u00e9 of Nithin Venkatramanan"></div>' +
      '<div class="popup-foot"><a download="Nithin_Venkatramanan_Resume.pdf">download</a></div>';
    document.body.appendChild(box);
    var sheet = $(".resume-sheet", box), save = $("a", box), foot = $(".popup-foot", box), picture = $("img", box);
    save.href = link.getAttribute("href");
    if (picture.decode) picture.decode().catch(function () {});

    // A folded square of paper sweeps out of the word on a crescent arc and opens along its creases:
    // first down the rows, then out to both sides. Closing plays the same folds backward.
    var stage = document.createElement("div"), root = document.createElement("div"), panels = {};
    var t = 0, dir = 0, raf = 0, last = 0, byKey = false, word, page, cellW, cellH;
    var OPEN = reduceMotion ? 1 : 1150, CLOSE = reduceMotion ? 1 : 800;
    var cast = document.createElement("div");
    cast.className = "fold-cast";
    stage.className = "fold";
    root.className = "fold-root";
    stage.hidden = true;
    stage.appendChild(cast);
    stage.appendChild(root);
    box.appendChild(stage);

    function unit(x) { return Math.max(0, Math.min(1, x)); }
    function curve(x) { x = unit(x); return x * x * x * (x * (6 * x - 15) + 10); }
    function span(a, b) { return unit((t - a) / (b - a)); }

    function layer(cls, parent) {
      var el = document.createElement("div");
      el.className = cls;
      parent.appendChild(el);
      return el;
    }
    function panel(col, row, parent, name) {
      var el = layer("fold-panel", parent), front = layer("fold-face", el), back = layer("fold-face fold-back", el);
      panels[name] = { el: el, front: front, ghost: layer("fold-ghost", back), shade: layer("fold-shade", front), backShade: layer("fold-shade", back), col: col, row: row };
      return el;
    }
    function row(r, parent) {
      var el = panel(1, r, parent, "m" + r);
      panel(0, r, el, "l" + r);
      panel(2, r, el, "r" + r);
      return el;
    }
    var middle = row(1, root);
    row(0, middle);
    row(3, row(2, middle));

    function prepare() {
      word = link.getBoundingClientRect();
      page = sheet.getBoundingClientRect();
      cellW = page.width / 3;
      cellH = page.height / 4;
      root.style.width = cellW + "px";
      root.style.height = cellH + "px";
      cast.style.width = page.width + "px";
      cast.style.height = page.height + "px";
      Object.keys(panels).forEach(function (name) {
        var p = panels[name], kind = name.charAt(0), e = p.el.style;
        e.width = cellW + "px"; e.height = cellH + "px";
        e.left = kind === "l" ? -cellW + "px" : kind === "r" ? cellW + "px" : "0";
        e.top = kind === "m" ? (name === "m0" ? -cellH + "px" : name === "m1" ? "0" : cellH + "px") : "0";
        e.transformOrigin = kind === "l" ? "100% 50%" : kind === "r" ? "0 50%" : name === "m0" ? "50% 100%" : "50% 0";
        [p.front, p.ghost].forEach(function (face) {
          face.style.backgroundImage = "url(" + picture.getAttribute("src") + ")";
          face.style.backgroundSize = page.width + "px " + page.height + "px";
          // Faces overlap by a pixel, so no hairline shows between panels.
          face.style.backgroundPosition = 1 - p.col * cellW + "px " + (1 - p.row * cellH) + "px";
        });
        // Each panel is darkest at its crease and lightens toward its free edge.
        var away = kind === "l" ? "to left" : kind === "r" ? "to right" : name === "m0" ? "to top" : "to bottom";
        p.backShade.style.background = "#000000";
        p.shade.style.backgroundImage = "linear-gradient(" + away + ", rgba(0,0,0,.3), rgba(0,0,0,0))";
      });
    }
    function bend(name, axis, sign, progress, fold) {
      var angle = fold * (1 - curve(progress)), p = panels[name];
      p.el.style.transform = "rotate" + axis + "(" + (sign * angle).toFixed(2) + "deg)";
      var light = Math.abs(Math.sin(angle * Math.PI / 180));
      p.shade.style.opacity = light.toFixed(3);
      p.backShade.style.opacity = (0.1 * light).toFixed(3);
    }
    function draw() {
      var fly = curve(span(0, 0.3)), rest = curve(span(0.3, 1));
      var wx = word.left + word.width / 2, wy = word.top + word.height / 2;
      var mx = page.left + page.width / 2, my = page.top + page.height / 2, ey = page.top + page.height * 0.375;
      // Crescent: a curve whose midpoint bows out to the side of the straight line from the word.
      var dx = mx - wx, dy = my - wy, reach = Math.max(180, Math.hypot(dx, dy) * 0.7), out = wx <= mx ? -1 : 1;
      var length = Math.hypot(dx, dy) || 1, cx = (wx + mx) / 2 + out * Math.abs(dy) / length * reach, cy = (wy + my) / 2 - out * dx / length * reach * (dy < 0 ? -1 : 1);
      var x = (1 - fly) * (1 - fly) * wx + 2 * fly * (1 - fly) * cx + fly * fly * mx;
      var y = (1 - fly) * (1 - fly) * wy + 2 * fly * (1 - fly) * cy + fly * fly * my + (ey - my) * rest;
      var spin = 26 * out * (1 - fly);
      var size = 0.12 + 0.88 * fly;
      root.style.opacity = unit(t * 40);
      // The packet leans back while it opens, so each fold reads in depth, then lies flat.
      var lean = fly * (1 - rest);
      root.style.transform = "translate(" + (x - cellW / 2).toFixed(2) + "px," + (y - cellH / 2).toFixed(2) + "px) rotate(" + spin.toFixed(2) + "deg) rotateX(" + (18 * lean).toFixed(2) + "deg) rotateY(" + (-14 * lean).toFixed(2) + "deg) scale(" + size.toFixed(4) + ")";

      var down = span(0.3, 0.46), up = span(0.42, 0.58), low = span(0.54, 0.7), left = span(0.66, 0.84), right = span(0.8, 1);
      var tall = curve(down) + curve(up) + curve(low), wide = curve(left) + curve(right);
      var shadowH = cellH * size * (1 + tall), shadowW = cellW * size * (1 + wide);
      var shadowTop = y - cellH * size / 2 - curve(up) * cellH * size + 26, shadowLeft = x - cellW * size / 2 - curve(left) * cellW * size;
      cast.style.transform = "translate(" + shadowLeft.toFixed(1) + "px," + shadowTop.toFixed(1) + "px) scale(" + (shadowW / page.width).toFixed(4) + "," + (shadowH / page.height).toFixed(4) + ")";
      cast.style.opacity = (fly * (1 - curve(span(0.86, 1)))).toFixed(3);

      bend("m1", "X", 1, 1, 0);
      bend("m2", "X", -1, down, 176);
      bend("m0", "X", 1, up, 170);
      bend("m3", "X", -1, low, 173);
      for (var r = 0; r < 4; r++) {
        bend("l" + r, "Y", 1, left, 175);
        bend("r" + r, "Y", -1, right, 178);
      }
    }
    function done() {
      raf = 0;
      if (dir > 0) {
        sheet.style.visibility = foot.style.visibility = "";
        requestAnimationFrame(function () { if (dir > 0 && t === 1) stage.hidden = true; });
        foot.animate([{ opacity: 0 }, { opacity: 1 }], { duration: reduceMotion ? 1 : 200, easing: "ease-out" });
      } else {
        stage.hidden = true;
        box.close();
        if (!byKey) link.blur();
        document.documentElement.style.overflow = "";
      }
    }
    function step(now) {
      var ms = Math.min(Math.max(0, now - last), 24);
      last = now;
      t = unit(t + dir * ms / (dir > 0 ? OPEN : CLOSE));
      draw();
      if ((dir > 0 && t === 1) || (dir < 0 && t === 0)) done();
      else raf = requestAnimationFrame(step);
    }
    function run(next) {
      dir = next;
      stage.hidden = false;
      draw();
      foot.style.visibility = "hidden";
      requestAnimationFrame(function () { if (dir && !(dir > 0 && t === 1)) sheet.style.visibility = "hidden"; });
      last = performance.now();
      if (!raf) raf = requestAnimationFrame(step);
    }
    function open(keyboard) {
      if (box.open) return;
      byKey = keyboard;
      document.documentElement.style.overflow = "hidden";
      box.showModal();
      sheet.focus({ preventScroll: true });
      prepare();
      t = 0;
      void box.offsetWidth;
      box.classList.add("on");
      run(1);
    }
    function close() {
      if (!box.open || dir < 0) return;
      box.classList.remove("on");
      run(-1);
    }
    link.addEventListener("click", function (e) {
      if (e.metaKey || e.ctrlKey || e.shiftKey || e.button !== 0) return;
      e.preventDefault();
      open(e.detail === 0);
    });
    box.addEventListener("click", function (e) { if (!sheet.contains(e.target) && e.target !== save) close(); });
    box.addEventListener("cancel", function (e) { e.preventDefault(); close(); });
  });

  // Escape on an inner page goes home, once anything open on top has been closed first.
  var homeLink = $("[data-back]");
  if (homeLink) document.addEventListener("keydown", function (e) {
    if (e.key === "Escape" && !e.defaultPrevented && !$("dialog[open]")) homeLink.click();
  });

  // Back link: return to the previous page (keeping its scroll position) when it was this site.
  $$("[data-back]").forEach(function (a) {
    a.addEventListener("click", function (e) {
      try {
        if (history.length > 1 && document.referrer && new URL(document.referrer).origin === location.origin) { e.preventDefault(); history.back(); }
      } catch (err) {}
    });
  });

  // footer clock, Providence time, ticking on the second
  var clock = $("[data-clock]");
  if (clock) {
    var clockFmt = new Intl.DateTimeFormat("en-US", { hour: "numeric", minute: "2-digit", second: "2-digit", timeZone: "America/New_York" });
    (function tickClock() {
      var now = new Date();
      clock.textContent = "Made with <3 by Nithin Venkatramanan in Providence, RI, where it\u2019s " + clockFmt.format(now) + ".";
      setTimeout(tickClock, 1000 - now.getMilliseconds());
    })();
  }

  // Assembled here so the address is not in the page source for scrapers.
  $$("[data-mail]").forEach(function (a) { a.href = "mailto:" + ["nithin", "venkatramanan"].join("_") + "@" + ["brown", "edu"].join("."); });

  // rotating contact link
  var PLATFORMS = [
    { label: "GitHub", href: "https://github.com/ninv29" },
    { label: "LinkedIn", href: "https://www.linkedin.com/in/nithin-venkatramanan" },
    { label: "Instagram", href: "https://www.instagram.com/nithinvenkatramanan/" }
  ];

  function Rotator(wrap) {
    var a = $("a", wrap), idx = 0, paused = false;
    var measure = document.createElement("span");
    measure.className = "rot-measure";
    measure.setAttribute("aria-hidden", "true");
    wrap.appendChild(measure);

    function widthOf(text) { measure.textContent = text; return Math.ceil(measure.getBoundingClientRect().width); }
    function letters(text, cls) {
      var box = document.createElement("span");
      box.setAttribute("aria-hidden", "true");
      text.split("").forEach(function (ch, i) {
        var l = document.createElement("span");
        l.className = cls ? "l " + cls : "l";
        l.style.setProperty("--i", i);
        l.textContent = ch;
        box.appendChild(l);
      });
      return box;
    }
    function fit() { a.style.width = widthOf(PLATFORMS[idx].label) + "px"; }

    // The outgoing word fades out on top while the new one drops in, so there is no blank gap.
    function swap() {
      var old = a.firstElementChild;
      var item = PLATFORMS[idx = (idx + 1) % PLATFORMS.length];
      var fresh = letters(item.label, "pre");
      old.className = "rot-old";
      a.appendChild(fresh);
      a.href = item.href;
      a.setAttribute("aria-label", item.label);
      fit();
      void a.offsetWidth;
      $$(".l", old).forEach(function (l) { l.classList.add("out"); });
      $$(".l", fresh).forEach(function (l) { l.classList.remove("pre"); });
      setTimeout(function () { old.remove(); }, 500);
    }
    function tick() {
      if (!paused && !document.hidden) swap();
      setTimeout(tick, 2000);
    }

    a.textContent = "";
    a.appendChild(letters(PLATFORMS[0].label));
    a.href = PLATFORMS[0].href;
    a.setAttribute("aria-label", PLATFORMS[0].label);
    fit();
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(fit);
    ["pointerenter", "focusin"].forEach(function (ev) { wrap.addEventListener(ev, function () { paused = true; }); });
    ["pointerleave", "focusout"].forEach(function (ev) { wrap.addEventListener(ev, function () { paused = false; }); });
    setTimeout(tick, 2000);
  }
  $$("[data-rotate]").forEach(Rotator);
})();
