// 20261111_alignview
// Alignment view for the results page.
//
// Draws the same alignment the CLI hands to IGV. make_visualization_sam() in
// CRISPRlungo_visualization.py does not pass minimap2's output through: it
// rebuilds one record per read from the reference and the mutations that
// survived filtering - matches take reference bases, a substitution takes the
// alternative base, an insertion its sequence, a deletion becomes a D, and an
// inversion becomes a second, reverse-strand segment tagged YC:Z:107,189,69.
// This module rebuilds the same picture in the browser from
// window.lungoReadMutations (the per-read filtered allele string) plus each
// read's span, so one read is one row, only real editing outcomes are drawn,
// and substitutions show per base as they do in IGV.
//
// Unticking "Filtered mutations only" falls back to the raw aligner output.
//
// Unlike GuehoLab/ExCas-Analyzer's viewer, whose canvas grows one row per read
// so the page just gets longer, the canvas here sits in a box of fixed height
// with its own scrollbar, and the ruler is a separate canvas pinned above it.

(function () {
  'use strict';

  var ROW_H = 11;
  var READ_H = 8;
  var RULER_H = 56;
  var MAX_ROWS = 200000;   // the box scrolls, so there is no need to cap the stack

  var PALETTE = {
    WT: '#9aa0a6', Sub: '#2a78d6', Ins: '#4a3aa7', Del: '#eda100',
    LargeIns: '#1baf7a', LargeDel: '#d03b3b', Complex: '#e87ba4', Inv: '#6bbd45',
    fwd: '#2a78d6', rev: '#eb6834'
  };
  // IGV's base colours, for the substitution ticks
  var BASE_COLOUR = { A: '#2a9d3a', C: '#2a78d6', G: '#d79a00', T: '#d03b3b', N: '#898781' };
  var INV_COLOUR = 'rgb(107,189,69)';   // the YC tag the CLI writes for inversions


  var state = {
    items: [], rows: [], refLen: 0, viewStart: 0, viewEnd: 0,
    filtered: true,
    canvas: null, ruler: null, scrollBox: null, scrollSpacer: null, drag: null, pick: '', total: 0,
    refKey: 'WT',
    colourBy: 'igv',
    largeDel: 50, largeIns: 50, cuts: [], spacers: [], hits: [],
    lines: null, mutations: null, spans: null
  };

  // ---------------------------------------------------------------- helpers

  function baseReadId(id) {
    id = String(id);
    var at = id.indexOf('_AI_');
    if (at >= 0) return id.slice(0, at);
    at = id.indexOf('_alignStrandInfo_');
    if (at >= 0) return id.slice(0, at);
    return id.replace(/_split_\d+$/, '');
  }

  function parseCigar(text) {
    var out = [], re = /(\d+)([MIDNSHP=X])/g, m;
    while ((m = re.exec(text)) !== null) out.push({ len: +m[1], op: m[2] });
    return out;
  }

  // Where each read sits on the reference, and which way round it is.
  function spansFromSam(lines) {
    var spans = {};
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i];
      if (!line || line.charAt(0) === '@') continue;
      var f = line.split('\t');
      if (f.length < 11) continue;
      var flag = +f[1];
      if (flag & 0x100) continue;
      var pos = +f[3] - 1;
      if (!(pos >= 0)) continue;
      var ref = pos;
      parseCigar(f[5]).forEach(function (c) {
        if (c.op === 'M' || c.op === 'D' || c.op === 'N' || c.op === '=' || c.op === 'X') ref += c.len;
      });
      var id = baseReadId(f[0]);
      var s = spans[id];
      if (!s) spans[id] = { start: pos, end: ref, reverse: !!(flag & 16) };
      else { if (pos < s.start) s.start = pos; if (ref > s.end) s.end = ref; }
    }
    return spans;
  }

  // "81_264:LargeDel_184,224_224:Sub_A>G,81_265:LargeIns_ACGT_inv"
  function parseAllele(text) {
    var out = [];
    if (!text || text === 'WT') return out;
    text.split(',').forEach(function (token) {
      token = token.trim();
      if (!token || token.indexOf(':') < 0) return;
      var span = token.split(':')[0];
      var detail = token.slice(span.length + 1);
      var se = span.split('_');
      var start = Number(se[0]), end = Number(se[1]);
      if (!isFinite(start) || !isFinite(end)) return;
      var inverted = /_inv$/.test(detail);
      if (inverted) detail = detail.replace(/_inv$/, '');
      var kind = detail.split('_')[0];
      var rest = detail.slice(kind.length + 1);
      var mut = { start: start, end: end, kind: kind, inverted: inverted };
      if (kind === 'Del' || kind === 'LargeDel') mut.len = Number(rest) || (end - start + 1);
      else if (kind === 'Ins' || kind === 'LargeIns') { mut.seq = rest; mut.len = rest.length; }
      else { mut.alt = (rest.split('>')[1] || 'N').charAt(0); mut.len = 1; }
      out.push(mut);
    });
    out.sort(function (a, b) { return a.start - b.start; });
    return out;
  }

  // One record per read, rebuilt the way make_visualization_sam() does.
  function buildFromMutations(spans, mutations) {
    var items = [];
    Object.keys(mutations).forEach(function (key) {
      var id = baseReadId(key);
      var span = spans[id];
      if (!span) return;
      var muts = parseAllele(mutations[key]);
      var blocks = [], dels = [], ins = [], subs = [], inversions = [];
      var pos = span.start;
      var maxDel = 0, maxIns = 0;

      muts.forEach(function (m) {
        if (m.start > pos) { blocks.push({ s: pos, e: m.start }); pos = m.start; }
        if (m.kind === 'Del' || m.kind === 'LargeDel') {
          if (m.inverted) {
            // the CLI draws the inverted stretch as its own reverse-strand
            // segment on the same read rather than as a gap
            inversions.push({ s: m.start, e: m.start + m.len });
          } else {
            dels.push({ s: m.start, e: m.start + m.len, len: m.len });
            if (m.len > maxDel) maxDel = m.len;
          }
          pos = m.start + m.len;
        } else if (m.kind === 'Ins' || m.kind === 'LargeIns') {
          if (!m.inverted) {
            ins.push({ at: m.start, len: m.len, seq: m.seq || '' });
            if (m.len > maxIns) maxIns = m.len;
          }
        } else {
          subs.push({ at: m.start, alt: m.alt });
          if (m.start >= pos) pos = m.start;
        }
      });
      if (span.end > pos) blocks.push({ s: pos, e: span.end });
      if (!blocks.length && !inversions.length) blocks.push({ s: span.start, e: span.end });

      items.push({
        name: id, start: span.start, end: span.end, reverse: span.reverse,
        blocks: blocks, dels: dels, ins: ins, subs: subs, inversions: inversions,
        maxDel: maxDel, maxIns: maxIns, inverted: inversions.length > 0
      });
    });
    return items;
  }

  // Fallback: the aligner's own records, one row per read, nothing filtered.
  function buildFromSam(lines) {
    var byId = {};
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i];
      if (!line || line.charAt(0) === '@') continue;
      var f = line.split('\t');
      if (f.length < 11) continue;
      var flag = +f[1];
      if (flag & 0x100) continue;
      var pos = +f[3] - 1;
      if (!(pos >= 0)) continue;
      var id = baseReadId(f[0]);
      var it = byId[id] || (byId[id] = {
        name: id, start: pos, end: pos, reverse: !!(flag & 16),
        blocks: [], dels: [], ins: [], subs: [], inversions: [],
        maxDel: 0, maxIns: 0, inverted: false, strands: []
      });
      it.strands.push(!!(flag & 16));
      var ref = pos, qpos = 0, query = f[9];
      parseCigar(f[5]).forEach(function (c) {
        if (c.op === 'M' || c.op === '=' || c.op === 'X') { it.blocks.push({ s: ref, e: ref + c.len }); ref += c.len; }
        else if (c.op === 'D' || c.op === 'N') {
          it.dels.push({ s: ref, e: ref + c.len, len: c.len });
          if (c.len > it.maxDel) it.maxDel = c.len;
          ref += c.len;
        } else if (c.op === 'I') {
          it.ins.push({ at: ref, len: c.len,
            seq: (query && query !== '*') ? query.substr(qpos, c.len) : '' });
          if (c.len > it.maxIns) it.maxIns = c.len;
        }
        if (c.op === 'M' || c.op === '=' || c.op === 'X' || c.op === 'I'
            || c.op === 'S') qpos += c.len;
      });
      if (pos < it.start) it.start = pos;
      if (ref > it.end) it.end = ref;
    }
    return Object.keys(byId).map(function (k) {
      var it = byId[k];
      it.inverted = it.strands.some(function (r) { return r !== it.strands[0]; });
      // a jump between two records of one read is a deletion the aligner split
      it.blocks.sort(function (a, b) { return a.s - b.s; });
      for (var i = 1; i < it.blocks.length; i++) {
        var gap = it.blocks[i].s - it.blocks[i - 1].e;
        if (gap > 0 && !it.dels.some(function (d) { return d.s === it.blocks[i - 1].e; })) {
          it.dels.push({ s: it.blocks[i - 1].e, e: it.blocks[i].s, len: gap });
          if (gap > it.maxDel) it.maxDel = gap;
        }
      }
      return it;
    });
  }

  function outcomeOf(it) {
    if (it.inverted) return 'Inv';
    if (it.maxDel > state.largeDel && it.maxIns > state.largeIns) return 'Complex';
    if (it.maxDel > state.largeDel) return 'LargeDel';
    if (it.maxIns > state.largeIns) return 'LargeIns';
    if (it.maxDel > 0 && it.maxIns > 0) return 'Complex';
    if (it.maxDel > 0) return 'Del';
    if (it.maxIns > 0) return 'Ins';
    if (it.subs && it.subs.length) return 'Sub';
    return 'WT';
  }

  var IGV_GREY = '#c8c7c1';
  var IGV_GREY_REV = '#d7d6d0';

  function colourOf(it) {
    // 'igv' is the default and matches the track the CLI writes: every read body
    // the same neutral grey, so the mutations are the only thing with colour
    if (state.colourBy === 'igv') return it.reverse ? IGV_GREY_REV : IGV_GREY;
    if (state.colourBy === 'strand') return it.reverse ? PALETTE.rev : PALETTE.fwd;
    return PALETTE[it.outcome] || PALETTE.WT;
  }

  function prepare(items) {
    var pa = window.paReadCategory || null;
    items.forEach(function (it) {
      it.outcome = outcomeOf(it);
      it.paCategory = pa ? (pa[it.name] || null) : null;
    });
    items.sort(function (a, b) { return (a.start - b.start) || (a.end - b.end); });
    return items;
  }

  // "outcome:LargeDel" or "pa:Precise_PA"; '' is every read.
  function matchesPick(it, pick) {
    if (!pick) return true;
    var sep = pick.indexOf(':');
    var kind = pick.slice(0, sep), name = pick.slice(sep + 1);
    return kind === 'pa' ? it.paCategory === name : it.outcome === name;
  }

  function fillPicker(items) {
    var sel = document.getElementById('alignViewPick');
    if (!sel) return;
    var byOutcome = {}, byPa = {};
    items.forEach(function (it) {
      byOutcome[it.outcome] = (byOutcome[it.outcome] || 0) + 1;
      if (it.paCategory) byPa[it.paCategory] = (byPa[it.paCategory] || 0) + 1;
    });
    function group(label, counts, kind, order) {
      var names = order.filter(function (n) { return counts[n]; })
        .concat(Object.keys(counts).filter(function (n) { return order.indexOf(n) < 0; }));
      if (!names.length) return '';
      return '<optgroup label="' + label + '">' + names.map(function (n) {
        return '<option value="' + kind + ':' + n + '">' + n.replace(/_/g, ' ')
          + ' (' + counts[n].toLocaleString() + ')</option>';
      }).join('') + '</optgroup>';
    }
    var keep = state.pick;
    sel.innerHTML = '<option value="">All reads (' + items.length.toLocaleString() + ')</option>'
      + group('PA category', byPa, 'pa', PA_ORDER)
      + group('Outcome', byOutcome, 'outcome', OUTCOME_ORDER);
    var exists = Array.prototype.some.call(sel.options, function (o) { return o.value === keep; });
    sel.value = exists ? keep : '';
    state.pick = sel.value;
  }

  var OUTCOME_ORDER = ['WT', 'Sub', 'Ins', 'Del', 'LargeIns', 'LargeDel', 'Inv', 'Complex'];
  var PA_ORDER = ['Unedited', 'Precise_PA', 'Flap_install_wo_indels', 'Flap_integration_wDel',
    'Flap_integration_wIns', 'Flap_install_with_Large_deletion', 'Large_deletion',
    'Indel_at_nick_site', 'Others', 'Not_analysed'];

  function pack(items, bpPerPx) {
    var rowEnds = [], rows = [], padBp = bpPerPx * 4;
    for (var i = 0; i < items.length; i++) {
      var it = items[i], placed = -1;
      for (var r = 0; r < rowEnds.length; r++) {
        if (it.start > rowEnds[r] + padBp) { placed = r; break; }
      }
      if (placed < 0) {
        if (rowEnds.length >= MAX_ROWS) continue;
        rowEnds.push(it.end); rows.push([it]);
      } else { rowEnds[placed] = it.end; rows[placed].push(it); }
    }
    return rows;
  }

  function formatBp(v) {
    if (Math.abs(v) < 1000) return String(Math.round(v));
    var kb = v / 1000;
    return (Math.abs(kb) < 10 ? kb.toFixed(1) : Math.round(kb)) + 'kb';
  }

  function niceStep(raw) {
    var p = Math.pow(10, Math.floor(Math.log10(Math.max(raw, 1))));
    var n = raw / p;
    return (n <= 1 ? 1 : n <= 2 ? 2 : n <= 5 ? 5 : 10) * p;
  }

  function cssWidth(canvas) {
    var w = canvas.parentNode.clientWidth;
    if (!w && canvas.parentNode.parentNode) w = canvas.parentNode.parentNode.clientWidth;
    return w || 800;
  }

  function setupCanvas(canvas, w, h) {
    var dpr = window.devicePixelRatio || 1;
    canvas.style.width = w + 'px';
    canvas.style.height = h + 'px';
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  // ---------------------------------------------------------------- drawing

  function drawRuler() {
    var canvas = state.ruler;
    if (!canvas) return;
    // the reads canvas lives inside the scrolling box, so it is narrower than the
    // ruler by the scrollbar; matching it keeps the ticks over the right bases
    var w = (state.canvas && state.scrollBox && state.scrollBox.clientWidth)
      ? state.scrollBox.clientWidth : cssWidth(canvas);
    var ctx = setupCanvas(canvas, w, RULER_H);
    ctx.clearRect(0, 0, w, RULER_H);
    ctx.fillStyle = '#fcfcfb';
    ctx.fillRect(0, 0, w, RULER_H);

    var span = state.viewEnd - state.viewStart;
    var bpPerPx = span / w;

    ctx.fillStyle = 'rgba(42,120,214,0.11)';
    state.spacers.forEach(function (sp) {
      var x0 = (sp.start - state.viewStart) / bpPerPx;
      var x1 = (sp.end - state.viewStart) / bpPerPx;
      ctx.fillRect(x0, 26, Math.max(1, x1 - x0), RULER_H - 26);
    });

    var step = niceStep(span / 8);
    ctx.strokeStyle = '#e1e0d9';
    ctx.fillStyle = '#5f6368';
    ctx.font = '10px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.lineWidth = 1;
    for (var p = Math.ceil(state.viewStart / step) * step; p <= state.viewEnd; p += step) {
      var x = Math.round((p - state.viewStart) / bpPerPx) + 0.5;
      ctx.textAlign = x < 18 ? 'left' : (x > w - 18 ? 'right' : 'center');
      ctx.fillText(formatBp(p), Math.min(Math.max(x, 1), w - 1), 20);
      ctx.textAlign = 'center';
      ctx.beginPath(); ctx.moveTo(x, 26); ctx.lineTo(x, RULER_H); ctx.stroke();
    }

    ctx.save();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = '#d03b3b';
    ctx.fillStyle = '#d03b3b';
    ctx.lineWidth = 1.5;
    // labels of marks that sit close together are staggered, because the PA
    // allele has four of them within a few dozen bases
    var lastX = -1e9, row = 0;
    (state.marks || []).forEach(function (m) {
      var cx = Math.round((m.pos - state.viewStart) / bpPerPx) + 0.5;
      ctx.beginPath(); ctx.moveTo(cx, 26); ctx.lineTo(cx, RULER_H); ctx.stroke();
      row = (cx - lastX < 56) ? (row + 1) % 2 : 0;
      lastX = cx;
      ctx.textAlign = cx < 30 ? 'left' : (cx > w - 30 ? 'right' : 'center');
      ctx.fillText(m.label, Math.min(Math.max(cx, 1), w - 1), 36 + row * 9);
      ctx.textAlign = 'center';
    });
    ctx.restore();

    ctx.strokeStyle = '#c3c2b7';
    ctx.beginPath(); ctx.moveTo(0, RULER_H - 0.5); ctx.lineTo(w, RULER_H - 0.5); ctx.stroke();
  }

  function drawReads() {
    var canvas = state.canvas;
    if (!canvas) return;
    var w = cssWidth(canvas);
    var bpPerPx = (state.viewEnd - state.viewStart) / w;
    var x = function (p) { return (p - state.viewStart) / bpPerPx; };

    state.rows = pack(state.items, bpPerPx);

    // The stack can be tens of thousands of rows deep, so only the slice the box
    // is showing is drawn. The spacer carries the full height and gives the
    // scrollbar its range; the canvas sticks to the top of the viewport and is
    // repainted as it scrolls, so nothing is ever left undrawn.
    var total = state.rows.length * ROW_H + 16;
    if (state.scrollSpacer) state.scrollSpacer.style.height = total + 'px';
    var boxH = state.scrollBox ? state.scrollBox.clientHeight : 460;
    var h = Math.max(60, Math.min(boxH, total));
    var ctx = setupCanvas(canvas, w, h);

    var scrollTop = state.scrollBox ? state.scrollBox.scrollTop : 0;
    var firstRow = Math.max(0, Math.floor((scrollTop - 8) / ROW_H));
    var lastRow = Math.min(state.rows.length, firstRow + Math.ceil(h / ROW_H) + 1);

    state.hits = [];
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = '#fcfcfb';
    ctx.fillRect(0, 0, w, h);

    ctx.save();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = 'rgba(208,59,59,0.45)';
    ctx.lineWidth = 1;
    state.cuts.forEach(function (c) {
      var cx = Math.round(x(c)) + 0.5;
      ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, h); ctx.stroke();
    });
    ctx.restore();

    for (var r = firstRow; r < lastRow; r++) {
      var y = 8 + r * ROW_H - scrollTop;
      if (y + READ_H < 0 || y > h) continue;
      var row = state.rows[r];
      for (var i = 0; i < row.length; i++) {
        var it = row[i];
        if (it.end < state.viewStart || it.start > state.viewEnd) continue;
        var colour = colourOf(it);

        // deletions: a join, never a blank, so the two sides stay one read
        for (var d = 0; d < it.dels.length; d++) {
          var g = it.dels[d];
          var gx0 = x(g.s), gx1 = x(g.e);
          var big = g.len > state.largeDel;
          var igv = state.colourBy === 'igv';
          var joinColour = igv ? '#4a4a46' : (big ? '#d03b3b' : colour);
          if (!igv) {
            ctx.save();
            ctx.globalAlpha = 0.18;
            ctx.fillStyle = joinColour;
            ctx.fillRect(gx0, y + 1, Math.max(1, gx1 - gx0), READ_H - 2);
            ctx.restore();
          }
          ctx.strokeStyle = joinColour;
          ctx.lineWidth = big ? 2.2 : 1.8;
          ctx.beginPath();
          ctx.moveTo(gx0, y + READ_H / 2); ctx.lineTo(gx1, y + READ_H / 2);
          ctx.stroke();
          ctx.lineWidth = 1.4;
          ctx.beginPath();
          ctx.moveTo(gx0 + 0.5, y); ctx.lineTo(gx0 + 0.5, y + READ_H);
          ctx.moveTo(gx1 - 0.5, y); ctx.lineTo(gx1 - 0.5, y + READ_H);
          ctx.stroke();
          if (gx1 - gx0 > 24) {
            // IGV puts the length in the gap itself, on the join line
            var label = formatBp(g.len);
            ctx.font = '9px system-ui, sans-serif';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            var lw = ctx.measureText(label).width + 4;
            var lx = (gx0 + gx1) / 2;
            ctx.fillStyle = '#fcfcfb';
            ctx.fillRect(lx - lw / 2, y + 1, lw, READ_H - 2);
            ctx.fillStyle = igv ? '#3f3e38' : (big ? '#d03b3b' : '#49483f');
            ctx.fillText(label, lx, y + READ_H / 2);
            ctx.textBaseline = 'alphabetic';
          }
        }

        ctx.fillStyle = colour;
        for (var b = 0; b < it.blocks.length; b++) {
          var bx0 = x(it.blocks[b].s), bx1 = x(it.blocks[b].e);
          if (bx1 < 0 || bx0 > w) continue;
          ctx.fillRect(bx0, y, Math.max(1, bx1 - bx0), READ_H);
        }

        for (var v = 0; v < it.inversions.length; v++) {
          var ix0 = x(it.inversions[v].s), ix1 = x(it.inversions[v].e);
          var iw = Math.max(1, ix1 - ix0);
          ctx.fillStyle = INV_COLOUR;
          ctx.fillRect(ix0, y, iw, READ_H);
          ctx.fillStyle = '#ffffff';
          for (var hx = ix0 + 2; hx < ix1 - 1; hx += 4) ctx.fillRect(hx, y + 1, 1, READ_H - 2);
          if (iw > 26) {
            ctx.fillStyle = '#2f6b1f';
            ctx.font = '8px system-ui, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText('inv', (ix0 + ix1) / 2, y - 1);
          }
        }

        for (var sIdx = 0; sIdx < it.subs.length; sIdx++) {
          var sx = x(it.subs[sIdx].at);
          if (sx < -2 || sx > w + 2) continue;
          ctx.fillStyle = BASE_COLOUR[it.subs[sIdx].alt] || BASE_COLOUR.N;
          ctx.fillRect(sx - 0.75, y, Math.max(1.5, 1 / bpPerPx), READ_H);
        }

        for (var n = 0; n < it.ins.length; n++) {
          var ins = it.ins[n];
          var nx = x(ins.at);
          if (nx < -2 || nx > w + 2) continue;
          var bigIns = ins.len > state.largeIns;
          ctx.fillStyle = '#4a3aa7';
          ctx.fillRect(nx - (bigIns ? 2 : 1.25), y - 1, bigIns ? 4 : 2.5, READ_H + 2);
          // a large insertion is worth a number, the way a deletion gets one
          if (bigIns) {
            var ilab = formatBp(ins.len);
            ctx.font = '9px system-ui, sans-serif';
            ctx.textAlign = 'left';
            var iw = ctx.measureText(ilab).width + 4;
            ctx.fillStyle = 'rgba(252,252,251,0.85)';
            ctx.fillRect(nx + 3, y, iw, READ_H);
            ctx.fillStyle = '#4a3aa7';
            ctx.textBaseline = 'middle';
            ctx.fillText(ilab, nx + 5, y + READ_H / 2);
            ctx.textBaseline = 'alphabetic';
            ctx.textAlign = 'center';
          }
          // click targets, in canvas coordinates, for the rows actually drawn
          state.hits.push({ x: nx, y: y, h: READ_H, read: it.name,
                            at: ins.at, len: ins.len, seq: ins.seq || '' });
        }
      }
    }

    var note = document.getElementById('alignViewNote');
    if (note) {
      var sel = document.getElementById('alignViewPick');
      var picked = state.pick && sel ? sel.options[sel.selectedIndex].textContent : '';
      note.textContent = state.items.length.toLocaleString() + ' reads in '
        + state.rows.length.toLocaleString() + ' rows'
        + (state.pick ? ' \u00b7 showing ' + picked.replace(/ \(\d[\d,]*\)$/, '')
            + ' only, of ' + state.total.toLocaleString() + ' reads' : '')
        + (state.filtered
          ? ' \u00b7 rebuilt from the filtered mutations, as in the IGV track'
          : ' \u00b7 raw aligner output, every record as minimap2 left it')
        + ' \u00b7 view ' + formatBp(state.viewStart) + '\u2013' + formatBp(state.viewEnd);
    }
  }

  function redraw() { hideInsertion(); drawRuler(); drawReads(); }

  function rebuild() {
    var all = prepare(
      (state.filtered && state.mutations && Object.keys(state.mutations).length)
        ? buildFromMutations(state.spans, state.mutations)
        : buildFromSam(state.lines)
    );
    fillPicker(all);
    state.total = all.length;
    state.items = state.pick
      ? all.filter(function (it) { return matchesPick(it, state.pick); })
      : all;
  }

  function observe() {
    if (state.observer || typeof ResizeObserver === 'undefined') return;
    state.observer = new ResizeObserver(function () {
      var w = state.scrollBox ? state.scrollBox.clientWidth : 0;
      if (w && w !== state.lastWidth) { state.lastWidth = w; redraw(); }
    });
    state.observer.observe(state.scrollBox);
  }

  function zoom(factor, anchorFrac) {
    var span = state.viewEnd - state.viewStart;
    var next = Math.max(60, Math.min(state.refLen, span * factor));
    var frac = anchorFrac === undefined ? 0.5 : anchorFrac;
    var start = state.viewStart + span * frac - next * frac;
    start = Math.max(0, Math.min(state.refLen - next, start));
    state.viewStart = start;
    state.viewEnd = start + next;
    redraw();
  }

  // Clicking an insertion opens what was actually inserted: the tick alone says
  // only that something went in there.
  function insertionPopup() {
    var box = document.getElementById('alignViewInsPopup');
    if (box) return box;
    box = document.createElement('div');
    box.id = 'alignViewInsPopup';
    box.style.cssText = 'position:fixed;z-index:1080;max-width:min(560px,90vw);display:none;'
      + 'background:#fff;border:1px solid #d7d6d0;border-radius:8px;'
      + 'box-shadow:0 8px 24px rgba(0,0,0,.14);padding:10px 12px;font-size:12px;';
    document.body.appendChild(box);
    return box;
  }

  function hideInsertion() {
    var box = document.getElementById('alignViewInsPopup');
    if (box) box.style.display = 'none';
  }

  function showInsertion(hit, clientX, clientY) {
    var box = insertionPopup();
    var seq = hit.seq || '';
    var pretty = seq ? seq.replace(/(.{60})/g, '$1\n') : '';
    box.innerHTML = '';
    var head = document.createElement('div');
    head.className = 'd-flex align-items-center justify-content-between gap-2 mb-1';
    head.innerHTML = '<strong>' + formatBp(hit.len) + ' inserted at ' + hit.at + '</strong>';
    var close = document.createElement('button');
    close.type = 'button';
    close.className = 'btn-close btn-sm';
    close.setAttribute('aria-label', 'Close');
    close.onclick = hideInsertion;
    head.appendChild(close);
    box.appendChild(head);

    var sub = document.createElement('div');
    sub.className = 'text-secondary mb-1';
    sub.textContent = 'read ' + hit.read;
    box.appendChild(sub);

    if (seq) {
      var pre = document.createElement('pre');
      pre.style.cssText = 'white-space:pre-wrap;word-break:break-all;margin:0 0 6px;'
        + 'font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;'
        + 'font-size:11px;max-height:180px;overflow:auto;';
      pre.textContent = pretty;
      box.appendChild(pre);
      var copy = document.createElement('button');
      copy.type = 'button';
      copy.className = 'btn btn-sm btn-outline-secondary';
      copy.textContent = 'Copy sequence';
      copy.onclick = function () {
        navigator.clipboard.writeText(seq);
        copy.textContent = 'Copied';
        setTimeout(function () { copy.textContent = 'Copy sequence'; }, 1200);
      };
      box.appendChild(copy);
    } else {
      var none = document.createElement('div');
      none.className = 'text-secondary';
      none.textContent = 'The inserted bases are not in this view\u2019s source.';
      box.appendChild(none);
    }

    box.style.display = '';
    var r = box.getBoundingClientRect();
    var left = Math.min(Math.max(8, clientX - r.width / 2), window.innerWidth - r.width - 8);
    var top = clientY + 14;
    if (top + r.height > window.innerHeight - 8) top = Math.max(8, clientY - r.height - 14);
    box.style.left = left + 'px';
    box.style.top = top + 'px';
  }

  function hitAt(canvas, e) {
    var rect = canvas.getBoundingClientRect();
    var px = e.clientX - rect.left, py = e.clientY - rect.top;
    var best = null, bestDist = 7;
    for (var i = 0; i < state.hits.length; i++) {
      var h = state.hits[i];
      if (py < h.y - 2 || py > h.y + h.h + 2) continue;
      var d = Math.abs(px - h.x);
      if (d <= bestDist) { bestDist = d; best = h; }
    }
    return best;
  }

  function wire() {
    var canvas = state.canvas;
    canvas.addEventListener('wheel', function (e) {
      if (!e.ctrlKey && !e.shiftKey) return;
      e.preventDefault();
      var rect = canvas.getBoundingClientRect();
      zoom(e.deltaY > 0 ? 1.25 : 0.8, (e.clientX - rect.left) / rect.width);
    }, { passive: false });

    canvas.addEventListener('mousedown', function (e) {
      state.drag = { x: e.clientX, start: state.viewStart, end: state.viewEnd, moved: false };
    });

    canvas.addEventListener('click', function (e) {
      if (state.drag && state.drag.moved) return;   // that was a pan, not a click
      var hit = hitAt(canvas, e);
      if (hit) showInsertion(hit, e.clientX, e.clientY);
      else hideInsertion();
    });

    canvas.addEventListener('mousemove', function (e) {
      if (state.drag) return;
      canvas.style.cursor = hitAt(canvas, e) ? 'pointer' : 'grab';
    });
    window.addEventListener('mousemove', function (e) {
      if (!state.drag) return;
      if (Math.abs(e.clientX - state.drag.x) > 3) state.drag.moved = true;
      var bpPerPx = (state.drag.end - state.drag.start) / cssWidth(canvas);
      var span = state.drag.end - state.drag.start;
      var start = Math.max(0, Math.min(state.refLen - span,
        state.drag.start + (state.drag.x - e.clientX) * bpPerPx));
      state.viewStart = start; state.viewEnd = start + span;
      redraw();
    });
    window.addEventListener('mouseup', function () {
      // the click handler runs after mouseup, so the pan flag has to outlive it
      setTimeout(function () { state.drag = null; }, 0);
    });

    var get = function (id) { return document.getElementById(id); };
    if (get('alignViewZoomIn')) get('alignViewZoomIn').onclick = function () { zoom(0.8); };
    if (get('alignViewZoomOut')) get('alignViewZoomOut').onclick = function () { zoom(1.25); };
    if (get('alignViewReset')) get('alignViewReset').onclick = function () {
      state.viewStart = 0; state.viewEnd = state.refLen; redraw();
    };
    if (get('alignViewColour')) get('alignViewColour').onchange = function () {
      state.colourBy = this.value; redraw();
    };
    if (get('alignViewHeight')) get('alignViewHeight').onchange = function () {
      state.scrollBox.style.maxHeight = this.value + 'px'; redraw();
    };
    if (get('alignViewFilter')) get('alignViewFilter').onchange = function () {
      state.filtered = this.checked; rebuild(); redraw();
    };
    if (get('alignViewPick')) get('alignViewPick').onchange = function () {  // 20261111_alignview
      state.pick = this.value;
      rebuild();
      if (state.scrollBox) state.scrollBox.scrollTop = 0;
      redraw();
    };
    if (get('alignViewRef')) get('alignViewRef').onchange = function () {  // 20261111_alignview
      window.renderAlignmentView({ refKey: this.value });
    };
    if (state.scrollBox) {
      state.scrollBox.addEventListener('scroll', function () {
        if (state.scrollPending) return;
        state.scrollPending = true;
        requestAnimationFrame(function () { state.scrollPending = false; drawReads(); });
      });
    }
    window.addEventListener('resize', redraw);
  }

  // 20261111_alignview: PA mode aligns the reads against the precise-PA allele as
  // well, and that view is the one that shows whether the donor sits in cleanly,
  // so the viewer can be pointed at either reference.
  function refChoices() {
    var out = [{ key: 'WT', label: 'Wild-type reference ('
      + ((window.reference || '').length) + ' bp)' }];
    var extra = window.paAlignRefs || {};
    Object.keys(extra).forEach(function (k) {
      out.push({ key: k, label: extra[k].label || k });
    });
    return out;
  }

  function refSource(key) {
    var extra = window.paAlignRefs || {};
    if (key && key !== 'WT' && extra[key]) {
      var r = extra[key];
      return {
        sam: r.sam, mutations: r.mutations, reference: r.seq,
        marks: r.marks || (r.cuts || []).map(function (p) { return { pos: p, label: 'cut' }; })
      };
    }
    return {
      sam: window.treatedAlignSam, mutations: window.lungoReadMutations,
      reference: window.reference,
      marks: [window.cleavagePos, window.cleavagePos2]
        .filter(function (c) { return typeof c === 'number'; })
        .map(function (p) { return { pos: p, label: 'cut ' + p }; })
    };
  }

  function fillRefSelector() {
    var sel = document.getElementById('alignViewRef');
    if (!sel) return;
    var choices = refChoices();
    var wrap = document.getElementById('alignViewRefWrap');
    if (wrap) wrap.style.display = choices.length > 1 ? '' : 'none';
    var keep = state.refKey;
    sel.innerHTML = choices.map(function (c) {
      return '<option value="' + c.key + '">' + c.label + '</option>';
    }).join('');
    sel.value = choices.some(function (c) { return c.key === keep; }) ? keep : 'WT';
    state.refKey = sel.value;
  }

  window.renderAlignmentView = function renderAlignmentView(options) {
    options = options || {};
    var section = document.getElementById('alignViewSection');
    if (!section) return;
    if (options.refKey) state.refKey = options.refKey;
    fillRefSelector();
    var src = refSource(state.refKey);
    var lines = options.sam || src.sam;
    if (!lines || !lines.length) { section.style.display = 'none'; return; }
    if (typeof lines === 'string') lines = lines.split('\n');

    state.lines = lines;
    state.spans = spansFromSam(lines);
    state.mutations = options.mutations || src.mutations || {};
    state.largeDel = Number(options.largeDel || window.largeDelLen || 50);
    state.largeIns = Number(options.largeIns || window.largeInsLen || 50);

    var filterBox = document.getElementById('alignViewFilter');
    var haveTable = Object.keys(state.mutations).length > 0;
    if (filterBox) {
      filterBox.disabled = !haveTable;
      filterBox.checked = haveTable;          // a fresh run always opens as IGV shows it
      filterBox.title = haveTable ? '' : 'No filtered mutation table for this run';
    }
    state.filtered = haveTable;

    var ref = options.reference || src.reference || '';
    state.refLen = ref.length || Object.keys(state.spans).reduce(function (m, k) {
      return Math.max(m, state.spans[k].end);
    }, 1000);
    state.viewStart = 0;
    state.viewEnd = state.refLen;

    state.marks = (options.marks || src.marks || []).filter(function (m) {
      return m && typeof m.pos === 'number' && m.pos >= 0;
    }).sort(function (a, b) { return a.pos - b.pos; });
    state.cuts = state.marks.map(function (m) { return m.pos; });

    state.spacers = [];
    [window.target, window.target2].forEach(function (t) {
      if (!t) return;
      var at = ref.indexOf(t);
      if (at < 0) {
        var rc = t.split('').reverse().map(function (b) {
          return { A: 'T', T: 'A', C: 'G', G: 'C' }[b] || b;
        }).join('');
        at = ref.indexOf(rc);
      }
      if (at >= 0) state.spacers.push({ start: at, end: at + t.length });
    });

    section.style.display = '';
    state.ruler = document.getElementById('alignViewRuler');
    state.canvas = document.getElementById('alignViewCanvas');
    state.scrollBox = document.getElementById('alignViewScroll');
    state.scrollSpacer = document.getElementById('alignViewSpacer');
    var pickBox = document.getElementById('alignViewPick');
    if (pickBox) state.pick = pickBox.value || '';
    var colourBox = document.getElementById('alignViewColour');
    if (colourBox) state.colourBy = colourBox.value;
    if (!state.canvas || !state.scrollBox) return;
    rebuild();
    if (state.scrollBox.scrollTop) state.scrollBox.scrollTop = 0;
    if (!state.wired) { wire(); state.wired = true; }
    observe();
    redraw();
    if (!state.scrollBox.clientWidth) {
      requestAnimationFrame(function () { if (state.scrollBox.clientWidth) redraw(); });
    }
  };
})();
