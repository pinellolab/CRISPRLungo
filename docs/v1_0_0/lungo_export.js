/* ============================================================================
 * lungo_export.js
 * ----------------------------------------------------------------------------
 * Builds CLI-style result files from a finished web analysis and downloads
 * them as one ZIP laid out like a CRISPRlungo CLI output directory:
 *
 *   results/read_classification.txt      results/mutation_patter_count.txt
 *   results/mutation_summary_count.txt   results/preprocess_count.txt
 *   results/input_summary.txt            results/Allele_table.txt
 *   results/mutation_pattern_p_values.txt (control runs)
 *   align/Treated_alignment.sam          align/Control_alignment.sam
 *   ref_seq/ref_wo_umi.fasta             plots/*.png, plots/allele_plot.svg
 *
 * Formats follow src/CRISPRlungo_mutation_analysis.py (classify_mut_mild,
 * process_mutations) and src/CRISPRlungo_visualization.py (write_read_count,
 * allele_table) so the folder can be read by CRISPRlungoAllele.
 *
 * Depends on globals from main.js / CRISPRlungo_ai.js (mainResult,
 * filteredSigMut, reference, cleavagePos, ...) and on JSZip.
 * ==========================================================================*/

var LUNGO_EXPORT_REF_NAME = 'Ref';
var LUNGO_PARTIAL_INDUCE_CUTOFF = 0.8;
var LX_RANGE_ALIGN_END = 100;

// src/possible_insertion.fasta, as CRISPRlungo_insert_analysis.confirm_insertion_seq reads it
// (names keep the '>' and sequences keep their case, as in the CLI).
var LX_POSSIBLE_INSERTION = [
  ['>nCas9', 'gacaagaagtacagcatcggcctggacatcggcaccaactctgtgggctgggccgtgatcaccgacgagtacaaggtgcccagcaagaaattcaaggtgctgggcaacaccgaccggcacagcatcaagaagaacctgatcggagccctgctgttcgacagcggcgaaacagccgaggccacccggctgaagagaaccgccagaagaagatacaccagacggaagaaccggatctgctatctgcaagagatcttcagcaacgagatggccaaggtggacgacagcttcttccacagactggaagagtccttcctggtggaagaggataagaagcacgagcggcaccccatcttcggcaacatcgtggacgaggtggcctaccacgagaagtaccccaccatctaccacctgagaaagaaactggtggacagcaccgacaaggccgacctgcggctgatctatctggccctggcccacatgatcaagttccggggccacttcctgatcgagggcgacctgaaccccgacaacagcgacgtggacaagctgttcatccagctggtgcagacctacaaccagctgttcgaggaaaaccccatcaacgccagcggcgtggacgccaaggccatcctgtctgccagactgagcaagagcagaaagctggaaaatctgatcgcccagctgcccggcgagaagaagaatggcctgttcggaaacctgattgccctgagcctgggcctgacccccaacttcaagagcaacttcgacctggccgaggatgccaaactgcagctgagcaaggacacctacgacgacgacctggacaacctgctggcccagatcggcgaccagtacgccgacctgtttctggccgccaagaacctgtccgacgccatcctgctgagcgacatcctgagagtgaacaccgagatcaccaaggcccccctgagcgcctctatgatcaagagatacgacgagcaccaccaggacctgaccctgctgaaagctctcgtgcggcagcagctgcctgagaagtacaaagagattttcttcgaccagagcaagaacggctacgccggctacattgacggcggagccagccaggaagagttctacaagttcatcaagcccatcctggaaaagatggacggcaccgaggaactgctcgtgaagctgaagagagaggacctgctgcggaagcagcggaccttcgacaacggcagcatcccccaccagatccacctgggagagctgcacgccattctgcggcggcaggaagatttttacccattcctgaaggacaaccgggaaaagatcgagaagatcctgaccttccgcatcccctactacgtgggccctctggccaggggaaacagcagattcgcctggatgaccagaaagagcgaggaaaccatcaccccctggaacttcgaggaagtggtggacaagggcgcttccgcccagagcttcatcgagcggatgaccaacttcgataagaacctgcccaacgagaaggtgctgcccaagcacagcctgctgtacgagtacttcaccgtgtataacgagctgaccaaagtgaaatacgtgaccgagggaatgagaaagcccgccttcctgagcggcgagcagaaaaaggccatcgtggacctgctgttcaagaccaaccggaaagtgaccgtgaagcagctgaaagaggactacttcaagaaaatcgagtgcttcgactccgtggaaatctccggcgtggaagatcggttcaacgcctccctgggcacataccacgatctgctgaaaattatcaaggacaaggacttcctggacaatgaggaaaacgaggacattctggaagatatcgtgctgaccctgacactgtttgaggacagagagatgatcgaggaacggctgaaaacctatgcccacctgttcgacgacaaagtgatgaagcagctgaagcggcggagatacaccggctggggcaggctgagccggaagctgatcaacggcatccgggacaagcagtccggcaagacaatcctggatttcctgaagtccgacggcttcgccaacagaaacttcatgcagctgatccacgacgacagcctgacctttaaagaggacatccagaaagcccaggtgtccggccagggcgatagcctgcacgagcacattgccaatctggccggcagccccgccattaagaagggcatcctgcagacagtgaaggtggtggacgagctcgtgaaagtgatgggccggcacaagcccgagaacatcgtgatcgaaatggccagagagaaccagaccacccagaagggacagaagaacagccgcgagagaatgaagcggatcgaagagggcatcaaagagctgggcagccagatcctgaaagaacaccccgtggaaaacacccagctgcagaacgagaagctgtacctgtactacctgcagaatgggcgggatatgtacgtggaccaggaactggacatcaaccggctgtccgactacgatgtggacgctatcgtgcctcagagctttctgaaggacgactccatcgacaacaaggtgctgaccagaagcgacaagaaccggggcaagagcgacaacgtgccctccgaagaggtcgtgaagaagatgaagaactactggcggcagctgctgaacgccaagctgattacccagagaaagttcgacaatctgaccaaggccgagagaggcggcctgagcgaactggataaggccggcttcatcaagagacagctggtggaaacccggcagatcacaaagcacgtggcacagatcctggactcccggatgaacactaagtacgacgagaatgacaagctgatccgggaagtgaaagtgatcaccctgaagtccaagctggtgtccgatttccggaaggatttccagttttacaaagtgcgcgagatcaacaactaccaccacgcccacgacgcctacctgaacgccgtcgtgggaaccgccctgatcaaaaagtaccctaagctggaaagcgagttcgtgtacggcgactacaaggtgtacgacgtgcggaagatgatcgccaagagcgagcaggaaatcggcaaggctaccgccaagtacttcttctacagcaacatcatgaactttttcaagaccgagattaccctggccaacggcgagatccggaagcggcctctgatcgagacaaacggcgaaaccggggagatcgtgtgggataagggccgggattttgccaccgtgcggaaagtgctgagcatgccccaagtgaatatcgtgaaaaagaccgaggtgcagacaggcggcttcagcaaagagtctatcctgcccaagaggaacagcgataagctgatcgccagaaagaaggactgggaccctaagaagtacggcggcttcgacagccccaccgtggcctattctgtgctggtggtggccaaagtggaaaagggcaagtccaagaaactgaagagtgtgaaagagctgctggggatcaccatcatggaaagaagcagcttcgagaagaatcccatcgactttctggaagccaagggctacaaagaagtgaaaaaggacctgatcatcaagctgcctaagtactccctgttcgagctggaaaacggccggaagagaatgctggcctctgccggcgaactgcagaagggaaacgaactggccctgccctccaaatatgtgaacttcctgtacctggccagccactatgagaagctgaagggctcccccgaggataatgagcagaaacagctgtttgtggaacagcacaagcactacctggacgagatcatcgagcagatcagcgagttctccaagagagtgatcctggccgacgctaatctggacaaagtgctgtccgcctacaacaagcaccgggataagcccatcagagagcaggccgagaatatcatccacctgtttaccctgaccaatctgggagcccctgccgccttcaagtactttgacaccaccatcgaccggaagaggtacaccagcaccaaagaggtgctggacgccaccctgatccaccagagcatcaccggcctgtacgagacacggatcgacctgtctcagctgggaggtgac'],
  ['>Scaffold', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'],
  ['>DC_cluister', 'cagagccagtgggctagcaccagaggcggacacaaggaacggcctaccgtgatgctcttcccattccgc'],
  ['>ITR', 'aggaacccctagtgatggagttggccactccctctctgcgcgctcgctcgctcactgaggccgggcgaccaaaggtcgcccgacgcccgggctttgcccgggcggcctcagtgagcgagcgagcgcgcag'],
  ['>ITR_2', 'aggaacccctagtgatggagttggccactccctctctgcgcgctcgctcgctcactgaggccgggcgaccaaaggtcgcccgacgcccgggcggcctcagtgagcgagcgagcgcgcagctgcctgcagg']
];

function lxRefName() {
  var el = document.getElementById('exportRefName');
  var v = el ? el.value.trim().replace(/^>/, '').split(/\s+/)[0] : '';
  return v || LUNGO_EXPORT_REF_NAME;
}

function lxNum(id, fallback) {
  var el = document.getElementById(id);
  if (!el) return fallback;
  var v = Number(el.value);
  return isNaN(v) ? fallback : v;
}

function lxRound2(x) {
  return Math.round(x * 100) / 100;
}

// Python-like str() for a float rounded to 2 decimals: 42.76, 0.5 -> 0.5, 100 -> 100.0
function lxPyFloat(x) {
  var r = lxRound2(x);
  var s = String(r);
  if (s.indexOf('.') === -1 && s.indexOf('e') === -1) s += '.0';
  return s;
}

// "1157_1158:Ins_AA" / "..:Del_4" / "..:Sub_T>C" (web, optional _inv, optional Large prefix)
function lxParseWebMut(m) {
  var inv = false;
  if (m.slice(-4) === '_inv') { inv = true; m = m.slice(0, -4); }
  var c = m.indexOf(':');
  var pos = m.slice(0, c).split('_');
  var rest = m.slice(c + 1);
  var u = rest.indexOf('_');
  var type = rest.slice(0, u).replace('Large', '');
  var detail = rest.slice(u + 1);
  var o = { start: Number(pos[0]), end: Number(pos[1]), type: type, inv: inv, raw: m };
  if (type === 'Del') { o.len = Number(detail); }
  else if (type === 'Ins') { o.seq = detail; o.len = detail.length; }
  else if (type === 'Sub') { var fr = detail.split('>'); o.from = fr[0]; o.to = fr[1]; o.len = fr[0].length; }
  return o;
}

// CLI token, as written by classify_mut_mild
function lxCliToken(p) {
  var span = p.start + '_' + p.end + ':';
  if (p.type === 'Del') return span + (p.len <= largeDelLen ? 'Del' : 'LargeDel') + '_' + p.len;
  if (p.type === 'Ins') return span + (p.len <= largeInsLen ? 'Ins' : 'LargeIns') + '_' + p.len + '_' + p.seq;
  return span + 'Sub_' + p.len + '_' + p.from + '>' + p.to;
}

function lxInducedKey(p) {
  if (p.type === 'Del') return p.start + '_' + p.end + ':Del_' + p.len;
  if (p.type === 'Ins') return p.start + '_' + p.end + ':Ins_' + p.seq;
  return p.start + '_' + p.end + ':Sub_' + p.from + '>' + p.to;
}

// classify_mut_mild() on parsed web mutations
function lxClassify(muts, insertInfo, inducedKeys) {
  var hasInduced = inducedKeys.size > 0;
  if (muts.length === 0) {
    return { cls: 'WT', info: 'None', integ: '', induce: hasInduced ? 'X' : '-' };
  }
  var kinds = [];
  var tokens = [];
  var matched = new Set();
  var nonInduced = 0;
  muts.forEach(function (p) {
    var k = lxInducedKey(p);
    if (inducedKeys.has(k)) { matched.add(k); } else { nonInduced += 1; }
    if (p.type === 'Del') kinds.push(p.len <= largeDelLen ? 'Del' : 'LargeDel');
    else if (p.type === 'Ins') kinds.push(p.len <= largeInsLen ? 'Ins' : 'LargeIns');
    else kinds.push('Sub');
    tokens.push(lxCliToken(p));
  });
  var induce = '-';
  if (hasInduced) {
    if (matched.size / inducedKeys.size > LUNGO_PARTIAL_INDUCE_CUTOFF) {
      if (matched.size === inducedKeys.size) induce = nonInduced === 0 ? 'Precise' : 'Precise_with_mutations';
      else induce = nonInduced === 0 ? 'Partial' : 'Partial_with_mutations';
    } else {
      induce = 'X';
    }
  }
  var has = function (k) { return kinds.indexOf(k) !== -1; };
  var cls;
  if (insertInfo.indexOf('inversion') !== -1) cls = 'Inv';
  else if (has('LargeDel') || has('LargeIns')) {
    cls = (has('LargeDel') && has('LargeIns')) ? 'ComplexWithLargeMut' : (has('LargeDel') ? 'LargeDel' : 'LargeIns');
  }
  else if (has('Del') && has('Ins')) cls = 'Complex';
  else if (has('Del')) cls = 'Del';
  else if (has('Ins')) cls = 'Ins';
  else if (has('Sub')) cls = 'Sub';
  else cls = 'Complex';
  return { cls: cls, info: tokens.join(','), integ: insertInfo, induce: induce };
}

/* ---------------- CLI-identical insertion annotation ---------------- */

function lxParseCigar(c) {
  var ops = [], re = /(\d+)([MIDNSHP=X])/g, m;
  while ((m = re.exec(c))) ops.push(['MIDNSHP=X'.indexOf(m[2]), Number(m[1])]);
  return ops;
}

// Alignment records per read, in file order (= FASTQ order).
function lxSamByRead(lines) {
  var map = new Map();
  (lines || []).forEach(function (line) {
    if (!line || line[0] === '@') return;
    var f = line.split('\t');
    var flag = Number(f[1]);
    if ((flag & 4) || (flag & 0x100)) return;
    var name = f[0];
    var cut = name.indexOf('_AI_');
    var base = cut === -1 ? name : name.slice(0, cut);
    var strands = cut === -1 ? [] : name.slice(cut + 4).split('_');
    var prod = 1;
    strands.slice(1).forEach(function (x) { prod *= Number(x); });
    var ops = lxParseCigar(f[5]);
    var pos = Number(f[3]) - 1;
    var refEnd = pos;
    ops.forEach(function (o) { if (o[0] === 0 || o[0] === 2 || o[0] === 3 || o[0] === 7 || o[0] === 8) refEnd += o[1]; });
    if (!map.has(base)) map.set(base, []);
    map.get(base).push({ flag: flag, pos: pos, refEnd: refEnd, ops: ops, seqLen: (f[9] && f[9] !== '*') ? f[9].length : 0, toSeq: prod });
  });
  return map;
}

// Insertions with their read coordinates (seq_pos), as quant_unique_indels /
// analyze_SA_reads in CRISPRlungo_mutation_analysis.py produce them.
function lxCliInsertions(recs, refLen) {
  var out = [];
  var prim = recs.filter(function (r) { return !(r.flag & 0x800); })[0];
  if (!prim) return out;
  var sups = recs.filter(function (r) { return r.flag & 0x800; });
  if (!sups.length) {
    var rp = prim.pos, qp = 0;
    prim.ops.forEach(function (o) {
      if (o[0] === 1) out.push({ pos: rp, len: o[1], qs: qp, qe: qp + o[1] });
      if (o[0] === 0 || o[0] === 2 || o[0] === 3) rp += o[1];
      if (o[0] === 0 || o[0] === 1 || o[0] === 4) qp += o[1];
    });
    return out;
  }
  var qlen = prim.seqLen, R = LX_RANGE_ALIGN_END;
  var segs = [];
  [prim].concat(sups).forEach(function (r) {
    var first = r.ops[0], last = r.ops[r.ops.length - 1];
    var ps = (first[0] === 4 || first[0] === 5) ? first[1] : 0;
    var pe = (last[0] === 4 || last[0] === 5) ? qlen - last[1] - 1 : qlen - 1;
    if (!(r.pos < R || r.refEnd > refLen - R)) return;
    var inRef = (r.pos < R && r.refEnd > refLen - R) ? 3 : ((r.flag & 16) ? -1 : 1);
    segs.push({ rs: r.pos, re: r.refEnd, ps: ps, pe: pe, inRef: inRef, toSeq: r.toSeq, ops: r.ops });
  });
  segs.sort(function (a, b) { return a.ps - b.ps; });
  segs.push({ rs: 2, re: 2, ps: 2, pe: 2, inRef: 2, toSeq: 2, ops: [] });
  var lastX = 0;
  var walk = function (sub, ost, reverse) {
    var rp2 = sub.rs, sp = 0;
    sub.ops.forEach(function (o) {
      if (o[0] === 0) { if (o[1] > 0) lastX = o[1] - 1; }
      else if (o[0] === 1 && sub.rs + ost < rp2) out.push({ pos: rp2, len: o[1], qs: sp, qe: reverse ? sp + lastX : sp + o[1] });
      if (o[0] === 0 || o[0] === 2 || o[0] === 3) rp2 += o[1];
      if (o[0] === 0 || o[0] === 1 || o[0] === 4) sp += o[1];
    });
  };
  for (var i = 0; i < segs.length - 1; i++) {
    var a = segs[i], b = segs[i + 1];
    if (b.inRef === 3) continue;
    if (a.inRef === 3) {
      var rp3 = a.rs, sp3 = 0;
      a.ops.forEach(function (o) {
        if (o[0] === 0) { if (o[1] > 0) lastX = o[1] - 1; }
        if (o[0] === 1) out.push({ pos: rp3, len: o[1], qs: sp3, qe: sp3 + o[1] });
        if (o[0] === 0 || o[0] === 2 || o[0] === 3) rp3 += o[1];
        if (o[0] === 0 || o[0] === 1 || o[0] === 4) sp3 += o[1];
      });
      continue;
    }
    var overlap;
    if (a.toSeq === b.toSeq && a.inRef === b.inRef && a.toSeq === 1 && a.re < b.re) {
      overlap = b.rs < a.re ? a.re - b.rs : 0;
      if (Math.abs(a.pe - b.ps) - 1 > 20) out.push({ pos: a.re, len: Math.abs(a.pe - b.ps) - 1 + overlap, qs: a.pe, qe: a.pe + Math.abs(a.pe - b.ps) });
      walk(a, 0, false); walk(b, overlap, false);
    } else if (a.toSeq === b.toSeq && a.inRef === b.inRef && a.toSeq === -1 && b.re < a.re) {
      overlap = b.re > a.rs ? b.re - a.rs : 0;
      if (Math.abs(b.ps - a.pe) - 1 > 20) out.push({ pos: b.re, len: Math.abs(b.ps - a.pe) - 1 + overlap, qs: a.pe, qe: a.pe + (b.ps - a.pe - 1) });
      walk(a, 0, true); walk(b, overlap, true);
    }
  }
  return out;
}

// edlib.align(q, t, mode="HW", task="locations"): distance and first location.
function lxEdlibHW(q, t) {
  var n = q.length, m = t.length;
  var prev = new Int32Array(n + 1), cur = new Int32Array(n + 1);
  for (var i = 0; i <= n; i++) prev[i] = i;
  var dist = Infinity, end = -1;
  for (var j = 1; j <= m; j++) {
    cur[0] = 0;
    var tc = t.charCodeAt(j - 1);
    for (var k = 1; k <= n; k++) {
      var d = prev[k - 1] + (q.charCodeAt(k - 1) === tc ? 0 : 1);
      var e = prev[k] + 1, f = cur[k - 1] + 1;
      cur[k] = d < e ? (d < f ? d : f) : (e < f ? e : f);
    }
    if (cur[n] < dist) { dist = cur[n]; end = j - 1; }
    var tmp = prev; prev = cur; cur = tmp;
  }
  return { dist: dist, end: end, start: end < 0 ? -1 : lxEdlibStart(q, t, end, dist) };
}

function lxEdlibStart(q, t, end, dist) {
  var n = q.length;
  var rq = q.split('').reverse().join('');
  var lo = Math.max(0, end - n - dist - 1);
  var rt = t.slice(lo, end + 1).split('').reverse().join('');
  var m = rt.length;
  var prev = new Int32Array(n + 1), cur = new Int32Array(n + 1);
  for (var i = 0; i <= n; i++) prev[i] = i;
  var best = Infinity, bestJ = 0;
  for (var j = 1; j <= m; j++) {
    cur[0] = j;
    var tc = rt.charCodeAt(j - 1);
    for (var k = 1; k <= n; k++) {
      var d = prev[k - 1] + (rq.charCodeAt(k - 1) === tc ? 0 : 1);
      var e = prev[k] + 1, f = cur[k - 1] + 1;
      cur[k] = d < e ? (d < f ? d : f) : (e < f ? e : f);
    }
    if (cur[n] <= best) { best = cur[n]; bestJ = j; }
    var tmp = prev; prev = cur; cur = tmp;
  }
  return end - bestJ + 1;
}

function lxRevCompUpper(s) {
  var c = { A: 'T', T: 'A', G: 'C', C: 'G' };
  var out = '';
  for (var n = s.length - 1; n >= 0; n--) out += c[s[n]] || s[n];
  return out;
}

// CRISPRlungo_insert_analysis.align_edlib()
function lxAlignEdlib(ins, qs, qe, refDict, oriName, cache) {
  var L = ins.length;
  var res;
  if (cache.has(ins)) res = cache.get(ins);
  else {
    var insR = lxRevCompUpper(ins);
    res = refDict.map(function (rd) {
      var name = rd[0], seq = rd[1], r = [false];
      if (L <= seq.length) {
        var a = lxEdlibHW(ins, seq);
        if (a.dist < L * 0.05) r = [a.dist, a.start, a.end, '+'];
        var b = lxEdlibHW(insR, seq);
        if (b.dist < L * 0.05 && (r[0] === false || r[0] === 0 || b.dist < r[0])) r = [b.dist, b.start, b.end, '-'];
      } else {
        var a2 = lxEdlibHW(seq, ins);
        if (a2.end - a2.start > 20) {
          var a3 = lxEdlibHW(ins.slice(a2.start, a2.end), seq);
          if (a3.dist < L * 0.05) r = [a3.dist, a3.start, a3.end, '+'];
        }
        var b2 = lxEdlibHW(seq, insR);
        if (b2.end - b2.start > 20) {
          var b3 = lxEdlibHW(insR.slice(b2.start, b2.end), seq);
          if (b3.dist < L * 0.05 && (r[0] === false || r[0] === 0 || b3.dist < r[0])) r = [b3.dist, b3.start, b3.end, '-'];
        }
      }
      return [name, r];
    });
    res = res.map(function (x, n) { return [x, n]; })
      .sort(function (p, q) { return ((p[0][1][0] || 0) - (q[0][1][0] || 0)) || (p[1] - q[1]); })
      .map(function (p) { return p[0]; });
    cache.set(ins, res);
  }
  var out = [];
  res.forEach(function (x) {
    var r = x[1];
    if (r[0] === false || r[0] === 0) return;
    var name = (x[0] === oriName && r[3] === '-') ? 'inversion' : x[0];
    out.push([name, qs, qe, r[1], r[2], r[3]]);
  });
  return out;
}

function lxRunMinimap2(reference, query) {
  return new Promise(function (resolve, reject) {
    var w;
    try { w = new Worker('lungo_export_worker.js'); } catch (e) { reject(e); return; }
    w.onmessage = function (ev) { w.terminate(); ev.data.ok ? resolve(ev.data.sam) : reject(new Error(ev.data.error)); };
    w.onerror = function (ev) { w.terminate(); reject(new Error(ev.message || 'minimap2 worker failed')); };
    w.postMessage({ reference: reference, query: query });
  });
}

// Long-insert branch of confirm_insertion_seq(): parse minimap2 SAM like the CLI.
function lxParseInsertSam(sam, oriName) {
  var lines = sam.split('\n');
  var result = {};
  for (var li = 0; li < lines.length; li++) {
    var line = lines[li];
    if (!line || line[0] === '@') continue;
    var sp = line.trim().split(/\s+/);
    if (sp[5] === '*' || Number(sp[4]) < 40) continue;
    if (sp[1] !== '0' && sp[1] !== '16') continue;
    var rid = sp[0].split('_').map(Number);
    var reads = [sp];
    if (line.indexOf('SA:Z:') !== -1) {
      var sa = sp.filter(function (x) { return x.slice(0, 5) === 'SA:Z:'; })[0];
      var need = sa.split(';').length;
      while (reads.length < need && li + 1 < lines.length) {
        li += 1;
        var sec = lines[li].split('\t');
        if (Number(sec[1]) & 0x800) reads.push(sec);
      }
    }
    var muts = [];
    reads.forEach(function (r) {
      var flag = Number(r[1]), alignLen = 0, s = '';
      for (var c = 0; c < r[5].length; c++) {
        var ch = r[5][c];
        if (ch >= '0' && ch <= '9') s += ch;
        else { var v = Number(s); if (ch === 'M' || ch === 'D') alignLen += v; s = ''; }
      }
      if (alignLen < 20) return;
      var refName = r[2], strand = (flag & 0x10) ? -1 : 1;
      if (refName === oriName && strand === -1) refName = 'inversion';
      var rs = Number(r[3]);
      muts.push([refName, rid[2], rid[3], rs, rs + alignLen, strand]);
    });
    var key = rid[0] + '_' + rid[1];
    result[key] = muts;
  }
  return result;
}

// annotations[rowIndex][keptIndex] = list of (name, qs, qe, refStart, refEnd, strand)
async function lxAnnotateInsertions(reads, refName) {
  var sam = lxSamByRead(treatedAlignSam);
  var withControl = !!controlFile;
  var refDict = [[refName, reference.toUpperCase()]];
  if (!withControl) LX_POSSIBLE_INSERTION.forEach(function (p) { refDict.push([p[0], p[1]]); });
  var cache = new Map();
  var ann = {}, longQ = '', stats = { short: 0, long: 0, noCoords: 0, minimap2: 'not needed' };
  reads.forEach(function (rd, ri) {
    var cands = null;
    rd.kept.forEach(function (p, mi) {
      if (p.type !== 'Ins' || p.seq.length <= 20) return;
      if (!cands) cands = lxCliInsertions(sam.get(rd.base) || [], reference.length);
      var hit = null;
      for (var c = 0; c < cands.length; c++) {
        if (!cands[c].used && cands[c].pos === p.start && cands[c].len === p.len) { hit = cands[c]; hit.used = true; break; }
      }
      if (!hit) { stats.noCoords += 1; hit = { qs: 'NA', qe: 'NA' }; }
      if (p.seq.length < 100) {
        stats.short += 1;
        var a = lxAlignEdlib(p.seq, hit.qs, hit.qe, refDict, refName, cache);
        if (a.length) { ann[ri] = ann[ri] || {}; ann[ri][mi] = a; }
      } else {
        stats.long += 1;
        longQ += '>' + ri + '_' + mi + '_' + hit.qs + '_' + hit.qe + '\n' + p.seq + '\n';
      }
    });
  });
  if (longQ) {
    var refFa = '>' + refName + '\n' + reference + '\n';
    if (!withControl) LX_POSSIBLE_INSERTION.forEach(function (p) { refFa += p[0] + '\n' + p[1] + '\n'; });
    try {
      var parsed = lxParseInsertSam(await lxRunMinimap2(refFa, longQ), refName);
      Object.keys(parsed).forEach(function (k) {
        var x = k.split('_');
        ann[x[0]] = ann[x[0]] || {};
        ann[x[0]][x[1]] = parsed[k];
      });
      stats.minimap2 = 'ok';
    } catch (err) {
      console.warn('minimap2 for insert annotation failed', err);
      stats.minimap2 = 'failed: ' + err.message;
    }
  }
  return { ann: ann, stats: stats };
}

function lxIntegrationString(tuples) {
  var s = '';
  (tuples || []).forEach(function (t) { s += t[0] + '_' + t[1] + 'to' + t[2] + 'inSeq_' + t[3] + 'to' + t[4] + 'inRef_' + t[5] + ','; });
  return s;
}

function lxInducedKeySet() {
  var keys = new Set();
  (mainResult.induced_mutations || []).forEach(function (m) {
    if (typeof mutationInWindow === 'function' &&
        !mutationInWindow(m, cleavagePos, cleavagePos2, windowRange, wholeWindow === 1)) return;
    keys.add(lxInducedKey(lxParseWebMut(m)));
  });
  return keys;
}

// Reads in alignment (= FASTQ) order with their kept / other mutations.
function lxCollectReads() {
  var sig = new Set(filteredSigMut || []);
  var order = new Map(), n = 0;
  (treatedAlignSam || []).forEach(function (line) {
    if (!line || line[0] === '@') return;
    var q = line.slice(0, line.indexOf('\t'));
    var cut = q.indexOf('_AI_');
    var b = cut === -1 ? q : q.slice(0, cut);
    if (!order.has(b)) order.set(b, n++);
  });
  var reads = Object.entries(mainResult.treated_read_to_mut).map(function (e) {
    var id = e[0], val = e[1];
    var cut = id.indexOf('_AI_');
    var base = cut === -1 ? id : id.slice(0, cut);
    var kept = [], other = [];
    (val[0] || []).forEach(function (m) { (sig.has(m) ? kept : other).push(m); });
    (val[1] || []).forEach(function (m) { other.push(m); });
    return {
      id: id, base: base, ord: order.has(base) ? order.get(base) : Infinity,
      kept: kept.map(lxParseWebMut),
      other: other.map(lxParseWebMut).sort(function (a, b) { return a.start - b.start; })
    };
  });
  var splitNo = function (id) { var m = id.match(/_s(\d+)$/); return m ? Number(m[1]) : -1; };
  reads.sort(function (a, b) { return (a.ord - b.ord) || (splitNo(a.id) - splitNo(b.id)); });
  return reads;
}

function lxCliReadId(id) {
  return id.replace('_AI_', '_alignStrandInfo_').replace(/_s(\d+)$/, '_split_$1');
}

// One row per used read, CLI read_classification.txt semantics.
function lxBuildReadRows(reads, ann) {
  var inducedKeys = lxInducedKeySet();
  return reads.map(function (rd, ri) {
    var integ = '';
    rd.kept.forEach(function (p, mi) {
      if (p.type === 'Ins' && ann[ri] && ann[ri][mi]) integ += lxIntegrationString(ann[ri][mi]);
    });
    var c = lxClassify(rd.kept, integ, inducedKeys);
    return {
      id: lxCliReadId(rd.id),
      cls: c.cls, info: c.info, integ: c.integ, induce: c.induce,
      whole: rd.other.length ? rd.other.map(lxCliToken).join(',') : '[]'
    };
  });
}

function lxReadClassificationTxt(rows) {
  var out = ['Read_id\tClassification\tMutation_info\tIntegration_info\tInduce_type\twhole_mutation'];
  rows.forEach(function (r) {
    out.push([r.id, r.cls, r.info, r.integ, r.induce, r.whole].join('\t'));
  });
  return out.join('\r\n') + '\r\n';
}

// write_read_count(): allele table + summary
function lxPatternCount(rows, minCnt, minFreq) {
  var dict = new Map();
  rows.forEach(function (r) {
    if (!dict.has(r.info)) dict.set(r.info, { cnt: 0, integ: r.integ, induce: r.induce, cls: r.cls });
    dict.get(r.info).cnt += 1;
  });
  var all = rows.length;
  var sorted = Array.from(dict.entries()).sort(function (a, b) { return b[1].cnt - a[1].cnt; });
  var kept = [], lowCnt = 0, filtAll = 0;
  sorted.forEach(function (kv) {
    if (kv[1].cnt < minCnt || lxRound2(kv[1].cnt * 100 / all) < minFreq) { lowCnt += kv[1].cnt; return; }
    kept.push(kv);
    filtAll += kv[1].cnt;
  });
  var mutCnt = {}, preciseCnt = {};
  ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv', 'Complex', 'ComplexWithLargeMut'].forEach(function (k) { mutCnt[k] = 0; });
  ['Precise', 'Precise_with_mutations', 'Partial', 'Partial_with_mutations', 'X'].forEach(function (k) { preciseCnt[k] = 0; });
  var lines = ['Mutation_patterns\tCount\t(%)\tIntegration_information\tInduced_mutation\tDetail_mutation_informations'];
  var table = [];
  kept.forEach(function (kv) {
    var v = kv[1];
    var pct = lxRound2(v.cnt * 100 / filtAll);
    lines.push([v.cls, v.cnt, lxPyFloat(pct), v.integ.slice(0, -1), v.induce, kv[0]].join('\t'));
    table.push({ cls: v.cls, cnt: v.cnt, pct: pct, integ: v.integ.slice(0, -1), induce: v.induce, info: kv[0] });
    mutCnt[v.cls] = (v.cls in mutCnt) ? mutCnt[v.cls] + v.cnt : 1;
    preciseCnt[v.induce] = (v.induce in preciseCnt) ? preciseCnt[v.induce] + v.cnt : 1;
  });
  return { text: lines.join('\n') + '\n', table: table, mutCnt: mutCnt, preciseCnt: preciseCnt, lowCnt: lowCnt };
}

// The wasm reports read counts only for the treated sample; running it once
// more with the control alignment as "treated" gives the control counts.
function lxWasmCounts(samLines) {
  return new Promise(function (resolve, reject) {
    var w;
    try { w = new Worker('lungo_worker.js', { type: 'module' }); } catch (e) { reject(e); return; }
    w.onmessage = function (ev) {
      w.terminate();
      if (ev.data && ev.data.type === 'success') resolve(ev.data.result.treated_align_cnt);
      else reject(new Error((ev.data && ev.data.message) || 'wasm failed'));
    };
    w.onerror = function (ev) { w.terminate(); reject(new Error(ev.message || 'wasm worker failed')); };
    var s = samLines.join('\n');
    w.postMessage({
      controlSamString: s, treatedSamString: s,
      inducedMutationString: (desiredAlignSam && desiredAlignSam.join) ? desiredAlignSam.join('\n') : '',
      reference: reference, cvPos: cleavagePos, cvPos2: cleavagePos2,
      window: windowRange, wholeWindow: wholeWindow, filter1: true, windowFilter: true
    });
  });
}

async function lxPreprocess() {
  var t = mainResult.treated_align_cnt;
  var tAll = (treatedInputReadCount == null) ? t.All_reads + t.Unmapped : treatedInputReadCount;
  var tUn = (treatedUnmappedCount == null) ? t.Unmapped : treatedUnmappedCount;
  if (controlFile) {
    var c = null;
    try { c = await lxWasmCounts(controlAlignSam); } catch (err) { console.warn('control counts failed', err); }
    var cUn = (controlUnmappedCount == null) ? (c ? c.Unmapped : 'NA') : controlUnmappedCount;
    var cAll = (controlInputReadCount == null) ? (c ? c.All_reads + cUn : 'NA') : controlInputReadCount;
    var head = ['Control_all_reads', 'Control_unmapped', 'Control_low_quality', 'Control_short', 'Control_used'];
    var vals = [cAll, cUn, c ? c.Low_quality : 'NA', c ? c.Short : 'NA', c ? c.Used : 'NA'];
    var head2 = ['Treated_all_reads', 'Treated_unmapped', 'Treated_low_quality', 'Treated_short', 'Treated_used'];
    var vals2 = [tAll, tUn, t.Low_quality, t.Short, t.Used];
    var text = head.join('\t') + '\t\t' + head2.join('\t') + '\t\t\n' + vals.join('\t') + '\t\t' + vals2.join('\t') + '\t';
    var dict = {};
    head.forEach(function (k, n) { dict[k] = vals[n]; });
    head2.forEach(function (k, n) { dict[k] = vals2[n]; });
    return { text: text, dict: dict };
  }
  var keys = ['all_reads', 'short', 'unmapped', 'low_quality', 'used'];
  var v = [tAll, t.Short, tUn, t.Low_quality, t.Used];
  var d = {};
  keys.forEach(function (k, n) { d[k] = v[n]; });
  return { text: keys.join('\t') + '\n' + v.join('\t') + '\t\n', dict: d };
}

function lxSummaryCount(pre, pc) {
  var s1 = '', s2 = '';
  ['Treated_all_reads', 'Treated_unmapped', 'Treated_low_quality', 'Treated_short', 'all_reads', 'unmapped', 'low_quality', 'short', 'used_with_supple'].forEach(function (k) {
    if (k in pre) { s1 += k + '\t'; s2 += pre[k] + '\t'; }
  });
  if ('Treated_used' in pre) { s1 += 'Treated_low_cnt_reads\tTreated_used\t'; s2 += pc.lowCnt + '\t' + (Number(pre.Treated_used) - pc.lowCnt) + '\t'; }
  else if ('used' in pre) { s1 += 'low_cnt_reads\tused\t'; s2 += pc.lowCnt + '\t' + (Number(pre.used) - pc.lowCnt) + '\t'; }
  s1 += '\t'; s2 += '\t';
  ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv', 'Complex', 'ComplexWithLargeMut'].forEach(function (k) {
    if (k in pc.mutCnt) { s1 += k + '\t'; s2 += pc.mutCnt[k] + '\t'; }
  });
  s1 += '\t'; s2 += '\t';
  ['Precise', 'Precise_with_mutations', 'Partial', 'Partial_with_mutations', 'X'].forEach(function (k) {
    if (k in pc.preciseCnt) { s1 += k + '\t'; s2 += pc.preciseCnt[k] + '\t'; }
  });
  return s1 + '\n' + s2 + '\n';
}

function lxPyBool(b) { return b ? 'True' : 'False'; }

function lxInducedMutStr() {
  var keys = Array.from(lxInducedKeySet());
  if (!keys.length) return 'False';
  return keys.map(function (k) { return lxCliToken(lxParseWebMut(k)); }).join(',');
}

function lxInputSummary(opts) {
  var t1 = target, t2 = target2 ? target2 : 'False';
  var cv2 = (cleavagePos2 === '' || cleavagePos2 == null) ? 'False' : cleavagePos2;
  var st2 = (cleavagePos2 === '' || cleavagePos2 == null) ? 'False' : strand2;
  var lines = [
    'Target_1 :' + t1,
    'Target_2 :' + t2,
    'Ref_seq :' + reference,
    'CleavagePos_1 :' + cleavagePos,
    'CleavageStrand_1 :' + strand,
    'CleavagePos_2 :' + cv2,
    'CleavageStrand_2 :' + st2,
    'Window : ' + windowRange,
    'Window_between_cleavage : ' + lxPyBool(wholeWindow === 1 && !!target2),
    'induced_mut :' + opts.inducedStr,
    'cut_pos_in_target :' + cleavagePosTarget,
    'original_target :1',
    'minimum_read_count : ' + opts.minCnt,
    'minimum_read_frequency : ' + opts.minFreq,
    'induced_sequence_path : ' + (desiredSeq ? 'web_desired_sequence.fasta' : 'None'),
    'largeins_cutlen : ' + largeInsLen,
    'largedel_cutlen : ' + largeDelLen,
    '',
    '# Run',
    'CRISPRlungo_version : web',
    'Started : ' + opts.date,
    'Treated_file : ' + (treatedFile ? treatedFile.name : ''),
    'Control_file : ' + (controlFile ? controlFile.name : 'False'),
    '',
    '# Effective values',
    'reference_name : ' + opts.refName,
    'reference_length : ' + reference.length,
    'allele_plot_window_used : ' + plotWindow,
    'minimap2_longjoin_bandwidth : ' + longjoinBandWidth,
    'minimap2_chaining_bandwidth : ' + chainingBandWidth,
    'induced_mutation_count : ' + lxInducedKeySet().size
  ];
  return lines.join('\n') + '\n';
}

// p-values the wasm reports for the mutations it kept (significant ones only).
function lxPValues() {
  var lines = ['mutation\tcount_edited\tp_value\tsignificant'];
  Object.entries(mainResult.significant_keys).forEach(function (kv) {
    var p = lxParseWebMut(kv[0]);
    var tup;
    if (p.type === 'Del') tup = "('deletion', " + p.start + ', ' + p.len + ')';
    else if (p.type === 'Ins') tup = "('insertion', " + p.start + ', ' + p.len + ", '" + p.seq + "')";
    else tup = "('substitution', " + p.start + ', ' + p.len + ", '" + p.from + "', '" + p.to + "')";
    lines.push([tup, kv[1][1], kv[1][0], 'True'].join('\t'));
  });
  return lines.join('\n') + '\n';
}

// allele_table() port (CRISPRlungo_visualization.py)
function lxAlleleTable(table, inducedStr) {
  var ref = reference;
  var pw = plotWindow;
  var cv = cleavagePos;
  var cv2 = (cleavagePos2 === '' || cleavagePos2 == null) ? false : Number(cleavagePos2);
  var out = ['ref_seq\tmut_seq\tRaw_count\t%\tCIGAR\tMutation_info'];
  var hasInducedStr = inducedStr !== 'False';
  table.forEach(function (row, n) {
    if (row.info === 'None') {
      var endWT = cv2 ? cv2 + pw : cv + pw;
      var seg = ref.slice(cv - pw, endWT);
      var m = cv2 ? (cv2 - cv + pw * 2) : pw * 2;
      out.push([seg, seg, row.cnt, lxPyFloat(row.pct), m + 'M', row.info].join('\t'));
      return;
    }
    var toks = row.info.split(',');
    var parse = function (t) {
      return t.replace(':', '_').replace('>', '_').split('_').map(function (x) { return /^\d+$/.test(x) ? Number(x) : x; });
    };
    var ex = parse(toks[0]);
    ex[1] = Number(ex[0]) - pw - 1;
    var alignRef = '', alignMut = '', cig = '';
    var mut = ex;
    toks.forEach(function (tok) {
      mut = parse(tok);
      if (ex[1] + 1 < mut[0]) {
        alignRef += ref.slice(ex[1] + 1, mut[0]);
        alignMut += ref.slice(ex[1] + 1, mut[0]);
      }
      var dist = mut[0] - ex[1] - 1;
      if (dist > 0) cig += dist + 'M';
      if (mut[2] === 'Sub') {
        alignRef += mut[4]; alignMut += mut[5]; cig += mut[3] + '=';
      } else if (String(mut[2]).indexOf('Del') !== -1) {
        alignRef += ref.slice(mut[0], mut[1] + 1);
        alignMut += '-'.repeat(mut[3]);
        cig += mut[3] + 'D';
      } else if (String(mut[2]).indexOf('Ins') !== -1) {
        alignRef += '-'.repeat(mut[3]);
        alignMut += mut[4];
        cig += mut[3] + 'I';
        if (String(ex[2]).indexOf('Del') !== -1 && ex[0] === mut[0]) mut[1] -= 1;
        else mut[1] -= 2;
      }
      ex = mut;
    });
    var cv2v = cv2 === false ? 0 : cv2;
    var tail = (mut[1] + pw > cv2v + pw) ? ref.slice(mut[1] + 1, mut[1] + 1 + pw) : ref.slice(mut[1] + 1, cv2v + pw);
    alignRef += tail; alignMut += tail; cig += tail.length + 'M';
    if (!hasInducedStr || n !== 0) {
      out.push([alignRef, alignMut, row.cnt, lxPyFloat(row.pct), cig, row.info].join('\t'));
    }
  });
  return out.join('\n') + '\n';
}

var LX_PLOTS = [
  ['alignSumPlot', 'Alignment_summary'], ['allPiePlot', 'Mutation_pie_chart'], ['mutPiePlot', 'Mutation_type_pie_chart'],
  ['combinedPlot', 'Deletion_and_insertions_per_position'], ['insertionsPlot', 'Insertions_per_position'],
  ['deletionsPlot', 'Deletions_per_position'], ['largeDeletionPlot', 'Large_deletions'],
  ['substitutionsPlot', 'Substitutions_per_position'], ['customPiePlot', 'Custom_classification_pie_chart']
];

async function lxPlotFiles(zip) {
  for (var n = 0; n < LX_PLOTS.length; n++) {
    var gd = document.getElementById(LX_PLOTS[n][0]);
    if (!gd || !gd.data || !gd.data.length || typeof Plotly === 'undefined') continue;
    try {
      var url = await Plotly.toImage(gd, { format: 'png', width: 1000, height: 600 });
      zip.file('plots/' + LX_PLOTS[n][1] + '.png', url.split(',')[1], { base64: true });
    } catch (err) { console.warn('plot export failed', LX_PLOTS[n][0], err); }
  }
  var svgs = [['allelePlot', 'allele_plot.svg'], ['customAllelePlot', 'custom_allele_plot.svg']];
  svgs.forEach(function (s) {
    var svg = document.querySelector('#' + s[0] + ' svg');
    if (!svg) return;
    var clone = svg.cloneNode(true);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    zip.file('plots/' + s[1], '<?xml version="1.0" encoding="UTF-8"?>\n' + clone.outerHTML);
  });
}

var LX_README = [
  'CRISPRLungo web tool - result files',
  '',
  'Folder layout and file formats follow the CRISPRlungo command-line tool, so the',
  'folder can be passed to CRISPRlungoAllele (e.g. CRISPRlungoAllele -i categories.txt <this folder>).',
  '',
  'Reads are called by the web engine (Rust/WebAssembly). The CLI (Python) engine can call a',
  'small number of reads differently, mostly reads split into supplementary alignments.',
  '',
  'Not produced by the web tool:',
  '- mutation_pattern_p_values.txt lists only mutations the web engine kept as significant',
  '  (columns: mutation, count_edited, p_value, significant; p_value -1 = absent from control).',
  '- No BAM / IGV session (align/*.sam hold the alignments made in the browser), no combined_graphs.html.',
  '- Plots are the web plots (PNG / SVG), not the CLI matplotlib figures.',
  '- input_summary.txt has no "# Options" block (web settings are listed instead).',
  ''
].join('\n');

// Base qualities per read name from the uploaded FASTQ (the browser aligner
// works on FASTA, so its SAM has QUAL '*').
async function lxReadQualities(file) {
  var quals = new Map();
  if (!file) return quals;
  var text;
  try {
    var head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
    if (head[0] === 0x1f && head[1] === 0x8b) {
      if (typeof DecompressionStream === 'undefined') return quals;
      text = await new Response(file.stream().pipeThrough(new DecompressionStream('gzip'))).text();
    } else {
      text = await file.text();
    }
  } catch (err) { console.warn('could not read qualities', err); return quals; }
  var lines = text.split(/\r?\n/);
  for (var i = 0; i + 3 < lines.length; i += 4) {
    if (lines[i][0] !== '@') { i -= 3; continue; }
    quals.set(lines[i].slice(1).split(/\s/)[0], lines[i + 3]);
  }
  return quals;
}

function lxSamText(lines, refName, quals) {
  var out = (lines || []).map(function (l) {
    if (l[0] === '@') return l.replace(/\tSN:Ref\t/, '\tSN:' + refName + '\t');
    var f = l.split('\t');
    var flag = Number(f[1]);
    f[0] = f[0].replace('_AI_', '_alignStrandInfo_');
    if (f[2] === 'Ref') f[2] = refName;
    if (quals && f[10] === '*' && !(flag & 0x800) && !(flag & 0x100) && !(flag & 4)) {
      var q = quals.get(f[0].split('_alignStrandInfo_')[0]);
      if (q && q.length === f[9].length) f[10] = (flag & 16) ? q.split('').reverse().join('') : q;
    }
    return f.join('\t');
  });
  return out.join('\n') + '\n';
}

async function lxBuildFiles(progress) {
  progress = progress || function () { };
  var date = new Date().toISOString();
  var refName = lxRefName();
  var minCnt = lxNum('readCountMin', 0);
  var minFreq = lxNum('readFreqMin', 0) * 100;
  var inducedStr = lxInducedMutStr();
  progress('Collecting reads ...');
  var reads = lxCollectReads();
  progress('Annotating inserted sequences ...');
  var an = await lxAnnotateInsertions(reads, refName);
  var rows = lxBuildReadRows(reads, an.ann);
  var pc = lxPatternCount(rows, minCnt, minFreq);
  progress('Counting reads ...');
  var pre = await lxPreprocess();
  var files = {};
  files['results/read_classification.txt'] = lxReadClassificationTxt(rows);
  files['results/mutation_patter_count.txt'] = pc.text;
  files['results/preprocess_count.txt'] = pre.text;
  files['results/mutation_summary_count.txt'] = lxSummaryCount(pre.dict, pc);
  files['results/input_summary.txt'] = lxInputSummary({ inducedStr: inducedStr, minCnt: minCnt, minFreq: minFreq, date: date, refName: refName });
  files['results/Allele_table.txt'] = lxAlleleTable(pc.table, inducedStr);
  if (controlFile) files['results/mutation_pattern_p_values.txt'] = lxPValues();
  files['ref_seq/ref_wo_umi.fasta'] = '>' + refName + '\n' + reference + '\n';
  if (desiredSeq) files['web_desired_sequence.fasta'] = '>Induced\n' + desiredSeq + '\n';
  if (treatedAlignSam && treatedAlignSam.length) {
    progress('Writing alignments ...');
    files['align/Treated_alignment.sam'] = lxSamText(treatedAlignSam, refName, await lxReadQualities(treatedFile));
    if (treatedInputReadCount != null) files['align/Treated_alignment.sam.align_stats'] = 'input_reads\t' + treatedInputReadCount + '\nunmapped\t' + treatedUnmappedCount + '\n';
  }
  if (controlAlignSam && controlAlignSam.length) {
    files['align/Control_alignment.sam'] = lxSamText(controlAlignSam, refName, await lxReadQualities(controlFile));
    if (controlInputReadCount != null) files['align/Control_alignment.sam.align_stats'] = 'input_reads\t' + controlInputReadCount + '\nunmapped\t' + controlUnmappedCount + '\n';
  }
  if (typeof customClassifierCode !== 'undefined' && customClassifierCode) {
    files['custom_results/custom_classifier.js'] = customClassifierCode + '\n';
    var lines = ['Category\tReads'];
    Object.entries(customCategoryCounts).sort(function (a, b) { return b[1] - a[1]; })
      .forEach(function (kv) { lines.push(kv[0] + '\t' + kv[1]); });
    files['custom_results/custom_category_count.txt'] = lines.join('\n') + '\n';
  }
  files['README_web_export.txt'] = LX_README + '\nInsertion annotation: ' + an.stats.short + ' inserts 21-99 bp (edlib-equivalent), ' +
    an.stats.long + ' inserts >= 100 bp (minimap2: ' + an.stats.minimap2 + '), ' + an.stats.noCoords + ' without read coordinates.\n';
  return { files: files, rows: rows, patterns: pc, date: date, annotation: an.stats };
}

async function downloadResultFiles() {
  var statusEl = document.getElementById('exportStatus');
  var btn = document.getElementById('exportResultsBtn');
  if (!mainResult) { if (statusEl) statusEl.textContent = 'Run an analysis first.'; return; }
  if (typeof JSZip === 'undefined') { if (statusEl) statusEl.textContent = 'ZIP library failed to load.'; return; }
  if (btn) btn.disabled = true;
  if (statusEl) statusEl.textContent = 'Preparing files ...';
  try {
    var built = await lxBuildFiles(function (msg) { if (statusEl) statusEl.textContent = msg; });
    var stamp = built.date.replace(/[:.]/g, '-').slice(0, 19);
    var base = 'CRISPRlungo_web_' + stamp + '/';
    var zip = new JSZip();
    var folder = zip.folder(base.slice(0, -1));
    Object.entries(built.files).forEach(function (kv) { folder.file(kv[0], kv[1]); });
    await lxPlotFiles(folder);
    var blob = await zip.generateAsync({ type: 'blob', compression: 'DEFLATE', compressionOptions: { level: 6 } });
    var link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = 'CRISPRlungo_web_' + stamp + '.zip';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    setTimeout(function () { URL.revokeObjectURL(link.href); }, 2000);
    if (statusEl) statusEl.textContent = 'Downloaded ' + Object.keys(built.files).length + ' files + plots.';
  } catch (err) {
    console.error(err);
    if (statusEl) statusEl.textContent = 'Error: ' + err.message;
  } finally {
    if (btn) btn.disabled = false;
  }
}
