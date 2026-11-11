// 20261111_pa
// Prime assembly (PA) outcome classification in the browser.
//
// Mirrors src/CRISPRlungo_pa.py. Five alleles are derived from the reference, the
// two cleavage sites and the donor - wild type, precise PA, dual flap without
// donor, 5' flap only, 3' flap only. Each is aligned and analysed with the same
// workers the normal run uses, then every read is assigned to the allele that
// explains the most of it, and classified from the mutations in the windows
// around that allele's cleavage, flap and fragment-junction positions.

(function () {
  'use strict';

  var REF_ORDER = ['WT', 'PA', 'DUAL', 'F5', 'F3'];
  var REF_LABEL = {
    WT: 'wild type', PA: 'precise PA', DUAL: 'dual flap, no donor',
    F5: "5' flap only", F3: "3' flap only"
  };
  var PATTERN_OF_REF = { WT: 'none', PA: 'full_donor', DUAL: 'dual_flap', F5: '5p_flap', F3: '3p_flap' };
  var CATEGORIES = ['Precise_PA', 'Flap_install_wo_indels', 'Flap_integration_wDel',
    'Flap_integration_wIns', 'Flap_install_with_Large_deletion', 'Large_deletion',
    'Unedited', 'Indel_at_nick_site', 'Others', 'Not_analysed'];

  // Shown in this order, with the unedited reads first and under the same name
  // the outcome table above uses for them.
  var DISPLAY_ORDER = ['Unedited', 'Precise_PA', 'Flap_install_wo_indels',
    'Flap_integration_wDel', 'Flap_integration_wIns', 'Flap_install_with_Large_deletion',
    'Large_deletion', 'Indel_at_nick_site', 'Others'];
  var DISPLAY_NAME = { Unedited: 'WT' };

  function parseFlapLen(spec) {
    var parts = String(spec || '32').split(',').map(function (x) { return x.trim(); })
      .filter(function (x) { return x !== ''; });
    var values;
    if (parts.length === 1) values = [Number(parts[0]), Number(parts[0])];
    else if (parts.length === 2) values = [Number(parts[0]), Number(parts[1])];
    else throw new Error('Flap length takes one value, or two separated by a comma.');
    if (!values.every(function (v) { return v > 0 && isFinite(v); })) {
      throw new Error('Flap length values must be positive whole numbers.');
    }
    // main.js sorts nothing, but the second protospacer can sit before the first
    // in the reference; the flap lengths follow the protospacers, not the order
    // the sites happen to fall in.
    if (window.cleavagePos > window.cleavagePos2 && window.cleavagePos2 !== '') values = values.reverse();
    return values;
  }

  function buildAlleles(ref, cv, cv2, fragments, f5, f3) {
    var lo = Math.min(cv, cv2), hi = Math.max(cv, cv2);
    var donor = fragments.map(function (r) { return r.seq; }).join('');
    var left = ref.slice(0, lo + 1), right = ref.slice(hi + 1);
    var f5seq = donor.slice(0, f5), f3seq = donor.slice(donor.length - f3);

    var junctions = [], offset = 0;
    for (var i = 0; i < fragments.length - 1; i++) {
      offset += fragments[i].seq.length;
      junctions.push(lo + offset);
    }

    var alleles = {};
    alleles.WT = { seq: ref, cleavage: [lo, hi], flap: [], junction: [] };
    alleles.PA = {
      seq: left + donor + right,
      cleavage: [lo, lo + donor.length],
      flap: [lo + f5, lo + donor.length - f3],
      junction: junctions
    };
    alleles.DUAL = {
      seq: left + f5seq + f3seq + right,
      cleavage: [lo, lo + f5 + f3], flap: [lo + f5], junction: []
    };
    alleles.F5 = { seq: left + f5seq + right, cleavage: [lo, lo + f5], flap: [], junction: [] };
    alleles.F3 = { seq: left + f3seq + right, cleavage: [lo, lo + f3], flap: [], junction: [] };
    return { alleles: alleles, donor: donor, flapSeqs: [f5seq, f3seq] };
  }

  // The web mutation strings are START_END:TYPE_DETAIL, with Del_N, Ins_SEQ,
  // Sub_A>G, LargeDel_N and LargeIns_SEQ. (The CLI spells insertions and
  // substitutions differently; only this form reaches us here.)
  function parseMut(token) {
    if (!token || token.indexOf(':') < 0) return null;
    var span = token.split(':')[0];
    var detail = token.slice(span.length + 1);
    var se = span.split('_');
    var start = Number(se[0]), end = Number(se[1]);
    if (!isFinite(start) || !isFinite(end)) return null;
    var kind = detail.split('_')[0];
    var rest = detail.slice(kind.length + 1);
    var len;
    if (kind === 'Del' || kind === 'LargeDel') len = Number(rest) || (end - start + 1);
    var seq = '';
    if (kind === 'Ins' || kind === 'LargeIns') {
      seq = (rest || '').replace(/_inv$/, '').toUpperCase();
      len = seq.length;
    } else if (kind !== 'Del' && kind !== 'LargeDel') len = 1;
    return { start: start, end: end, kind: kind, len: len, seq: seq };
  }

  function baseReadId(id) {
    return String(id).replace(/_AI_.*$/, '').replace(/_alignStrandInfo_.*$/, '').replace(/_s\d+$/, '').replace(/_split_\d+$/, '');
  }

  function overlaps(m, site, half) {
    return !(m.end < site - half || m.start > site + 1 + half);
  }

  function windowIdentity(muts, cleavage) {
    var lo = Math.min.apply(null, cleavage), hi = Math.max.apply(null, cleavage) + 1;
    var span = Math.max(1, hi - lo + 1), sub = 0, ins = 0, del = 0;
    muts.forEach(function (m) {
      if (m.end < lo || m.start > hi) return;
      if (m.kind === 'Sub') sub += m.len;
      else if (m.kind === 'Ins' || m.kind === 'LargeIns') ins += m.len;
      else del += Math.min(m.len, span);
    });
    var matches = Math.max(0, span - sub - del);
    var denom = matches + sub + ins;
    return denom ? matches / denom : 0;
  }

  // Per read: how much of it this allele accounts for. Without this an allele far
  // shorter than the read wins by default, because whatever it cannot explain is
  // clipped away and costs nothing.
  function coverageFromSam(lines) {
    var aligned = {}, total = {};
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i];
      if (!line || line.charAt(0) === '@') continue;
      var f = line.split('\t');
      if (f.length < 11) continue;
      var flag = Number(f[1]);
      if (flag & 0x100) continue;
      var id = baseReadId(f[0]);
      var re = /(\d+)([MIDNSHP=X])/g, m, a = 0, t = 0;
      while ((m = re.exec(f[5])) !== null) {
        var len = Number(m[1]), op = m[2];
        if (op === 'M' || op === 'I' || op === '=' || op === 'X') { a += len; t += len; }
        else if (op === 'S' || op === 'H') { t += len; }
      }
      aligned[id] = (aligned[id] || 0) + a;
      total[id] = Math.max(total[id] || 0, t);
    }
    var out = {};
    Object.keys(aligned).forEach(function (id) {
      if (total[id]) out[id] = aligned[id] / total[id];
    });
    return out;
  }

  // minimap2's -r bandwidths are derived from the reference length, so each
  // allele needs its own: the precise-PA allele is several kb longer than the
  // wild-type one the main run sized them for.
  function bandWidths(seq) {
    var longjoin = parseInt(seq.length * 0.3, 10);
    return { longjoin: longjoin, chaining: longjoin < 500 ? longjoin : 500 };
  }

  function align(fastq, reference, fileType) {
    var band = bandWidths(reference);
    return new Promise(function (resolve, reject) {
      var worker = new Worker('minimap2_worker_prog.js');
      worker.onmessage = function (e) {
        if (e.data.error) { worker.terminate(); reject(new Error(e.data.error)); return; }
        if (e.data.type === 2) { worker.terminate(); resolve(e.data.result); }
      };
      worker.onerror = function (e) { worker.terminate(); reject(e); };
      worker.postMessage({
        fastqFile: fastq,
        reference: reference,
        longjoinBandWidth: band.longjoin,
        chainingBandWidth: band.chaining,
        fileType: fileType,
        extraOpt: window.minimap2Opt,
        saLenThreshold: window.alignSaLenThreshold
      });
    });
  }

  function analyse(controlSam, treatedSam, allele) {
    return new Promise(function (resolve, reject) {
      var worker = new Worker('lungo_worker.js', { type: 'module' });
      worker.onmessage = function (e) {
        worker.terminate();
        // 20261111_worker: a wasm panic comes back as an error message, and used to
        // resolve with undefined, which only showed up far downstream
        if (!e.data || e.data.type === 'error' || !e.data.result) {
          reject(new Error((e.data && e.data.message) || 'the analysis returned nothing'));
        } else resolve(e.data.result);
      };
      worker.onerror = function (e) { worker.terminate(); reject(e); };
      worker.postMessage({
        controlSamString: controlSam.join('\n'),
        treatedSamString: treatedSam.join('\n'),
        inducedMutationString: '',
        reference: allele.seq,
        cvPos: allele.cleavage[0],
        cvPos2: allele.cleavage[1],
        window: window.paWindow,
        wholeWindow: 1,
        filter1: true,
        windowFilter: true,
        alpha: (typeof window.mutAllowance === 'number' ? window.mutAllowance : 0.05),
        length_min: (typeof window.lengthMin === 'number' ? window.lengthMin : 10)
      });
    });
  }

  // significant_keys drops nanopore noise, but on an allele reference the
  // difference between that allele and the control is itself a multi-kb event
  // present in both, so anything longer than the large-deletion cutoff is kept
  // whatever the test said. Same guard as --background_max_len in the CLI.
  function readsFromResult(result, largeCut, tokensOut) {
    var significant = new Set(Object.keys(result.significant_keys || {}));
    var table = {};
    Object.entries(result.treated_read_to_mut || {}).forEach(function (entry) {
      var id = baseReadId(entry[0]);
      var muts = table[id] || (table[id] = []);
      if (tokensOut && !(id in tokensOut)) tokensOut[id] = '';  // 20261111_alignview
      (entry[1][0] || []).forEach(function (token) {
        var m = parseMut(token);
        if (!m) return;
        if (!significant.has(token) && m.len <= largeCut) return;
        muts.push(m);
        // the alignment viewer parses the raw token string, so keep it as well
        if (tokensOut) tokensOut[id] += token + ',';  // 20261111_alignview
      });
    });
    return table;
  }

  // Fraction of other's k-mers present in seq: robust to the indels an ONT read
  // puts inside an inserted stretch, which a base-by-base comparison is not.
  function kmerShare(seq, other, k) {
    k = k || 6;
    if (!seq || !other || other.length < k) return 0;
    var pool = {}, i;
    for (i = 0; i + k <= seq.length; i++) pool[seq.substr(i, k)] = 1;
    var hits = 0, total = other.length - k + 1;
    for (i = 0; i < total; i++) if (pool[other.substr(i, k)]) hits++;
    return hits / total;
  }

  // "Flap integration without indels" is the flap written at a nick with the rest
  // of the amplicon intact - not a flap with the sequence between the nicks
  // removed, which is an integration with a deletion. Against the wild-type allele
  // that is an insertion at a nick carrying the flap's own sequence, or, when the
  // flap replaces an equal length, a run of substitutions about as long as a flap.
  function flapWrittenIn(muts, inWindow, flapSeqs, sites) {
    var found = false;
    inWindow.forEach(function (m) {
      if (found || !(m.kind === 'Ins' || m.kind === 'LargeIns') || !m.seq) return;
      flapSeqs.forEach(function (f) {
        // the share is of the FLAP's k-mers, so a short insertion cannot reach the
        // threshold however lucky it is; the length guard only makes that explicit
        if (!found && f && m.seq.length >= 0.5 * f.length
            && kmerShare(m.seq, f) >= 0.35) found = true;
      });
    });
    if (found) return true;

    // the flap replacing an equal length shows up as a dense run of substitutions
    // over the flap's footprint, which is far wider than the nick window
    if (sites.length < 2 || flapSeqs.length < 2) return false;
    var lo = Math.min(sites[0], sites[1]), hi = Math.max(sites[0], sites[1]);
    var spans = [];
    if (flapSeqs[0]) spans.push([lo + 1, lo + flapSeqs[0].length]);
    if (flapSeqs[1]) spans.push([hi - flapSeqs[1].length + 1, hi]);
    return spans.some(function (sp) {
      var len = sp[1] - sp[0] + 1;
      if (len < 8) return false;
      var subs = 0;
      muts.forEach(function (m) {
        if (m.kind === 'Sub' && m.end >= sp[0] && m.start <= sp[1]) subs += m.len;
      });
      return subs >= 0.4 * len;
    });
  }

  // How much of the donor an allele is missing, so a read assigned to it reports
  // the deletion its allele already implies rather than 0.
  function donorMissing(name, donorLen, f5, f3) {
    if (name === 'DUAL') return Math.max(0, donorLen - f5 - f3);
    if (name === 'F5') return Math.max(0, donorLen - f5);
    if (name === 'F3') return Math.max(0, donorLen - f3);
    return 0;
  }

  function classify(perRef, alleles, half, largeCut, covered, ctx) {
    ctx = ctx || { flapSeqs: [], donorLen: 0, f5: 0, f3: 0 };
    var names = Object.keys(perRef);
    if (!names.length) {
      return { category: 'Not_analysed', best: '.', pattern: '.', del: 0, ins: 0, identity: {} };
    }
    var identity = {}, burden = {};
    names.forEach(function (n) {
      identity[n] = windowIdentity(perRef[n], alleles[n].cleavage);
      burden[n] = perRef[n].reduce(function (a, m) { return a + m.len; }, 0);
    });
    function rank(n) {
      return [-Math.round((covered[n] || 0) * 100) / 100, burden[n], -identity[n]];
    }
    function less(a, b) {
      var ra = rank(a), rb = rank(b);
      for (var i = 0; i < 3; i++) { if (ra[i] !== rb[i]) return ra[i] < rb[i]; }
      return false;
    }
    var best = null;
    REF_ORDER.forEach(function (n) {
      if (!(n in perRef)) return;
      if (best === null || less(n, best)) best = n;
    });

    var a = alleles[best];
    var sites = a.cleavage.concat(a.flap, a.junction);
    var inWindow = perRef[best].filter(function (m) {
      return sites.some(function (s) { return overlaps(m, s, half); });
    });
    var del = 0, ins = 0, largeDel = false;
    inWindow.forEach(function (m) {
      if (m.kind === 'Del' || m.kind === 'LargeDel') { del += m.len; if (m.len > largeCut) largeDel = true; }
      else if (m.kind === 'Ins' || m.kind === 'LargeIns') ins += m.len;
    });

    // 20261111_pa_cat: the dual-flap, 5'- and 3'-flap alleles are the donor with its
    // body missing, so a read assigned to one of them has lost material even when it
    // matches that allele perfectly. Both flaps with everything between them gone is
    // an integration WITH a deletion, not an integration without indels.
    var missing = donorMissing(best, ctx.donorLen, ctx.f5, ctx.f3);
    var category;
    if (best === 'WT') {
      if (largeDel) category = 'Large_deletion';
      else if (flapWrittenIn(perRef[best], inWindow, ctx.flapSeqs, alleles.WT.cleavage)) {
        if (inWindow.some(function (m) { return m.kind === 'Del' || m.kind === 'LargeDel'; })) {
          category = 'Flap_integration_wDel';
        } else if (inWindow.filter(function (m) {
          return m.kind === 'Ins' || m.kind === 'LargeIns';
        }).length > 1) {
          category = 'Flap_integration_wIns';
        } else category = 'Flap_install_wo_indels';
      } else if (del || ins) category = 'Indel_at_nick_site';
      else if (!inWindow.length) category = 'Unedited';
      else category = 'Others';
    } else if (largeDel) category = 'Flap_install_with_Large_deletion';
    // material missing wins over material added: the paper's "integration +
    // insertion" rows are the complete product plus extra, never a gapped one
    else if (missing > 0 || del > ins) category = 'Flap_integration_wDel';
    else if (ins) category = 'Flap_integration_wIns';
    else if (best === 'PA') category = 'Precise_PA';
    else category = 'Others';

    return { category: category, best: best, pattern: PATTERN_OF_REF[best],
             del: del + missing, ins: ins, identity: identity };
  }

  function status(text) {
    var el = document.getElementById('paStatus');
    if (el) el.textContent = text;
    // 20261111_progress: while PA runs the results page is not up yet, so the
    // message belongs on the progress page
    if (window.lungoProgress && text) {
      window.lungoProgress.say(text);
      window.lungoProgress.detail('pa', text.replace(/^PA /, '').split(' \u2014 ')[0]);
    }
  }

  // Exposed so the pieces can be checked without running a whole analysis.
  window.__paInternals = {
    parseFlapLen: parseFlapLen, buildAlleles: buildAlleles, parseMut: parseMut,
    windowIdentity: windowIdentity, coverageFromSam: coverageFromSam,
    classify: classify, baseReadId: baseReadId, REF_ORDER: REF_ORDER,
    render: function () { return render.apply(null, arguments); }, CATEGORIES: CATEGORIES
  };

  window.runPaAnalysis = async function runPaAnalysis(opts) {
    var defer = !!(opts && opts.defer);
    status('Preparing the five PA alleles ...');

    var fragments = (window.paDonorFragments && window.paDonorFragments.length)
      ? window.paDonorFragments
      : [{ name: 'donor', seq: window.paDonorSeq }];
    var flaps = parseFlapLen(window.paFlapLen);
    var f5 = flaps[0], f3 = flaps[1];
    var built = buildAlleles(window.reference, window.cleavagePos, window.cleavagePos2,
      fragments, f5, f3);
    var alleles = built.alleles;
    if (f5 + f3 >= built.donor.length) {
      status('The two flaps together are longer than the donor.');
      return;
    }

    var paView = null;  // 20261111_alignview
    var largeCut = Number(window.largeDelLen || 50);
    var half = Number(window.paSiteWindow || 1);
    var tables = {}, coverage = {};

    for (var i = 0; i < REF_ORDER.length; i++) {
      var name = REF_ORDER[i];
      status('PA ' + (i + 1) + '/5: ' + REF_LABEL[name] + ' — aligning ...');
      if (name === 'WT') {
        if (!window.mainResult) {
          throw new Error('the main analysis did not finish, so PA has nothing to build on');
        }
        // The run already analysed this allele. Reusing its result instead of
        // repeating it with the PA window is what keeps "Unedited" here equal to
        // "WT" in the outcome table above.
        tables.WT = readsFromResult(window.mainResult, largeCut);
        coverage.WT = coverageFromSam(window.treatedAlignSam);
      } else {
        var treatedSam = await align(window.treatedFile, alleles[name].seq, 'file');
        var controlSam = (window.controlFile && window.controlFile !== false)
          ? await align(window.controlFile, alleles[name].seq, 'file') : [];
        status('PA ' + (i + 1) + '/5: ' + REF_LABEL[name] + ' — analysing ...');
        var result = await analyse(controlSam, treatedSam, alleles[name]);
        // 20261111_alignview: the precise-PA allele is worth looking at in the
        // alignment viewer - it is the one that shows whether the donor sits in
        // cleanly - so its alignment and mutation table are kept
        var tokens = (name === 'PA') ? {} : null;
        tables[name] = readsFromResult(result, largeCut, tokens);
        coverage[name] = coverageFromSam(treatedSam);
        if (name === 'PA') {
          paView = { sam: treatedSam, mutations: tokens };
        }
      }
      console.log('PA allele ' + name + ': ' + Object.keys(tables[name]).length
        + ' reads, ' + Object.keys(coverage[name]).length + ' with coverage, '
        + alleles[name].seq.length + ' bp');
    }

    status('Merging the five views ...');
    var ids = [], seen = new Set();
    REF_ORDER.forEach(function (n) {
      Object.keys(tables[n]).forEach(function (id) {
        if (!seen.has(id)) { seen.add(id); ids.push(id); }
      });
    });

    var paCtx = { flapSeqs: built.flapSeqs, donorLen: built.donor.length,
                  f5: f5, f3: f3 };  // 20261111_pa_cat
    var counts = {}, patterns = {}, chosen = {};
    var readCategory = {};  // 20261111_pa: read id -> category, for the allele plot
    CATEGORIES.forEach(function (c) { counts[c] = 0; });
    var rows = [['Read_id', 'Category', 'Best_reference', 'Flap_pattern', 'Deleted_bp', 'Inserted_bp']
      .concat(REF_ORDER.map(function (n) { return 'Identity_' + n; })).join('\t')];

    ids.forEach(function (id) {
      var perRef = {}, covered = {};
      REF_ORDER.forEach(function (n) {
        if (id in tables[n]) { perRef[n] = tables[n][id]; covered[n] = coverage[n][id] || 0; }
      });
      var r = classify(perRef, alleles, half, largeCut, covered, paCtx);  // 20261111_pa_cat
      readCategory[id] = r.category;  // 20261111_pa
      counts[r.category] = (counts[r.category] || 0) + 1;
      patterns[r.pattern] = (patterns[r.pattern] || 0) + 1;
      chosen[r.best] = (chosen[r.best] || 0) + 1;
      rows.push([id, r.category, r.best, r.pattern, r.del, r.ins]
        .concat(REF_ORDER.map(function (n) {
          return n in r.identity ? r.identity[n].toFixed(4) : '.';
        })).join('\t'));
    });

    var analysed = ids.length - counts.Not_analysed;
    window.paReadCategory = readCategory;  // 20261111_pa

    // 20261111_alignview: what the alignment viewer can be pointed at
    if (paView) {
      window.paAlignRefs = {
        PA: {
          label: 'Precise PA allele (' + alleles.PA.seq.length + ' bp)',
          seq: alleles.PA.seq,
          sam: paView.sam,
          mutations: paView.mutations,
          // named marks, so the ruler can say what each dashed line is
          marks: [].concat(
            (alleles.PA.cleavage || []).map(function (p, i) {
              return { pos: p, label: 'nick ' + (i + 1) };
            }),
            (alleles.PA.flap || []).map(function (p, i) {
              return { pos: p, label: (i === 0 ? "5' flap" : "3' flap") };
            }),
            (alleles.PA.junction || []).map(function (p, i) {
              return { pos: p, label: 'fragment ' + (i + 1) };
            }))
        }
      };
    }
    window.paResult = {
      counts: counts, patterns: patterns, chosen: chosen, analysed: analysed,
      alleles: alleles, flap: [f5, f3], rows: rows
    };
    if (!defer) render(counts, patterns, chosen, analysed, alleles, f5, f3);
    if (window.lungoProgress) {
      window.lungoProgress.done('pa', analysed.toLocaleString() + ' reads classified');
    }
    status('');
  };

  // Called once the results page is up: see showResultsPage() in CRISPRlungo.js.
  window.renderPaResult = function renderPaResult() {  // 20261111_pa
    var r = window.paResult;
    if (!r) return;
    render(r.counts, r.patterns, r.chosen, r.analysed, r.alleles, r.flap[0], r.flap[1]);
  };

  // Only builds the row -> category map, without redrawing: the allele plot has
  // not been drawn yet at this point.
  window.applyPaAlleleLabels = function applyPaAlleleLabels() {  // 20261111_pa
    labelAllelePlot(false);
  };


  // The run draws two pies from its own classification: one over every outcome
  // and one edited / not edited. In PA mode both are redrawn over the PA
  // categories, so the picture and the table tell the same story.
  function redrawPies(counts, analysed) {  // 20261111_pa
    if (typeof Plotly === 'undefined') return;
    var shown = DISPLAY_ORDER.filter(function (c) { return (counts[c] || 0) > 0; });
    var colours = ['#CECFD1', '#4357A7', '#6BBD45', '#FABF37', '#98D9ED',
                   '#BF71AE', '#ED2024', '#8c6d4a', '#000000'];
    try {
      Plotly.newPlot('allPiePlot', [{
        labels: shown.map(function (c) { return DISPLAY_NAME[c] || c.replace(/_/g, ' '); }),
        values: shown.map(function (c) { return counts[c] || 0; }),
        type: 'pie',
        marker: { colors: shown.map(function (c, i) { return c === 'Unedited' ? '#CECFD1' : colours[(i % colours.length)]; }) }
      }], { title: 'Prime assembly outcomes' });
    } catch (e) { console.error('PA outcome pie:', e); }

    var unedited = counts.Unedited || 0;
    try {
      Plotly.newPlot('mutPiePlot', [{
        labels: ['Unedited', 'Edited'],
        values: [unedited, Math.max(0, analysed - unedited)],
        type: 'pie',
        marker: { colors: ['#CECFD1', '#4357A7'] }
      }], { title: 'Edited reads (PA)' });
    } catch (e) { console.error('PA edited pie:', e); }
  }

  // Each allele-plot row is one mutation pattern, and the reads behind it may not
  // all land in the same PA category, so the row is labelled with the category
  // that holds most of them.
  var labelsApplied = false;  // 20261111_pa
  function labelAllelePlot(redraw) {  // 20261111_pa
    var byRead = window.paReadCategory, keyOf = window.lungoReadAlleleKey;
    if (!byRead || !keyOf) return;
    var tally = {};
    Object.keys(byRead).forEach(function (id) {
      var key = keyOf[id];
      if (key === undefined) return;
      (tally[key] || (tally[key] = {}))[byRead[id]] =
        ((tally[key] || {})[byRead[id]] || 0) + 1;
    });
    var label = {};
    Object.keys(tally).forEach(function (key) {
      var cats = tally[key], best = null;
      Object.keys(cats).forEach(function (c) {
        if (!best || cats[c] > cats[best]) best = c;
      });
      if (!best) return;
      label[key] = DISPLAY_NAME[best] || best.replace(/_/g, ' ');
    });
    window.paAlleleCategory = label;
    if (redraw === false) { labelsApplied = true; return; }
    // the deferred flow applies the labels before the plot is first drawn, so
    // there is nothing to redraw
    if (labelsApplied) return;
    if (typeof drawAllelePlot === 'function' && typeof d3 !== 'undefined') {
      try {
        d3.select('#allelePlot').selectAll('*').remove();
        drawAllelePlot();
      } catch (e) { console.error('allele plot relabel:', e); }
    }
  }

  function render(counts, patterns, chosen, analysed, alleles, f5, f3) {
    var pct = function (v) { return v + ' (' + (analysed ? (v * 100 / analysed).toFixed(1) : '0.0') + ' %)'; };
    var shown = DISPLAY_ORDER;

    // 20261111_pa: PA mode replaces the run's own outcome table rather than
    // adding a second one. Two count tables over the same reads, classified two
    // different ways, are only confusing, and the pie charts above belonged to
    // the one that is no longer the headline. The heading and the explanation
    // sit directly above the table, so the reader knows what they are looking at
    // before they read the numbers.
    var head = document.getElementById('paHeading');
    if (head) head.style.display = '';
    var note = document.getElementById('paHeadingNote');
    if (note) {
      note.innerHTML = 'Every read is compared with the wild-type, precise-PA, dual-flap, '
        + '5&prime;-flap and 3&prime;-flap alleles and assigned to the one that explains the most '
        + 'of it. The table, the pie charts and the allele-plot labels below are that '
        + 'classification; outcome names follow Levesque et al., Nature 2026. Flaps 5&prime; '
        + f5 + ' nt and 3&prime; ' + f3 + ' nt; on the precise-PA allele the cleavage sites are at '
        + alleles.PA.cleavage.join(' and ') + ' and the flap boundaries at '
        + alleles.PA.flap.join(' and ') + '. Percentages are of the ' + analysed.toLocaleString()
        + ' reads at least one of the five alleles could place \u2014 another allele places reads '
        + 'the wild-type analysis drops, so this is not the same total as a run without PA. '
        + 'The per-position insertion, deletion and substitution plots further down stay on the '
        + 'wild-type reference, since PA has no equivalent of them.';
    }

    var main = document.getElementById('cnt-table-container');
    if (main) {
      main.innerHTML = '';
      renderCountTable(main, [
        {
          headers: window.lungoReadCounts
            ? window.lungoReadCounts.headers
            : ['Total reads', 'Unmapped', 'Reads used in analysis'],
          values: [window.lungoReadCounts ? window.lungoReadCounts.values[0] : '-',
                   window.lungoReadCounts ? window.lungoReadCounts.values[1] : '-',
                   analysed]
        },
        {
          headers: shown.map(function (c) { return DISPLAY_NAME[c] || c.replace(/_/g, ' '); }),
          values: shown.map(function (c) { return pct(counts[c] || 0); })
        }
      ]);
    }

    redrawPies(counts, analysed);
    var tab = document.getElementById('mutPiePlotBtn');  // 20261111_pa
    if (tab) tab.textContent = 'Unedited vs edited';
    labelAllelePlot();

    var btn = document.getElementById('paDownload');
    if (btn) btn.onclick = function () {
      var blob = new Blob([window.paResult.rows.join('\n')], { type: 'text/plain' });
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'pa_read_classification.txt';
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
    };
  }
})();
