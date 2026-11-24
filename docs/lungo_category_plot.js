// 20261124_category_plot
// Top alleles per category, chosen from a menu on the allele plot.
//
// The allele plot shows the most frequent patterns of the whole run, so a rare
// class - a handful of partial edits, the few large insertions - never reaches
// it. This draws the same plot once per outcome category, from only the reads
// of that category, so each class shows its own top alleles: what the precise
// edits look like, what the partial ones carry besides the edit, which
// deletions and insertions occur.
//
// The categories are the ones the outcome table counts. A read can sit in two
// of them, as it does in the table: "Partially induced editing" is a flag on top
// of the read's mutation type (Ins, Del, ...), not a class of its own. With PA
// mode on, the PA outcomes come first; their alleles are still drawn on the
// wild-type reference, so reads only another allele could place are counted but
// cannot be drawn.

(function () {
  'use strict';

  var STANDARD = [
    ['Precise', 'Precisely induced editing'],
    ['Partial', 'Partially induced editing'],
    ['Ins', 'Ins'], ['Del', 'Del'], ['Sub', 'Sub'],
    ['LargeIns', 'LargeIns'], ['LargeDel', 'LargeDel'],
    ['Inv', 'Inversion'], ['Complex', 'Complex']
  ];
  var PA_FALLBACK_ORDER = ['Unedited', 'Precise_PA', 'Flap_install_wo_indels',
    'Flap_integration_wDel', 'Flap_integration_wIns', 'Flap_install_with_Large_deletion',
    'Large_deletion', 'Indel_at_nick_site', 'Others'];

  var state = { cats: [] };

  function paName(c) {
    var pi = window.__paInternals || {};
    var names = pi.DISPLAY_NAME || { Unedited: 'WT' };
    return names[c] || c.replace(/_/g, ' ');
  }

  function add(map, cat, key) {
    var b = map[cat] || (map[cat] = { reads: 0, alleles: {} });
    b.reads++;
    if (key !== undefined) b.alleles[key] = (b.alleles[key] || 0) + 1;
  }

  // Most frequent first; ties keep the order of the run's own allele table.
  function sortAlleles(alleles) {
    var order = {};
    Object.keys(window.filteredTreatedReads || {}).forEach(function (k, i) { order[k] = i; });
    return Object.keys(alleles).sort(function (a, b) {
      return (alleles[b] - alleles[a]) || ((order[a] || 0) - (order[b] || 0));
    }).map(function (k) { return [k, alleles[k]]; });
  }

  function collect() {
    var keyOf = window.lungoReadAlleleKey || {};
    var typeOf = window.lungoReadMutType || {};
    var cats = [];

    var paCat = window.paReadCategory;
    if (paCat && window.paResult) {
      var pa = {};
      Object.keys(paCat).forEach(function (id) { add(pa, paCat[id], keyOf[id]); });
      var order = (window.__paInternals && window.__paInternals.DISPLAY_ORDER) || PA_FALLBACK_ORDER;
      order.forEach(function (c) {
        if (!pa[c] || !pa[c].reads) return;
        cats.push({ id: 'pa:' + c, group: 'pa', name: paName(c), reads: pa[c].reads,
                    total: window.paResult.analysed, alleles: sortAlleles(pa[c].alleles) });
      });
    }

    var std = {};
    Object.keys(keyOf).forEach(function (id) {
      var key = keyOf[id];
      if (/_PreciseEdit$/.test(key)) add(std, 'Precise', key);
      else if (/_PartialEdit$/.test(key)) add(std, 'Partial', key);
      var t = typeOf[id];
      if (t && t !== 'WT') add(std, t, key);
    });
    var used = window.mainResult && window.mainResult.treated_align_cnt
      ? window.mainResult.treated_align_cnt.Used : 0;
    STANDARD.forEach(function (c) {
      if (!std[c[0]] || !std[c[0]].reads) return;
      cats.push({ id: 'std:' + c[0], group: 'std', name: c[1], reads: std[c[0]].reads,
                  total: used, alleles: sortAlleles(std[c[0]].alleles) });
    });
    return cats;
  }

  function plotted(c) {
    return c.alleles.reduce(function (s, a) { return s + a[1]; }, 0);
  }

  function count(n, word) {
    return n.toLocaleString() + ' ' + word + (n === 1 ? '' : 's');
  }

  function caption(c, rows) {
    var pct = c.total ? ' (' + (c.reads * 100 / c.total).toFixed(1) + ' %)' : '';
    var n = c.alleles.length;
    var txt = c.name + ' \u2014 ' + count(c.reads, 'read') + pct + ', ' + count(n, 'allele');
    if (n > rows) txt += ', top ' + rows + ' shown';
    var missing = c.reads - plotted(c);
    if (missing > 0) {
      txt += '; ' + count(missing, 'read') + (missing === 1 ? ' has' : ' have')
        + ' no wild-type alignment and cannot be drawn';
    }
    return txt;
  }

  function plotData(c) {
    var out = {};
    c.alleles.forEach(function (a) { out[a[0]] = a[1]; });
    return out;
  }

  // PA rows carry the PA label in the main plot; inside a PA category every read
  // of the row is already that category, so the label would only repeat it.
  function plotOpts(c, rows) {
    return c.group === 'pa' ? { rows: rows, labels: {} } : { rows: rows };
  }

  // drawAllelePlot leaves the viewBox much taller than the rows, which is fine for
  // one plot but makes a stack of them mostly white space; keep the width and
  // trim the height to what was drawn.
  function trimSvg(box, top) {
    var svg = box.querySelector('svg');
    var g = svg && svg.querySelector('g.main-container');
    if (!g || typeof d3 === 'undefined') return svg;
    var bb;
    try { bb = g.getBBox(); } catch (e) { return svg; }
    if (!bb || !bb.height) return svg;
    var vb = svg.getAttribute('viewBox').split(/\s+/).map(Number);
    var t = d3.zoomTransform(svg);
    var pad = 8 + (top || 0);
    var y0 = t.y + t.k * bb.y - pad;
    var h = t.k * bb.height + pad + 8;
    // the PA label after the read count can run past the right edge
    var right = Math.max(vb[0] + vb[2], t.x + t.k * (bb.x + bb.width) + 8);
    svg.setAttribute('viewBox', vb[0] + ' ' + y0 + ' ' + (right - vb[0]) + ' ' + h);
    return svg;
  }

  // d3 sizes the plot from its container's width, which is not there while the
  // container or one of its parents is hidden; give the box a width of its own.
  function boxWidth() {
    var main = document.getElementById('allelePlot');
    var w = main && main.parentNode ? main.parentNode.clientWidth : 0;
    return w > 0 ? w : 960;
  }

  function rowsWanted() {
    var el = document.getElementById('categoryPlotRows');
    return Math.max(1, parseInt(el && el.value, 10) || 5);
  }

  // "All reads" is the run's own allele plot, untouched; a category swaps it for
  // that category's plot in the same place, and back.
  function draw() {
    var host = document.getElementById('categoryAllelePlot');
    var main = document.getElementById('allelePlot');
    var pick = document.getElementById('categoryPlotPick');
    var rowsEl = document.getElementById('categoryPlotRows');
    var note = document.getElementById('categoryPlotNote');
    if (!host || !main || !pick || typeof drawAllelePlot !== 'function') return;
    host.innerHTML = '';
    var c = null;
    state.cats.forEach(function (x) { if (x.id === pick.value) c = x; });
    if (!c) {
      host.style.display = 'none';
      main.style.display = '';
      if (rowsEl) rowsEl.style.display = 'none';
      if (note) note.style.display = 'none';
      return;
    }
    var rows = rowsWanted();
    main.style.display = 'none';
    host.style.display = '';
    if (rowsEl) rowsEl.style.display = '';
    if (note) {
      note.innerHTML = '';
      var b = document.createElement('span');
      b.className = 'fw-semibold text-body';
      b.textContent = caption(c, rows) + '. ';
      note.appendChild(b);
      note.appendChild(document.createTextNode('Only the reads of this category; row percentages '
        + 'are of the category. A read can count in two categories, as in the table: partially '
        + 'induced editing is counted on top of the mutation type.'
        + (c.group === 'pa' ? ' Drawn on the wild-type reference.' : '')
        + ' Choose “All reads” to return to the full plot.'));
      note.style.display = '';
    }
    if (!c.alleles.length) return;
    var box = document.createElement('div');
    box.id = 'categoryAllelePlot_0';
    box.style.width = boxWidth() + 'px';
    host.appendChild(box);
    try {
      drawAllelePlot(plotData(c), '#' + box.id, plotOpts(c, rows));
      trimSvg(box);
    } catch (e) { console.error('category allele plot', c.id, e); }
  }

  function fillPicker() {
    var pick = document.getElementById('categoryPlotPick');
    pick.innerHTML = '';
    var all = document.createElement('option');
    all.value = '';
    all.textContent = 'All reads';
    pick.appendChild(all);
    var groups = [['pa', 'PA outcome'], ['std', 'Category']];
    groups.forEach(function (g) {
      var members = state.cats.filter(function (c) { return c.group === g[0]; });
      if (!members.length) return;
      var og = document.createElement('optgroup');
      og.label = g[1];
      members.forEach(function (c) {
        var o = document.createElement('option');
        o.value = c.id;
        o.textContent = c.name + ' (' + c.reads.toLocaleString() + ')';
        og.appendChild(o);
      });
      pick.appendChild(og);
    });
    pick.value = '';
  }

  // Each run opens on "All reads", the plot as it always was.
  window.renderCategoryAllelePlots = function renderCategoryAllelePlots() {
    var sec = document.getElementById('categoryPlotSection');
    var pick = document.getElementById('categoryPlotPick');
    if (!sec || !pick) return;
    state.cats = collect();
    fillPicker();
    draw();
    sec.style.display = state.cats.length ? '' : 'none';
    if (!pick.dataset.bound) {
      pick.addEventListener('change', draw);
      var rowsEl = document.getElementById('categoryPlotRows');
      if (rowsEl) rowsEl.addEventListener('change', draw);
      pick.dataset.bound = '1';
    }
  };

  function safeName(s) {
    return String(s).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
  }

  // For the result ZIP: one SVG per category, drawn off-screen at the width of
  // the main plot, and a table of the same top alleles.
  window.lungoCategoryPlotExport = function lungoCategoryPlotExport() {
    var cats = collect();
    if (!cats.length || typeof drawAllelePlot !== 'function') return [];
    var rows = rowsWanted();
    var width = boxWidth();
    var tmp = document.createElement('div');
    tmp.style.cssText = 'position:absolute;left:-100000px;top:0;visibility:hidden;width:' + width + 'px;';
    document.body.appendChild(tmp);
    var files = [];
    var table = ['Group\tCategory\tCategory_reads\tRank\tReads\tPct_of_category\tAllele'];
    try {
      cats.forEach(function (c, i) {
        var group = c.group === 'pa' ? 'PA' : 'Standard';
        c.alleles.slice(0, rows).forEach(function (a, r) {
          table.push([group, c.name, c.reads, r + 1, a[1], (a[1] * 100 / c.reads).toFixed(2), a[0]].join('\t'));
        });
        if (!c.alleles.length) return;
        var box = document.createElement('div');
        box.id = 'lungoCategoryExport_' + i;
        box.style.width = width + 'px';
        tmp.appendChild(box);
        drawAllelePlot(plotData(c), '#' + box.id, plotOpts(c, rows));
        var svg = trimSvg(box, 24);
        if (!svg) return;
        var clone = svg.cloneNode(true);
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        var t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        t.setAttribute('x', '10');
        t.setAttribute('y', String(Number(svg.getAttribute('viewBox').split(/\s+/)[1]) + 18));
        t.setAttribute('font-size', '14');
        t.setAttribute('font-family', 'sans-serif');
        t.setAttribute('font-weight', 'bold');
        t.textContent = caption(c, rows);
        clone.insertBefore(t, clone.firstChild);
        var num = (i + 1 < 10 ? '0' : '') + (i + 1);
        files.push([num + '_' + group + '_' + safeName(c.name) + '.svg',
          '<?xml version="1.0" encoding="UTF-8"?>\n' + clone.outerHTML]);
      });
    } finally {
      document.body.removeChild(tmp);
    }
    files.push(['top_alleles_by_category.txt', table.join('\n') + '\n']);
    return files;
  };

  window.__categoryPlotInternals = { collect: collect, caption: caption };
})();
