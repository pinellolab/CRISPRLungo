// 20261111_options
// Reads the Detailed options panel and the PA panel, and applies them to the
// globals the analysis already uses. Every control here is wired to something
// real: an option the wasm engine takes, a value the JS classification reads, or
// a plotting global. Options the browser cannot honour (UMI, threads, Badread
// simulation, visualization-only reruns) are deliberately absent rather than
// present and ignored.

(function () {
  'use strict';

  function num(id, fallback) {
    var el = document.getElementById(id);
    if (!el || el.value === '') return fallback;
    var v = Number(el.value);
    return isNaN(v) ? fallback : v;
  }

  function str(id, fallback) {
    var el = document.getElementById(id);
    if (!el) return fallback;
    var v = String(el.value || '').trim();
    return v === '' ? fallback : v;
  }

  function bool(id, fallback) {
    var el = document.getElementById(id);
    return el ? !!el.checked : fallback;
  }

  // Applied before the run starts, from main().
  window.applyAdvancedOptions = function applyAdvancedOptions() {
    // Mutation classification
    window.partialEditThreshold = num('inducedPartialSimilarity', 0.8);
    window.mixTag = bool('mixTag', false);

    // Background filtering against the control
    window.useAllMutations = bool('useAllMutations', false);
    window.lengthMin = num('lengthMin', 10);
    window.mutAllowance = num('mutAllowance', 0.05);

    // Allele table filtering
    window.readCountMin = num('readCountMin', 0);
    window.readFreqMin = num('readFreqMin', 0);

    // Alignment
    window.alignSaLenThreshold = num('alignSaLenThreshold', 100);
    window.minimap2Opt = str('minimap2Opt', '');

    // Allele plot
    window.plotWindow = num('allelePlotWindow', 20);
    window.numRows = num('allelePlotLines', 10);
    window.showAllBetweenAllele = bool('showAllBetweenAllele', false);
    window.removeLargeMutationsInPlot = num('removeLargeMutationsInPlot', -1);

    // Prime assembly
    window.paEnabled = bool('paEnable', false);
    window.paDonorSeq = (document.getElementById('paDonor')
      ? document.getElementById('paDonor').value.replace(/^>.*$/gm, '').replace(/\s+/g, '').toUpperCase()
      : '');
    window.paDonorFragments = (document.getElementById('paDonor')
      ? splitFasta(document.getElementById('paDonor').value)
      : []);
    window.paFlapLen = str('paFlapLen', '32');
    window.paWindow = num('paWindow', 2);
    window.paSiteWindow = num('paSiteWindow', 1);

    if (window.paEnabled && !window.paDonorSeq) {
      alert('PA mode is on but no donor sequence was given.');
      return false;
    }
    if (window.paEnabled && !document.getElementById('targetSequence2').value.trim()) {
      alert('PA mode needs two protospacers. Fill in the additional protospacer.');
      return false;
    }
    return true;
  };

  // One table shape, used by the normal outcome table and by the PA one, so the
  // two read the same. Split in two so the read counts do not have to share a
  // column grid with ten outcome columns.
  window.renderCountTable = function renderCountTable(container, blocks) {
    if (typeof container === 'string') container = document.getElementById(container);
    if (!container) return;
    blocks.forEach(function (block) {
      var table = document.createElement('table');
      table.className = 'table table-bordered table-striped w-100 result-table mb-3';
      var colgroup = document.createElement('colgroup');
      block.headers.forEach(function () {
        var col = document.createElement('col');
        col.style.width = (100 / block.headers.length) + '%';
        colgroup.appendChild(col);
      });
      table.appendChild(colgroup);
      var thead = document.createElement('thead');
      thead.className = 'thead-light';
      var hr = document.createElement('tr');
      block.headers.forEach(function (h) {
        var th = document.createElement('th');
        th.textContent = h;
        hr.appendChild(th);
      });
      thead.appendChild(hr);
      table.appendChild(thead);
      var tbody = document.createElement('tbody');
      var row = document.createElement('tr');
      block.values.forEach(function (v) {
        var td = document.createElement('td');
        td.textContent = v;
        row.appendChild(td);
      });
      tbody.appendChild(row);
      table.appendChild(tbody);
      container.appendChild(table);
    });
  };

  // 20261111_cli: the same run as a command line. Only options that differ from
  // the CLI's own defaults are written out, so what is printed is exactly what
  // changes the result. Pasted sequences become the FASTA files the button saves.
  var CLI_DEFAULT = {
    window: 5, cleavage_pos: 16, largeins_cutlen: 50, largedel_cutlen: 50,
    induced_paritial_similiarity: 0.8, length_min: 10, min_read_cnt: 0, min_read_freq: 0,
    align_sa_len_threshold: 100, allele_plot_window: 20, allele_plot_lines: 20,
    remove_large_mutations_in_plot: -1, pa_flap_len: '32', pa_window: 2, pa_site_window: 1
  };

  function cliFiles() {
    var files = [{ name: 'reference.fasta', text: '>Ref\n' + (window.reference || '') + '\n' }];
    var desired = (document.getElementById('desiredSequence') || {}).value || '';
    desired = desired.replace(/^>.*$/gm, '').replace(/\s+/g, '').toUpperCase();
    if (desired) files.push({ name: 'desired.fasta', text: '>desired\n' + desired + '\n' });
    if (window.paEnabled && window.paDonorFragments && window.paDonorFragments.length) {
      files.push({
        name: 'donor.fasta',
        text: window.paDonorFragments.map(function (r, i) {
          return '>' + (r.name || ('fragment' + (i + 1))) + '\n' + r.seq + '\n';
        }).join('')
      });
    }
    return files;
  }
  window.cliInputFiles = cliFiles;

  window.buildCliCommand = function buildCliCommand() {
    var treated = (document.getElementById('treatedFile').files[0] || {}).name || 'treated.fastq';
    var controlEl = document.getElementById('controlFile');
    var control = controlEl.files.length ? controlEl.files[0].name : null;
    var target = (document.getElementById('targetSequence').value || '').trim().toUpperCase();
    var target2 = (document.getElementById('targetSequence2').value || '').trim().toUpperCase();
    var desired = ((document.getElementById('desiredSequence') || {}).value || '')
      .replace(/^>.*$/gm, '').replace(/\s+/g, '');

    var parts = ['CRISPRlungo reference.fasta ' + treated + ' crispr_output ' + target];
    var add = function (flag, value, def) {
      if (value === undefined || value === null || value === '' || value === def) return;
      parts.push(flag + ' ' + value);
    };
    if (target2) parts.push('--additional_target ' + target2);
    if (control) parts.push('--control ' + control);
    if (window.wholeWindow === 1) parts.push('--whole_window_between_targets');
    add('--window', Number(document.getElementById('analysisWindow').value), CLI_DEFAULT.window);
    add('--cleavage_pos', Number(document.getElementById('cleavageSite').value), CLI_DEFAULT.cleavage_pos);
    add('--largeins_cutlen', Number(document.getElementById('largeInsertionLen').value), CLI_DEFAULT.largeins_cutlen);
    add('--largedel_cutlen', Number(document.getElementById('largeDeletionLen').value), CLI_DEFAULT.largedel_cutlen);
    if (desired) parts.push('--induced_sequence_path desired.fasta');
    add('--induced_paritial_similiarity', window.partialEditThreshold, CLI_DEFAULT.induced_paritial_similiarity);
    if (window.mixTag) parts.push('--mix_tag');
    if (window.useAllMutations) parts.push('--using_all_mutations');
    add('--length_min', window.lengthMin, CLI_DEFAULT.length_min);
    add('--min_read_cnt', window.readCountMin, CLI_DEFAULT.min_read_cnt);
    add('--min_read_freq', window.readFreqMin, CLI_DEFAULT.min_read_freq);
    add('--align_sa_len_threshold', window.alignSaLenThreshold, CLI_DEFAULT.align_sa_len_threshold);
    if (window.minimap2Opt) parts.push("--minimap2_opt '" + window.minimap2Opt + "'");
    add('--allele_plot_window', window.plotWindow, CLI_DEFAULT.allele_plot_window);
    add('--allele_plot_lines', window.numRows, CLI_DEFAULT.allele_plot_lines);
    add('--remove_large_mutations_in_plot', window.removeLargeMutationsInPlot, CLI_DEFAULT.remove_large_mutations_in_plot);
    if (window.showAllBetweenAllele) parts.push('--show_all_between_allele');
    if (window.paEnabled) {
      parts.push('--pa_donor donor.fasta');
      add('--pa_flap_len', window.paFlapLen, CLI_DEFAULT.pa_flap_len);
      add('--pa_window', window.paWindow, CLI_DEFAULT.pa_window);
      add('--pa_site_window', window.paSiteWindow, CLI_DEFAULT.pa_site_window);
    }

    var notes = [];
    if (typeof window.mutAllowance === 'number' && window.mutAllowance !== 0.05) {
      notes.push('# BH significance alpha ' + window.mutAllowance
        + ' has no command-line option; the CLI always uses 0.05.');
    }
    return notes.concat([parts.join(' \\\n    ')]).join('\n');
  };

  window.renderCliCommand = function renderCliCommand() {
    var section = document.getElementById('cliCommandSection');
    var pre = document.getElementById('cliCommandText');
    if (!section || !pre) return;
    section.style.display = '';
    pre.textContent = window.buildCliCommand();
    var files = cliFiles();
    var list = document.getElementById('cliCommandFiles');
    if (list) {
      list.textContent = 'Run it beside ' + files.map(function (f) { return f.name; }).join(', ')
        + ' and the FASTQ file' + (document.getElementById('controlFile').files.length ? 's' : '') + '.';
    }
    var copy = document.getElementById('cliCommandCopy');
    if (copy) copy.onclick = function () {
      navigator.clipboard.writeText(pre.textContent);
      copy.textContent = 'Copied';
      setTimeout(function () { copy.textContent = 'Copy'; }, 1200);
    };
    var dl = document.getElementById('cliCommandDownload');
    if (dl) dl.onclick = function () {
      files.forEach(function (f, i) {
        setTimeout(function () {
          var a = document.createElement('a');
          a.href = URL.createObjectURL(new Blob([f.text], { type: 'text/plain' }));
          a.download = f.name;
          document.body.appendChild(a); a.click(); document.body.removeChild(a);
        }, i * 300);
      });
    };
  };

  // A donor pasted as multi-record FASTA is an assembly of fragments, in order.
  function splitFasta(text) {
    var out = [], cur = null;
    String(text || '').split(/\r?\n/).forEach(function (line) {
      line = line.trim();
      if (!line) return;
      if (line.charAt(0) === '>') { cur = { name: line.slice(1).split(/\s+/)[0], seq: '' }; out.push(cur); }
      else {
        if (!cur) { cur = { name: 'donor', seq: '' }; out.push(cur); }
        cur.seq += line.replace(/\s+/g, '').toUpperCase();
      }
    });
    return out.filter(function (r) { return r.seq.length; });
  }
  window.splitFastaRecords = splitFasta;

  // The PA card is a button, not part of the Detailed options collapse, because
  // it is reached far more often than the rest.
  document.addEventListener('DOMContentLoaded', function () {
    var toggle = document.getElementById('paEnable');
    var panel = document.getElementById('paPanel');
    var badge = document.getElementById('paStateBadge');
    if (!toggle || !panel) return;
    function sync() {
      panel.style.display = toggle.checked ? '' : 'none';
      if (badge) {
        badge.textContent = toggle.checked ? 'on' : 'off';
        badge.className = 'badge rounded-pill ' + (toggle.checked ? 'text-bg-primary' : 'text-bg-secondary');
      }
    }
    toggle.addEventListener('change', sync);
    sync();
  });
})();
