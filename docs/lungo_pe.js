// 20261128_paired_end
// Paired-end input. With the switch on, each sample takes an R1 and an R2 file;
// before alignment they are merged into one read per pair the way CRISPResso2
// does it, with fastp's overlap merge (lungo_pe_worker.js)
// and the merged FASTQ takes the place of the single file everywhere after
// that - alignment, PA, the CLI-format export - so nothing downstream has to
// know the reads were paired.

(function () {
  'use strict';

  function el(id) { return document.getElementById(id); }

  function num(id, fallback) {
    var e = el(id);
    if (!e || e.value === '') return fallback;
    var v = Number(e.value);
    return isNaN(v) ? fallback : v;
  }

  // Read from main(), before the run starts. false stops the run.
  window.applyPairedOptions = function applyPairedOptions() {
    window.lungoPairedEnabled = !!(el('peEnable') && el('peEnable').checked);
    window.lungoPeStats = null;
    window.lungoPeMerged = null;
    if (!window.lungoPairedEnabled) return true;
    window.peMinOverlap = Math.max(1, Math.round(num('peMinOverlap', 10)));
    window.peDiffLimit = Math.max(0, Math.round(num('peDiffLimit', 5)));
    window.peDiffPercent = Math.min(100, Math.max(0, num('peDiffPercent', 20)));
    if (!el('treatedFile2').files.length) {
      alert('Paired-end is on: choose the treated R2 file as well.');
      return false;
    }
    var c1 = el('controlFile').files.length, c2 = el('controlFile2').files.length;
    if (c1 !== c2) {
      alert(c1 ? 'Paired-end is on: choose the control R2 file as well.'
        : 'A control R2 was given without a control R1.');
      return false;
    }
    return true;
  };

  function baseName(file) {
    var n = String(file.name || 'reads').replace(/\.gz$/i, '').replace(/\.(fastq|fq)$/i, '');
    n = /[._-]R1(?=$|[._-])/i.test(n) ? n.replace(/[._-]R1(?=$|[._-])/i, '') : n.replace(/[._-]1$/, '');
    return n || 'reads';
  }

  function mergeOne(r1, r2, label) {
    return new Promise(function (resolve, reject) {
      var w = new Worker('lungo_pe_worker.js?v=20261128b');
      w.onmessage = function (ev) {
        var d = ev.data || {};
        if (d.type === 'progress') {
          if (window.lungoProgress) {
            window.lungoProgress.detail('merge', label + ': ' + d.done.toLocaleString()
              + ' / ' + d.total.toLocaleString() + ' pairs');
            window.lungoProgress.bar(Math.round(d.done * 100 / Math.max(1, d.total)));
          }
        } else if (d.type === 'done') {
          w.terminate();
          resolve(d);
        } else if (d.type === 'error') {
          w.terminate();
          reject(new Error(label + ': ' + d.message));
        }
      };
      w.onerror = function (e) {
        w.terminate();
        reject(new Error(label + ': ' + ((e && e.message) || 'the merge worker stopped')));
      };
      w.postMessage({
        r1: r1, r2: r2, overlapRequire: window.peMinOverlap,
        diffLimit: window.peDiffLimit, diffPercentLimit: window.peDiffPercent
      });
    });
  }

  function fail(message) {
    if (window.lungoProgress) {
      window.lungoProgress.fail('merge', 'failed');
      window.lungoProgress.say('<span class="text-danger">' + String(message).replace(/[<>]/g, '')
        + '</span><br><span class="small text-secondary">Use &ldquo;Back to input&rdquo; to change '
        + 'the files or the merge settings and try again.</span>');
    }
    var back = el('progressBackBtn');
    if (back) back.style.display = '';
    console.error('paired-end merge:', message);
  }

  // Merges every sample, then swaps the merged files in. Resolves false when the
  // run cannot go on; the progress page then already says why.
  window.lungoMergePairs = async function lungoMergePairs() {
    var jobs = [['treated', el('treatedFile').files[0], el('treatedFile2').files[0]]];
    if (el('controlFile').files.length) {
      jobs.push(['control', el('controlFile').files[0], el('controlFile2').files[0]]);
    }
    var stats = {}, merged = {};
    try {
      for (var j = 0; j < jobs.length; j++) {
        var label = jobs[j][0] === 'treated' ? 'Treated' : 'Control';
        if (window.lungoProgress) window.lungoProgress.say('Merging ' + label.toLowerCase() + ' R1 and R2 ...');
        var res = await mergeOne(jobs[j][1], jobs[j][2], label);
        if (!res.stats.merged) {
          throw new Error(label + ': none of the ' + res.stats.pairs.toLocaleString()
            + ' read pairs overlap. R1 and R2 must overlap by more than ' + window.peMinOverlap
            + ' nt, so the amplicon has to be shorter than R1 + R2 minus that.');
        }
        stats[jobs[j][0]] = res.stats;
        merged[jobs[j][0]] = new File([res.fastq], baseName(jobs[j][1]) + '_merged.fastq',
          { type: 'text/plain' });
      }
    } catch (e) {
      fail((e && e.message) || String(e));
      return false;
    }
    window.treatedFile = merged.treated;
    if (merged.control) window.controlFile = merged.control;
    window.lungoPeStats = stats;
    window.lungoPeMerged = merged;
    if (window.lungoProgress) {
      var t = stats.treated;
      window.lungoProgress.done('merge', t.merged.toLocaleString() + ' of ' + t.pairs.toLocaleString()
        + ' treated pairs merged');
      window.lungoProgress.bar(0);
    }
    return true;
  };

  function pct(a, b) { return b ? (a * 100 / b).toFixed(1) : '0.0'; }

  // One line above the outcome table: how many pairs became reads, and the
  // merged files, which are what the analysis (and the CLI command below) used.
  window.renderPairedSummary = function renderPairedSummary() {
    var box = el('peSummary');
    if (!box) return;
    box.innerHTML = '';
    var st = window.lungoPeStats;
    if (!window.lungoPairedEnabled || !st) { box.style.display = 'none'; return; }
    ['treated', 'control'].forEach(function (k) {
      var s = st[k];
      if (!s) return;
      var line = document.createElement('div');
      line.textContent = (k === 'treated' ? 'Treated' : 'Control') + ': '
        + s.pairs.toLocaleString() + ' read pairs, ' + s.merged.toLocaleString() + ' merged ('
        + pct(s.merged, s.pairs) + ' %), ' + s.notMerged.toLocaleString()
        + ' not merged and left out; merged reads ' + Math.round(s.meanLength) + ' nt on average, overlap '
        + Math.round(s.meanOverlap) + ' nt' + (s.readThrough ? ', ' + s.readThrough.toLocaleString()
          + ' shorter than the reads (adapter cut off)' : '') + '. ';
      var f = window.lungoPeMerged && window.lungoPeMerged[k];
      if (f) {
        var a = document.createElement('a');
        a.href = URL.createObjectURL(f);
        a.download = f.name;
        a.textContent = 'Download ' + f.name;
        line.appendChild(a);
      }
      box.appendChild(line);
    });
    var note = document.createElement('div');
    note.className = 'text-secondary';
    note.textContent = 'Paired-end: R1 and R2 were merged as CRISPResso2 does it (fastp --merge, '
      + 'overlap over ' + window.peMinOverlap + ' nt, at most ' + window.peDiffLimit + ' mismatches and '
      + window.peDiffPercent + ' %; pairs that do not overlap are left out). Total reads below are merged reads.';
    box.appendChild(note);
    box.style.display = '';
  };

  // Example data for paired-end mode: simulated Illumina 2 x 250 reads over a
  // 300-bp amplicon of the PD1 example locus (see data/Example_PE_*). Loaded by
  // the Examples buttons when the switch is on.
  var PE_EXAMPLE = {
    target: 'GGCGCCCTGGCCAGTCGTCT',
    reference: 'data/Example_PE_reference.fasta',
    desired: 'data/Example_PE_desired.fasta',
    files: [
      ['treatedFile', 'data/Example_PE_Treated_R1.fastq.gz'],
      ['treatedFile2', 'data/Example_PE_Treated_R2.fastq.gz'],
      ['controlFile', 'data/Example_PE_Control_R1.fastq.gz'],
      ['controlFile2', 'data/Example_PE_Control_R2.fastq.gz']
    ]
  };

  async function fetchOk(path) {
    var res = await fetch(path);
    if (!res.ok) throw new Error('Can not find ' + path.split('/').pop());
    return res;
  }

  function fastaSeq(text) {
    return text.replace(/^>.*$/gm, '').replace(/\s+/g, '').toUpperCase();
  }

  function setValue(id, value) {
    var e = el(id);
    if (!e) return;
    e.value = value;
    e.dispatchEvent(new Event('input', { bubbles: true }));
  }

  window.lungoLoadPairedExample = async function lungoLoadPairedExample() {
    try {
      setValue('targetSequence', PE_EXAMPLE.target);
      setValue('targetSequence2', '');
      setValue('reference', fastaSeq(await (await fetchOk(PE_EXAMPLE.reference)).text()));
      setValue('desiredSequence', fastaSeq(await (await fetchOk(PE_EXAMPLE.desired)).text()));
      window.lungoPeExampleLoaded = true;
      await Promise.all(PE_EXAMPLE.files.map(async function (f) {
        var bytes = new Uint8Array(await (await fetchOk(f[1])).arrayBuffer());
        var name = f[1].split('/').pop();
        // a server that decompresses .gz on the way must not leave a plain
        // FASTQ named .gz behind
        if (!(bytes[0] === 0x1f && bytes[1] === 0x8b)) name = name.replace(/\.gz$/, '');
        var dt = new DataTransfer();
        dt.items.add(new File([bytes], name, { type: 'application/octet-stream' }));
        var input = el(f[0]);
        input.files = dt.files;
        input.dispatchEvent(new Event('change', { bubbles: true }));
      }));
    } catch (e) {
      console.error('Error loading the paired-end example:', e);
      alert('Error: ' + ((e && e.message) || e));
    }
  };

  window.lungoDownloadPairedExample = async function lungoDownloadPairedExample() {
    var paths = PE_EXAMPLE.files.map(function (f) { return f[1]; })
      .concat([PE_EXAMPLE.reference, PE_EXAMPLE.desired]);
    try {
      for (var i = 0; i < paths.length; i++) {
        var blob = await (await fetchOk(paths[i])).blob();
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = paths[i].split('/').pop();
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        await new Promise(function (r) { setTimeout(r, 500); });
      }
    } catch (e) {
      console.error('Download error:', e);
      alert('Error: ' + ((e && e.message) || e));
    }
  };

  // R1 gets the same tag as R2 while the switch is on. The input is moved into
  // an input group, not replaced, so a file already chosen stays chosen.
  function tagR1(id, on) {
    var input = el(id);
    if (!input) return;
    var group = input.parentNode && input.parentNode.classList.contains('pe-r1') ? input.parentNode : null;
    if (on && !group) {
      group = document.createElement('div');
      group.className = 'input-group input-group-sm pe-r1';
      var tag = document.createElement('span');
      tag.className = 'input-group-text justify-content-center';
      tag.style.width = '2.6em';
      tag.textContent = 'R1';
      input.parentNode.insertBefore(group, input);
      group.appendChild(tag);
      group.appendChild(input);
    } else if (!on && group) {
      group.parentNode.insertBefore(input, group);
      group.parentNode.removeChild(group);
    }
  }

  // The switch sits on the Inputs heading; turning it on adds an R2 row inside
  // the treated and control file boxes, so the form does not grow a new card.
  document.addEventListener('DOMContentLoaded', function () {
    var toggle = el('peEnable');
    if (!toggle) return;
    var labels = [['treatedFile', 'Treated FASTQ'], ['controlFile', 'Control FASTQ']];
    function sync() {
      Array.prototype.forEach.call(document.querySelectorAll('.pe-r2'), function (row) {
        row.style.display = toggle.checked ? '' : 'none';
      });
      labels.forEach(function (l) {
        tagR1(l[0], toggle.checked);
        var lab = document.querySelector('label[for="' + l[0] + '"]');
        if (lab) lab.textContent = l[1] + (toggle.checked ? ' (R1 / R2)' : '');
      });
    }
    toggle.addEventListener('change', sync);
    sync();
  });
})();
