/* ============================================================================
 * lungo_agent.js
 * ----------------------------------------------------------------------------
 * "Ask about your results": a Claude agent that answers questions about the
 * current CRISPRLungo web analysis (e.g. "a large insertion occurred - what was
 * inserted?"). Claude calls the tools below, which run here in the browser on
 * the in-memory results, and answers from what they return.
 *
 * Tools: get_run_overview, list_alleles, get_insertions, map_sequence,
 *        get_reference_region, list_reads, get_custom_classification
 *
 * Depends on lungo_export.js (CLI-style read rows) and globals from main.js /
 * CRISPRlungo_ai.js / custom_classify.js.
 * ==========================================================================*/

var agentHistory = [];
var agentCache = { key: null, built: null };
var agentKmerIndex = {};
var AGENT_MAX_STEPS = 12;
var AGENT_K = 11;

/* ------------------------------ data access ------------------------------ */

async function agentData() {
  if (agentCache.key !== mainResult || !agentCache.built) {
    agentCache.built = await lxBuildFiles();
    agentCache.key = mainResult;
  }
  return agentCache.built;
}

function agentRevComp(s) {
  var c = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
  var out = '';
  for (var n = s.length - 1; n >= 0; n--) out += c[s[n]] || s[n];
  return out;
}

function agentShortSeq(s, max) {
  max = max || 80;
  return s.length <= max ? s : s.slice(0, max) + '...(' + s.length + ' bp)';
}

function agentShortInfo(info) {
  return info.replace(/(Ins_\d+_)([ACGTN]{40,})/g, function (m, p, seq) { return p + agentShortSeq(seq, 40); });
}

// Known elements searched by the web insertion step (parsed from confirmInsertion).
function agentKnownElements() {
  var out = [{ name: 'Reference', seq: reference.toUpperCase() }];
  try {
    var src = String(confirmInsertion);
    var m = src.match(/insReference\s*=\s*'[^']*'\s*\+\s*reference\s*\+\s*'([^']*)'/);
    if (m) {
      var parts = m[1].split('\\n');
      for (var n = 0; n + 1 < parts.length; n++) {
        if (parts[n][0] !== '>') continue;
        var name = parts[n].slice(1).trim();
        if (!name || name === 'ref' || !parts[n + 1]) continue;
        out.push({ name: name, seq: parts[n + 1].toUpperCase() });
      }
    }
  } catch (e) { }
  return out;
}

function agentIndex(t) {
  if (agentKmerIndex[t.name] && agentKmerIndex[t.name].len === t.seq.length) return agentKmerIndex[t.name].idx;
  var idx = new Map();
  for (var p = 0; p + AGENT_K <= t.seq.length; p++) {
    var k = t.seq.substr(p, AGENT_K);
    if (k.indexOf('N') !== -1) continue;
    var arr = idx.get(k);
    if (!arr) { arr = []; idx.set(k, arr); }
    arr.push(p);
  }
  agentKmerIndex[t.name] = { len: t.seq.length, idx: idx };
  return idx;
}

// k-mer seed + diagonal clustering; coordinates are 0-based, end exclusive.
function agentMapSequence(seq, maxHits) {
  seq = String(seq || '').toUpperCase().replace(/[^ACGTN]/g, '');
  maxHits = maxHits || 5;
  var res = { query_length: seq.length, gc_percent: 0, longest_homopolymer: '', hits: [] };
  if (!seq.length) return res;
  var gc = (seq.match(/[GC]/g) || []).length;
  res.gc_percent = Math.round(gc * 1000 / seq.length) / 10;
  var hp = seq.match(/(A+|C+|G+|T+)/g) || [];
  var longest = hp.reduce(function (a, b) { return b.length > a.length ? b : a; }, '');
  res.longest_homopolymer = longest[0] + ' x ' + longest.length;
  if (seq.length < AGENT_K) return res;
  var nK = seq.length - AGENT_K + 1;
  var targets = agentKnownElements();
  var clusters = [];
  targets.forEach(function (t) {
    var idx = agentIndex(t);
    [['+', seq], ['-', agentRevComp(seq)]].forEach(function (sq) {
      var bins = new Map();
      for (var i = 0; i < nK; i++) {
        var pos = idx.get(sq[1].substr(i, AGENT_K));
        if (!pos || pos.length > 20) continue;
        pos.forEach(function (p) {
          var b = Math.round((p - i) / 25);
          var c = bins.get(b);
          if (!c) { c = { qs: i, qe: i + AGENT_K, ts: p, te: p + AGENT_K, n: 0 }; bins.set(b, c); }
          c.n += 1;
          c.qs = Math.min(c.qs, i); c.qe = Math.max(c.qe, i + AGENT_K);
          c.ts = Math.min(c.ts, p); c.te = Math.max(c.te, p + AGENT_K);
        });
      }
      bins.forEach(function (c) {
        if (c.n < 3 || c.qe - c.qs < 25) return;
        var qs = c.qs, qe = c.qe;
        if (sq[0] === '-') { qs = seq.length - c.qe; qe = seq.length - c.qs; }
        clusters.push({
          target: t.name, strand: sq[0], query_start: qs, query_end: qe,
          target_start: c.ts, target_end: c.te,
          query_coverage_percent: Math.round((qe - qs) * 1000 / seq.length) / 10,
          seed_identity_percent: Math.min(100, Math.round(c.n * 1000 / Math.max(1, (c.qe - c.qs - AGENT_K + 1))) / 10)
        });
      });
    });
  });
  clusters.sort(function (a, b) { return (b.query_end - b.query_start) - (a.query_end - a.query_start); });
  res.hits = clusters.slice(0, maxHits);
  if (!res.hits.length) res.note = 'No ' + AGENT_K + '-mer match to the reference or the known elements (' + targets.map(function (t) { return t.name; }).join(', ') + ').';
  return res;
}

/* --------------------------------- tools --------------------------------- */

var AGENT_TOOLS = [
  {
    name: 'get_run_overview',
    description: 'Run settings and read counts of the current analysis: targets, cut sites, window, reference length, desired (induced) mutations, preprocessing counts, and read counts per classification (WT, Del, Ins, Sub, LargeDel, LargeIns, Inv, Complex, ComplexWithLargeMut) and per induced-edit status (Precise, Precise_with_mutations, Partial, Partial_with_mutations, X). Call this first.',
    input_schema: { type: 'object', properties: {} }
  },
  {
    name: 'list_alleles',
    description: 'Allele table (reads grouped by identical in-window mutation pattern), most frequent first. Mutation tokens are START_END:TYPE_DETAIL, e.g. 1157_1158:Ins_1_A, 1157_1160:Del_4, 1160_1160:Sub_1_T>C, LargeDel/LargeIns past the size cutoff. Long inserted sequences are shortened; use get_insertions for full sequences.',
    input_schema: {
      type: 'object',
      properties: {
        classification: { type: 'string', description: 'Only this classification (e.g. LargeIns, Del, Inv, WT).' },
        induced_status: { type: 'string', description: 'Only this induced-edit status (Precise, Precise_with_mutations, Partial, Partial_with_mutations, X).' },
        contains: { type: 'string', description: 'Only patterns whose mutation string contains this text.' },
        min_reads: { type: 'integer' },
        limit: { type: 'integer', description: 'Default 30, max 200.' }
      }
    }
  },
  {
    name: 'get_insertions',
    description: 'Inserted sequences found in reads (from the in-window mutations), grouped by position + sequence with read counts, the full inserted sequence, the web insertion-search annotation, and an automatic k-mer mapping of the insert against the reference and known elements (nCas9, sgRNA scaffold, DC_cluster, AAV ITR).',
    input_schema: {
      type: 'object',
      properties: {
        min_length: { type: 'integer', description: 'Minimum insertion length in bp. Default 20.' },
        limit: { type: 'integer', description: 'Default 15, max 100.' }
      }
    }
  },
  {
    name: 'map_sequence',
    description: 'Map any DNA sequence against the reference amplicon and the known elements with ' + AGENT_K + '-mer seeds on both strands. Returns hit regions (0-based, end exclusive), strand, query coverage, seed identity, GC% and the longest homopolymer. No hit means the sequence is not from the amplicon or those elements.',
    input_schema: { type: 'object', properties: { sequence: { type: 'string' } }, required: ['sequence'] }
  },
  {
    name: 'get_reference_region',
    description: 'Reference sequence between 0-based start (inclusive) and end (exclusive), max 3000 bp, with cut-site positions marked in the response.',
    input_schema: { type: 'object', properties: { start: { type: 'integer' }, end: { type: 'integer' } }, required: ['start', 'end'] }
  },
  {
    name: 'list_reads',
    description: 'Individual reads (id, classification, mutation pattern, insertion annotation, induced status, other mutations outside the window or not significant). Filter by text in the mutation pattern or by classification.',
    input_schema: {
      type: 'object',
      properties: {
        contains: { type: 'string' },
        classification: { type: 'string' },
        limit: { type: 'integer', description: 'Default 10, max 50.' }
      }
    }
  },
  {
    name: 'get_custom_classification',
    description: 'Category counts from the AI custom classification, if one has been run on this page.',
    input_schema: { type: 'object', properties: {} }
  }
];

async function agentRunTool(name, input) {
  input = input || {};
  var d = await agentData();
  if (name === 'get_run_overview') {
    var cv2 = (cleavagePos2 === '' || cleavagePos2 == null) ? null : Number(cleavagePos2);
    return {
      target_1: target, target_2: target2 || null,
      cleavage_position_1: cleavagePos, cleavage_position_2: cv2,
      window: windowRange, whole_window_between_targets: cv2 === null ? null : wholeWindow === 1,
      reference_length: reference.length,
      large_insertion_cutoff_bp: largeInsLen, large_deletion_cutoff_bp: largeDelLen,
      desired_mutations: lxInducedMutStr(),
      control_sample: !!controlFile,
      preprocess_count: d.files['results/preprocess_count.txt'],
      used_reads: d.rows.length,
      classification_counts: d.patterns.mutCnt,
      induced_status_counts: d.patterns.preciseCnt,
      distinct_alleles: d.patterns.table.length,
      coordinate_note: 'Positions are 0-based reference coordinates as written by CRISPRLungo.'
    };
  }
  if (name === 'list_alleles') {
    var lim = Math.min(Number(input.limit) || 30, 200);
    var t = d.patterns.table.filter(function (r) {
      if (input.classification && r.cls !== input.classification) return false;
      if (input.induced_status && r.induce !== input.induced_status) return false;
      if (input.contains && r.info.indexOf(input.contains) === -1) return false;
      if (input.min_reads && r.cnt < input.min_reads) return false;
      return true;
    });
    return {
      total_matching_patterns: t.length,
      total_matching_reads: t.reduce(function (a, r) { return a + r.cnt; }, 0),
      alleles: t.slice(0, lim).map(function (r) {
        return { classification: r.cls, reads: r.cnt, percent: r.pct, induced_status: r.induce, insertion_annotation: r.integ || null, mutations: agentShortInfo(r.info) };
      })
    };
  }
  if (name === 'get_insertions') {
    var minLen = Number(input.min_length) || 20;
    var lim2 = Math.min(Number(input.limit) || 15, 100);
    var groups = new Map();
    d.rows.forEach(function (r) {
      if (r.info === 'None') return;
      r.info.split(',').forEach(function (tok) {
        var m = tok.match(/^(\d+)_(\d+):(Large)?Ins_(\d+)_([A-Z]+)/);
        if (!m || Number(m[4]) < minLen) return;
        var key = m[1] + '|' + m[5];
        var g = groups.get(key);
        if (!g) { g = { position: Number(m[1]), length: Number(m[4]), sequence: m[5], reads: 0, read_classes: {}, annotations: {}, example_pattern: r.info }; groups.set(key, g); }
        g.reads += 1;
        g.read_classes[r.cls] = (g.read_classes[r.cls] || 0) + 1;
        if (r.integ) g.annotations[r.integ.replace(/,$/, '')] = 1;
      });
    });
    var list = Array.from(groups.values()).sort(function (a, b) { return b.reads - a.reads || b.length - a.length; });
    return {
      distinct_insertions: list.length,
      reads_with_insertion: list.reduce(function (a, g) { return a + g.reads; }, 0),
      insertions: list.slice(0, lim2).map(function (g) {
        return {
          position: g.position, length: g.length, reads: g.reads, read_classifications: g.read_classes,
          sequence: g.sequence.length > 3000 ? g.sequence.slice(0, 3000) + '...(truncated)' : g.sequence,
          web_insertion_search: Object.keys(g.annotations),
          kmer_mapping: agentMapSequence(g.sequence, 4),
          example_pattern: agentShortInfo(g.example_pattern)
        };
      })
    };
  }
  if (name === 'map_sequence') return agentMapSequence(input.sequence, 8);
  if (name === 'get_reference_region') {
    var st = Math.max(0, Number(input.start) || 0);
    var ed = Math.min(reference.length, Number(input.end) || st + 100, st + 3000);
    var cuts = [cleavagePos];
    if (!(cleavagePos2 === '' || cleavagePos2 == null)) cuts.push(Number(cleavagePos2));
    return { start: st, end: ed, sequence: reference.slice(st, ed), cut_sites_in_region: cuts.filter(function (c) { return c >= st && c < ed; }) };
  }
  if (name === 'list_reads') {
    var lim3 = Math.min(Number(input.limit) || 10, 50);
    var rs = d.rows.filter(function (r) {
      if (input.classification && r.cls !== input.classification) return false;
      if (input.contains && r.info.indexOf(input.contains) === -1) return false;
      return true;
    });
    return {
      matching_reads: rs.length,
      reads: rs.slice(0, lim3).map(function (r) {
        return { id: r.id, classification: r.cls, mutations: agentShortInfo(r.info), insertion_annotation: r.integ || null, induced_status: r.induce, other_mutations: agentShortSeq(r.whole, 400) };
      })
    };
  }
  if (name === 'get_custom_classification') {
    if (typeof customCategoryCounts === 'undefined' || !Object.keys(customCategoryCounts).length) return { available: false };
    return { available: true, request: (customClassifierMeta || {}).request || null, category_reads: customCategoryCounts };
  }
  return { error: 'Unknown tool ' + name };
}

/* ------------------------------- prompting ------------------------------- */

function agentSystemPrompt() {
  return [
    'You are a genome-editing data analyst embedded in the CRISPRLungo web tool, which analyses long-read',
    '(Nanopore/PacBio) amplicon sequencing of CRISPR-edited samples. The user is looking at the results of one',
    'analysis and asks you questions about it.',
    '',
    'How to work:',
    '- Get facts from the tools; do not guess numbers. Start with get_run_overview when you need context.',
    '- For insertion questions, call get_insertions and read its kmer_mapping and web_insertion_search fields.',
    '  A hit to Reference on + strand from elsewhere in the amplicon suggests a duplication / template switch;',
    '  a - strand hit suggests an inverted segment; hits to nCas9, Scaffold (sgRNA) or ITR point to vector or',
    '  plasmid DNA; a long homopolymer or low complexity suggests a sequencing or PCR artefact; no hit means the',
    '  insert is foreign to the amplicon - say what it resembles if you recognise it, state the confidence, and',
    '  suggest BLAST (show the sequence) for confirmation.',
    '- Positions are 0-based reference coordinates. Relate events to the cut site(s).',
    '- Separate what the data shows from your interpretation, and say how many reads support each point.',
    '- Answer in the language the user writes in. Be concise; use short lists and small tables when helpful.'
  ].join('\n');
}

/* -------------------------------- API loop ------------------------------- */

async function agentApi(apiKey, model, messages) {
  var response = await fetch('https://api.anthropic.com/v1/messages', {
    method: 'POST',
    headers: {
      'content-type': 'application/json',
      'x-api-key': apiKey,
      'anthropic-version': '2023-06-01',
      'anthropic-dangerous-direct-browser-access': 'true'
    },
    body: JSON.stringify({
      model: model,
      max_tokens: 8192,
      system: agentSystemPrompt(),
      tools: AGENT_TOOLS,
      messages: messages
    })
  });
  if (!response.ok) {
    var detail = '';
    try { detail = (await response.json()).error.message; } catch (e) { detail = await response.text(); }
    throw new Error('API error ' + response.status + ': ' + detail);
  }
  return await response.json();
}

async function askAgent() {
  var statusEl = document.getElementById('agentStatus');
  var btn = document.getElementById('agentAskBtn');
  var qEl = document.getElementById('agentQuestion');
  var apiKey = document.getElementById('anthropicApiKey').value.trim();
  var model = document.getElementById('claudeModel').value.trim();
  var question = qEl.value.trim();

  if (!mainResult) { statusEl.textContent = 'Run an analysis first.'; return; }
  if (!apiKey) { statusEl.textContent = 'Please enter your Anthropic API key.'; return; }
  if (!question) { statusEl.textContent = 'Please type a question.'; return; }

  btn.disabled = true;
  agentRender('user', question);
  qEl.value = '';
  agentHistory.push({ role: 'user', content: question });
  var turnStart = agentHistory.length - 1;

  try {
    for (var step = 0; step < AGENT_MAX_STEPS; step++) {
      statusEl.textContent = step === 0 ? 'Thinking ...' : 'Working with the data (step ' + (step + 1) + ') ...';
      var data = await agentApi(apiKey, model, agentHistory);
      var content = data.content || [];
      agentHistory.push({ role: 'assistant', content: content });

      var texts = content.filter(function (c) { return c.type === 'text' && c.text.trim(); });
      var calls = content.filter(function (c) { return c.type === 'tool_use'; });

      if (calls.length === 0 || data.stop_reason !== 'tool_use') {
        agentRender('assistant', texts.map(function (c) { return c.text; }).join('\n\n') || '(no answer)');
        statusEl.textContent = '';
        return;
      }
      var results = [];
      for (var ci = 0; ci < calls.length; ci++) {
        var call = calls[ci], out;
        agentRenderTool(call.name, call.input);
        try { out = await agentRunTool(call.name, call.input); }
        catch (err) { out = { error: String(err && err.message || err) }; }
        results.push({ type: 'tool_result', tool_use_id: call.id, content: JSON.stringify(out) });
      }
      agentHistory.push({ role: 'user', content: results });
    }
    agentRender('assistant', '(Stopped after ' + AGENT_MAX_STEPS + ' steps without a final answer.)');
    statusEl.textContent = '';
  } catch (err) {
    console.error(err);
    agentHistory = agentHistory.slice(0, turnStart);
    statusEl.textContent = 'Error: ' + err.message;
  } finally {
    btn.disabled = false;
  }
}

function clearAgent() {
  agentHistory = [];
  document.getElementById('agentLog').innerHTML = '';
  document.getElementById('agentStatus').textContent = '';
}

/* -------------------------------- render --------------------------------- */

function agentEscape(s) {
  return String(s).replace(/[&<>"']/g, function (ch) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch];
  });
}

function agentInline(s) {
  return s.replace(/`([^`]+)`/g, '<code>$1</code>')
    .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>');
}

// Small markdown subset: headings, lists, tables, code blocks, bold, inline code.
function agentMarkdown(md) {
  var lines = agentEscape(md).split('\n');
  var html = '', inCode = false, list = null, table = [];
  var flushList = function () { if (list) { html += '</' + list + '>'; list = null; } };
  var flushTable = function () {
    if (!table.length) return;
    var rows = table.filter(function (r) { return !/^\|?\s*:?-{2,}/.test(r); });
    html += '<div style="overflow-x:auto"><table class="table table-sm table-bordered small my-2">';
    rows.forEach(function (r, n) {
      var cells = r.replace(/^\|/, '').replace(/\|$/, '').split('|');
      var tag = n === 0 ? 'th' : 'td';
      html += '<tr>' + cells.map(function (c) { return '<' + tag + '>' + agentInline(c.trim()) + '</' + tag + '>'; }).join('') + '</tr>';
    });
    html += '</table></div>';
    table = [];
  };
  lines.forEach(function (line) {
    if (/^```/.test(line)) { flushList(); flushTable(); html += inCode ? '</pre>' : '<pre class="p-2 bg-light border rounded small" style="white-space:pre-wrap;word-break:break-all;">'; inCode = !inCode; return; }
    if (inCode) { html += line + '\n'; return; }
    if (/^\s*\|.*\|\s*$/.test(line)) { flushList(); table.push(line.trim()); return; }
    flushTable();
    var h = line.match(/^(#{1,4})\s+(.*)$/);
    if (h) { flushList(); html += '<div class="fw-semibold mt-2">' + agentInline(h[2]) + '</div>'; return; }
    var ul = line.match(/^\s*[-*]\s+(.*)$/);
    var ol = line.match(/^\s*\d+[.)]\s+(.*)$/);
    if (ul || ol) {
      var want = ul ? 'ul' : 'ol';
      if (list !== want) { flushList(); html += '<' + want + ' class="mb-1">'; list = want; }
      html += '<li>' + agentInline((ul || ol)[1]) + '</li>';
      return;
    }
    flushList();
    if (line.trim() === '') { html += '<div style="height:.4rem"></div>'; return; }
    html += '<div>' + agentInline(line) + '</div>';
  });
  flushList(); flushTable();
  if (inCode) html += '</pre>';
  return html;
}

function agentRender(role, text) {
  var log = document.getElementById('agentLog');
  var div = document.createElement('div');
  div.className = 'p-2 rounded mb-2 ' + (role === 'user' ? 'bg-light border' : 'border');
  div.style.wordBreak = 'break-word';
  div.innerHTML = '<div class="small text-secondary mb-1">' + (role === 'user' ? 'You' : 'Claude') + '</div>' +
    (role === 'user' ? '<div>' + agentEscape(text).replace(/\n/g, '<br>') + '</div>' : agentMarkdown(text));
  div.setAttribute('data-role', role);
  log.appendChild(div);
  log.scrollTop = log.scrollHeight;
}

function agentRenderTool(name, input) {
  var log = document.getElementById('agentLog');
  var div = document.createElement('div');
  div.className = 'small text-secondary mb-1';
  div.style.wordBreak = 'break-all';
  div.setAttribute('data-role', 'tool');
  var arg = JSON.stringify(input || {});
  if (arg.length > 120) arg = arg.slice(0, 120) + '...';
  div.innerHTML = '<i class="bi bi-tools me-1"></i>' + agentEscape(name) + ' ' + agentEscape(arg === '{}' ? '' : arg);
  log.appendChild(div);
}
