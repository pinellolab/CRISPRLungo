/* ============================================================================
 * custom_classify.js
 * ----------------------------------------------------------------------------
 * AI-assisted custom re-classification for the CRISPRLungo web tool.
 *
 * Flow:
 *   1. On the results page the user types a classification request + their
 *      Anthropic API key and picks a model.
 *   2. We send the request, a description of the allele data model, and a
 *      sample of the current alleles to the Claude Messages API.
 *   3. Claude returns a pure JS function `classifyAllele(allele, parsed, ctx)`.
 *   4. We auto-run that function over `filteredTreatedReads` to bucket every
 *      allele into a user-defined category, then render a pie chart, a count
 *      table, and a per-category allele plot (reusing drawAllelePlot).
 *
 * Depends on globals defined by main.js / CRISPRlungo_ai.js:
 *   filteredTreatedReads, mainResult, reference, cleavagePos, cleavagePos2,
 *   strand, strand2, largeInsLen, largeDelLen, windowRange, drawAllelePlot().
 * ==========================================================================*/

var customAlleleByCat = {};   // { category: { alleleString: count } }
var customCategoryCounts = {}; // { category: totalReadCount }
var customClassifierCode = null; // the JS source of the last-used classifier
var customClassifierMeta = {};   // { request, model, date } for the last generated classifier

/* ----------------------------- allele parsing ---------------------------- */

// Parse a single mutation token, e.g. "312_312:Sub_A>G", "310_315:Del_5",
// "311_311:Ins_ATG", with an optional trailing "_inv".
function parseMutationToken(tok) {
  if (!tok) return null;
  var inversion = false;
  if (tok.slice(-4) === '_inv') { inversion = true; tok = tok.slice(0, -4); }

  var colon = tok.indexOf(':');
  if (colon === -1) return null;

  var posParts = tok.slice(0, colon).split('_');
  var start = Number(posParts[0]);
  var end = Number(posParts[1]);

  var rest = tok.slice(colon + 1).split('_'); // [TYPE, DETAIL...]
  var type = rest[0];
  var detail = rest.slice(1).join('_');

  var m = {
    raw: tok, start: start, end: end, type: type,
    len: 0, seq: '', from: '', to: '', inversion: inversion
  };

  if (type === 'Sub') {
    var fr = detail.split('>');
    m.from = fr[0] || '';
    m.to = fr[1] || '';
    m.len = m.from.length || 1;
  } else if (type === 'Del' || type === 'LargeDel') {
    m.len = Number(detail) || (end - start) || 0;
  } else if (type === 'Ins' || type === 'LargeIns') {
    m.seq = detail;
    m.len = detail.length;
  }
  return m;
}

// Parse a full allele string into { allele, isWT, editStatus, muts:[...] }.
function parseAllele(allele) {
  var editStatus = 'none';
  var s = allele;
  if (s.slice(-'_PreciseEdit'.length) === '_PreciseEdit') {
    editStatus = 'precise';
    s = s.slice(0, -'_PreciseEdit'.length);
  } else if (s.slice(-'_PartialEdit'.length) === '_PartialEdit') {
    editStatus = 'partial';
    s = s.slice(0, -'_PartialEdit'.length);
  }

  var isWT = (s === 'WT' || s === '');
  var muts = [];
  if (!isWT) {
    var toks = s.split(',');
    for (var i = 0; i < toks.length; i++) {
      var m = parseMutationToken(toks[i]);
      if (m) muts.push(m);
    }
  }
  return { allele: allele, isWT: isWT, editStatus: editStatus, muts: muts };
}

// Analysis context handed to the generated classifier.
function buildClassifierContext() {
  return {
    cleavagePos: cleavagePos,
    cleavagePos2: (cleavagePos2 === '' || cleavagePos2 == null) ? null : cleavagePos2,
    reference: reference,
    referenceLength: reference.length,
    inducedMutations: (mainResult && mainResult.induced_mutations) ? mainResult.induced_mutations : [],
    largeInsLen: largeInsLen,
    largeDelLen: largeDelLen,
    window: windowRange,
    strand: strand,
    strand2: strand2
  };
}

/* ------------------------------ prompt build ----------------------------- */

function buildSystemPrompt() {
  return [
    "You are a bioinformatics assistant that writes JavaScript to classify CRISPR editing",
    "alleles detected by CRISPRLungo from long-read sequencing data.",
    "",
    "Respond with EXACTLY ONE JavaScript code block and NOTHING else (no prose before or after).",
    "The code must define a single pure function with this EXACT signature:",
    "",
    "    function classifyAllele(allele, parsed, ctx) { /* ... */ return 'CategoryName'; }",
    "",
    "Rules for the function:",
    "- Pure: no network, no DOM, no async, no external libraries, no console output.",
    "- Returns a short human-readable category name (string) for the given allele.",
    "- Any allele that does not match a requested category must return 'Other'.",
    "- Make categories mutually exclusive; keep names short.",
    "- The function will be saved and re-applied to OTHER experiments (different reference, cut site, edit).", // 20261110_ai_classify
    "  Never hard-code absolute coordinates, sequences or allele strings seen in the examples;", // 20261110_ai_classify
    "  express positions relative to ctx.cleavagePos / ctx.cleavagePos2 and sizes via ctx values.", // 20261110_ai_classify
    "  If the request itself names fixed positions, define them as named constants at the top of the function.", // 20261110_ai_classify
    "",
    "Inputs:",
    "- allele: the raw allele string (format below).",
    "- parsed: { allele, isWT (bool), editStatus ('precise'|'partial'|'none'), muts: [mutation,...] }",
    "  mutation = { raw, start (int, 0-based ref position), end (int), type",
    "    ('Sub'|'Del'|'LargeDel'|'Ins'|'LargeIns'), len (int), seq (inserted bases for Ins),",
    "    from/to (reference>alternate bases for Sub), inversion (bool) }",
    "- ctx: { cleavagePos (int), cleavagePos2 (int or null for single-cut), reference (string),",
    "    referenceLength (int), inducedMutations (array of raw mutation strings),",
    "    largeInsLen (int), largeDelLen (int), window (int), strand, strand2 }",
    "",
    "Allele string format: comma-separated mutation tokens, or 'WT' for wild-type.",
    "Each token is START_END:TYPE_DETAIL. Examples:",
    "  '312_312:Sub_A>G'  = substitution A to G at position 312",
    "  '310_315:Del_5'    = 5 bp deletion starting at 310",
    "  '311_311:Ins_ATG'  = insertion of ATG after position 311",
    "  'LargeDel'/'LargeIns' tokens mark large events; a trailing '_inv' marks an inversion.",
    "",
    "Helpful facts:",
    "- Deletion length and insertion length are in mutation.len (insertion len == seq.length).",
    "- 'in-frame' means (len % 3 === 0); 'frameshift' means (len % 3 !== 0).",
    "- Distance to the cut site: Math.abs(mutation.start - ctx.cleavagePos).",
    "- parsed.isWT is true for unedited reads; use parsed.editStatus for intended-edit status.",
    "",
    "Return only the code block."
  ].join("\n");
}

function buildUserPrompt(request, ctx) {
  // Provide a sample of the most frequent alleles so the model can match logic.
  var entries = Object.entries(filteredTreatedReads); // already sorted desc
  var sample = entries.slice(0, 40)
    .map(function (e) { return "  " + e[0] + "   (" + e[1] + " reads)"; })
    .join("\n");

  return [
    "Classification request:",
    request,
    "",
    "Experiment context:",
    "- cleavagePos: " + ctx.cleavagePos,
    "- cleavagePos2: " + (ctx.cleavagePos2 === null ? "none (single cut)" : ctx.cleavagePos2),
    "- reference length: " + ctx.referenceLength,
    "- analysis window: " + ctx.window,
    "- large insertion cutoff: " + ctx.largeInsLen + " bp, large deletion cutoff: " + ctx.largeDelLen + " bp",
    "- intended/induced mutations: " + (ctx.inducedMutations.length ? ctx.inducedMutations.join(", ") : "none"),
    "- number of distinct alleles: " + entries.length,
    "",
    "Example allele strings from the current results (most frequent first):",
    sample || "  (none)",
    "",
    "Write classifyAllele(allele, parsed, ctx) to satisfy the request."
  ].join("\n");
}

/* ------------------------------- API call -------------------------------- */

async function callAnthropic(apiKey, model, systemPrompt, userPrompt) {
  const response = await fetch("https://api.anthropic.com/v1/messages", {
    method: "POST",
    headers: {
      "content-type": "application/json",
      "x-api-key": apiKey,
      "anthropic-version": "2023-06-01",
      "anthropic-dangerous-direct-browser-access": "true"
    },
    body: JSON.stringify({
      model: model,
      max_tokens: 4096,
      system: systemPrompt,
      messages: [{ role: "user", content: userPrompt }]
    })
  });

  if (!response.ok) {
    let detail = "";
    try { detail = (await response.json()).error.message; } catch (e) { detail = await response.text(); }
    throw new Error("API error " + response.status + ": " + detail);
  }

  const data = await response.json();
  const parts = (data.content || []).filter(function (c) { return c.type === "text"; });
  return parts.map(function (c) { return c.text; }).join("\n");
}

// Pull the JS out of a fenced code block; fall back to the whole text.
function extractJsCode(text) {
  var fence = text.match(/```(?:javascript|js)?\s*([\s\S]*?)```/i);
  return (fence ? fence[1] : text).trim();
}

/* --------------------------- run the classifier -------------------------- */

// Build the classifier function from code and bucket every allele.
function applyClassifier(code, ctx) {
  // Isolate the generated code and hand back the function reference.
  var factory = new Function(code + "\n; if (typeof classifyAllele !== 'function') { throw new Error('classifyAllele was not defined'); } return classifyAllele;");
  var classifyAllele = factory();

  var counts = {};
  var byCat = {};

  for (var allele in filteredTreatedReads) {
    var cnt = filteredTreatedReads[allele];
    var parsed = parseAllele(allele);
    var cat;
    try {
      cat = classifyAllele(allele, parsed, ctx);
    } catch (e) {
      cat = "Error";
    }
    if (cat == null || cat === "") cat = "Unclassified";
    cat = String(cat);

    counts[cat] = (counts[cat] || 0) + cnt;
    if (!byCat[cat]) byCat[cat] = {};
    byCat[cat][allele] = cnt;
  }
  return { counts: counts, byCat: byCat };
}

/* ------------------------------ rendering -------------------------------- */

function renderCustomResults() {
  var used = (mainResult && mainResult.treated_align_cnt) ? mainResult.treated_align_cnt.Used : 0;

  // Sort categories by count desc.
  var sortedCats = Object.entries(customCategoryCounts).sort(function (a, b) { return b[1] - a[1]; });

  // ---- Pie chart ----
  var pieColors = ['#4357A7', '#ED2024', '#6BBD45', '#FABF37', '#98D9ED', '#BF71AE',
    '#CECFD1', '#000000', '#F39C12', '#16A085', '#8E44AD', '#2C3E50'];
  var trace = {
    labels: sortedCats.map(function (c) { return c[0]; }),
    values: sortedCats.map(function (c) { return c[1]; }),
    type: 'pie',
    marker: { colors: pieColors }
  };
  Plotly.newPlot('customPiePlot', [trace], { title: 'Custom classification' });

  // ---- Count table (horizontal: categories laid out across columns) ----
  var container = document.getElementById('customTableContainer');
  container.innerHTML = "";
  var table = document.createElement("table");
  table.className = "table table-bordered table-striped result-table w-100";

  // Header row: leading label cell + one column per category.
  var thead = document.createElement("thead");
  thead.className = "thead-light";
  var headHtml = "<tr><th>Category</th>";
  sortedCats.forEach(function (c) { headHtml += "<th>" + escapeHtml(c[0]) + "</th>"; });
  headHtml += "</tr>";
  thead.innerHTML = headHtml;
  table.appendChild(thead);

  // Two data rows: read counts and percentage of used reads.
  var tbody = document.createElement("tbody");
  var readsRow = "<tr><th>Reads</th>";
  var pctRow = "<tr><th>% of used</th>";
  sortedCats.forEach(function (c) {
    var pct = used ? (c[1] * 100 / used).toFixed(2) + " %" : "-";
    readsRow += "<td>" + c[1] + "</td>";
    pctRow += "<td>" + pct + "</td>";
  });
  readsRow += "</tr>";
  pctRow += "</tr>";
  tbody.innerHTML = readsRow + pctRow;
  table.appendChild(tbody);
  container.appendChild(table);

  // ---- Category dropdown for allele plot ----
  var select = document.getElementById('customCategorySelect');
  select.innerHTML = "";
  sortedCats.forEach(function (c) {
    var opt = document.createElement("option");
    opt.value = c[0];
    opt.textContent = c[0] + " (" + c[1] + " reads)";
    select.appendChild(opt);
  });

  document.getElementById('customResults').style.display = "block";
  drawCustomAllelePlot();
}

// Redraw the allele plot for the currently-selected category.
function drawCustomAllelePlot() {
  var select = document.getElementById('customCategorySelect');
  if (!select || !select.value) return;
  var subset = customAlleleByCat[select.value] || {};

  // Sort the subset by read count (drawAllelePlot shows the top rows).
  var sorted = Object.fromEntries(
    Object.entries(subset).sort(function (a, b) { return b[1] - a[1]; })
  );

  d3.select('#customAllelePlot').selectAll('*').remove();
  drawAllelePlot(sorted, '#customAllelePlot');
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, function (ch) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch];
  });
}

/* -------------------------------- driver --------------------------------- */

async function runCustomClassification() {
  var statusEl = document.getElementById('classifyStatus');
  var btn = document.getElementById('classifyRunBtn');

  var apiKey = document.getElementById('anthropicApiKey').value.trim();
  var model = document.getElementById('claudeModel').value;
  var request = document.getElementById('classifyPrompt').value.trim();

  if (!filteredTreatedReads || Object.keys(filteredTreatedReads).length === 0) {
    statusEl.textContent = "No results to classify yet — run an analysis first.";
    return;
  }
  if (!apiKey) { statusEl.textContent = "Please enter your Anthropic API key."; return; }
  if (!request) { statusEl.textContent = "Please describe how you want the reads classified."; return; }

  btn.disabled = true;
  statusEl.textContent = "Asking Claude to write a classifier ...";

  try {
    var ctx = buildClassifierContext();
    var raw = await callAnthropic(apiKey, model, buildSystemPrompt(), buildUserPrompt(request, ctx));
    var code = extractJsCode(raw);

    customClassifierMeta = { request: request, model: model, date: new Date().toISOString() };

    statusEl.textContent = "Applying the classifier to " + Object.keys(filteredTreatedReads).length + " alleles ...";
    runClassifierCode(code, ctx);

    // Enable download of the freshly generated algorithm.
    document.getElementById('classifyDownloadBtn').disabled = false;
    statusEl.textContent = "Done. " + Object.keys(customCategoryCounts).length + " categories.";
  } catch (err) {
    console.error(err);
    statusEl.textContent = "Error: " + err.message;
  } finally {
    btn.disabled = false;
  }
}

// Shared: run a classifier code string over the current results and render.
function runClassifierCode(code, ctx) {
  ctx = ctx || buildClassifierContext();

  // Show the algorithm for transparency.
  document.getElementById('classifyCode').textContent = code;
  document.getElementById('classifyCodeWrap').style.display = "block";

  var result = applyClassifier(code, ctx);
  customClassifierCode = code;
  customCategoryCounts = result.counts;
  customAlleleByCat = result.byCat;

  renderCustomResults();
}

// Save the current classifier to a .js file (with a metadata header).
function downloadAlgorithm() {
  if (!customClassifierCode) return;
  var m = customClassifierMeta || {};
  var header = [
    "/* CRISPRLungo custom classification algorithm",
    " * Generated: " + (m.date || new Date().toISOString()),
    " * Model: " + (m.model || "unknown"),
    " * Built on: reference length " + (typeof reference === "string" ? reference.length : "?") + ", cleavagePos " + cleavagePos + (cleavagePos2 === "" || cleavagePos2 == null ? "" : ", cleavagePos2 " + cleavagePos2), // 20261110_ai_classify
    " * Request: " + String(m.request || "").replace(/\*\//g, "* /").replace(/\n/g, "\n *          "),
    " *",
    " * Defines classifyAllele(allele, parsed, ctx). Re-load it in the",
    " * CRISPRLungo web tool via 'Apply a saved algorithm file'.",
    " */",
    ""
  ].join("\n");

  var blob = new Blob([header + customClassifierCode + "\n"], { type: "text/javascript" });
  var stamp = (m.date || new Date().toISOString()).replace(/[:.]/g, "-").slice(0, 19);
  var link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "crisprlungo_classifier_" + stamp + ".js";
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(link.href);
}

// Load a saved algorithm file and apply it without calling the API.
async function applyUploadedAlgorithm() {
  var statusEl = document.getElementById('applyAlgoStatus');
  var fileInput = document.getElementById('algoFile');
  var btn = document.getElementById('applyAlgoBtn');

  if (!filteredTreatedReads || Object.keys(filteredTreatedReads).length === 0) {
    statusEl.textContent = "No results to classify yet — run an analysis first.";
    return;
  }
  if (!fileInput.files || fileInput.files.length === 0) {
    statusEl.textContent = "Please choose an algorithm .js file.";
    return;
  }

  btn.disabled = true;
  statusEl.textContent = "Reading algorithm file ...";
  try {
    var code = await fileInput.files[0].text();
    customClassifierMeta = { request: "(loaded from file: " + fileInput.files[0].name + ")", model: "-", date: new Date().toISOString() };

    statusEl.textContent = "Applying to " + Object.keys(filteredTreatedReads).length + " alleles ...";
    runClassifierCode(code);

    // Allow re-downloading the loaded algorithm too.
    document.getElementById('classifyDownloadBtn').disabled = false;
    statusEl.textContent = "Done. " + Object.keys(customCategoryCounts).length + " categories.";
  } catch (err) {
    console.error(err);
    statusEl.textContent = "Error: " + err.message;
  } finally {
    btn.disabled = false;
  }
}
