// 20261111_progress
// The progress page used to be one bar, one status line and one caret line that
// said "FASTQ file alignment" whatever was actually happening. A long-read run
// takes minutes, so it is worth saying which stage it is in, which stages are
// still to come, and how long each took.
(function () {
  'use strict';

  var state = { steps: [], index: {}, started: 0, timer: null };

  function el(id) { return document.getElementById(id); }

  function fmt(ms) {
    if (ms < 1000) return '<1 s';
    var s = ms / 1000;
    if (s < 60) return s.toFixed(1) + ' s';
    var m = Math.floor(s / 60);
    return m + ' min ' + Math.round(s - m * 60) + ' s';
  }

  function icon(status) {
    if (status === 'done') return '<i class="bi bi-check-circle-fill text-success"></i>';
    if (status === 'running') return '<span class="spinner-border spinner-border-sm text-primary" role="status"></span>';
    if (status === 'failed') return '<i class="bi bi-x-circle-fill text-danger"></i>';
    if (status === 'skipped') return '<i class="bi bi-dash-circle text-secondary"></i>';
    return '<i class="bi bi-circle text-secondary opacity-50"></i>';
  }

  function paint() {
    var list = el('progressSteps');
    if (!list) return;
    list.innerHTML = state.steps.map(function (s) {
      var right = s.detail || '';
      if (s.status === 'done' && s.elapsed) right = right ? right + ' · ' + fmt(s.elapsed) : fmt(s.elapsed);
      else if (s.status === 'running' && s.startedAt) {
        var run = fmt(Date.now() - s.startedAt);
        right = right ? right + ' · ' + run : run;
      }
      var cls = s.status === 'pending' ? 'text-secondary' : '';
      return '<div class="d-flex align-items-center gap-2 py-1 ' + cls + '">'
        + '<span style="width:1.2rem;text-align:center;">' + icon(s.status) + '</span>'
        + '<span class="flex-grow-1' + (s.status === 'running' ? ' fw-semibold' : '') + '">'
        + s.label + '</span>'
        + '<span class="small text-secondary text-end" style="min-width:9rem;">' + right + '</span>'
        + '</div>';
    }).join('');

    var total = el('progressElapsed');
    if (total && state.started) total.textContent = 'elapsed ' + fmt(Date.now() - state.started);
  }

  function tick() {
    if (!state.timer) state.timer = setInterval(paint, 500);
  }

  var api = {
    // steps: [{id, label}] - the ones that will not run are passed as skipped
    init: function (steps) {
      state.steps = steps.map(function (s) {
        return { id: s.id, label: s.label, status: s.status || 'pending', detail: s.detail || '' };
      });
      state.index = {};
      state.steps.forEach(function (s, i) { state.index[s.id] = i; });
      state.started = Date.now();
      api.bar(0);
      paint();
      tick();
    },

    start: function (id, detail) {
      var s = state.steps[state.index[id]];
      if (!s) return;
      s.status = 'running';
      s.startedAt = Date.now();
      if (detail !== undefined) s.detail = detail;
      paint();
    },

    detail: function (id, text) {
      var s = state.steps[state.index[id]];
      if (!s) return;
      s.detail = text;
      paint();
    },

    done: function (id, detail) {
      var s = state.steps[state.index[id]];
      if (!s) return;
      s.elapsed = s.startedAt ? Date.now() - s.startedAt : 0;
      s.status = 'done';
      if (detail !== undefined) s.detail = detail;
      paint();
    },

    fail: function (id, message) {
      var s = state.steps[state.index[id]];
      if (s) { s.status = 'failed'; s.detail = message || 'failed'; }
      paint();
    },

    // The one line under the bar: what is happening right now.
    say: function (text) {
      var st = el('status');
      if (st) { st.innerHTML = text; st.style.color = ''; }
    },

    bar: function (percent) {
      var b = el('progress-bar');
      if (!b) return;
      var v = Math.max(0, Math.min(100, percent));
      b.style.width = v.toFixed(1) + '%';
      b.textContent = v.toFixed(0) + '%';
      var wrap = b.parentNode;
      if (wrap && wrap.setAttribute) wrap.setAttribute('aria-valuenow', v.toFixed(0));
    },

    finish: function () {
      if (state.timer) { clearInterval(state.timer); state.timer = null; }
      api.bar(100);
      paint();
    },

    totalElapsed: function () { return state.started ? Date.now() - state.started : 0; }
  };

  window.lungoProgress = api;
})();
