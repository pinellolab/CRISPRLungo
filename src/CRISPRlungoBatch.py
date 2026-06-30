#!/usr/bin/env python
"""
CRISPRlungoBatch
================
Run CRISPRlungo on many samples described in a single batch file (like
CRISPRessoBatch), then build a combined summary table and an interactive
summary HTML in which each row links to that sample's full report.

Batch file
----------
A tab- (or comma-) separated file with a header row. One row per sample.
Column names map to CRISPRlungo arguments. Recognised columns:

  name                 (required)  sample name -> output sub-folder
  ref                  (required)  reference FASTA (path)
  treated              (required)  treated FASTQ(.gz) (path)
  target               (required)  target / protospacer sequence
  control                          control FASTQ(.gz) for background filtering
  additional_target                second target sequence
  induced_sequence_path            desired/induced sequence FASTA
  induced_mutation_patterns        desired mutation pattern file
  integration_file                 FASTA of possible integration sequences
  umi                              TRUE/FALSE -> use --umi
  cleavage_pos, window, clust_cutoff, min_read_cnt, min_read_freq,
  largeins_cutlen, largedel_cutlen, allele_plot_window, allele_plot_lines,
  minimap2_opt                     per-sample overrides (optional)

Any of the above (except name) can instead be supplied once on the command
line with --<col> VALUE to act as a default for every sample. A value in the
batch file always overrides the command-line default.

Usage
-----
  CRISPRlungoBatch -b batch.txt -o batch_output [-t 8] [--ref ref.fasta] ...

Outputs (in batch_output/)
  <name>/                              full CRISPRlungo run per sample
  CRISPRlungoBatch_summary.txt         tab-separated combined summary
  CRISPRlungoBatch_summary.html        interactive summary (click a sample)
"""

import argparse
import os
import sys
import subprocess
import time
import html as _html

# Editing categories as written in results/mutation_summary_count.txt
CATEGORIES = ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv',
              'Complex', 'ComplexWithLargeMut']
PRECISE_COLS = ['Precise', 'Precise_with_mutations', 'Partial',
                'Partial_with_mutations']

# Reserved columns handled specially (sample name + positional CRISPRlungo args)
RESERVED = {'name', 'ref', 'treated', 'target'}
POSITIONAL = ['ref', 'treated', 'target']
# CRISPRlungo boolean flags (store_true): a column value of TRUE/1/YES emits --<name>
FLAG_OPTS = {'umi', 'no_multipass_align', 'whole_window_between_targets',
             'using_all_mutations', 'only_consider_desired_mutation', 'mix_tag',
             'just_visualization', 'show_all_between_allele'}
# Option columns whose value is a file path (resolved relative to the batch file)
PATH_OPTS = {'control', 'induced_sequence_path', 'induced_mutation_patterns',
             'integration_file', 'large_induced_insertion'}
# ANY other column is passed straight through as --<column> <value> to CRISPRlungo,
# so new CRISPRlungo options can be used as batch columns with no code change.


def parse_batch_file(path):
    """Read a batch file -> list of dict rows.

    Columns may be separated by tabs, commas, or spaces (or a mix). Because
    sample names, file paths and sequences do not contain spaces, any run of
    whitespace is treated as a single delimiter (CSV is used only when commas
    are present and there are no tabs).
    """
    import re
    with open(path) as f:
        raw = [ln.rstrip('\n') for ln in f]
    lines = [ln for ln in raw if ln.strip() != '' and not ln.lstrip().startswith('#')]
    if not lines:
        sys.exit('ERROR: batch file is empty: ' + path)

    def split(ln):
        if ',' in ln and '\t' not in ln:
            return [c.strip() for c in ln.split(',')]
        return re.split(r'\s+', ln.strip())

    header = split(lines[0])
    rows = []
    for ln in lines[1:]:
        vals = split(ln)
        row = {header[i]: (vals[i].strip() if i < len(vals) else '')
               for i in range(len(header))}
        rows.append(row)
    return rows


def resolve_path(p, base):
    """Resolve a possibly relative path against the batch file's directory."""
    if not p:
        return p
    if os.path.isabs(p) or os.path.exists(p):
        return p
    cand = os.path.join(base, p)
    return cand if os.path.exists(cand) else p


def build_command(crispr_cmd, row, defaults, sample_out, threads, batch_dir):
    """Build the CRISPRlungo command list for one sample."""
    def get(key):
        v = row.get(key)
        if v is None or v == '':
            v = defaults.get(key)
        return v

    ref = resolve_path(get('ref'), batch_dir)
    treated = resolve_path(get('treated'), batch_dir)
    target = get('target')
    missing = [k for k, v in (('ref', ref), ('treated', treated), ('target', target))
               if not v]
    if missing:
        raise ValueError('missing required field(s): ' + ', '.join(missing))

    cmd = list(crispr_cmd) + [ref, treated, sample_out, target]

    # Every column (and CLI default) that is not a reserved/positional field is
    # passed straight through to CRISPRlungo as --<name> (flag) or --<name> VALUE.
    opt_keys = (set(row) | set(defaults)) - RESERVED
    opt_keys.discard('')
    for opt in sorted(opt_keys):
        v = get(opt)
        if v is None or str(v).strip() == '':
            continue
        if opt in FLAG_OPTS:
            if str(v).strip().upper() in ('TRUE', '1', 'YES', 'Y'):
                cmd.append('--' + opt)
        else:
            val = resolve_path(str(v), batch_dir) if opt in PATH_OPTS else str(v)
            cmd += ['--' + opt, val]
    cmd += ['--threads', str(threads)]
    return cmd


def read_summary_counts(sample_out):
    """Parse results/mutation_summary_count.txt -> {col: int}."""
    path = os.path.join(sample_out, 'results', 'mutation_summary_count.txt')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip() != '']
    if len(lines) < 2:
        return None
    header = lines[0].split('\t')
    values = lines[1].split('\t')
    d = {}
    for i, h in enumerate(header):
        h = h.strip()
        if h == '':
            continue
        v = values[i].strip() if i < len(values) else ''
        try:
            d[h] = int(v)
        except (ValueError, TypeError):
            d[h] = 0
    return d


def read_input_summary(sample_out):
    """Parse results/input_summary.txt -> {key: value}."""
    path = os.path.join(sample_out, 'results', 'input_summary.txt')
    info = {}
    if not os.path.exists(path):
        return info
    for ln in open(path):
        if ':' in ln:
            k, v = ln.split(':', 1)
            info[k.strip()] = v.strip()
    return info


def pct(n, d):
    return (100.0 * n / d) if d else 0.0


def write_summary_txt(out_dir, results):
    """results: list of dicts with name, status, counts(dict), info(dict)."""
    path = os.path.join(out_dir, 'CRISPRlungoBatch_summary.txt')
    cols = (['Name', 'Status', 'Target_1', 'Target_2', 'Total_reads',
             'Reads_used', 'Modified_pct'] +
            CATEGORIES + [c + '_pct' for c in CATEGORIES] + ['Precise', 'Partial'])
    with open(path, 'w') as fw:
        fw.write('\t'.join(cols) + '\n')
        for r in results:
            c = r.get('counts') or {}
            info = r.get('info') or {}
            used = c.get('used', 0)
            wt = c.get('WT', 0)
            modified = used - wt
            line = [r['name'], r['status'],
                    info.get('Target_1', ''), info.get('Target_2', ''),
                    str(c.get('all_reads', '')), str(used),
                    '%.2f' % pct(modified, used)]
            line += [str(c.get(cat, 0)) for cat in CATEGORIES]
            line += ['%.2f' % pct(c.get(cat, 0), used) for cat in CATEGORIES]
            line += [str(c.get('Precise', 0)), str(c.get('Partial', 0))]
            fw.write('\t'.join(line) + '\n')
    return path


def write_summary_html(out_dir, results, batch_name='CRISPRlungo Batch'):
    path = os.path.join(out_dir, 'CRISPRlungoBatch_summary.html')
    show_precise = any((r.get('counts') or {}).get('Precise', 0) or
                       (r.get('counts') or {}).get('Partial', 0) for r in results)

    head_cols = ['Sample', 'Total reads', 'Reads used', 'Modified %'] + \
                [c + ' %' for c in CATEGORIES]
    if show_precise:
        head_cols += ['Precise %', 'Partial %']
    head_cols += ['Report']

    def th(cols):
        return ''.join('<th>%s</th>' % _html.escape(c) for c in cols)

    body_rows = []
    for r in results:
        name = r['name']
        c = r.get('counts') or {}
        if 'used' not in c:
            used = c.get('Treated_used', 0)
        else:
            used = c.get('used', 0)
        ok = r['status'] == 'OK' and c
        if not ok:
            body_rows.append(
                '<tr class="failed"><td>%s</td><td colspan="%d">%s</td>'
                '<td>-</td></tr>' % (
                    _html.escape(name), len(head_cols) - 2,
                    _html.escape('FAILED: ' + r.get('error', r['status']))))
            continue
        wt = c.get('WT', 0)
        modified = used - wt
        cells = ['<td class="name">%s</td>' % _html.escape(name),
                 '<td>%s</td>' % format(c.get('all_reads', 0), ','),
                 '<td>%s</td>' % format(used, ','),
                 '<td class="mod">%.1f</td>' % pct(modified, used)]
        for cat in CATEGORIES:
            cells.append('<td>%.1f</td>' % pct(c.get(cat, 0), used))
        if show_precise:
            cells.append('<td>%.1f</td>' % pct(c.get('Precise', 0), used))
            cells.append('<td>%.1f</td>' % pct(c.get('Partial', 0), used))
        rel = '%s/combined_graphs.html' % name
        cells.append('<td><a class="btn" href="%s" target="_blank">View &#8599;</a></td>'
                     % _html.escape(rel))
        body_rows.append('<tr onclick="window.open(\'%s\',\'_blank\')">%s</tr>'
                         % (_html.escape(rel), ''.join(cells)))

    html_doc = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%(title)s</title>
<style>
:root{--bd:#e3e6ea;--hd:#1f2d3d;--accent:#2563eb;--mod:#b91c1c;}
*{box-sizing:border-box}
body{font-family:'Inter',-apple-system,Segoe UI,Roboto,sans-serif;margin:0;
 background:#f6f8fa;color:#1f2d3d;padding:28px}
h1{font-size:20px;margin:0 0 4px}
.sub{color:#6b7280;font-size:13px;margin-bottom:18px}
.wrap{overflow-x:auto;background:#fff;border:1px solid var(--bd);border-radius:10px}
table{border-collapse:collapse;width:100%%;font-size:13px;white-space:nowrap}
thead th{position:sticky;top:0;background:var(--hd);color:#fff;font-weight:600;
 padding:10px 12px;text-align:right;font-size:12px}
thead th:first-child{text-align:left}
tbody td{padding:9px 12px;text-align:right;border-top:1px solid var(--bd)}
tbody td.name{text-align:left;font-weight:600}
tbody td.mod{color:var(--mod);font-weight:600}
tbody tr{cursor:pointer}
tbody tr:hover{background:#eef4ff}
tr.failed{background:#fff5f5;color:#b91c1c;cursor:default}
a.btn{display:inline-block;padding:4px 10px;background:var(--accent);color:#fff;
 border-radius:6px;text-decoration:none;font-size:12px}
a.btn:hover{background:#1d4fd8}
.legend{margin-top:12px;color:#6b7280;font-size:12px}
</style></head>
<body>
<h1>%(title)s</h1>
<div class="sub">%(n)d sample(s) &middot; generated %(date)s &middot; click a row to open its full report</div>
<div class="wrap"><table>
<thead><tr>%(head)s</tr></thead>
<tbody>
%(rows)s
</tbody></table></div>
<div class="legend">Percentages are of reads used in analysis. Modified %% = (used &minus; WT) / used.</div>
</body></html>
""" % {
        'title': _html.escape(batch_name),
        'n': len(results),
        'date': time.strftime('%Y-%m-%d %H:%M'),
        'head': th(head_cols),
        'rows': '\n'.join(body_rows),
    }
    with open(path, 'w') as fw:
        fw.write(html_doc)
    return path


def main():
    parser = argparse.ArgumentParser(
        description='Run CRISPRlungo on a batch of samples and build a combined summary.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--batch_file', required=True,
                        help='Tab/comma separated batch file (one row per sample)')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Batch output directory')
    parser.add_argument('-t', '--threads', type=int, default=8,
                        help='Threads per CRISPRlungo run')
    parser.add_argument('-n', '--batch_name', default='CRISPRlungo Batch',
                        help='Title shown in the summary HTML')
    parser.add_argument('--crisprlungo_cmd', default=None,
                        help="How to invoke CRISPRlungo. Default: 'python <this_dir>/CRISPRlungo.py'")
    parser.add_argument('--skip_failed', action='store_true',
                        help='Continue the batch even if a sample fails (default: also continue, kept for clarity)')
    parser.add_argument('--rerun', action='store_true',
                        help='Re-run samples even if an output already exists')
    # Any other --<crisprlungo_option> [value] on the command line becomes a
    # default applied to every sample (a batch-file column of the same name overrides it).
    args, extra = parser.parse_known_args()

    # assemble CRISPRlungo invocation
    if args.crisprlungo_cmd:
        crispr_cmd = args.crisprlungo_cmd.split()
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        crispr_cmd = [sys.executable, os.path.join(here, 'CRISPRlungo.py')]

    # Parse leftover CLI tokens: '--opt value' or '--flag' -> global defaults
    defaults = {}
    i = 0
    while i < len(extra):
        tok = extra[i]
        if tok.startswith('--'):
            key = tok[2:]
            if key in FLAG_OPTS:
                defaults[key] = 'TRUE'; i += 1
            elif i + 1 < len(extra) and not extra[i + 1].startswith('--'):
                defaults[key] = extra[i + 1]; i += 2
            else:
                defaults[key] = 'TRUE'; i += 1
        else:
            print('  (ignoring unrecognised argument: %s)' % tok); i += 1

    batch_dir = os.path.dirname(os.path.abspath(args.batch_file))
    rows = parse_batch_file(args.batch_file)
    if not rows:
        sys.exit('ERROR: no samples found in batch file')
    if any('name' not in r or not r['name'] for r in rows):
        sys.exit("ERROR: batch file must have a non-empty 'name' column for every row")

    os.makedirs(args.output_dir, exist_ok=True)

    results = []
    t0 = time.time()
    for i, row in enumerate(rows, 1):
        name = row['name']
        sample_out = os.path.join(args.output_dir, name)
        print('\n[%d/%d] %s' % (i, len(rows), name))
        rec = {'name': name, 'status': 'OK'}

        done_marker = os.path.join(sample_out, 'combined_graphs.html')
        if os.path.exists(done_marker) and not args.rerun:
            print('  already done (use --rerun to force) -> reusing results')
        else:
            try:
                cmd = build_command(crispr_cmd, row, defaults, sample_out,
                                    args.threads, batch_dir)
            except ValueError as e:
                print('  SKIP: ' + str(e))
                rec.update(status='ERROR', error=str(e))
                results.append(rec)
                continue
            print('  $ ' + ' '.join(cmd))
            st = time.time()
            proc = subprocess.run(cmd)
            if proc.returncode != 0:
                print('  FAILED (exit %d)' % proc.returncode)
                rec.update(status='ERROR', error='exit code %d' % proc.returncode)
                results.append(rec)
                if not args.skip_failed:
                    # keep going through the batch regardless, but record it
                    pass
                continue
            print('  done in %.1f s' % (time.time() - st))

        counts = read_summary_counts(sample_out)
        if counts is None:
            rec.update(status='ERROR', error='no summary output produced')
        else:
            rec['counts'] = counts
            rec['info'] = read_input_summary(sample_out)
        results.append(rec)

    txt = write_summary_txt(args.output_dir, results)
    html = write_summary_html(args.output_dir, results, args.batch_name)

    ok = sum(1 for r in results if r['status'] == 'OK' and r.get('counts'))
    print('\n==================================================')
    print('Batch finished: %d/%d OK in %.1f s' % (ok, len(results), time.time() - t0))
    print('Summary table : ' + txt)
    print('Summary HTML  : ' + html)


if __name__ == '__main__':
    main()
