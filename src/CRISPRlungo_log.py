#!/usr/bin/env python
"""Console output helpers for CRISPRlungo.

Plain text only, no colours, no boxes.  Three things are provided:

  * ``banner()`` / ``finish()``   - run header (version + command) and footer
  * ``plan()`` / ``step()``       - numbered section headers, one per process
  * ``task()``                    - a one-line progress entry that updates in
                                    place while it runs and is finalised with
                                    its elapsed time when it is done

When stdout is not a terminal (redirected to a log file, piped, run inside a
batch) the in-place updates are suppressed and only the final line of each
task is written, so log files stay readable.

Nothing in this module touches analysis behaviour.
"""

import os
import shlex
import sys
import time
from datetime import datetime

VERSION = '0.2'

_RULE_WIDTH = 62
_INDENT = '  '

_state = {
	'steps': [],
	'current': 0,
	'start': None,
	'start_text': None,
	'open_len': 0,
}


# ---------------------------------------------------------------- internals

def _is_tty():
	try:
		return sys.stdout.isatty()
	except Exception:
		return False


def _out(text):
	sys.stdout.write(text)
	sys.stdout.flush()


def _clear_open_line():
	"""Erase the transient line currently on screen, if there is one."""
	if _state['open_len']:
		_out('\r' + ' ' * _state['open_len'] + '\r')
		_state['open_len'] = 0


def _draw_transient(text):
	"""Draw text on the current line; it will be overwritten later."""
	if not _is_tty():
		return
	pad = ' ' * max(0, _state['open_len'] - len(text))
	_out('\r' + text + pad)
	_state['open_len'] = len(text)


def _draw_final(text):
	"""Close off the current line for good."""
	if _is_tty() and _state['open_len']:
		pad = ' ' * max(0, _state['open_len'] - len(text))
		_out('\r' + text + pad + '\n')
	else:
		_out(text + '\n')
	_state['open_len'] = 0


def _fmt_elapsed(seconds):
	if seconds < 60:
		return f'{seconds:.2f} s'
	if seconds < 3600:
		return f'{int(seconds // 60)} m {int(seconds % 60):02d} s'
	return f'{int(seconds // 3600)} h {int((seconds % 3600) // 60):02d} m'


def started_at():
	"""Formatted start time recorded by banner(), for the run summary files."""
	return _state['start_text'] or datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def command_line():
	"""The command the user typed, with the script path shortened."""
	argv = list(sys.argv)
	if argv:
		argv[0] = os.path.basename(argv[0])
	try:
		return shlex.join(argv)
	except AttributeError:				# python < 3.8
		return ' '.join(shlex.quote(a) for a in argv)


# ------------------------------------------------------------------- header

def banner(tool='CRISPRlungo', **fields):
	"""Run header: tool version, the command line, start time, extra fields."""
	_state['start'] = time.time()
	_state['start_text'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	rows = [('Command', command_line()),
			('Started', _state['start_text'])]
	rows += [(k.replace('_', ' '), v) for k, v in fields.items() if v not in (None, False, '')]
	width = max(len(k) for k, _ in rows)

	_out('=' * _RULE_WIDTH + '\n')
	_out(f'{tool} v{VERSION}\n')
	_out('=' * _RULE_WIDTH + '\n')
	for key, value in rows:
		_out(f'{key.ljust(width)} : {value}\n')
	_out('=' * _RULE_WIDTH + '\n')


def finish(**fields):
	"""Run footer: total elapsed time plus any closing fields."""
	_clear_open_line()
	total = time.time() - _state['start'] if _state['start'] else 0
	_out('\n' + '=' * _RULE_WIDTH + '\n')
	_out(f'Finished in {_fmt_elapsed(total)}\n')
	for key, value in fields.items():
		if value in (None, False, ''):
			continue
		_out(f'{key.replace("_", " ")} : {value}\n')
	_out('=' * _RULE_WIDTH + '\n')


# -------------------------------------------------------------------- steps

def plan(steps):
	"""Register the section titles for this run so steps can be numbered."""
	_state['steps'] = list(steps)
	_state['current'] = 0


def step(title=None):
	"""Start a new section: blank line, '[n/total] Title', underline."""
	_clear_open_line()
	_state['current'] += 1
	n = _state['current']
	total = len(_state['steps'])
	if title is None:
		title = _state['steps'][n - 1] if n <= total else ''
	head = f'[{n}/{total}] {title}' if total else f'[{n}] {title}'
	_out('\n' + head + '\n')
	_out('-' * max(_RULE_WIDTH, len(head)) + '\n')


# ----------------------------------------------------------------- messages

def info(message, indent=1):
	_clear_open_line()
	_out(_INDENT * indent + str(message) + '\n')


def blank():
	_clear_open_line()
	_out('\n')


def warn(message):
	_clear_open_line()
	_out(_INDENT + 'WARNING: ' + str(message) + '\n')


def error(message, detail=None, exit_code=1):
	"""Report a fatal problem and (by default) stop the run."""
	_clear_open_line()
	_out('\nERROR: ' + str(message) + '\n')
	if detail:
		for line in str(detail).rstrip().splitlines():
			_out(_INDENT + line + '\n')
	if exit_code is not None:
		sys.exit(exit_code)


# -------------------------------------------------------------------- tasks

class Task:
	"""One line of output that updates in place, then finalises with a time.

	        my_task = log.task('Aligning treated.fastq')
	        my_task.update('pass 2')
	        my_task.done('3 passes')

	prints, while running:   ``  Aligning treated.fastq ... pass 2``
	and when finished:       ``  Aligning treated.fastq ... done 2.31 s  (3 passes)``
	"""

	def __init__(self, label, indent=1):
		self.label = str(label)
		self.prefix = _INDENT * indent
		self.start = time.time()
		self.closed = False
		_clear_open_line()
		_draw_transient(f'{self.prefix}{self.label} ...')

	@property
	def elapsed(self):
		return time.time() - self.start

	def update(self, detail):
		"""Redraw the line with a new progress detail (terminal only)."""
		if not self.closed:
			_draw_transient(f'{self.prefix}{self.label} ... {detail}')

	def progress(self, current, total, noun=''):
		suffix = f' {noun}' if noun else ''
		self.update(f'{current}/{total}{suffix}')

	def done(self, detail=None):
		"""Finalise the line as successful."""
		if self.closed:
			return
		line = f'{self.prefix}{self.label} ... done {_fmt_elapsed(self.elapsed)}'
		if detail:
			line += f'  ({detail})'
		_draw_final(line)
		self.closed = True

	def skip(self, reason='skipped'):
		if self.closed:
			return
		_draw_final(f'{self.prefix}{self.label} ... {reason}')
		self.closed = True

	def fail(self, reason='failed'):
		if self.closed:
			return
		_draw_final(f'{self.prefix}{self.label} ... {reason}')
		self.closed = True


def task(label, indent=1):
	return Task(label, indent=indent)


# -------------------------------------------------------------------- misc

def count(value):
	"""Thousands separator for read counts etc."""
	try:
		return f'{int(value):,}'
	except (TypeError, ValueError):
		return str(value)


def summary_table(title, counts, total=None, indent=1):
	"""A small aligned count/percentage block, e.g. the read classification."""
	rows = [(k, int(v)) for k, v in counts.items() if int(v)]
	if not rows:
		return
	if total is None:
		total = sum(v for _, v in rows)
	_clear_open_line()
	pad = _INDENT * indent
	name_w = max(len(k) for k, _ in rows)
	value_w = max(len(count(v)) for _, v in rows)
	_out('\n' + pad + title + '\n')
	for key, value in rows:
		share = f'{value * 100 / total:5.1f} %' if total else ''
		_out(f'{pad}  {key.ljust(name_w)}  {count(value).rjust(value_w)}  {share}\n')


def stream_stderr(process, label, indent=1):
	"""Show a subprocess' stderr as a single updating line (badread, etc)."""
	handle = task(label, indent=indent)
	last = ''
	for line in process.stderr:
		line = line.strip()
		if line:
			last = line
			handle.update(line[:70])
	process.wait()
	return handle, last
