// 20261128_paired_end
// Merges paired-end reads (Illumina R1/R2) into one read per pair, so the rest
// of the analysis sees what it gets from a long-read run: one read over the
// whole amplicon.
//
// This is CRISPResso2's merge: CRISPResso2 (v2.3+) runs
//   fastp -i R1 -I R2 --merge --overlap_len_require 10 --disable_adapter_trimming
//         --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering
// and analyses only the merged reads. With adapter trimming, quality, length and
// poly-G filtering all off, fastp's merge comes down to OverlapAnalysis::analyze
// and ::merge (src/overlapanalysis.cpp, fastp 1.4.0), ported here line for line:
//   - R2 is reverse-complemented (anything but ACGT becomes N);
//   - offsets are tried from 0 upward while the overlap is longer than
//     overlap_len_require, then from 0 downward (the insert is shorter than the
//     reads and adapter was read); the FIRST offset accepted wins;
//   - an offset is accepted when the mismatches in the first 50 bases of the
//     overlap are at most min(overlap_diff_limit, int(overlap * percent / 100));
//     mismatches past base 50 are counted but do not reject it;
//   - the merged read is R1 up to the end of the overlap, then the rest of the
//     reverse-complemented R2; R1's bases and qualities are kept in the overlap;
//   - pairs that do not overlap are left out.

// The gzip helpers are a copy of the ones in minimap2_worker_prog.js, which
// cannot be imported here without also loading its aligner.

function findGzipMemberStarts(bytes) {
  const starts = [];
  for (var i = 0; i + 2 < bytes.length; i++) {
    if (bytes[i] === 0x1f && bytes[i + 1] === 0x8b && bytes[i + 2] === 0x08) {
      starts.push(i);
    }
  }
  return starts;
}

// Decompress bytes[start, end) as exactly one complete gzip member, or return
// null. This doubles as the boundary test: a gzip member ends with a CRC32 and
// the uncompressed length, both verified by the decompressor, so a wrong end
// offset fails here rather than silently producing wrong data.
async function inflateGzipRange(bytes, start, end) {
  const ds = new DecompressionStream('gzip');
  const writer = ds.writable.getWriter();
  // Not awaited: these settle only once the reader below drains the stream,
  // and they reject if the range turns out not to be a whole member.
  writer.write(bytes.subarray(start, end)).catch(() => { });
  writer.close().catch(() => { });

  const reader = ds.readable.getReader();
  const parts = [];
  try {
    for (; ;) {
      const chunk = await reader.read();
      if (chunk.done) { break; }
      parts.push(chunk.value);
    }
  } catch (e) {
    return null;
  }
  return parts;
}

// Decode every member of a (possibly concatenated) gzip into one string.
async function gunzipToText(bytes) {
  const starts = findGzipMemberStarts(bytes);
  if (starts.length === 0 || starts[0] !== 0) {
    throw new Error('This file does not start with a gzip header, so it cannot be read as .gz. '
      + 'Upload an uncompressed .fastq instead.');
  }

  // Candidate end offsets, in order, plus the end of the file for the last member.
  const bounds = starts.slice(1);
  bounds.push(bytes.length);

  const decoder = new TextDecoder();
  var text = '';
  var cursor = 0;
  var members = 0;

  while (cursor < bytes.length) {
    var decoded = null;
    var memberEnd = -1;
    for (var k = 0; k < bounds.length; k++) {
      if (bounds[k] <= cursor) { continue; }
      decoded = await inflateGzipRange(bytes, cursor, bounds[k]);
      if (decoded) { memberEnd = bounds[k]; break; }
    }
    if (!decoded) {
      throw new Error('This .gz file could not be fully decompressed after ' + members + ' gzip member(s). '
        + 'Decompress it first (gunzip -c in.fastq.gz > in.fastq) and upload the FASTQ.');
    }
    for (var part of decoded) {
      // Streaming decode so a multi-byte character split across chunks -- or
      // across members -- still decodes correctly.
      text += decoder.decode(part, { stream: true });
    }
    members++;
    cursor = memberEnd;
  }
  text += decoder.decode();
  return text;
}

async function readText(file) {
  let isGz = !!(file.name && /\.gz$/i.test(file.name));
  if (!isGz) {
    try {
      const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
      isGz = head[0] === 0x1f && head[1] === 0x8b;
    } catch (e) { /* plain text */ }
  }
  if (!isGz) return await file.text();
  if (typeof DecompressionStream === 'undefined') {
    throw new Error('This browser cannot decompress .gz files. Use a recent Chrome/Edge/Firefox, or upload an unzipped .fastq.');
  }
  try {
    const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
    return await new Response(stream).text();
  } catch (e) {
    return await gunzipToText(new Uint8Array(await file.arrayBuffer()));
  }
}

// fastp's reverseComplement: lower case is complemented to upper case and
// anything that is not A/C/G/T becomes N.
var RC = { A: 'T', a: 'T', T: 'A', t: 'A', C: 'G', c: 'G', G: 'C', g: 'C' };

function revComp(s) {
  var out = new Array(s.length);
  for (var i = 0, j = s.length - 1; j >= 0; i++, j--) out[i] = RC[s[j]] || 'N';
  return out.join('');
}

function reverse(s) {
  return s.split('').reverse().join('');
}

// Read name without the leading '@', anything after the first blank, and a
// trailing /1 or /2 - only used to check that R1 and R2 are in the same order.
function pairName(header) {
  return header.slice(1).split(/\s/)[0].replace(/\/[12]$/, '');
}

var COMPLETE_COMPARE_REQUIRE = 50;

function countMismatches(a, ai, b, bi, len) {
  var d = 0;
  for (var k = 0; k < len; k++) if (a.charCodeAt(ai + k) !== b.charCodeAt(bi + k)) d++;
  return d;
}

// acceptNoGapOverlap() in fastp
function accept(a, ai, b, bi, len, limit) {
  var prefix = Math.min(len, COMPLETE_COMPARE_REQUIRE);
  var d = countMismatches(a, ai, b, bi, prefix);
  if (d > limit) return -1;
  if (len > COMPLETE_COMPARE_REQUIRE) d = countMismatches(a, ai, b, bi, len);
  return d;
}

// OverlapAnalysis::analyze(r1, r2, diffLimit, overlapRequire, diffPercentLimit), no gap
function analyze(str1, str2, diffLimit, overlapRequire, diffPercentLimit) {
  var len1 = str1.length, len2 = str2.length;
  var offset = 0, ol, lim, diff;
  while (offset < len1 - overlapRequire) {
    ol = Math.min(len1 - offset, len2);
    lim = Math.min(diffLimit, Math.trunc(ol * diffPercentLimit));
    diff = accept(str1, offset, str2, 0, ol, lim);
    if (diff >= 0) return { overlapped: true, offset: offset, overlapLen: ol, diff: diff };
    offset += 1;
  }
  offset = 0;
  while (offset > -(len2 - overlapRequire)) {
    ol = Math.min(len1, len2 - Math.abs(offset));
    lim = Math.min(diffLimit, Math.trunc(ol * diffPercentLimit));
    diff = accept(str1, 0, str2, -offset, ol, lim);
    if (diff >= 0) return { overlapped: true, offset: offset, overlapLen: ol, diff: diff };
    offset -= 1;
  }
  return { overlapped: false, offset: 0, overlapLen: 0, diff: 0 };
}

// OverlapAnalysis::merge
function merge(s1, q1, rcSeq2, rQual2, ov) {
  var ol = ov.overlapLen;
  var len1 = ol + Math.max(0, ov.offset);
  var len2 = ov.offset > 0 ? rcSeq2.length - ol : 0;
  var seq = s1.substr(0, len1), qual = q1.substr(0, len1);
  if (ov.offset > 0) {
    seq += rcSeq2.substr(ol, len2);
    qual += rQual2.substr(ol, len2);
  }
  return { seq: seq, qual: qual, len1: len1, len2: len2 };
}

function mergePair(s1, q1, s2, q2, p) {
  var rc2 = revComp(s2);
  var ov = analyze(s1, rc2, p.diffLimit, p.overlapRequire, p.diffPercentLimit / 100);
  if (!ov.overlapped) return null;
  var m = merge(s1, q1, rc2, reverse(q2), ov);
  m.offset = ov.offset;
  m.overlap = ov.overlapLen;
  m.diff = ov.diff;
  return m;
}

self.onmessage = async function (event) {
  var d = event.data;
  if (d.test) { self.postMessage({ type: 'test', result: runSelfTest() }); return; }
  try {
    var t1 = await readText(d.r1);
    var t2 = await readText(d.r2);
    var l1 = t1.split(/\r?\n/), l2 = t2.split(/\r?\n/);
    t1 = t2 = null;
    while (l1.length && l1[l1.length - 1] === '') l1.pop();
    while (l2.length && l2[l2.length - 1] === '') l2.pop();
    if (l1.length % 4 !== 0 || l2.length % 4 !== 0) {
      throw new Error('R1 or R2 is not a well-formed FASTQ (4 lines per read).');
    }
    if (l1.length !== l2.length) {
      throw new Error('R1 has ' + (l1.length / 4) + ' reads and R2 has ' + (l2.length / 4)
        + '. The two files of a pair must hold the same reads in the same order.');
    }
    var pairs = l1.length / 4;
    var out = [];
    var merged = 0, readThrough = 0, lenSum = 0, ovSum = 0;
    var step = Math.max(1, Math.floor(pairs / 50));
    for (var p = 0; p < pairs; p++) {
      var i = p * 4;
      if (l1[i][0] !== '@' || l2[i][0] !== '@') {
        throw new Error('Read ' + (p + 1) + ' does not start with @ in ' + (l1[i][0] !== '@' ? 'R1' : 'R2') + '.');
      }
      var name = pairName(l1[i]);
      if (name !== pairName(l2[i])) {
        throw new Error('Read ' + (p + 1) + ' is ' + name + ' in R1 but ' + pairName(l2[i])
          + ' in R2. The two files must hold the same reads in the same order.');
      }
      var m = mergePair(l1[i + 1], l1[i + 3], l2[i + 1], l2[i + 3], d);
      if (m) {
        merged++;
        lenSum += m.seq.length;
        ovSum += m.overlap;
        if (m.offset < 0) readThrough++;
        // fastp names the merged read after R1, with merged_<len1>_<len2> appended
        out.push(l1[i] + ' merged_' + m.len1 + '_' + m.len2 + '\n' + m.seq + '\n+\n' + m.qual + '\n');
      }
      if (p % step === 0) self.postMessage({ type: 'progress', done: p, total: pairs });
    }
    self.postMessage({
      type: 'done',
      fastq: out.join(''),
      stats: {
        pairs: pairs, merged: merged, notMerged: pairs - merged, readThrough: readThrough,
        meanLength: merged ? lenSum / merged : 0, meanOverlap: merged ? ovSum / merged : 0
      }
    });
  } catch (e) {
    self.postMessage({ type: 'error', message: (e && e.message) || String(e) });
  }
};

// fastp's own OverlapAnalysis::test() cases, to check the port.
function runSelfTest() {
  var r1 = 'CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC';
  var r2 = 'ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC';
  var ov = analyze(r1, revComp(r2), 2, 30, 0.2);
  var late1 = new Array(51).join('A') + new Array(31).join('C');
  var lateRc2 = new Array(51).join('A') + new Array(31).join('G');
  var ov2 = analyze(late1, revComp(revComp(lateRc2)), 0, 30, 0.0);
  return {
    basic: ov, basicPass: ov.overlapped && ov.offset === 10 && ov.overlapLen === 79 && ov.diff === 1,
    late: ov2, latePass: ov2.overlapped && ov2.offset === 0 && ov2.overlapLen === 80 && ov2.diff === 30
  };
}
