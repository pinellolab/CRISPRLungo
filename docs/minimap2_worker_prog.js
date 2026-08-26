
importScripts('https://biowasm.com/cdn/v3/aioli.js');

// Read fastqfile to fasta

function reverseComplement(seq) {
  const complement = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' };
  let revComp = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    const base = seq[i];
    revComp += complement[base] || base;
  }
  return revComp;
}

// A .gz produced by a basecaller/demultiplexer, or by `cat a.gz b.gz > c.gz`,
// is usually a *concatenated* gzip: several complete members in one file.
// DecompressionStream('gzip') decodes only the first member and then throws
// "Junk found after end of compressed data" -- and when the stream is consumed
// through a Response that surfaces as a bare "TypeError: Failed to fetch", so
// the worker died before its first postMessage and the UI sat at 0% forever.
//
// Note that we cannot simply decode the whole file and keep whatever arrived
// before the error: erroring the readable side *discards output still queued*,
// so a large member loses most of its data (measured in Chrome: a 860 KB
// member yielded only the first 131072 bytes). Each member therefore has to be
// handed to the decompressor on its own, with no trailing bytes.

// Every offset where a gzip member could begin. The 3-byte signature also
// occurs by chance inside compressed data, so these are candidates, not
// boundaries; inflateGzipRange below is what confirms one.
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

async function readFastqtoFasta(file) {
  // Transparently decompress gzip (.gz) input. Detect by extension OR by the
  // gzip magic bytes (0x1f 0x8b) so a renamed .gz is still handled.
  let isGz = !!(file.name && /\.gz$/i.test(file.name));
  if (!isGz) {
    try {
      const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
      isGz = head[0] === 0x1f && head[1] === 0x8b;
    } catch (e) { /* ignore, treat as plain text */ }
  }

  let text;
  if (isGz) {
    if (typeof DecompressionStream === "undefined") {
      throw new Error("This browser cannot decompress .gz files. Use a recent Chrome/Edge/Firefox, or upload an unzipped .fastq.");
    }
    try {
      // Fast path: a single-member .gz streams straight through without ever
      // holding the compressed bytes in memory.
      const stream = file.stream().pipeThrough(new DecompressionStream("gzip"));
      text = await new Response(stream).text();
    } catch (e) {
      // Concatenated members land here (as a generic "Failed to fetch"), and so
      // does a genuinely damaged file -- gunzipToText tells the two apart.
      text = await gunzipToText(new Uint8Array(await file.arrayBuffer()));
    }
  } else {
    text = await file.text();
  }

  const lines = text.split(/\r?\n/);
  let fasta = [];
  for (let i = 0; i < lines.length; i += 4) {
    if (!lines[i]) continue;
    let header = lines[i];
    if (header.startsWith('@')) {
      header = '>' + header.slice(1);
    } else {
      header = '>' + header;
    }
    let sequence = lines[i + 1] || "";
    fasta.push(header + "\n" + sequence + "\n");
  }
  return fasta;
}

function softClipped(alignOutput, n) {

  checkSA = 0;
  var readLines = [];
  var writeLines = [];
  var unmappedCnt = 0;

  for (var line of alignOutput.stdout.split('\n').filter(item => item !== '')) {

    if (line[0] == '@') {
      readLines.push(line);
      continue;
    }

    line = line.split('\t');

    var flag = parseInt(line[1]);
    if ((flag & 0x100) !== 0 || (flag & 0x800) !== 0) { continue; }
    // Unmapped reads are dropped here and never reach the analysis, which is
    // why "Total reads" used to be the aligned-read count and the unmapped
    // number was never available at all. Tally them before dropping them.
    if ((flag & 0x4) !== 0) { unmappedCnt++; continue; }
    var mapQ = parseInt(line[4]);
    if (mapQ < 30 && n != 1) { continue; }

    var clippedSeq = [[], []]
    var cigar = line[5].match(/(\d+)([MIDNSHP=X])/g);

    var queryName = line[0];

    if (!queryName.includes('_AI')) {
      queryName += '_AI';
    }
    var isReverse = (flag & 0x10) ? -1 : 1;
    queryName += '_' + isReverse;

    line[0] = queryName;

    var cigarFirstLen = parseInt(cigar[0].slice(0, -1));
    var cigarLastLen = parseInt(cigar[cigar.length - 1].slice(0, -1));

    var check = false;
    if (cigar[0].slice(-1) == 'S' && cigarFirstLen > 100) {
      writeLines.push('>' + queryName + '\n' + line[9].slice(0, cigarFirstLen) + '\n');
      check = true;
    }
    if (cigar[cigar.length - 1].slice(-1) == 'S' && cigarLastLen > 100) {
      writeLines.push('>' + queryName + '\n' + line[9].slice(-cigarLastLen,) + '\n');
      check = true;
    }

    readLines.push(line.join('\t'));

    if (check) { checkSA++; }

  }
  return [checkSA, readLines, writeLines, unmappedCnt];
}

function getQueryAlignmentSequence(querySequence, cigar) {

  var st = 0;
  var ed = querySequence.length;

  if (cigar[0].slice(-1) == 'S') {
    st = cigar[0].slice(0, -1);
  }
  if (cigar[cigar.length - 1].slice(-1) == 'S') {
    ed -= cigar[cigar.length - 1].slice(0, -1);
  }

  return querySequence.slice(st, ed);
}

self.onmessage = async function (event) {

  const { fastqFile, reference, longjoinBandWidth, chainingBandWidth, fileType } = event.data;


  if (fileType == "file") {
    var fileLines = await readFastqtoFasta(fastqFile);
  } else if (fileType == "String") {
    var fileLines = [fastqFile];
  }


  var readLineCnt = fileLines.length;
  self.postMessage({ type: 0, fileLen: readLineCnt });

  // One FASTA record per input read, so this is the input read count that the
  // result table reports as "Total reads".
  var inputReadCount = (fileType == "file")
    ? readLineCnt
    : (String(fastqFile).match(/^>/gm) || []).length;
  var unmappedTotal = 0;

  // Run Minimap2

  const CLI = await new Aioli(["minimap2/2.22"], {
    printInterleaved: false
  });

  var alignResList = [];
  var get_header = false;
  var totalIndex = Math.floor(readLineCnt / 500) + 1;

  for (var partIndex = 0; partIndex < totalIndex; partIndex++) {

    var partFileLines = fileLines.splice(0, 500).join('');
    var alignResDict = {};
    var n = 0;
    var checkSA = 1;
    var readOut, writeLines, oriSeq, oriStrand, partStart, partEnd, unmappedInPass;

    while (checkSA > 0) {

      n++;

      await CLI.mount([{
        name: "reference.fa",
        data: ">Ref\n" + reference
      }]);

      if (n == 1) {
        await CLI.mount([{
          name: "input.fa",
          data: partFileLines
        }]);
        partFileLines = null;
      } else {
        await CLI.mount([{
          name: "input.fa",
          data: writeLines.join('')
        }]);
      }

      const command = `minimap2 -ax map-ont -z 100 -p 0.5 reference.fa input.fa -r ` + chainingBandWidth + ',' + longjoinBandWidth;
      var output = await CLI.exec(command);

      [checkSA, readOut, writeLines, unmappedInPass] = softClipped(output, n);

      // Only the first pass sees whole reads; later passes realign
      // soft-clipped fragments of reads that were already counted.
      if (n == 1) { unmappedTotal += unmappedInPass; }

      if (!get_header) {
        get_header = true;
        for (var read of readOut) {
          if (read[0] != '@') {
            break;
          }
          alignResList.push(read.replace(/\n$/, ""))
        }
      }

      for (var read of readOut) {

        if (read[0] == '@') {
          continue;
        }

        read = read.split('\t');

        var queryName = read[0];
        var strandInfo = queryName.split('_AI_')[1].split('_');
        var strand = 1;
        for (var x of strandInfo) {
          strand *= x;
        }
        if (strand == -1) {
          read[1] = (read[1] | 0x10).toString();
        } else {
          read[1] = (read[1] & ~0x10).toString();
        }
        queryName = queryName.slice(0, queryName.indexOf('_AI_'));
        var suppleCheck = false;
        if (n != 1) {
          read[1] = (read[1] | 0x800).toString();
        }

        read_tags = [];
        for (i of read.slice(11,)) {
          var ii = i.slice(0, 2);
          if (ii == 'AS' || ii == 'NM') {
            read_tags.push(i);
          }
        }
        read = read.slice(0, 11).concat(read_tags);

        if (!alignResDict[queryName]) {
          alignResDict[queryName] = [[read, read[9], read[1]]];
        } else {
          oriSeq = alignResDict[queryName][0][1];
          oriStrand = (alignResDict[queryName][0][2] & 0x10) ? -1 : 1;
          cigar = read[5].match(/(\d+)([MIDNSHP=X])/g);
          alignedSeq = getQueryAlignmentSequence(read[9], cigar);
          if (oriStrand * strand == -1) {
            partStart = oriSeq.indexOf(reverseComplement(alignedSeq));
            partEnd = partStart + alignedSeq.length;
          } else {
            partStart = oriSeq.indexOf(alignedSeq);
            partEnd = partStart + alignedSeq.length;
          }
          if (cigar[0].slice(-1) == 'S') {
            read[9] = alignedSeq;
            cigar = cigar.slice(1,);
            cigar.unshift(partStart + 'H');
          } else {
            cigar.unshift(partStart + 'H');
          }

          if (cigar[cigar.length - 1].slice(-1) == 'S') {
            read[9] = alignedSeq;
            cigar = cigar.slice(0, -1);
            cigar.push((oriSeq.length - partEnd) + 'H');
          } else {
            cigar.push((oriSeq.length - partEnd) + 'H');
          }
          read[5] = cigar.join('');

          alignResDict[queryName].push(read);
        }
      }
      self.postMessage({ type: 1, progress: partIndex * 100 / totalIndex, alignedCount: alignResDict.length, n: n });
    }


    for (var [x, reads] of Object.entries(alignResDict)) {
      var read_n = 0;
      reads[0] = reads[0][0];
      if (reads.length == 1) {
        var tags = reads[0].slice(11,);
        try {
          var filteredTags = tags.filter(tags => tags.slice(0, 2) !== 'SA');
        } catch (error) {
          console.log(error);
        }
        reads[0] = reads[0].slice(0, 11).concat(filteredTags);
      } else {
        reads[0].push("SA:Z:" + 'N;'.repeat(reads.length).slice(0, -1));
      }
      for (read of reads) {
        read_n += 1;
        alignResList.push(read.join('\t').replace(/\n$/, ""));
      }
    }
  }

  // Make Bam format

  /*const CLI_samtools = await new Aioli(["samtools/1.17"], {
    printInterleaved: false
  });

  await CLI_samtools.mount([{
    name: "align.sam",
    data: alignResList.join('\n')
  }]);

  const command = "samtools view -bS align.sam -o align.bam";
  var output = await CLI_samtools.exec(command);
  var res = await CLI_samtools.FileReader('align.sam');
  var binaryOutput = await CLI_samtools.FileReader('align.bam', 'binary');*/

  self.postMessage({ type: 2, result: alignResList, inputReadCount: inputReadCount, unmappedCount: unmappedTotal });
  //} catch (error) {
  //  self.postMessage({ error: error.message });
  //}

};