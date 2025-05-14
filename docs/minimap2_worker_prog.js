
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

function readFastqtoFasta(file) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = function (e) {
      const text = e.target.result;
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

      resolve(fasta);
    };
    reader.onerror = function (err) {
      reject(err);
    };
    reader.readAsText(file);
  });
}

function softClipped(alignOutput) {

  checkSA = 0;
  var readLines = [];
  var writeLines = [];

  for (var line of alignOutput.stdout.split('\n').filter(item => item !== '')) {

    if (line[0] == '@') {
      readLines.push(line);
      continue;
    }

    line = line.split('\t');

    var flag = parseInt(line[1]);
    if ((flag & 0x4) !== 0 || (flag & 0x100) !== 0 || (flag & 0x800) !== 0) { continue; }

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
  return [checkSA, readLines, writeLines];
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
    var readOut, writeLines, oriSeq, oriStrand, partStart, partEnd;

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

      const command = `minimap2 -ax map-ont -p 0.5 reference.fa input.fa -r ` + chainingBandWidth + ',' + longjoinBandWidth;
      var output = await CLI.exec(command);

      [checkSA, readOut, writeLines] = softClipped(output);

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
        queryName = queryName.slice(0, queryName.indexOf('_AI_'));
        var suppleCheck = false;
        if (n != 1) {
          read[1] = (read[1] | 0x800).toString();
        }
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
      self.postMessage({ type: 1, progress: partIndex*100/totalIndex, alignedCount: alignResDict.length, n: n });
    }
  

    for (var [x, reads] of Object.entries(alignResDict)) {
      var read_n = 0;
      reads[0] = reads[0][0];
      if (reads.length == 1) {
        var tags = reads[0].slice(11,);
        try{
          var filteredTags = tags.filter(tags => tags.slice(0, 2) !== 'SA');
        } catch(error) {
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

  self.postMessage({ type: 2, result: alignResList });
  //} catch (error) {
  //  self.postMessage({ error: error.message });
  //}

};