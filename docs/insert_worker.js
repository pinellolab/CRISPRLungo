importScripts('https://biowasm.com/cdn/v3/aioli.js');

self.onmessage = async function (event) {

    var { treated_read_to_mut, reference } = event.data;

    var fileLines = [];
    for (i of Object.entries(treated_read_to_mut)) {
        var readName = '>' + i[0];
        var m_n = 0;
        for (m of i[1][0]) {
            x = m.split(':')[1].split('_');
            if (x[0] == 'Ins' && x[1].length > 20) {
                fileLines.push(readName + '_' + m_n + '\n' + x[1] + '\n');
                m_n += 1
            }
            
        }
        treated_read_to_mut[i[0]].push([false]);
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

        await CLI.mount([{
        name: "reference.fa",
        data: reference
        }]);

        await CLI.mount([{
            name: "input.fa",
            data: partFileLines
        }]);
        partFileLines = null;

        const command = `minimap2 -ax map-ont reference.fa input.fa `;
        var output = await CLI.exec(command);

        for (var line of output.stdout.split('\n').filter(item => item !== '')) {
            if (line[0] == '@') {
                continue;
            }
            var line_sp = line.split('\t');
            if (line_sp[5] == '*' || parseInt(line_sp[4],10) < 40) {
                continue;
            }
            if (line_sp[1] !== '0' && line_sp[1] !== '16' && line_sp[1] !== '2048' && line_sp[1] !== '2064') {
                continue;
            }
            var read_id = [];
            for (i of line_sp[0].split('_')) {
                read_id.push(i);
            }
            var reads = [line_sp];
            var invCheck = false;
            var querySeq = line_sp[9];
            var query_seq_len = querySeq.length;

            var flag = parseInt(line_sp[1], 10);
            var align_len = 0;
            var clip_len = [0,0];
            var front_clip = true;
            var s = '';

            for (i of line_sp[5]) {
                if (!isNaN(parseInt(i, 10))) {
                    s += i;
                } else {
                    s = parseInt(s, 10);
                    if (i == "M" || i == "D") {
                        align_len += s;
                    } else if (i == 'S' || i == 'H') {
                        if (front_clip) {
                            clip_len[0] = s;
                        } else {
                            clip_len[1] = s;
                        }        
                    }
                    front_clip = false;
                    s = '';
                }
            }

            if (align_len < 20) {
                continue
            }

            ref_name = line_sp[2];
            strand = 1;
            if (flag & 0x10) {
                clip_len = clip_len.slice().reverse();
                strand = -1;
            }
            read_id = read_id.slice(0, read_id.length-1).join('_');
            if (ref_name == 'ref' && strand == -1) {
                treated_read_to_mut[read_id][2][0] = true;
                ref_name = 'inversion';
            }
            var ref_start = parseInt(line_sp[3], 10);
            var ref_end = ref_start + align_len;
            treated_read_to_mut[read_id][2].push([ref_name, ref_start, ref_end, strand].join('_'));
        }
        self.postMessage({ type: 1, progress: partIndex*100/totalIndex, alignedCount: alignResDict.length, n: n });
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

  self.postMessage({ type: 2, result: treated_read_to_mut });
  //} catch (error) {
  //  self.postMessage({ error: error.message });
  //}

};