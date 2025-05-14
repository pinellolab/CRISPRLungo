importScripts("https://biowasm.com/cdn/v3/aioli.js");

var mutationArray = {};
var mutationList = [];

function alignLength(cigar) {

    let refLength = 0;

    for (const match of cigar) {
        const [, length, type] = match.match(/(\d+)([MIDNSHP=X])/);
        const num = parseInt(length, 10);

        if (type === "M" || type === "D") {
            refLength += num;
        }
    }

    return refLength;
}

function getReverseComplement(dnaSequence) {

    const complementMap = {
        A: 'T',
        T: 'A',
        C: 'G',
        G: 'C',
        a: 't',
        t: 'a',
        c: 'g',
        g: 'c'
    };

    let complementSequence = '';
    for (let i = dnaSequence.length - 1; i >= 0; i--) {
        const nucleotide = dnaSequence[i];
        if (complementMap[nucleotide]) {
            complementSequence += complementMap[nucleotide];
        } else {
            complementSequence += nucleotide;
        }
    }

    return complementSequence;
}


function analyzeSAReads(partialReads, querySeq, refSeq, refLen) {

    var partialReadInfo = [];

    for (let i = 0; i < partialReads.length; i++) {

        let read = partialReads[i];
        
        let cigar = read[0];
        let flag = read[1] * 1;
        let alignLen = alignLength(cigar);

        let posRefst = read[2] * 1;
        let posRefed = read[3];

        let posInRef = false;
        let partialSeq = read[4];
        let alignStrandInRefseq = 2;
        let alignStrandToSeq = 2;

        let posInSeq, posInSeqEd;

        if (cigar[0].slice(-1,) == 'S') {
            posInSeq = cigar[0].slice(0,-1) * 1;
            posInSeqEd = posInSeq;
            for (const x of cigar) {
                if (x.slice(-1,) in ['M', 'I']) {
                    posInSeqEd += x.slice(0,-1);
                }
            }
            alignStrandtoSeq = 1;
        } else {
            posInSeq = querySeq.indexOf(partialSeq);
            alignStrandToSeq = 1;
            if (posInSeq == -1) {
                posInSeq = querySeq.indexOf(getReverseComplement(partialSeq));
                alignStrandToSeq = -1;
            }
            posInSeqEd = posInSeq + partialSeq.length;
        }
        
        if (posInSeq ==-1) continue;

        if (posRefst < 100 && posRefed > refLen - 100) {
            alignStrandInRefseq = 3;
        } else if ((flag & 16) == 0) {
            alignStrandInRefseq = 1;
        } else if ((flag & 16) == 1) {
            alignStrandInRefseq = -1;
        } else continue;

        if (posRefst < 100) {
            posInRef = 1;
        } else if (posRefed > refLen - 100) {
            posInRef = -1;
        } else continue;

        partialReadInfo.push((posRefst, posRefed, posInSeq, posInSeqEd, alignStrandInRefseq, alignStrandToSeq, cigar));

    }

    partialReadInfo = partialReadInfo.sort((a, b) => a[2] - b[2]);
    partialReadInfo.push((2, 2, 2, 2, 2, 2, ''));

    var mutationsInReads = [];

    for (let i = 0; i < partialReadInfo.length - 1; i++) {
        let mutationsInRead = [];
        let info = partialReadInfo.slice(i, i+2);

        let refPos, seqPos, usingQuerySeq, cigar, insEndPos, rcQuerySeq, querySeqLen;

        if (info[0][4] == 3) {
            refPos = info[0][0];
            seqPos = info[0][2];
            usingQuerySeq = querySeq;

            if (info[0][5] == -1) {
                usingQuerySeq = getReverseComplement(querySeq);
                seqPos = querySeq.length - info[0][3];
            } 

            cigar = info[0][6]

            let operation, length 

            for (const c of cigar) {
                operation = c.slice(-1,);
                length = c.slice(0,-1) * 1;
                if (operation == 'M') {
                    for (let x = 0; x < length; x++) {
                        if (usingQuerySeq[seqPos + x] != refSeq[refPos + x]) {
                            mutationsInRead.push('0_' + (refPos + x) + '_1_' + refSeq[refPos + x] + '_' + usingQuerySeq[seqPos + x]);
                        }
                    }
                } else if (operation == 'I') {
                    mutationsInRead.push('1_' + refPos +"_"+ length+"_"+ usingQuerySeq.slice(seqPos, seqPos + length));
                } else if (operation == 'D') {
                    mutationsInRead.push('2_' + refPos +"_"+ length);
                }
                if (operation == 'M' || operation == 'D') refPos += length;
                if (operation == 'M' || operation == 'I') seqPos += length;
            }

            mutationsInReads.push(mutationsInRead);
            continue
        }

        if (info[0][4] * info[1][4] != 1) continue;
        if (info[0][5] * info[1][5] != 1) continue;

        if (info[0][5] == 1) {
            if (info[0][1] < info[1][0] - 100) mutationsInRead.push('2_' + info[0][1] +'_'+ (info[1][0] - info[0][1]));
            if (Math.abs(info[0][3] - info[1][2]) > 20) {
                insEndPos = info[1][0];
                if (insEndPos <= info[0][1]) {
                    insEndPos = info[0][1] + 1;
                }
                mutationsInRead.push('1_' + info[0][1]+'_'+ Math.abs(info[0][3] - info[1][2])+'_'+ querySeq.slice(info[0][3], info[0][3] + Math.abs(info[0][3] - info[1][2]))+'_'+ insEndPos);
            }
            for (const subInfo of info) {
                refPos = subInfo[0];
                seqPos = subInfo[2];
                for (const c of subInfo[6]) {
                    operation = c.slice(-1,);
                    length = c.slice(0,-1) * 1;
                    if (operation == 'M') {
                        for (let x = 0; x < length; x++) {
                            if (querySeq[seqPos + x] != refSeq[refPos + x]) {
                                mutationsInRead.push(('0_' + refPos + x+'_'+ 1+'_'+ refSeq[refPos + x]+'_'+ querySeq[seqPos + x]));
                            }
                        } 
                    } else if (operation == 'I') {
                        mutationsInRead.push(('1_' + refPos+'_'+ length+'_'+ querySeq.slice(seqPos, seqPos + length)+'_'+ refPos + 1));
                    } else if (operation == 'D') {
                        mutationsInRead.push(('2_' + refPos+'_'+ length));
                    }
                    if (operation == 'M' || operation == 'D') refPos += length;
                    if (operation == 'M' || operation == 'I') seqPos += length;
                }
            }
        } else if (info[0][5] == -1) {
            if (info[1][1] < info[0][0] - 100) {
                mutationsInRead.push((2, info[1][1], info[0][0] - info[1][1] - 1));
            }
            if (info[1][2] - info[0][3] > 20) {
                insEndPos = info[0][0];
                if (insEndPos <= info[1][1]) {
                    insEndPos = info[1][1] + 1
                }
                mutationsInRead.push((1, info[1][1], info[1][2] - info[0][3] - 1, getReverseComplement(querySeq.slice(info[0][3], info[0][3] + (info[1][2] - info[0][3] - 1)), insEndPos)));
            }
            rcQuerySeq = getReverseComplement(querySeq);
            querySeqLen = querySeq.length;
            for (const subInfo of info) {
                refPos = subInfo[0];
                seqPos = subInfo[3];
                seqPos = querySeqLen - seqPos;
                for (const c of subInfo[6]) {
                    operation = c.slice(-1,);
                    length = c.slice(0,-1) * 1;
                    if (operation == 'M') {
                        for (let x = 0; x < length; x++) {
                            if (rcQuerySeq[seqPos + x] != refSeq[refPos + x]) {
                                mutationsInRead.push('0_' + (refPos + x)+'_'+ 1+'_'+ refSeq[refPos + x]+'_'+ rcQuerySeq[seqPos + x]);
                            }
                        }
                    } else if (operation == 'I') {
                        mutationsInRead.push('1_' + refPos+'_'+ length+'_'+ rcQuerySeq.slice(seqPos, seqPos + x)+'_'+ (refPos + 1));
                    } else if (operation == 'D') {
                        mutationsInRead.push('2_' + refPos+'_'+ length);
                    }
                    if (operation == 'M' || operation == 'D') refPos += length;
                    if (operation == 'M' || operation == 'I') seqPos += length;
                }
            }
        }
        mutationsInReads.push(mutationsInRead);
    }

    return mutationsInReads;
}


function quantUniqueIndels(alignReads, refSeq, refLen) {

    var mutations = {};
    var totalReads = 0;
    var counter = 0;


    for (let i = 0; i < alignReads.length; i++) {
        let read = alignReads[i];
        if (read[0] == '@' || read == "") continue;
        let readSplit = read.split('\t');
        let flag = readSplit[1] * 1;
        if ((flag & (0x4 | 0x100 | 0x800)) !== 0) continue;

        let mutationsInRead = [];
        let cigar = readSplit[5].match(/(\d+)([MIDNSHP=X])/g);
        let refSt = readSplit[3] * 1;
        let refEd = refSt + alignLength(cigar);
        let refPos = refSt - 1;
        let SAReads, SAn, oriQuerySeq, nextciagr, nextRead, x, querySeq, queryPos, operation, length, tags, key;

        if (refSt > 100 && refEd < refLen - 100) continue;

        tags = {};
        for (const tag of readSplit.slice(11,)) {
            tags[tag.slice(0,2)] = tag.slice(5,);
        }

        if ('SA' in tags) {
            SAReads = [[cigar, flag, refSt, refEd, readSplit[9]]];
            SAn = tags['SA'].split(';').length;
            oriQuerySeq = readSplit[9];
            x = 0;
            nextRead = '';
            nextcigar = '';
            while (SAReads.length < SAn) {
                x += 1;
                nextRead = alignReads[i+x].split('\t');
                if ((nextRead[1] & (0x800)) != 0) {
                    nextciagr = nextRead[5].match(/(\d+)([MIDNSHP=X])/g);
                    SAReads.push([nextciagr, nextRead[1], nextRead[3], nextRead[3] + alignLength(nextcigar), nextRead[9]]);
                }
            }
            if (SAReads.length != SAn) continue;
            for (const SAMutations of analyzeSAReads(SAReads, oriQuerySeq, refSeq, refLen)) {
                totalReads += 1;
                for (const mut of SAMutations) {
                    if (!mutationArray[mut]){
                        mutationArray[mut] = 1
                    } else {
                        mutationArray[mut] += 1
                    }
                }
            }
        } else {
            if (readSplit[4] <= 30) continue;
            totalReads += 1;
            querySeq = readSplit[9];
            queryPos = 0;
            for (const c of cigar) {
                operation = c.slice(-1,);
                length = c.slice(0,-1) * 1;
                if (operation == 'M') {
                    for (x = 0; x < length; x++) {
                        if (querySeq[queryPos + x] != refSeq[refPos + x]) {
                            key = '0_' + (refPos + x)+'_'+ 1+'_'+ refSeq[refPos + x]+'_'+ querySeq[queryPos + x];
                            if (!mutationArray[key]) mutationArray[key] = 1;
                            else mutationArray[key] += 1;
                        }
                    } 
                } else if (operation == 'I') {
                    key = '1_' + refPos+'_'+ length+'_'+ querySeq.slice(refPos, refPos + length);
                    if (!mutationArray[key]) mutationArray[key] = 1;
                    else mutationArray[key] += 1;
                } else if (operation == 'D') {
                    key = '2_' + refPos +'_'+ length;
                    if (!mutationArray[key]) mutationArray[key] = 1;
                    else mutationArray[key] += 1;
                } 
                if (operation == 'M' || operation == 'D') refPos += length;
                if (operation == 'M' || operation == 'I' || operation == 'S') queryPos += length;
            counter += 1
            }
        }
    }
    return totalReads;
} 

self.addEventListener("message", async (event) => {

    const { referenceFile, readsFile } = event.data;
    const maxReads = 1000;

    try {
        self.postMessage({ type: "status", message: "Loading Minimap2..." });
        const CLI = await new Aioli(["minimap2/2.22"], { printInterleaved: false });

        self.postMessage({ type: "status", message: "Reading reference file..." });
        const referenceData = await referenceFile.text();
        const refSeq = referenceData.split('\n').slice(1,).join('').replaceAll('\n','').toUpperCase();
        const refLen = refSeq.length;
        await CLI.fs.writeFile("/reference.fasta", referenceData)

        self.postMessage({ type: "status", message: "Aligning reads..." });
        const reader = readsFile.stream().getReader();
        const decoder = new TextDecoder();

        let partialLine = "";
        let reads = "";
        let readCount = 0;
        let analyzedCount = 0;
        let maxReads = 1000;

        var result = '';
        var start_time = new Date();
        var end_time = '';
        var result_input = '';

        while (true) {
            const { done, value } = await reader.read();
            if (done) {
                await CLI.fs.writeFile("/reads.fastq", reads);
                result = await CLI.exec("minimap2 -a /reference.fasta /reads.fastq");
                result_input = result.stdout.split('\n');
                analyzedCount += quantUniqueIndels(result_input, refSeq, refLen);
                end_time = new Date() - start_time;
                self.postMessage({ type: "status", message: readCount + " reads  " + end_time + 'ms ... Done!' });
                break
            };

            const chunk = decoder.decode(value, { stream: true });
            const lines = (partialLine + chunk).split("\n");
            partialLine = lines.pop(); 
            
            for (let i = 0; i < lines.length; i += 1) {
                if (readCount % maxReads == 0) {
                    await CLI.fs.writeFile("/reads.fastq", reads);
                    result = await CLI.exec("minimap2 -a /reference.fasta /reads.fastq");
                    result_input = result.stdout.split('\n');
                    analyzedCount += quantUniqueIndels(result_input, refSeq, refLen);
                    end_time = new Date() - start_time;
                    self.postMessage({ type: "status", message: readCount + " reads  " + end_time + 'ms' });
                    start_time = new Date();
                    reads = '';
                }
                reads += lines[i] + "\n";
                readCount++;
            }
        }

        self.postMessage({ type: "output", message: result.stdout || "No alignment results." });

    } catch (error) {
        self.postMessage({ type: "error", message: error.message });
    }
});


