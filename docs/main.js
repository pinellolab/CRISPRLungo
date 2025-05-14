
const form = document.getElementById("main-form");
const statusDiv = document.getElementById("status");
const outputDiv = document.getElementById("output");

var target, target2, reference, desiredSeq, treatedFile, controlFile, cleavagePos, cleavagePos2, strand, strand2, longjoinBandWidth, chainingBandWidth;
var statusElem, treatedAlignSam, controlAlignSam, mainResult, filteredSigMut, filteredTreatedReads;
var desiredAlignSam = "";
var mutTypeCnt = {"WT":0, "Ins": 0, "Del": 0, "Sub": 0, "LargeIns": 0, "LargeDel": 0, "Inv": 0, "Complex": 0, "Precise": 0, "Partial": 0};

var subPlotList = [];
var insPlotList = [[],[],[]]; //Pos, Len, Freq.
var delPlotList = [[],[],[]];

var subPosList = [];
var insPosList = []; 
var delPosList = [];

var subProportionList = [[],[],[],[]]



var windowRange = 5;
var cleavagePosTarget = 17;

function reverseComplement(seq) {
    const complement = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'};
    let revComp = "";
    for (let i = seq.length - 1; i >= 0; i--) {
        const base = seq[i];
        revComp += complement[base] || base; 
    }
    return revComp;
}

function main() {

    target = document.getElementById("targetSequence").value.trim().toUpperCase();
    target2 = document.getElementById("targetSequence2").value.trim().toUpperCase();
    reference = document.getElementById("reference").value.trim().toUpperCase();
    desiredSeq = document.getElementById("desiredSequence").value.trim().toUpperCase();
    treatedFile =   document.getElementById("treatedFile").files[0];
    controlFile =   document.getElementById("controlFile").files[0];

    statusElem = document.getElementById("status");
    
    if (target === "" || reference === "") {
        statusElem.innerHTML = "Please enter both Target and Reference sequences";
        statusElem.style.color = "red";
        return;
    }

    if  (!treatedFile || !controlFile) {
        statusElem.innerHTML = "No fastq file";
        statusElem.style.color = "red";
        return;
    }

    if (reference.includes(target)) {
        cleavagePos = reference.indexOf(target) + 16;
        strand = 1;
    } else if (reference.includes(reverseComplement(target))) {
        cleavagePos = reference.indexOf(reverseComplement(target)) + 3;
        strand = -1;
    } else {
        statusElem.innerHTML = "Target sequence is not found in the Reference";
        statusElem.style.color = "red";
        return ;
    }
    
    
    if (target2 !== '') {
        if (reference.includes(target2)) {
            cleavagePos2 = reference.indexOf(target2) + 16;
            strand2 = 1;
        } else if (reference.includes(reverseComplement(target2))) {
            cleavagePos2 = reference.indexOf(reverseComplement(target2)) + 3;
            strand2 = -1;
        } else {
            statusElem.innerHTML = "Second target sequence is not found in the Reference";
            statusElem.style.color = "red";
            return ;
        }
    }

    let status_str = 'Cleavage site is ' + cleavagePos;
    if (target2 !== "") {
        status_str += ' and second cleavage site is ' + cleavagePos2;
    }
    statusElem.innerHTML = status_str;
    statusElem.style.color = "green";

    longjoinBandWidth = parseInt(reference.length * 0.3);
    chainingBandWidth = 500;
    if (longjoinBandWidth < 500) {
        chainingBandWidth = longjoinBandWidth;
    }

    for (i = 0; i < reference.length; i++) {
        subPosList.push(0);
        insPosList.push(0);
        delPosList.push(0);
        subProportionList.push([0,0,0,0,0,0,0]); //ATGCN, Del, Ins
    }

    
    runAlignDesired();

}