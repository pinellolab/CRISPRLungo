
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



var windowRange = Number(document.getElementById("analysisWindow").value);
var cleavagePosTarget = Number(document.getElementById("cleavageSite").value);
const largeInsLen = Number(document.getElementById("largeInsertionLen").value);
const largeDelLen = Number(document.getElementById("largeDeletionLen").value);
const wholeWindow = Number(document.querySelector('input[name="wholeWindow"]:checked').value)


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
    if (document.getElementById("controlFile").files.length == 0) {
        controlFile = false;
    } else {
        controlFile =   document.getElementById("controlFile").files[0];
    }

    statusElem = document.getElementById("status");
    
    if (target === "" || reference === "") {
        alert("Please enter both Target and Reference sequences");
        return;
    }

    if  (!treatedFile) {
        alert("No Treate fastq file");
        return;
    }

    if (reference.includes(target)) {
        cleavagePos = reference.indexOf(target) + cleavagePosTarget;
        strand = 1;
    } else if (reference.includes(reverseComplement(target))) {
        cleavagePos = reference.indexOf(reverseComplement(target)) + target.length - cleavagePosTarget - 2;
        strand = -1;
    } else {
         alert("Target sequence is not found in the Reference");
        return ;
    }
    
    
    if (target2 !== '') {
        if (target2.length < 10) {
            alert('Additional target is too short!')
            return ;
        }
        if (reference.includes(target2)) {
            cleavagePos2 = reference.indexOf(target2) + cleavagePosTarget;
            strand2 = 1;
        } else if (reference.includes(reverseComplement(target2))) {
            cleavagePos2 = reference.indexOf(reverseComplement(target2)) + target2.length - cleavagePosTarget - 2;
            strand2 = -1;
        } else {
             alert("Second target sequence is not found in the Reference");
            return ;
        }
    } else {
        cleavagePos2 = '';
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