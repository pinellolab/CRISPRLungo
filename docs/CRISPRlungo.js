

var margin = { top: 40, right: 20, bottom: 20, left: 100 };
var cellSize = 16; // base 하나당 폭
var cellSizeHeight = 17;
var cellSizeWidth = 13;
var numRows = 10;
var plotWindow = 20;

var preciseEditThreshold = 1.0;
var partialEditThreshold = 0.8;

var showAllBetweenAllele = false;
var minaResult;


async function runAlignDesired() {

  document.getElementById("inputPage").classList.remove("active");
  document.getElementById("progressPage").classList.add("active");
  document.body.style.backgroundColor = "#fdfdfd";


    if (desiredSeq !== "") {

      const minimap2Worker_desired = new Worker('minimap2_worker_prog.js');
      minimap2Worker_desired.postMessage({
        fastqFile: ">desired\n" + desiredSeq,
        reference: reference,
        longjoinBandWidth: longjoinBandWidth,
        chainingBandWidth: chainingBandWidth,
        fileType: "String"
      });
  
      const desiredPromise = new Promise((resolve, reject) => {
        minimap2Worker_desired.onmessage = function(event) {
          if (event.data.error) {
            reject(event.data.error);
          } else {
            if (event.data.type === 0) {
              a = 1;
            }  else if (event.data.type === 1) {
              b = 1;
            } else if (event.data.type === 2) {
              resolve(event.data.result);
            }
          }
        };
      });
  
      try {
        desiredAlignSam = await desiredPromise;
        console.log("Worker 완료됨:", desiredAlignSam);
      } catch (error) {
        console.error("Worker 처리 중 에러 발생:", error);
      } finally {
        minimap2Worker_desired.terminate();
        runAlignFiles();
      }
    } else {
        desiredAlignSam = [];
        runAlignFiles();
    }
  }

async function runAlignFiles() {

  let minimap2Worker_treated = null;
  let minimap2Worker_control = null;

  minimap2Worker_treated = new Worker('minimap2_worker_prog.js');
  minimap2Worker_treated.postMessage({
      fastqFile : treatedFile,
      reference: reference,
      longjoinBandWidth : longjoinBandWidth,
      chainingBandWidth : chainingBandWidth,
      fileType : "file"
  });

  var progress_max = 1;

  if (controlFile !== false) {
    minimap2Worker_control = new Worker('minimap2_worker_prog.js');
    minimap2Worker_control.postMessage({
      fastqFile : controlFile,
      reference: reference,
      longjoinBandWidth : longjoinBandWidth,
      chainingBandWidth : chainingBandWidth,
      fileType : "file"
    });
    progress_max = 0.5;
  }

  const progressBarUpdate = document.getElementById('progress-bar');
  var alignProgressTreat = 0;
  var alignProgressControl = 0;

  const treatedPromise = new Promise((resolve, reject) => {
      minimap2Worker_treated.onmessage = function(event) {
      if (event.data.error) {
          console.error("Error:", event.data.error);
      } else {
          if (event.data.type === 0) {
          statusElem.innerHTML = "Treated file has " + event.data.fileLen + " lines";
          }  else if (event.data.type === 1) {
            alignProgressTreat = event.data.progress;
            progressBarUpdate.style.width = ((alignProgressControl + alignProgressTreat)*progress_max).toFixed(1) + '%';
            progressBarUpdate.textContent = ((alignProgressControl + alignProgressTreat)*progress_max).toFixed(1) + '%';
          console.log("Mapping result:", event.data.stdout);
          } else if (event.data.type === 2) {
              resolve(event.data.result);
          }
        }
      };
  });

  let controlPromise
  if (minimap2Worker_control) {
    controlPromise = new Promise((resolve, reject) => {
        minimap2Worker_control.onmessage = function(event) {
        if (event.data.error) {
            console.error("Error:", event.data.error);
        } else {
            if (event.data.type === 0) {
            statusElem.innerHTML = "Control file has " + event.data.fileLen + " lines";
            }  else if (event.data.type === 1) {
              alignProgressControl = event.data.progress;
              progressBarUpdate.style.width = ((alignProgressControl + alignProgressTreat)*progress_max).toFixed(1) + '%';
              progressBarUpdate.textContent = ((alignProgressControl + alignProgressTreat)*progress_max).toFixed(1) + '%';
            } else if (event.data.type === 2) {
            resolve(event.data.result);
            }
        }
      };
    });
  } else {
    controlPromise = Promise.resolve(null);
  }

  Promise.all([treatedPromise, controlPromise])
    .then(([x, y]) => {
    treatedAlignSam = x;
    controlAlignSam = y;
    console.log("두 워커 모두 완료됨");
    console.log("Treated result:", treatedAlignSam);
    console.log("Control result:", controlAlignSam);
    minimap2Worker_treated.terminate();
    if (controlFile == false) { 
      //minimap2Worker_control.terminate();
      controlAlignSam = [];
    }
    runMainAlgorithm();
    })
    .catch(error => {
    console.error("워커 처리 중 에러 발생:", error);
    });

}

function downloadText(variable, filename = "data.txt") {
  const blob = new Blob([String(variable)], { type: "text/plain" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

async function runMainAlgorithm() {

  const mainWorker = new Worker('lungo_worker.js', { type: 'module' });
  document.getElementById("progress_1").innerHTML =
  '<i class="bi bi-caret-right-fill"></i> Analyzing mutations ... <i class="bi bi-caret-left-fill"></i>';

  mainWorker.postMessage({
    controlSamString : controlAlignSam.join('\n'),
    treatedSamString : treatedAlignSam.join('\n'),
    inducedMutationString : desiredAlignSam.join('\n'),
    reference: reference,
    cvPos : cleavagePos,
    cvPos2 : cleavagePos2,
    window : windowRange,
    wholeWindow : wholeWindow,
    filter1 : true,
    windowFilter : true,
    alpha: 0.05,
    length_min: 10
  });

  mainWorker.onmessage = (event) => {
    if (event.data && event.data.type === 'error') {
      console.error("Worker에서 에러 발생:", event.data.message);
      console.error("Stack:", event.data.stack);
    } else {
      mainResult = event.data.result
      console.log("Worker 결과:", event.data);
      confirmInsertion();
    }
  };
  
  mainWorker.onerror = (event) => {
    console.error("Worker 에러 발생 (onerror):", event);
  };
}

async function confirmInsertion() {

  const insertWorker = new Worker('insert_worker.js');
  const insReference = '>ref\n' + reference + '\n>nCas9\ngacaagaagtacagcatcggcctggacatcggcaccaactctgtgggctgggccgtgatcaccgacgagtacaaggtgcccagcaagaaattcaaggtgctgggcaacaccgaccggcacagcatcaagaagaacctgatcggagccctgctgttcgacagcggcgaaacagccgaggccacccggctgaagagaaccgccagaagaagatacaccagacggaagaaccggatctgctatctgcaagagatcttcagcaacgagatggccaaggtggacgacagcttcttccacagactggaagagtccttcctggtggaagaggataagaagcacgagcggcaccccatcttcggcaacatcgtggacgaggtggcctaccacgagaagtaccccaccatctaccacctgagaaagaaactggtggacagcaccgacaaggccgacctgcggctgatctatctggccctggcccacatgatcaagttccggggccacttcctgatcgagggcgacctgaaccccgacaacagcgacgtggacaagctgttcatccagctggtgcagacctacaaccagctgttcgaggaaaaccccatcaacgccagcggcgtggacgccaaggccatcctgtctgccagactgagcaagagcagaaagctggaaaatctgatcgcccagctgcccggcgagaagaagaatggcctgttcggaaacctgattgccctgagcctgggcctgacccccaacttcaagagcaacttcgacctggccgaggatgccaaactgcagctgagcaaggacacctacgacgacgacctggacaacctgctggcccagatcggcgaccagtacgccgacctgtttctggccgccaagaacctgtccgacgccatcctgctgagcgacatcctgagagtgaacaccgagatcaccaaggcccccctgagcgcctctatgatcaagagatacgacgagcaccaccaggacctgaccctgctgaaagctctcgtgcggcagcagctgcctgagaagtacaaagagattttcttcgaccagagcaagaacggctacgccggctacattgacggcggagccagccaggaagagttctacaagttcatcaagcccatcctggaaaagatggacggcaccgaggaactgctcgtgaagctgaagagagaggacctgctgcggaagcagcggaccttcgacaacggcagcatcccccaccagatccacctgggagagctgcacgccattctgcggcggcaggaagatttttacccattcctgaaggacaaccgggaaaagatcgagaagatcctgaccttccgcatcccctactacgtgggccctctggccaggggaaacagcagattcgcctggatgaccagaaagagcgaggaaaccatcaccccctggaacttcgaggaagtggtggacaagggcgcttccgcccagagcttcatcgagcggatgaccaacttcgataagaacctgcccaacgagaaggtgctgcccaagcacagcctgctgtacgagtacttcaccgtgtataacgagctgaccaaagtgaaatacgtgaccgagggaatgagaaagcccgccttcctgagcggcgagcagaaaaaggccatcgtggacctgctgttcaagaccaaccggaaagtgaccgtgaagcagctgaaagaggactacttcaagaaaatcgagtgcttcgactccgtggaaatctccggcgtggaagatcggttcaacgcctccctgggcacataccacgatctgctgaaaattatcaaggacaaggacttcctggacaatgaggaaaacgaggacattctggaagatatcgtgctgaccctgacactgtttgaggacagagagatgatcgaggaacggctgaaaacctatgcccacctgttcgacgacaaagtgatgaagcagctgaagcggcggagatacaccggctggggcaggctgagccggaagctgatcaacggcatccgggacaagcagtccggcaagacaatcctggatttcctgaagtccgacggcttcgccaacagaaacttcatgcagctgatccacgacgacagcctgacctttaaagaggacatccagaaagcccaggtgtccggccagggcgatagcctgcacgagcacattgccaatctggccggcagccccgccattaagaagggcatcctgcagacagtgaaggtggtggacgagctcgtgaaagtgatgggccggcacaagcccgagaacatcgtgatcgaaatggccagagagaaccagaccacccagaagggacagaagaacagccgcgagagaatgaagcggatcgaagagggcatcaaagagctgggcagccagatcctgaaagaacaccccgtggaaaacacccagctgcagaacgagaagctgtacctgtactacctgcagaatgggcgggatatgtacgtggaccaggaactggacatcaaccggctgtccgactacgatgtggacgctatcgtgcctcagagctttctgaaggacgactccatcgacaacaaggtgctgaccagaagcgacaagaaccggggcaagagcgacaacgtgccctccgaagaggtcgtgaagaagatgaagaactactggcggcagctgctgaacgccaagctgattacccagagaaagttcgacaatctgaccaaggccgagagaggcggcctgagcgaactggataaggccggcttcatcaagagacagctggtggaaacccggcagatcacaaagcacgtggcacagatcctggactcccggatgaacactaagtacgacgagaatgacaagctgatccgggaagtgaaagtgatcaccctgaagtccaagctggtgtccgatttccggaaggatttccagttttacaaagtgcgcgagatcaacaactaccaccacgcccacgacgcctacctgaacgccgtcgtgggaaccgccctgatcaaaaagtaccctaagctggaaagcgagttcgtgtacggcgactacaaggtgtacgacgtgcggaagatgatcgccaagagcgagcaggaaatcggcaaggctaccgccaagtacttcttctacagcaacatcatgaactttttcaagaccgagattaccctggccaacggcgagatccggaagcggcctctgatcgagacaaacggcgaaaccggggagatcgtgtgggataagggccgggattttgccaccgtgcggaaagtgctgagcatgccccaagtgaatatcgtgaaaaagaccgaggtgcagacaggcggcttcagcaaagagtctatcctgcccaagaggaacagcgataagctgatcgccagaaagaaggactgggaccctaagaagtacggcggcttcgacagccccaccgtggcctattctgtgctggtggtggccaaagtggaaaagggcaagtccaagaaactgaagagtgtgaaagagctgctggggatcaccatcatggaaagaagcagcttcgagaagaatcccatcgactttctggaagccaagggctacaaagaagtgaaaaaggacctgatcatcaagctgcctaagtactccctgttcgagctggaaaacggccggaagagaatgctggcctctgccggcgaactgcagaagggaaacgaactggccctgccctccaaatatgtgaacttcctgtacctggccagccactatgagaagctgaagggctcccccgaggataatgagcagaaacagctgtttgtggaacagcacaagcactacctggacgagatcatcgagcagatcagcgagttctccaagagagtgatcctggccgacgctaatctggacaaagtgctgtccgcctacaacaagcaccgggataagcccatcagagagcaggccgagaatatcatccacctgtttaccctgaccaatctgggagcccctgccgccttcaagtactttgacaccaccatcgaccggaagaggtacaccagcaccaaagaggtgctggacgccaccctgatccaccagagcatcaccggcctgtacgagacacggatcgacctgtctcagctgggaggtgac\n>Scaffold\nGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC\n>DC_cluster\ncagagccagtgggctagcaccagaggcggacacaaggaacggcctaccgtgatgctcttcccattccgc\n>ITR\naggaacccctagtgatggagttggccactccctctctgcgcgctcgctcgctcactgaggccgggcgaccaaaggtcgcccgacgcccgggctttgcccgggcggcctcagtgagcgagcgagcgcgcag\n';

  insertWorker.postMessage({
    treated_read_to_mut: mainResult.treated_read_to_mut,
    reference: insReference,
  });

  const progressBarUpdate = document.getElementById('progress-bar');
  insertWorker.onmessage = (event) => {
    if (event.data && event.data.type === 'error') {
      console.error("Worker에서 에러 발생:", event.data.message);
      console.error("Stack:", event.data.stack);
    } else if (event.data.type == 1) {
      alignProgressControl = event.data.progress;
      progressBarUpdate.style.width = (alignProgressControl).toFixed(1) + '%' ;
      progressBarUpdate.textContent = (alignProgressControl).toFixed(1) + '%';
    } else if (event.data.type == 2) {
      mainResult.treated_read_to_mut = event.data.result
      filterMutations();
    }
  };
  
  insertWorker.onerror = (event) => {
    console.error("Worker 에러 발생 (onerror):", event);
  };
}

function filterMutations() {

  filteredSigMut = [];
  var key, val, cnt, info;
  for (val of Object.entries(mainResult.significant_keys)) {
    //if (val[1][0] < 0.002) {
    filteredSigMut.push(val[0]);
    //}
  }

  filteredTreatedReads = {};
  var liens_res = [];
  var filteredMutStr = '';
  var inducedMutPatCnt = mainResult.induced_mutations.length;

  for (val of Object.entries(mainResult.treated_read_to_mut)) {
    var filteredMutStr = '';
    cnt = 1

    for (i of val[1][0]) {
      if (filteredSigMut.includes(i)) {
        if (i.includes('Del')) {
          if (Number(i.split('_')[2]) > largeDelLen) {
            i = i.replace('Del', 'LargeDel');
          }
        } else if (i.includes('Ins')) {
          if (i.split('_')[2].length > largeInsLen) {
            i = i.replace('Ins', 'LargeIns');
          }
        }
        if (val[1][2][0] == true) {
          var st = Number(i.split('_')[0]);
          var ed = Number(i.split(':')[0].split('_')[1]);
          for (x of val[1][2].slice(1,)) {
            if (x.includes('inversion')) {
              var inv_st = Number(x.split('_')[1]);
              var inv_ed = Number(x.split('_')[2]);
              if (Math.abs(st - inv_st) < 5 && Math.abs(ed - inv_ed) < 5) {
                i += '_inv';
              }
            }
          }
        }
        filteredMutStr += i + ',';
      } else {

      }
    }


    
    liens_res.push(val[0] + ': ' + filteredMutStr);

    if (filteredMutStr == "") {
      filteredMutStr = "WT";
      mutTypeCnt["WT"] += cnt;
    } else {
      var mutType = "";
      if (val[1][2][0] === true && val[1][2][1].includes('inversion')) {
        mutType = 'Inv';
      } else if (filteredMutStr.indexOf("Large") != -1) {
        if (filteredMutStr.indexOf("LargeIns") != -1 && filteredMutStr.indexOf("LargeDel") != -1) {
          mutType = "Complex";
        } else if (filteredMutStr.indexOf("LargeDel") != -1) {
          mutType = "LargeDel";
        } else if (filteredMutStr.indexOf("LargeIns") != -1) {
          mutType = "LargeIns";
        }
      } else if (filteredMutStr.indexOf("Ins") != -1 && filteredMutStr.indexOf("Del") != -1) {
        mutType = "Complex";
      } else if (filteredMutStr.indexOf("Del") != -1) {
        mutType = "Del";
      } else if (filteredMutStr.indexOf("Ins") != -1) {
        mutType = "Ins";
      } else {
        if (filteredMutStr.indexOf("Sub") != -1) {
          mutType = "Sub";
        } else {
          mutType = "Complex";
        }
      }
      mutTypeCnt[mutType] += cnt;
      mainResult.treated_read_to_mut[val[0]].push(filteredMutStr);
      mainResult.treated_read_to_mut[val[0]].push(mutType);
      filteredMutStr = filteredMutStr.slice(0,-1); //+ "_" + mutType;
    }
    preciseMutCnt = 0
    for (mut of filteredMutStr.split(',')) {
      if (mainResult.induced_mutations.includes(mut)) {
        preciseMutCnt += 1
      }
    }
    if (preciseMutCnt / inducedMutPatCnt >= preciseEditThreshold) {
      filteredMutStr += "_PreciseEdit";
      mutTypeCnt["Precise"] += cnt;
    } else if (preciseMutCnt / inducedMutPatCnt >= partialEditThreshold) {
      filteredMutStr += "_PartialEdit";
      mutTypeCnt["Partial"] += cnt;
    }
    console.log()

    filteredTreatedReads[filteredMutStr] = (filteredTreatedReads[filteredMutStr] || 0) + cnt;

  }

  var sortedEntries = Object.entries(filteredTreatedReads).sort((a, b) => b[1] - a[1]);
  filteredTreatedReads = Object.fromEntries(sortedEntries);

  var all_cnt = 0
  for (key in filteredTreatedReads) {
    cnt = filteredTreatedReads[key];
    all_cnt += cnt
    for (i of key.split(',')) {
      i = i.replace(':','_').replace('>','_').split('_');
      if (i[0] == 'WT') {
        continue
      } else if (i[2] == 'Sub') {
        for (x = 0; x < i[3].length; x ++) {
          subPosList[Number(i[0])+x] += cnt;
          subProportionList[Number(i[0])+x]['ATGCN'.indexOf(i[4][x])] += cnt;
        }
      } else if (i[2].indexOf('Del') != -1) {
        delPlotList[0].push(Number(i[0]));
        delPlotList[1].push(Number(i[3]));
        delPlotList[2].push((cnt/mainResult.treated_align_cnt.Used).toFixed(2))
        for (x = 0; x < i[3]; x ++) {
          delPosList[Number(i[0])+x] += cnt;
          subProportionList[Number(i[0])+x][5] += cnt;
        }
      } else if (i[2].indexOf('Ins') != -1) {
        insPlotList[0].push(Number(i[0]));
        insPlotList[1].push(i[3].length);
        insPlotList[2].push((cnt/mainResult.treated_align_cnt.Used).toFixed(2))
        insPosList[Number(i[0])] += cnt;
        subProportionList[Number(i[0])][6] += cnt;
      }
    }
  }

  for (i = 0; i < reference.length; i++) {
    subPosList[i] = (subPosList[i]*100/all_cnt).toFixed(2);
    insPosList[i] = (insPosList[i]*100/all_cnt).toFixed(2);
    delPosList[i] = (delPosList[i]*100/all_cnt).toFixed(2);
  }

  document.getElementById("progressPage").classList.remove("active");
  document.getElementById("resultsPage").classList.add("active"); 
  document.body.style.backgroundColor = "#ffffff";

  runVisualization()

}


function drawAllelePlot() {

  var countSum = Object.values(filteredTreatedReads).reduce((acc, cur) => acc + cur, 0);
  var xlimMin = 0;
  var width, height, maxLength;

  function drawD3(drawInfo, xPos, section) {
    var subMarker;

    if (drawInfo[0].length != 0) {
      for (i of drawInfo[0][1][0]) {
        base = i[0];
        rowG.append("rect")
          .attr("x", xPos)
          .attr("y", 0)
          .attr("width", 13)
          .attr("height", 17.5)
          .attr("fill", backgroundColor[base]);

        rowG.append("text")
          .attr("x", xPos + cellSizeWidth/2)
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "middle")
          .attr("dominant-baseline", "middle")
          .style("fill", 'black')
          .text(base);

        xPos += cellSizeWidth;
      }
      if (drawInfo[0][2] = true) {
        rowG.append("rect")
          .attr("x", xPos)
          .attr("y", 0)
          .attr("width", 2)
          .attr("height", 17.5)
          .attr("class", "out-ins")
          .attr("stroke", "red")      
          .attr("fill", "none")
          .attr("stroke-width", 2); 
      }
    } else {
      xPos += 10 * cellSizeWidth;
    }
    for (x=0; x<drawInfo[1].length; x++) {
      i = drawInfo[1][x];
      base = i[0];
      var largeDelSt, largeDelEd;
      if (i[0] == '-' && i.split('_')[1] * 1 > 100) {
        largeDelSt = 0;
       while(i[0] == '-' && x < plotWindow*2) {
          i = drawInfo[1][x];
          x += 1
        }
        x-=1
        
        largeDelEd = x;
        if (section == 0) {
          largeDelEd -= plotWindow + 1;
        }

        rowG.append("line")
          .attr("x1", xPos + cellSizeWidth * (largeDelSt))   
          .attr("y1", 0 + 0.5*cellSizeHeight)   
          .attr("x2", xPos + cellSizeWidth * (largeDelEd+5))   
          .attr("y2", 0 + 0.5*cellSizeHeight)  
          .attr("stroke", "black")              
          .attr("stroke-width", 1);
        xPos += cellSizeWidth * largeDelEd;
        i = drawInfo[1][x];
        base = i[0];
        if (x == plotWindow*2 - 1) {
          continue
        }
      } 
      
      if (i[1] == 'S'){
        subMarker = 'bold';
      } else {
        subMarker = 'normal';
      }

      rowG.append("rect")
        .attr("x", xPos)
        .attr("y", 0)
        .attr("width", 13)
        .attr("height", 17.5)
        .attr("fill", backgroundColor[base]);

      rowG.append("text")
        .attr("x", xPos + cellSizeWidth/2)
        .attr("y", cellSizeHeight*0.6)
        .attr("font-size", 12)
        .attr("class", "cell-text")
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "middle")
        .style("font-weight", subMarker)
        .style("fill", 'black')
        .text(base);
      
      if (i[1] == 'I' && ( x == 0 || drawInfo[1][x-1][1] !='I')) {
        var InsSt = x;
        var InsEd;
        for (var InsEd = x; InsEd < drawInfo[1].length; InsEd++) {
          if (drawInfo[1][InsEd][1] != 'I') {
            break
          }
        }
        
        rowG.append("rect")
          .attr("x", xPos+1)
          .attr("y", 0)
          .attr("width", 13 * (InsEd - InsSt)+2)
          .attr("height", 17.5)
          .attr("stroke", "red")      
          .attr("fill", "none")
          .attr("stroke-width", 2); 
      }

      xPos += cellSizeWidth;
      
    }  
    if (drawInfo[2].length != 0) {
      if (drawInfo[2][2] = true) {
        rowG.append("rect")
          .attr("x", xPos-2)
          .attr("y", 0)
          .attr("width", 2)
          .attr("height", 17.5)
          .attr("class", "out-ins")
          .attr("stroke", "red")      
          .attr("fill", "none")
          .attr("stroke-width", 2); 
      }
      for (i of drawInfo[2][1][0]) {
        base = i[0];
        rowG.append("rect")
          .attr("x", xPos)
          .attr("y", 0)
          .attr("width", 13)
          .attr("height", 17.5)
          .attr("fill", backgroundColor[base]);

        rowG.append("text")
          .attr("x", xPos + cellSizeWidth/2)
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "middle")
          .attr("dominant-baseline", "middle")
          .style("fill", 'black')
          .text(base);

        xPos += cellSizeWidth;
      }
      
    }
  }


  function drawLine(mutList, read, per) {
    var aligned = reference.split('');
    var checkMiddleLD = false;
    var checkMiddleInv = false;
    for (i = 0; i < mutList.length; i++) {
      mutInfo = mutList[i].replace(';Complex','').replace(':', '_').split('_');
      if (mutList[i].includes('inv') == true) {
        checkMiddleInv = true;
      }
      m = mutInfo[2];
      var mutSeq = '';
      var pos = 0;
      if (m == 'Sub') {
        mutSeq = mutInfo[3].split('>')[1];
        pos = +mutInfo[0];
        for (x = 0; x < mutSeq.length; x++) {
          aligned[pos + x] = mutSeq[x] + 'S';
        }
      } else if (m == 'Del' || m == 'LargeDel') {
        mutSeq = +mutInfo[3];
        pos = +mutInfo[0];
        if (pos < cleavagePos + plotWindow && mutInfo[1]*1 > cleavagePos2 - plotWindow) {
          checkMiddleLD = mutSeq;
        }
        for (x = 0; x < mutSeq; x++) {
          aligned[pos+x] = '-_' + mutSeq + '_' + (x);
          if (mutList[i].includes('inv') == true) {
            aligned[pos+x] += '_inv';
          } 
        }
      } else if (m == 'Ins' || m == 'LargeIns' || mutList[i].includes('inv') == false) {
        mutSeq = mutInfo[3];
        pos = +mutInfo[0];
        aligned[pos] = aligned[pos]+'I'+mutSeq;
      }
    }
    
    var drawInfo_1 = [[], [], [], [false, false]]; //[front_info, middle, back_info, LD]
    var m, mDel;
    for (i = 0; i < plotWindow*2; i ++) {
      pos = cleavagePos + i - plotWindow+1;
      m = aligned[pos];  
      if (i == 0 && m[0] == '-') {
        var mDel = m.split('_');
        drawInfo_1[0] = [[],[],[]]
        drawInfo_1[0][0] = mDel[2];
        drawInfo_1[0][1].push(aligned.slice(pos - drawInfo_1[0][0] - 10 ,pos - drawInfo_1[0][0]));
        drawInfo_1[0][2] = false;
        if (m.indexOf('I') != -1) {
          drawInfo_1[0][2] = true;
        }
        drawInfo_1[3][0] = true;
      }
      if (m[0] != '-') {
        drawInfo_1[1].push(m.slice(0,2));
        if (m[1] == 'I') {
          for (x of m.slice(2,)){
            drawInfo_1[1].push(x + 'I');
          }
        }
      } else {
        if (m.indexOf('inv') == -1 && m.indexOf('I') != -1) {
          for (x of m.slice(m.indexOf('I')+1,)) {
            drawInfo_1[1].push(x + 'I');
          }
        }
        drawInfo_1[1].push(m);
      }
      drawInfo_1[1] = drawInfo_1[1].slice(0,plotWindow*2);
      if (i == plotWindow*2-1 && m[0] == '-') {
        var mDel = m.split('_');
        drawInfo_1[2] = [[],[],[]]
        drawInfo_1[2][0] = mDel[2] - mDel[1];
        drawInfo_1[2][1].push(aligned.slice(pos + drawInfo_1[2][0], pos - drawInfo_1[2][0] + 10));
        drawInfo_1[2][2] = false;
        if (m.indexOf('I') != -1) {
          drawInfo_1[2][2] = true;
        }
        drawInfo_1[3][1] = true;
      }
    }

    if (cleavagePos2 != "") {
      var drawInfo_2 = [[], [], [], [false, false]];
      for (i = 0; i < plotWindow*2; i ++) {
        pos = cleavagePos2 + i - plotWindow;
        m = aligned[pos];  
        if (i == 0 && m[0] == '-') {
          var mDel = m.split('_');
          drawInfo_2[0] = [[],[],[]];
          drawInfo_2[0][0] = mDel[2];
          drawInfo_2[0][1].push(aligned.slice(pos-drawInfo_2[0][0] - 10 ,pos - drawInfo_2[0][0]));
          drawInfo_2[0][2] = false;
          if (m.indexOf('I') != -1) {
            drawInfo_1[0][2] = true;
          }
          drawInfo_2[3][0] = true;
        }
        if (m[0] != '-') {
          drawInfo_2[1].push(m.slice(0,2));
          if (m[1] == 'I') {
            for (x of m.slice(2,)){
              drawInfo_2[1].push(x + 'I');
            }
          }
        } else {
          if (m.indexOf('I') != -1) {
            for (x of m.slice(m.indexOf('I')+1,)) {
              drawInfo_2[1].push(x + 'I');
            }
          }
          drawInfo_2[1].push(m);
        }
        drawInfo_2[1] = drawInfo_2[1].slice(0,plotWindow*2);
        if (i == plotWindow*2-1 && m[0] == '-') {
          var mDel = m.split('_');
          drawInfo_2[2] = [[],[],[]];
          drawInfo_2[2][0] = mDel[1] - mDel[2];
          drawInfo_2[2][1].push(aligned.slice(pos + drawInfo_2[2][0], pos + drawInfo_2[2][0] + 10));
          drawInfo_2[2][2] = false;
          if (m.indexOf('I') != -1) {
            drawInfo_2[2][2] = true;
          }
          drawInfo_2[3][1] = true;
        }
      }
      
      if (cleavagePos2 != "") {
        var middleCigar = [];
        var middleCigarSum = {"S": 0, "D": 0, "I": 0};
        for (pos = cleavagePos-plotWindow; pos < cleavagePos2+plotWindow; pos++) {
          mutInfo = aligned[pos];
          if (mutInfo[0] != '-') {
            if (mutInfo[1] == 'S') {
              if (middleCigar.slice(-1) != 'S') {
                middleCigar.push(0);
                middleCigar.push('S');
              }
              middleCigar[middleCigar.length-2] += 1;
              middleCigarSum['S'] += 1;
            } else {
              if (middleCigar.slice(-1) != 'M') {
                middleCigar.push(0);
                middleCigar.push('M');
              }
              middleCigar[middleCigar.length-2] += 1
            }
            if (mutInfo[1] == 'I') {
              middleCigar.push(mutInfo.slice(2,).length);
              middleCigar.push('I');
              middleCigarSum['I'] += mutInfo.slice(2,).length;
            }
          } else {
            if (mutInfo[1] == 'I') {
              middleCigar.push(mutInfo.slice(2,).length);
              middleCigar.push('I');
              middleCigarSum['I'] += mutInfo.slice(2,).length;
            }
            if (middleCigar.slice(-1) != 'D') {
              middleCigar.push(0);
              middleCigar.push('D');
            }
            middleCigar[middleCigar.length-2] += 1;
            middleCigarSum['D'] += 1;
          }
        }
        if (checkMiddleLD === false) {
          var middleCigarStr = '';
          //middleCigarStr = middleCigar.join('');
          for (x of 'SID') {
            if (middleCigarSum[x] != 0) {
              middleCigarStr += middleCigarSum[x] + x;
            }
          }
          rowG.append("rect")
            .attr("x", cellSizeWidth * (plotWindow*2 + 10))
            .attr("y", 2)
            .attr("width", cellSizeWidth*10)
            .attr("class", "middle-cell-text")
            .attr("height", 13.5)
            .attr("fill", '#F5F5F5');

          rowG.append("text")
            .attr("x", cellSizeWidth * (plotWindow*2 + 10) + cellSize*(10-2)/2)
            .attr("y", cellSizeHeight*0.6)
            .attr("font-size", 12)
            .attr("class", "middle-cell-text")
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .text(middleCigarStr);
        } else {
          var inputStr = checkMiddleLD + ' bp';
          if (checkMiddleInv) {
            inputStr += ' (inv)'
          }
          rowG.append("line")
            .attr("x1", cellSizeWidth * (10 + plotWindow*2))   
            .attr("y1", 0 + 0.5*cellSizeHeight)   
            .attr("x2", cellSizeWidth * (10 + plotWindow*2 + 10))   
            .attr("y2", 0 + 0.5*cellSizeHeight)  
            .attr("class", "middle-cell-text")
            .attr("stroke", "black")              
            .attr("stroke-width", 1);
          var textElement = rowG.append("text")
            .attr("x", cellSizeWidth * (plotWindow*2 + 10) + cellSize*(10-2)/2)
            .attr("y", cellSizeHeight*0.6)
            .attr("font-size", 12)
            .attr("class", "middle-cell-text")
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .text(inputStr);
          bbox = textElement.node().getBBox();
          textWidth = bbox.width + 20;
          rowG.append("rect")
              .attr("x", cellSizeWidth * (plotWindow*2 + 10) + cellSize*(10-2)/2 - textWidth/2)
              .attr("y", 0)
              .attr("width", textWidth)
              .attr("class", "middle-cell-text")
              .attr("height", 17.5)
              .attr("fill", 'white');
          rowG.append("text")
              .attr("x", cellSizeWidth * (plotWindow*2 + 10) + cellSize*(10-2)/2)
              .attr("y", cellSizeHeight*0.6)
              .attr("font-size", 12)
              .attr("class", "middle-cell-text")
              .attr("text-anchor", "middle")
              .attr("dominant-baseline", "middle")
              .text(inputStr);
        }
        drawInfo_1[2] = [];
        drawInfo_2[0] = [];
      }
      xPos = drawD3(drawInfo_1, 0, 0);
      xPos = drawD3(drawInfo_2, cellSizeWidth * (plotWindow*2 + 10), 1);
      
      if (read != '') {
        rowG.append("text")
          .attr("x", cellSizeWidth * (plotWindow*4 + 10*3 + 1))
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "left")
          .attr("dominant-baseline", "left")
          .text(per + ' % (' + read + ' reads)');
      }
    } else {
      xPos = drawD3(drawInfo_1, 0, 1);
      
      if (read != '') {
        rowG.append("text")
          .attr("x", cellSizeWidth * (plotWindow*2 + 10 + 1))
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "left")
          .attr("dominant-baseline", "left")
          .text(per + ' % (' + read + ' reads)');
      }
    }
  }
   
  var plotWidth
  windowStart = cleavagePos - plotWindow;
  if (numRows > Object.values(filteredTreatedReads).length) {
    numRows = Object.values(filteredTreatedReads).length;
  }
  if (showAllBetweenAllele) {
    windowEnd = cleavagePos2 + plotWindow;
    maxLength = cleavagePos2 - cleavagePos + 2*plotWindow;
    plotWidth = cellSizeWidth * maxLength + margin.left + margin.right;
  } else {
    windowEnd = cleavagePos + plotWindow;
    if (cleavagePos2) {
      maxLength = 4*plotWindow + 10;
      windowStart2 = cleavagePos2 - plotWindow;
      windowEnd2 = cleavagePos2 + plotWindow;
    } else {
      maxLength = 2*plotRef;
    }
    plotWidth = cellSizeWidth * maxLength + margin.left + margin.right;
  }

  var plotRef = reference.slice(cleavagePos - plotWindow + 1, windowEnd+1);
  if (cleavagePos2 != "" && showAllBetweenAllele == false) {
    var plotRef2 = reference.slice(cleavagePos2 - plotWindow, cleavagePos2 + plotWindow);
  }

  var backgroundColor = {
      "A": "#F0FFF0",
      "T": "#FFE4E1",
      "G": "#FFFFE0",
      "C": "#F0F8FF",
      "-": "#F5F5F5",
      "N": "#F5F5F5",
      ' ': "#FFFFFF"
  };


  
  var width = parseInt(d3.select("#allelePlot").style("width"));
  height = cellSizeHeight * (numRows+5) + margin.top  + margin.bottom ; 

  var initialScale = width/plotWidth * 0.8;
  var translateX = (width - width * initialScale);
  var translateY = (height - height * initialScale) / 2;


  const svg = d3.select("#allelePlot")
    .append("svg")
    .attr("viewBox", `0 0 ${width} ${height+100}`)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .style("width", "100%") 
    .style("height", "auto");

  const mainG = svg.append("g")
      .attr("class", "main-container");

 var zoom = d3.zoom()
    .scaleExtent([0.1, 10]) 
    .on("zoom", (event) => {
        mainG.attr("transform", event.transform);
        d3.select("#zoomSlider").property("value", event.transform.k);
    });

  svg.call(zoom);

  svg.call(zoom.transform, d3.zoomIdentity.translate(0, translateY).scale(initialScale));
  
  d3.select("#zoomSlider").property("value", initialScale);
  

  d3.select("#zoomSlider").on("input", function() {
    const newScale = +this.value;
    var currentTransform = d3.zoomTransform(svg.node());
    var newTransform = d3.zoomIdentity
                            .translate(currentTransform.x, currentTransform.y)
                            .scale(newScale);
    svg.call(zoom.transform, newTransform);
  })
 

  var xPos = 10 *cellSizeWidth;
  var betweenX = 0;
  var seq = plotRef;

  if (cleavagePos2 != "" && showAllBetweenAllele == false) {
    betweenX = 10;
  }

  var rowG = mainG.append("g")
  .attr("transform", `translate(${margin.left}, ${margin.top + 1 * cellSizeHeight})`);

  // Reference line

  rowG.append("text")
      .attr("x", 8 * cellSizeWidth)
      .attr("y", cellSizeHeight*0.6)
      .attr("font-size", 16)
      .attr("font-weight", "bold")
      .attr("font-family", "sans-serif")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "middle")
      .text("Ref.");

  for (i = 0; i < plotRef.length; i ++ ) {

      let base = plotRef[i];

      rowG.append("rect")
      .attr("x", xPos)
      .attr("y", 0)
      .attr("width", 13)
      .attr("height", 17.5)
      .attr("fill", backgroundColor[base]);

      rowG.append("text")
          .attr("x", xPos + cellSizeWidth/2)
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "middle")
          .attr("dominant-baseline", "middle")
          .style("fill", 'black')
          .text(base);

      xPos += cellSizeWidth;
  }

  if (cleavagePos2 != "" && showAllBetweenAllele == false) {

    rowG.append("rect")
      .attr("x", xPos)
      .attr("y", 2)
      .attr("width", cellSizeWidth*10)
      .attr("height", 13.5)
      .attr("fill", '#F5F5F5');
    
    rowG.append("text")
      .attr("x", xPos + cellSize*4)
      .attr("y", cellSizeHeight*0.6)
      .attr("font-size", 12)
      .attr("class", "cell-text")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "middle")
      .text((cleavagePos2 - cleavagePos - plotWindow*2) + " bp");

    xPos += cellSizeWidth * 10;

    for (i = 0; i < plotRef2.length; i ++ ) {

      let base = plotRef2[i];

      rowG.append("rect")
      .attr("x", xPos)
      .attr("y", 0)
      .attr("width", 13)
      .attr("height", 17.5)
      .attr("fill", backgroundColor[base]);

      rowG.append("text")
          .attr("x", xPos + cellSizeWidth/2)
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 12)
          .attr("class", "cell-text")
          .attr("text-anchor", "middle")
          .attr("dominant-baseline", "middle")
          .style("fill", 'black')
          .text(base);

      xPos += cellSizeWidth;
    }
  }
  
  // Target visual line
  rowG = mainG.append("g").attr("transform", `translate(${margin.left}, ${margin.top + 2 * cellSizeHeight})`);

  xPos = 0;

  if (strand == 1) {
    rowG.append("rect")
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(target.length - cleavagePosTarget - 1))
      .attr("y", 2)
      .attr("width", cellSizeWidth * 20)
      .attr("height", 13.5)
      .attr("fill", '#C9C9C9');

    rowG.append("text")
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(target.length - cleavagePosTarget - 2))
      .attr("y", cellSizeHeight*0.6)
      .attr("font-size", 14)
      .attr("font-family", "sans-serif")
      .attr("font-weight", "bold")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "middle")
      .text("sg1");
  } else {
    rowG.append("rect")
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(cleavagePosTarget))
      .attr("y", 2)
      .attr("width", cellSizeWidth * target.length)
      .attr("height", 13.5)
      .attr("fill", '#C9C9C9');

    rowG.append("text")
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(cleavagePosTarget - 2))
      .attr("y", cellSizeHeight*0.6)
      .attr("font-size", 14)
      .attr("font-family", "sans-serif")
      .attr("font-weight", "bold")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "middle")
      .text("sg1");
  }
  if (cleavagePos2 != "") {
    if (strand2 == 1) {
      rowG.append("rect")
        .attr("x", 10 *cellSizeWidth + cellSizeWidth*(plotWindow*3 + 10 - cleavagePosTarget))
        .attr("y", 2)
        .attr("width", cellSizeWidth * target.length)
        .attr("height", 13.5)
        .attr("fill", '#C9C9C9');

      rowG.append("text")
        .attr("x", 10 *cellSizeWidth + cellSizeWidth*(plotWindow*3 + 10 - cleavagePosTarget - 2))
        .attr("y", cellSizeHeight*0.6)
        .attr("font-size", 16)
        .attr("font-family", "sans-serif")
        .attr("font-weight", "bold")
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "middle")
        .text("sg2");
    } else {
      rowG.append("rect")
        .attr("x", 10 *cellSizeWidth + cellSizeWidth*(plotWindow*2 + 10 + cleavagePosTarget))
        .attr("y", 2)
        .attr("width", cellSizeWidth * 20)
        .attr("height", 13.5)
        .attr("fill", '#C9C9C9');

      rowG.append("text")
        .attr("x", 10 *cellSizeWidth + cellSizeWidth*(plotWindow*2 + 10 + cleavagePosTarget - 3))
        .attr("y", cellSizeHeight*0.6)
        .attr("font-size", 16)
        .attr("font-family", "sans-serif")
        .attr("font-weight", "bold")
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "middle")
        .text("sg2");
    }
  }

  yPos = 3;
  rowG = mainG.append("g").attr("transform", `translate(${margin.left}, ${margin.top + yPos * cellSizeHeight})`);

  drawLine(mainResult.induced_mutations, '', '');
  yPos += 2

  const totalRead = Object.values(filteredTreatedReads).reduce((acc, cur) => acc + cur, 0);
  var n = 0;
  for (mutStr in filteredTreatedReads) {
    if (n == numRows) {
      break
    }
    n ++;
    read = filteredTreatedReads[mutStr];
    per = read/totalRead*100;
    per = per.toFixed(2);
    rowG = mainG.append("g").attr("transform", `translate(${margin.left}, ${margin.top + yPos * cellSizeHeight})`);
    drawLine(mutStr.split(','), read, per);
    yPos += 1;
  }

  mainG.append("line")
    .attr("x1", cellSizeWidth*(+plotWindow + 10)+margin.left)   
    .attr("y1", 0 + margin.top)   
    .attr("x2", cellSizeWidth*(plotWindow + 10)+margin.left)   
    .attr("y2", cellSizeHeight*(numRows+6) + margin.top)  
    .attr("stroke", "black")              
    .attr("stroke-width", 2)              
    .attr("stroke-dasharray", "5,5");

  if (cleavagePos2 != "") {
    mainG.append("line")
      .attr("x1", cellSizeWidth*(3*plotWindow + 10+10+1)+margin.left)   
      .attr("y1", 0 + margin.top)   
      .attr("x2", cellSizeWidth*(3*plotWindow + 10+10+1)+margin.left)   
      .attr("y2", cellSizeHeight*(numRows+6) + margin.top)  
      .attr("stroke", "red")              
      .attr("stroke-width", 2)              
      .attr("stroke-dasharray", "5,5");
  }
  d3.selectAll(".middle-cell-text").raise();
  d3.selectAll(".out-ins").raise();  
}


function runVisualization() {

  //Draw Table

  var table = document.createElement("table");
  table.className = "table table-bordered table-striped w-100 result-table";

  var headers = ["Total reads", "Reads used in analysis", "WT", "Ins", "Del", "Sub", "LargeDel", "LargeIns", "Inversion", "Complex",  "Precisely induced editing", "Partially induced editing"];

  var colgroup = document.createElement("colgroup");
  for (let i = 0; i < headers.length; i++) {
    var col = document.createElement("col");
    col.style.width = (100 / headers.length) + "%"; // 각 열 너비를 균등하게 분배
    colgroup.appendChild(col);
  }
  table.appendChild(colgroup);

  var thead = document.createElement("thead");
  thead.className = "thead-light";
  var headerRow = document.createElement("tr");

  headers.forEach(headerText => {
    var th = document.createElement("th");
    th.textContent = headerText;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);
  table.appendChild(thead);

  var tbody = document.createElement("tbody");
  var dataRow = document.createElement("tr");

  var data = [
    mainResult.treated_align_cnt.All_reads,
    mainResult.treated_align_cnt.Used,
    mutTypeCnt.WT + ' (' + (mutTypeCnt.WT*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Ins + ' (' + (mutTypeCnt.Ins*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Del + ' (' + (mutTypeCnt.Del*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Sub + ' (' + (mutTypeCnt.Sub*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.LargeDel + ' (' + (mutTypeCnt.LargeDel*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.LargeIns + ' (' + (mutTypeCnt.LargeIns*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Inv + ' (' + (mutTypeCnt.Inv*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Complex + ' (' + (mutTypeCnt.Complex*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Precise + ' (' + (mutTypeCnt.Precise*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)',
    mutTypeCnt.Partial + ' (' + (mutTypeCnt.Partial*100/mainResult.treated_align_cnt.Used).toFixed(1) + ' %)'
  ];

  data.forEach(value => {
    const td = document.createElement("td");
    td.textContent = value;
    dataRow.appendChild(td);
  });
  tbody.appendChild(dataRow);
  table.appendChild(tbody);

  document.getElementById("cnt-table-container").appendChild(table);

  // pvalue distribution
  var threshold = Math.log10(0.002);
  var x = [];
  for (i of Object.entries(mainResult.significant_keys)) {
    i = i[2];
    i += 0.0000000001;
    x.push(Math.log10(i));
  }
  var histogram = d3.histogram()
  .domain([d3.min(x), d3.max(x)])
  .thresholds(30);
  var bins = histogram(x);
  var binCenters = bins.map(function(bin) {
    return (bin.x0 + bin.x1) / 2;
  });
  var counts = bins.map(function(bin) {
      return bin.length;
  });
  var trace = {
      x: binCenters,
      y: counts,
      mode: 'lines',
      type: 'scatter',
      name: 'Frequency Polygon'
  };
  var layout = {
    title: 'p-value distribution',
    xaxis: { title: 'p-value (log)' },
    yaxis: { title: 'Accumulation count'},
    shapes: [
      {
        type: 'line',
        x0: threshold,
        x1: threshold,
        y0: 0,
        y1: 1,
        xref: 'x',
        yref: 'paper',  
        line: {
          color: 'red',
          width: 2,
          dash: 'dash'
        }
      }
    ],
    annotations: [
      {
        x: threshold,
        y: 1,  
        xref: 'x',
        yref: 'paper',
        text: 'Threshold',
        showarrow: false,
        xanchor: 'left',
        font: {
          color: 'red',
          size: 12
        }
      }
    ]
  };
  //Plotly.newPlot('pvalDistPlot', [trace], layout);


  // freq distribution
  var threshold = 0;
  var x = [];
  for (i of Object.entries(mainResult.significant_keys)) {
    i = i[2];
    x.push(i);
  }
  var histogram = d3.histogram()
  .domain([d3.min(x), d3.max(x)])
  .thresholds(30);
  var bins = histogram(x);
  var binCenters = bins.map(function(bin) {
    return (bin.x0 + bin.x1) / 2;
  });
  var counts = bins.map(function(bin) {
      return bin.length;
  });
  var trace = {
      x: binCenters,
      y: counts,
      mode: 'lines',
      type: 'scatter',
      name: 'Frequency Polygon'
  };
  var layout = {
    title: 'frequency distribution',
    xaxis: { title: 'frequency' },
    yaxis: { title: 'Accumulation count'},
    shapes: [
      {
        type: 'line',
        x0: threshold,
        x1: threshold,
        y0: 0,
        y1: 1,
        xref: 'x',
        yref: 'paper',  
        line: {
          color: 'red',
          width: 2,
          dash: 'dash'
        }
      }
    ],
      }
  //Plotly.newPlot('freqDistPlot', [trace], layout);


  // Align summary plot

  var alignSum = mainResult.treated_align_cnt;
  var data = {'All reads':alignSum.All_reads, 'Aligned':alignSum.All_reads-alignSum.Unmapped, 'Used':alignSum.Used};
  var keys = Object.keys(data);
  var values = Object.values(data);

  var trace = {
    x: keys,
    y: values,
    type: 'bar',
    marker: {color: 'gray'}
  };
  var layout = {
    title: 'Alignment Summary',
    xaxis: { title: 'Reads' },
    yaxis: { title: '' }
  };
  Plotly.newPlot('alignSumPlot', [trace], layout);

 
  // Mutation count pie plot

  var labels = Object.keys(mutTypeCnt);
  var values = Object.values(mutTypeCnt);
  var customColors = ['#CECFD1', '#4357A7', '#ED2024', '#BF71AE', '#98D9ED', '#FABF37', '#6BBD45', '#000000'];
  var trace = {
    labels: labels,
    values: values,
    type: 'pie',
    marker: { colors: customColors }
  };

  var layout = {
    title: 'Mutation pattern pie chart'
  };

  Plotly.newPlot('allPiePlot', [trace], layout);


  var editedCnt = 0;
  for (i of ["Ins", "Del", "Sub", "LargeIns", "LargeDel", "Inv", "Complex"]) {
    editedCnt += mutTypeCnt[i];
  }
  var trace = {
    labels: ['Non-edited', 'Edited'],
    values: [mutTypeCnt.WT, editedCnt],
    type: 'pie',
    marker: { colors: customColors }
  };

  var layout = {
    title: 'Mutation Pie Chart'
  };

  Plotly.newPlot('mutPiePlot', [trace], layout);

  // Mutation distribution line plot

  var x = Array.from({length: reference.length}, (_, i) => i);
  var trace1 = { x, y: insPosList, mode: 'lines', name: 'Ins' };
  var trace2 = { x, y: delPosList, mode: 'lines', name: 'Del' };
  var trace3 = { x, y: subPosList, mode: 'lines', name: 'Sub' };
  
  var data = [trace1, trace2, trace3];
  
  var layout = {
    title: 'Mutation distributin line plot',
    xaxis: {
      title: '',
      rangeslider: { visible: true },
      type: 'linear'
    },
    yaxis: {
      title: 'Frequency (%)',
      range: [0, 100]
    }
  };
  
  Plotly.newPlot('combinedPlot', data, layout);

  // CV line

  shapes = [
    {
      type: 'line',
      xref: 'x',         
      yref: 'paper',     
      x0: cleavagePos,
      y0: 0,             
      x1: cleavagePos,
      y1: 1,             
      line: {
        color: 'black',  
        width: 2,        
        dash: 'dash'     
      }
    }
  ]

  if (cleavagePos2 != '') {
    shapes.push({
        type: 'line',
        xref: 'x',         
        yref: 'paper',     
        x0: cleavagePos2,
        y0: 0,             
        x1: cleavagePos2,
        y1: 1,             
        line: {
          color: 'red',  
          width: 2,        
          dash: 'dash'     
        }
    })
  }


  //Insertion distribution plot

  var trace = {
    x: insPlotList[0],
    y: insPlotList[1],
    mode: 'markers',
    type: 'scatter',
    marker: {            
      color: insPlotList[2],        
      colorscale: 'Viridis',         
      showscale: true               
    }
  };
  
  var maxy = insPlotList[1].reduce((max, current) => {
    return current > max ? current : max;
  }, -Infinity);


  var layout = {
    xaxis: { title: 'Insertion length' },
    yaxis: { title: 'Inseriton position' },
    shapes: shapes
  };
  
  Plotly.newPlot('insertionsPlot', [trace], layout);


  //Insertion distribution plot

  var trace = {
    x: delPlotList[0],
    y: delPlotList[1],
    mode: 'markers',
    type: 'scatter',
    marker: {            
      color: delPlotList[2],        
      colorscale: 'Viridis',         
      showscale: true               
    }
  };
  
  var maxy = delPlotList[1].reduce((max, current) => {
    return current > max ? current : max;
  }, -Infinity);


  var layout = {
    xaxis: { title: 'Deletion position' },
    yaxis: { title: 'Deletion length' },
    shapes: shapes
  };
  
  Plotly.newPlot('deletionsPlot', [trace], layout);


  //Substitution plot

  var windowStart = cleavagePos - plotWindow
	var windowEnd = cleavagePos + plotWindow

	if (cleavagePos2 != '') { windowEnd = cleavagePos2 + plotWindow; } 

	var baseProportions = {'A': [], 'T': [], 'G': [], 'C': [], 'N': [], 'Del': []};
  var refProportions = {'A': [], 'T': [], 'G': [], 'C': [], 'N': []}
  var positions = [];

	
	for (i=0; i<reference.length; i++) {
		if (windowStart > i || windowEnd <= i) {
			continue
    }

    var sumCnt = 0;

    var refnt = reference[i];
    for (x of 'ATGCN') {
      if (refnt != x) {
        baseProportions[x].push(subProportionList[i]['ATGCN'.indexOf(x)] * 100 / mainResult.treated_align_cnt.Used);
        refProportions[x].push(0);
        sumCnt += subProportionList[i]['ATGCN'.indexOf(x)] * 100 / mainResult.treated_align_cnt.Used
      }
    }
    positions.push(i);
    baseProportions['Del'].push(subProportionList[i][5] * 100 / mainResult.treated_align_cnt.Used);
		baseProportions[reference[i]].push(100 - sumCnt);
    refProportions[reference[i]].push(-1);
  
  }
	
	var baseColors = {'A': 'honeydew',
		'T': 'mistyrose',
		'G': 'lightyellow',
		'C': 'aliceblue',
    'N': 'lightgray',
    'Del': 'black'
  };

  var data = [];
  var minx =  Math.min(...positions);
  var maxx =  Math.max(...positions);

  for (i of ['A', 'T', 'G','C', 'N']){
    data.push({
      x: positions,
      y: baseProportions[i],
      type: 'bar',
      name: i,
      color: baseColors[i],
      offset: 0,
    });
  }

	var layout = {
    barmode: 'stack',
    bargap: 0, 
    xaxis: {
      title: 'Position',
      rangeslider: {
        visible: true
      },
      autorange: false,
      range: [minx, maxx],
    },
    yaxis: {
      title: '%'
    }
  };

  Plotly.newPlot('substitutionsPlot', data, layout)

  drawAllelePlot();
}
