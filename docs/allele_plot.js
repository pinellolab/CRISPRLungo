
var cleavagePos = 1979;
var cleavagePos2 = 2086;
var strand = 1;
var strand2 = -1;
var windowRange = 5;
var target = "cacttcaatctggcccctct";
var target2 = "ggcttagatatgacctgggt";
var cleavagePosTarget = 17


const reference = 'CTCTTACCGCCCTTTTCCGGGGCAAGGGAAGCTAGTAGCGGAGCCGGAAGTGAGGCACCCTCGGGCTCGAGACAGCGGCGACGTTTAAAGCTGAGCGACCCAGTGCCACTGGAGACGGTCAGCTTCTCCACTCAGGCTCCTCCAGCCCGAGCCAGAAGACCCCCTCCCCCAGAATTCTGGGGGCCGATGGAAGGGAGCCGAGTCAGATCGCGAGGTACCCAGAGCCGACAGACCGGAGCGACAGGGAGTTGCCAGAAGCCCCGCCCCTAGGAGTGATCGGAAAGCCTCACCCATCCGGGTGAGGAACCCGGAGGGACCGCCTCCGGGCGGAGCCCGCCGACCATGGCTACGCCCCTGGTGGCGGGTCCCGCAGCTCTACGCTTCGCCGCCGCGGCTAGCTGGCAGGTTGTGCGCGGACGCTGCGTGGAACATTTTCCGCGAGTACTGGAGTTTCTGCGATCTCTGCGCGCTGTTGCCCCTGGCTTGGTTCGCTACCGGCACCACGAACGCCTTTGTATGGGCCTAAAGGCCAAGGTATTGGAGCAAGTAGGACCTGGAAGGGGAAAGAAAGAAAGGGGCGGGATGCGAGCACCTGACAAGGGCCCACAGCTTGGGGATGAGACTAGTTGGAAAAGGTCAGCTTGCCCCGAAGTGGGGATTCTGAATTCAGTACCGTACAGGTGGTGGTGGAGCTGATCCTGCAGGGCCGGCCTTGGGCCCAAGTCCTGAAAGCCCTGAATCACCACTTTCCAGAATCTGGACCTATAGTGCGGGATCCCAAGGCTGTGAGTAATCCCCGGAACAAGCCCTGACCCCAGTTACACTTGGTGCAGCAAAGTTGCTCCTCCTGTCTACTGGATGTTGGTGCTAACTTCTCTGTCTTCTGTATCCTCAACAGACAAAGCAGGATCTGAGGAAGATTTTGGAGGCACAGGAAACTTTTTACCAGCAGGTGAAGCAGCTGTCAGAGGCTCCTGTGGATTTGGCCTCGAAGCTGCAGGTGAGACTGGTTTGAAGGCTATTATGTGGCTATTTTCTCTAACCCATTTTTTTTTTTTTTTGAGGAATCTTGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCACAATCTCGGCTCACTGCAACCTCTGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACAAATGCCCAGCTAATTTTTGTATTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAAGTGATCCGCCCACCTCGGCCTCCCAAATGGCCAGGATTACAGGCGTGAGCCACCGTGCCTGGCCACTAACCCACTTTTGACCTAATCTCCCCTGCTGTAGGAACTTGAACAAGAGTATGGGGAACCCTTTCTGGCTGCCATGGAAAAGCTGCTTTTTGAGTACTTGTGTCAGCTGGAGAAAGCACTGCCTACACCGCAGGCACAGCAGGTTTCGCTGGGGCAAAACATGTAAGGGCAGAGTGTGGGGTGGCTGTTTCAAGGATTGTGGTCTTCATACCTCTGTCCCTGCCCACAGCTTCAGGATGTGCTGAGTTGGATGCAGCCTGGAGTCTCTATCACCTCTTCTCTTGCCTGGAGACAATATGGTGTGGACATGGGGTGGCTGCTTCCAGGTACTAGGAATTTGGAGGTGTAGTGTTTAGCCTGAGACCTTTTGAGGCAGTCCACTGGAATAGTTCCATGGGCTCCGGGCATAAGAAACCAGTTATTCTGTTACTATCTGTTCTCTTTATAAGGGGCCTCATTTTCTTTTTCAGAGTGCTCTGTTACTGACTCAGTGAACCTGGCTGAGCCCATGGAACAGAATCCTCCTCAGCAACAAAGACTAGCACTCCACAATCCCCTGCCAAAAGCCAAGCCTGGCACACATCTTCCTCAGGGACCATCTTCAAGGACGCACCCAGAACCTCTAGCTGGCCGACACTTCAATCTGGCCCCTCTAGGCCGACGAAGAGTTCAGTCCCAATGGGCCTCCACTAGGGGAGGCCATAAGGAGCGCCCCACAGTCATGCTGTTTCCCTTTAGGAATCTCGGCTCACCAACCCAGGTCATATCTAAGCCTGAGAGCAAGGAAGAACATGCGATATACACAGCAGACCTAGCCATGGGCACAAGAGCAGCCTCCACTGGGAAGTCTAAGAGTCCATGCCAGACCCTGGGGGGAAGGGCTCTGAAGGAGAACCCAGTTGACTTGCCTGCCACAGAGCAAAAGGAGTGAGTGGAACAGAGTTGCTTCTTACTAGGAGCACATTCTTTGCCTGCCTTCCCTTCATCCTATCCTCTTTGCTTGCTCTCACCTCAGGAATTGCTTGGATTGCTACATGGACCCCCTGAGACTATCATTATTACCTCCTAGGGCCAGGAAGCCAGGTAGGTAGTCTGAGTCAGGATTGGATCAACAGCCTCCTCTCTTGGGGACTCTCAAGAGCCTGTGTTCATCTAGAAGTAGTAGTTTGATTCTGGTTTCCCTCCTACAGTGTGTCCTCCGTCTCTGTGCAGCTCCGTCATTACCATAGGGGACTTGGTTTTAGACTCTGATGAGGAAGAAAATGGCCAGGGGGAAGGAAAGGTGAGTGGGAAGGAGCAGAAAGCTGGGAAAGGGGATGGGTAGAACAAGACTGAGAAATCCACATGCTTCAGAATTCAGAGGGTTCAGGGAATGGTTTCGGATAGTAGGCTCTCCCTGCTCCCTTCTCTACAGGAATCTCTGGAAAACTATCAGAAGACAAAGTTTGACACCTTGATACCCACTCTCTGTGAATACCTACCCCCTTCTGGCCACGGTGCCATACCTGTTTCTTCCTGTGACTGTAGAGACAGTTCTAGACCTTTGTGATAGAACTAAAATGCTCTCTGTACTCTAGTCTCCTGCCTCCTCAGCTCTGCAAGTAGTTTAGTAGGAATGAAGTGGAAGTCCAGGCTTGGATTGCCTAACTACACTGCTAAAAATATTTGTAATCCTTAATAATTAAACTTTGGATTTGTTAAAATACTGCCTTGATGTGAAGGAGGG';
const filteredTreatedReads = {'114_2085:Del_1972,114_2085:Ins_GTGCAACCTGGGAAGCCCTACA;Complex': 35, 
'290_2471:Del_2182,290_2471:Ins_CTCCTTTTGCTGTTTACCATCATCGTCTGGTCATTTGATTTTACTGTGCCCAGCTAAATAAGACTCATTTCATGTGTCCCAAGAGCAAAGATCCCCAACTTTAATGGGAGTAATGATCTCCTGGGAAGCTGTGCAAAATGAAAAATGCAGATGTGCTGTCCCATTCCTCAGAGCCTATTTTTAGGTATGGCTGGGGCTCAGGAGTCTGCACTGTGCAGCAATCTAGGAGATGGTCAGATGGCCCTGGCCCACACTTGGAGTCTCAGA_Complex':2,
'949_2085:Del_1137,949_2085:Ins_GGTTGCGCAACCTGGGAAGCCCTACA,2087_2089:Sub_AGG>TCT,2091_2095:Sub_CATAT>ACCGC_Complex': 2
};
const desiredSeq = 'ctcttaccgcccttttccggggcaagggaagctagtagcggagccggaagtgaggcaccctcgggctcgagacagcggcgacgtttaaagctgagcgacccagtgccactggagacggtcagcttctccactcaggctcctccagcccgagccagaagaccccctcccccagaattctgggggccgatggaagggagccgagtcagatcgcgaggtacccagagccgacagaccggagcgacagggagttgccagaagccccgcccctaggagtgatcggaaagcctcacccatccgggtgaggaacccggagggaccgcctccgggcggagcccgccgaccatggctacgcccctggtggcgggtcccgcagctctacgcttcgccgccgcggctagctggcaggttgtgcgcggacgctgcgtggaacattttccgcgagtactggagtttctgcgatctctgcgcgctgttgcccctggcttggttcgctaccggcaccacgaacgcctttgtatgggcctaaaggccaaggtattggagcaagtaggacctggaaggggaaagaaagaaaggggcgggatgcgagcacctgacaagggcccacagcttggggatgagactagttggaaaaggtcagcttgccccgaagtggggattctgaattcagtaccgtacaggtggtggtggagctgatcctgcagggccggccttgggcccaagtcctgaaagccctgaatcaccactttccagaatctggacctatagtgcgggatcccaaggctgtgagtaatccccggaacaagccctgaccccagttacacttggtgcagcaaagttgctcctcctgtctactggatgttggtgctaacttctctgtcttctgtatcctcaacagacaaagcaggatctgaggaagattttggaggcacaggaaactttttaccagcaggtgaagcagctgtcagaggctcctgtggatttggcctcgaagctgcaggtgagactggtttgaaggctattatgtggctattttctctaacccatttttttttttttttgaggaatcttgctctgtcgcccaggctggagtgcagtggcacaatctcggctcactgcaacctctgcctcccgggttcaagcgattctcctgcctcagcctcccgagtagctgggattacaaatgcccagctaatttttgtattttagtagagatggggtttcaccatgttggccaggctggtctcgaactcctgacctcaagtgatccgcccacctcggcctcccaaatggccaggattacaggcgtgagccaccgtgcctggccactaacccacttttgacctaatctcccctgctgtaggaacttgaacaagagtatggggaaccctttctggctgccatggaaaagctgctttttgagtacttgtgtcagctggagaaagcactgcctacaccgcaggcacagcaggtttcgctggggcaaaacatgtaagggcagagtgtggggtggctgtttcaaggattgtggtcttcatacctctgtccctgcccacagcttcaggatgtgctgagttggatgcagcctggagtctctatcacctcttctcttgcctggagacaatatggtgtggacatggggtggctgcttccaggtactaggaatttggaggtgtagtgtttagcctgagaccttttgaggcagtccactggaatagttccatgggctccgggcataagaaaccagttattctgttactatctgttctctttataaggggcctcattttctttttcagagtgctctgttactgactcagtgaacctggctgagcccatggaacagaatcctcctcagcaacaaagactagcactccacaatcccctgccaaaagccaagcctggcacacatcttcctcagggaccatcttcaaggacgcacccagaacctctagctggccgacacttcaatctggccccattgggaagacggagagtgcagagccagtgggctagcaccagaggcggacacaaggaacggcctaccgtgatgctcttcccattccgcaacctgggaagccctacacaggtcatatctaagcctgagagcaaggaagaacatgcgatatacacagcagacctagccatgggcacaagagcagcctccactgggaagtctaagagtccatgccagaccctggggggaagggctctgaaggagaacccagttgacttgcctgccacagagcaaaaggagtgagtggaacagagttgcttcttactaggagcacattctttgcctgccttcccttcatcctatcctctttgcttgctctcacctcaggaattgcttggattgctacatggaccccctgagactatcattattacctcctagggccaggaagccaggtaggtagtctgagtcaggattggatcaacagcctcctctcttggggactctcaagagcctgtgttcatctagaagtagtagtttgattctggtttccctcctacagtgtgtcctccgtctctgtgcagctccgtcattaccataggggacttggttttagactctgatgaggaagaaaatggccagggggaaggaaaggtgagtgggaaggagcagaaagctgggaaaggggatgggtagaacaagactgagaaatccacatgcttcagaattcagagggttcagggaatggtttcggatagtaggctctccctgctcccttctctacaggaatctctggaaaactatcagaagacaaagtttgacaccttgatacccactctctgtgaatacctacccccttctggccacggtgccatacctgtttcttcctgtgactgtagagacagttctagacctttgtgatagaactaaaatgctctctgtactctagtctcctgcctcctcagctctgcaagtagtttagtaggaatgaagtggaagtccaggcttggattgcctaactacactgctaaaaatatttgtaatccttaataattaaactttggatttgttaaaatactgccttgatgtgaaggaggg';
const showAllBetweenAllele = false;
const induced_mutations = ['2019_2019:Sub_T>C', '2076_2079:Sub_CTCA>AAGC', '1986_1987:Sub_CC>AA', '2085_2085:Sub_C>A', '1980_1981:Sub_TC>AT', '2073_2073:Sub_C>G', '2040_2040:Sub_C>G', '2002_2003:Sub_TC>AG', '2013_2015:Sub_CTC>TAG', '2037_2037:Sub_G>A', '2049_2049:Sub_C>G', '2082_2082:Sub_A>T', '2043_2043:Sub_C>T', '2055_2055:Sub_G>C', '2028_2028:Sub_C>A', '2070_2070:Sub_T>C', '2007_2007:Sub_A>G', '2022_2022:Sub_G>A', '1998_1998:Sub_T>G', '1983_1983:Sub_A>G', '2025_2025:Sub_A>C', '2031_2031:Sub_T>C', '2046_2046:Sub_A>C', '1992_1992:Sub_A>G', '2061_2061:Sub_C>A', '2058_2058:Sub_T>C', '2064_2065:Sub_TA>CC', '2067_2067:Sub_G>C'];
const inducedMutationStr = induced_mutations.join(',');
const plotWindow = 20;


var windowStart, windowEnd, windowStart2, windowEnd2;

const margin = { top: 40, right: 20, bottom: 20, left: 100 };
const cellSize = 16; // base 하나당 폭
const cellSizeHeight = 17;
const cellSizeWidth = 13;
var numRows = 10;

const countSum = Object.values(filteredTreatedReads).reduce((acc, cur) => acc + cur, 0);

var xlimMin = 0;
var width, height, maxLength;

function drawAllelePlot() {

  function drawD3(drawInfo, xPos) {
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
        largeDelEd = 0;
       while(i[0] == '-' && x < plotWindow*2) {
          i = drawInfo[1][x];
          x += 1;
          largeDelEd += 1;
        }
        x-=1
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
      
      if (i[1] == 'I') {
        rowG.append("rect")
          .attr("x", xPos)
          .attr("y", 0)
          .attr("width", 15)
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
    for (i = 0; i < mutList.length; i++) {
      mutInfo = mutList[i].replace(';Complex','').replace(':', '_').split('_');
      m = mutInfo[2];
      var mutSeq = '';
      var pos = 0;
      if (m == 'Sub') {
        mutSeq = mutInfo[3].split('>')[1];
        pos = +mutInfo[0];
        for (x = 0; x < mutSeq.length; x++) {
          aligned[pos + x] = mutSeq[x] + 'S';
        }
      } else if (m == 'Del') {
        mutSeq = +mutInfo[3];
        pos = +mutInfo[0];
        if (pos < cleavagePos + plotWindow && mutInfo[1]*1 > cleavagePos2 - plotWindow) {
          checkMiddleLD = mutSeq;
        }
        for (x = 0; x < mutSeq; x++) {
          aligned[pos+x] = '-_' + mutSeq + '_' + (x);
        }
      } else if (m == 'Ins') {
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
        if (m.indexOf('I') != -1) {
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
        drawInfo_1[2][1].push(aligned.slice(pos - drawInfo_1[2][0], pos - drawInfo_1[2][0] + 10));
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
          drawInfo_2[2][1].push(aligned.slice(pos + drawInfo_2[2][0]+1, pos + drawInfo_2[2][0]+1 + 10));
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
            .attr("font-size", 8)
            .attr("class", "middle-cell-text")
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .text(middleCigarStr);
        } else {
          var inputStr = checkMiddleLD + ' bp';
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
            .attr("font-size", 8)
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
              .attr("font-size", 8)
              .attr("class", "middle-cell-text")
              .attr("text-anchor", "middle")
              .attr("dominant-baseline", "middle")
              .text(inputStr);
        }
        drawInfo_1[2] = [];
        drawInfo_2[0] = [];
      }
      xPos = drawD3(drawInfo_1, 0);
      xPos = drawD3(drawInfo_2, cellSizeWidth * (plotWindow*2 + 10));
      
      if (read != '') {
        rowG.append("text")
          .attr("x", cellSizeWidth * (plotWindow*4 + 10*3 + 1))
          .attr("y", cellSizeHeight*0.6)
          .attr("font-size", 8)
          .attr("class", "cell-text")
          .attr("text-anchor", "left")
          .attr("dominant-baseline", "left")
          .text(per + ' % (' + read + ' reads)');
      }

    }
  }
   
  windowStart = cleavagePos - plotWindow;
  if (numRows > Object.values(filteredTreatedReads).length) {
    numRows = Object.values(filteredTreatedReads).length;
  }
  if (showAllBetweenAllele) {
    windowEnd = cleavagePos2 + plotWindow;
    maxLength = cleavagePos2 - cleavagePos + 2*plotWindow;
    width = cellSizeWidth * maxLength + margin.left + margin.right;
  } else {
    windowEnd = cleavagePos + plotWindow;
    if (cleavagePos2) {
      maxLength = 4*plotWindow + 10;
      windowStart2 = cleavagePos2 - plotWindow;
      windowEnd2 = cleavagePos2 + plotWindow;
    } else {
      maxLength = 2*plotRef;
    }
    width = cellSizeWidth * maxLength + margin.left + margin.right;
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


  

  height = cellSizeHeight * (numRows+5) + margin.top  + margin.bottom ; 

  var initialScale = 0.7;
  var translateX = (width - width * initialScale) / 2;
  var translateY = (height - height * initialScale) / 2;


  const svg = d3.select("body")
    .append("svg")
    .attr("width", width)
    .attr("height", height+300);

  const mainG = svg.append("g")
      .attr("class", "main-container");

 var zoom = d3.zoom()
    .scaleExtent([0.1, 10]) 
    .on("zoom", (event) => {
        mainG.attr("transform", event.transform);
        d3.select("#zoomSlider").property("value", event.transform.k);
    });

  svg.call(zoom);

  svg.call(zoom.transform, d3.zoomIdentity.translate(translateX, translateY).scale(initialScale));
  
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
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(target.length - cleavagePosTarget))
      .attr("y", 2)
      .attr("width", cellSizeWidth * 20)
      .attr("height", 13.5)
      .attr("fill", '#C9C9C9');

    rowG.append("text")
      .attr("x", 10 *cellSizeWidth + cellSizeWidth*(target.length - cleavagePosTarget - 2))
      .attr("y", cellSizeHeight*0.6)
      .attr("font-size", 16)
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
      .attr("font-size", 16)
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

  drawLine(inducedMutationStr.split(','), '', '');
  yPos += 2

  const totalRead = Object.values(filteredTreatedReads).reduce((acc, cur) => acc + cur, 0);
  for (mutStr in filteredTreatedReads) {
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
      .attr("x1", cellSizeWidth*(3*plotWindow + 10+10)+margin.left)   
      .attr("y1", 0 + margin.top)   
      .attr("x2", cellSizeWidth*(3*plotWindow + 10+10)+margin.left)   
      .attr("y2", cellSizeHeight*(numRows+6) + margin.top)  
      .attr("stroke", "red")              
      .attr("stroke-width", 2)              
      .attr("stroke-dasharray", "5,5");
  }
  d3.selectAll(".middle-cell-text").raise();
  d3.selectAll(".out-ins").raise();


  
}
  
drawAllelePlot();
