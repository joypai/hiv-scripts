# USAGE: python position_analysis.py aa_alignment.fasta binding_positions.txt

import sys
import os.path
import operator
from collections import Counter, defaultdict
import xlsxwriter
from Bio import AlignIO

COLORS = ["#57a927","#7436c8","#ffe981","#cf45c4","#86953b","#7562da","#da682d","#7a82d8","#b48539",
          "#565160","#d5383b","#3f997e","#da4281","#4c6324","#c373b8","#894624","#518fc7","#d67373",
          "#863376","#95364b","#c4bfd6","#c4bfd6"]
AA = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','-','X']

def make_chart(workbook, worksheet, start, colors, timepoint, offset, end, total):
    chart = workbook.add_chart({'type':'doughnut'})
    chart.set_title({'name':timepoint})
    chart.add_series({
        'categories': [worksheet.get_name(),start,0,end,0],
        'values': [worksheet.get_name(),start,1,end,1],
        'points': colors
    })

    chart.set_size({'width':350, 'height':288})

    worksheet.insert_chart(start, 4+offset, chart)
    worksheet.insert_textbox(start+8,6+offset, str(total),
                             {'font': {'size':12}, 'width': 60, 'height':30, 'align':{'horizontal':'left'}, 'border':{'none':True}})


dir, sample_name = os.path.split(sys.argv[1])
if dir == "": dir = '.'
sample_name = sample_name.split('_')[0]
outf = open(dir+"/"+sample_name+"_aa_table.txt","w")
workbook = xlsxwriter.Workbook(dir+"/"+sample_name+"_position_analysis.xlsx")
bold = workbook.add_format({'bold':True})


# read alignment
alignment = AlignIO.read(sys.argv[1],"fasta")
print("num sequences:", len(alignment))

# get important positions
pos = [ int(x) for x in open(sys.argv[2]).read().strip().split(',') ]

# split sequences by time point
sample_ids = [ sample.id.replace('-','_') for sample in alignment ]
outf.write("pos\t"+"\t".join(sample_ids[1:])+"\n")
tp_idx = defaultdict(list)
for i, s in enumerate(sample_ids[1:]): # start from 1 to skip HXB2
    tp = s.split('_')[1]
    tp_idx[tp].append(i)
timepoints = sorted(tp_idx.keys())

# convert alignment column number to HXB2 residue number
res_num = []
idx = 1
for x in alignment[0].seq:
    if x != '-':
        res_num.append(idx)
        idx += 1
    else:
        res_num.append(-1)


for p in pos:
    row = 0

    # convert residue number to alignment column
    col = res_num.index(p)
    all_aa = alignment[1:,col]      # 1 to skip HXB2
    outf.write(str(p)+"\t"+"\t".join(all_aa) + "\n")

    # amino acid breakdown for all samples
    counts = Counter(all_aa)

    if len(counts) > 1 or p in (332,333,334,325):
        worksheet = workbook.add_worksheet("^"+str(p))
    else:
        continue

    #worksheet = workbook.add_worksheet("^"+str(p)) if len(counts) > 1 else workbook.add_worksheet(str(p))
    worksheet.write_row(row,0,['all timepoints'], bold)
    worksheet.write_row(row+1,0,['amino_acid','count','frequency'])       # header
    offset = 0
    row += 2
    top_row = row
    all_colors, aa_list = [], []
    for i in sorted(counts.items(), key=operator.itemgetter(1), reverse=True):
        worksheet.write_row(row, 0, list(i)+[i[1]/len(all_aa)])
        aa_list.append(i[0])
        all_colors.append({'fill': {'color': COLORS[AA.index(i[0])]}, 'border':{'color':'#606060'}})
        row += 1
    make_chart(workbook, worksheet, top_row, all_colors, 'all timepoints', offset, row-1, len(all_aa))
    tp_charts_needed = True if len(counts) > 1 else False       # no need to make chart for each timepoint if position has no variation

    row += 4
    # amino acid breakdown by timepoint
    for t in timepoints:
        colors = []
        offset = 7 if offset == 0 else 0
        tp_aa = [ all_aa[x] for x in tp_idx[t] ]
        counts = Counter(tp_aa)
        #print(t,tp_aa)

        worksheet.write_row(row,0,[t], bold)
        worksheet.write_row(row+1,0,['amino_acid','count','frequency'])       # header
        row += 2
        top_row = row
        for i in sorted(counts.items(), key=operator.itemgetter(1), reverse=True):
            worksheet.write_row(row, 0, list(i)+[i[1]/len(tp_aa)])
            colors.append(all_colors[aa_list.index(i[0])])
            row += 1
        if tp_charts_needed:
            make_chart(workbook, worksheet, top_row, colors, t, offset, row-1, len(tp_aa))

        row += 4

workbook.close()
