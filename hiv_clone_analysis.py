import re
import sys
import os.path
import xlsxwriter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CLONE_COLORS = ['#FF0000','#FF6666','#FF6600','#FFCC99','#FFFF00','#FFFFCC', # red, orange, yellow
          '#008000','#00FF00', '#CCFFCC', '#00FFFF','#99CCFF','#0000FF', # green, blue
          '#800080','#660099','#FF00FF','#CC99FF','#FF99CC','#E0E0E0']*4                      # purple, pink, grey

#SOURCES = ['D-14','D14','W23','W31']
#SOURCES = ['D14','W23']
SOURCES = sys.argv[2:]
rebound_idx = 2
REBOUND = SOURCES[rebound_idx:]
seq_dict = {}
timepoint_counts = { x:0 for x in SOURCES }
unique_seqs = []

# output files
filename, _ = os.path.splitext(sys.argv[1])
outf_counts = filename + "_counts.xlsx"
outf_unique_fasta = filename + "_unique.fasta"

workbook = xlsxwriter.Workbook(outf_counts)
worksheet = workbook.add_worksheet(filename)
chart_worksheet = workbook.add_worksheet('charts')


def make_chart(workbook, worksheet, start, colors, timepoint, offset, end, total):
    chart = workbook.add_chart({'type':'doughnut'})
    chart.set_title({'name':timepoint})
    chart.add_series({
        'categories': [worksheet.get_name(),start,0,end,0],
        'values': [worksheet.get_name(),start,1,end,1],
        'points': colors
    })

    chart.set_size({'width':350, 'height':288})
    #chart.set_legend({'layout': {'x':0.8, 'y': 0.2, 'width': 0.3, 'height':0.7}})

    worksheet.insert_chart(start, 4+offset, chart)
    worksheet.insert_textbox(start+8,6+offset, str(total),
                             {'font': {'size':12}, 'width': 60, 'height':30, 'align':{'horizontal':'left'}, 'border':{'none':True}})

# read sequences from fasta and identify identical sequences
pat_seq_counts = {}
with open(sys.argv[1]) as fh:
    for record in SeqIO.parse(fh, "fasta"):
        if record.id.find('HXB2') != -1:    # disregard HXB2 from identical sequences
            unique_seqs.append(record)
            continue

        # extract timepoint from sequence id
        #source = re.search(".*[-_]({src})[-_].*".format(src="|".join(SOURCES)), record.id).group(1)
        match = re.search("(?P<patient>\d*)[-_](?P<source>{src})[-_].*".format(src="|".join(SOURCES)), record.id)
        source = match.group('source')
        patient = match.group('patient')
        print(record.id, "\t", source, "\t", patient)

        if (patient, source) in pat_seq_counts:
            pat_seq_counts[(patient, source)] += 1
        else:
            pat_seq_counts[(patient, source)] = 1

        # count number of occurrences for each unique sequence, by timepoint
        seq = str(record.seq).replace('-','')       # get rid of gaps
        if seq not in seq_dict:
            seq_dict[seq] = { x:0 for x in SOURCES }
            seq_dict[seq]['names'] = []             # keep track of original sequence names

        seq_dict[seq][source] += 1
        seq_dict[seq]['names'].append(record.id)

        """
        seq = record.seq.replace('-','')
        if str(record.seq) not in seq_dict:
            seq_dict[str(record.seq)] = { x:0 for x in SOURCES }
            seq_dict[str(record.seq)]['names'] = []             # keep track of original sequence names

        seq_dict[str(record.seq)][source] += 1
        seq_dict[str(record.seq)]['names'].append(record.id)
        """

        timepoint_counts[source] += 1

print('number of unique env sequences:', len(seq_dict))


print(pat_seq_counts)

# report clusters and timepoint breakdown
worksheet.write_row(0, 0, ["seq_id","total_num_seqs"]+SOURCES+["sequences_in_clone","sequence"])

count = 1
chart_count = 0
clone_counts = { x:0 for x in SOURCES }
for s in sorted(seq_dict, key=lambda x: len(seq_dict[x]['names']), reverse=True):
    new_sid = "seq"+str(count)+"."+".".join([k+"_"+str(v) for k,v in seq_dict[s].items() if v != 0 and k != "names"])
    unique_seqs.append(SeqRecord(Seq(s),id=new_sid,description=""))

    data = [new_sid,len(seq_dict[s]['names'])]+[seq_dict[s][x] for x in SOURCES]+[",".join(seq_dict[s]['names']), s]
    worksheet.write_row(count, 0, data)
    count += 1

    if len(seq_dict[s]['names']) > 1: # clone
        for t in SOURCES:
            clone_counts[t] += seq_dict[s][t]

count += 1

# report percent of sequences in clones by timepoint
for t in SOURCES:
    worksheet.write_row(count, 0, ["# in clones",t,clone_counts[t]])
    worksheet.write_row(count+1, 0, ["{0} in clones".format("%"),t,(clone_counts[t]/timepoint_counts[t])])
    count+=2
worksheet.write_row(count, 0, ["# in clones", "all timepoints", sum([clone_counts[x] for x in clone_counts])])
worksheet.write_row(count+1, 0, ["{0} in clones".format("%"), "all timepoints", sum([clone_counts[x] for x in clone_counts])/sum([timepoint_counts[x] for x in timepoint_counts])])


# write chart
start_row, offset = 0, 0
for t in SOURCES:
    chart_worksheet.write_row(start_row,0,[t,'# sequences', '% of timepoint'])
    start_row += 1
    clones, counts, colors, ppts, patients = [], [], [], [], []
    singles = 0
    for i, s in enumerate(sorted(seq_dict, key=lambda x: len(seq_dict[x]['names']), reverse=True)):
        #if seq_dict[s][t] > 1:
        if len(seq_dict[s]['names']) > 1 and seq_dict[s][t] != 0:
            clones.append('clone'+str(i+1))
            counts.append(seq_dict[s][t])
            colors.append({'fill': {'color': CLONE_COLORS[i]}, 'border':{'color':'#606060'}})

            patient = seq_dict[s]['names'][0].split('_')[0]
            ppt = seq_dict[s][t] / pat_seq_counts[(patient,t)]
            ppts.append(ppt)
            patients.append(patient)
        elif seq_dict[s][t] == 1:
            singles += 1

    clones.append('singles')
    counts.append(singles)
    colors.append({'fill': {'color': 'white'}, 'border':{'color':'#606060'}})
    chart_worksheet.write_column(start_row, 0, clones)
    chart_worksheet.write_column(start_row, 1, counts)
    chart_worksheet.write_column(start_row, 2, ppts)
    chart_worksheet.write_column(start_row, 3, patients)

    make_chart(workbook, chart_worksheet, start_row, colors, t, offset, start_row+len(clones)-1, timepoint_counts[t])

    offset = 7 if offset == 0 else 0

    start_row += len(clones) + 5

# combined chart with all pre-rebound timepoints
pre_clones, pre_counts, pre_colors = [], [], []
pre_singles = 0

# combined chart with all timepoints
chart_worksheet.write_row(start_row,0,["all timepoints", '# sequences'])
start_row += 1
clones, counts, colors = [], [], []
singles = 0
for i, s in enumerate(sorted(seq_dict, key=lambda x: len(seq_dict[x]['names']), reverse=True)):
    if len(seq_dict[s]['names']) > 1:
        clones.append('clone'+str(i+1))
        counts.append(len(seq_dict[s]['names']))
        colors.append({'fill': {'color': CLONE_COLORS[i]}, 'border':{'color':'#606060'}})
    else:
        singles += 1

    # pre rebound
    pre_rebound_sum = sum([seq_dict[s][t] for t in SOURCES[:rebound_idx]])
    #if pre_rebound_sum > 1:
    if len(seq_dict[s]['names']) > 1 and pre_rebound_sum != 0:
        pre_clones.append('clone'+str(i+1))
        pre_counts.append(pre_rebound_sum)
        pre_colors.append({'fill': {'color': CLONE_COLORS[i]}, 'border':{'color':'#606060'}})
    elif len(seq_dict[s]['names']) == 1 and pre_rebound_sum == 1:
        pre_singles += 1

clones.append('singles')
counts.append(singles)
colors.append({'fill': {'color': 'white'}, 'border':{'color':'#606060'}})
chart_worksheet.write_column(start_row, 0, clones)
chart_worksheet.write_column(start_row, 1, counts)

make_chart(workbook, chart_worksheet, start_row, colors, 'all timepoints', offset, start_row+len(clones)-1, sum([timepoint_counts[x] for x in timepoint_counts]))
start_row += len(clones) + 5
offset = 7 if offset == 0 else 0

chart_worksheet.write_row(start_row,0,["before rebound timepoints", '# sequences'])
start_row += 1
pre_clones.append('singles')
pre_counts.append(pre_singles)
pre_colors.append({'fill': {'color': 'white'}, 'border':{'color':'#606060'}})
chart_worksheet.write_column(start_row, 0, pre_clones)
chart_worksheet.write_column(start_row, 1, pre_counts)

make_chart(workbook, chart_worksheet, start_row, pre_colors, 'before rebound timepoints', offset, start_row+len(pre_clones)-1, sum([timepoint_counts[x] for x in timepoint_counts if x not in REBOUND]))


workbook.close()

# output fasta containing unique sequences
with open(outf_unique_fasta,"w") as outf:
    SeqIO.write(unique_seqs,outf,"fasta")

