# python script to perform BLAST on target sequences of Xenium and convert those results to a BED file
# for graphical representation in genome browser
#

import sys, os
import numpy as np
from bisect import bisect
from operator import itemgetter

fasta = sys.argv[1]
db = sys.argv[2]
gtf = sys.argv[3]

outbed = fasta.replace(".fa", "_IGV.bed")
fout = open(outbed, "w")

## function to read fasta file and storing each sequence into dictionary
def read_fasta(fastafile):
    from Bio import SeqIO
    fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    return fasta_dict
fasta_sequence_db = read_fasta(db)

# parse the GTF file to store the name of the transcript and gene name (to later calculate gene length and convert co-ordinates)
print("Creating GTF dictionary")
def parseGTF(gtffile):
    gtf_dict_f = {}       # dictionary with key as transcript ID and value as its chr, start, stop coordinates, strand and gene name
    exon_dict_f = {}      # dictionary with key as transcript ID and value as list of exon sorted based on coordinates
    for lgtf in open(gtffile).readlines():
        if lgtf.startswith("#"):    continue
        egtf = lgtf.strip().split("\t")
        if (egtf[2] == "transcript"):
            temp = egtf[8].split(';')
            gene_name = [s for s in temp if "gene_name" in s][0].replace("\"", "").replace(" ", "").replace("gene_name", "")
            transcript_id = [s for s in temp if "transcript_id" in s][0].replace("\"", "").replace(" ", "").replace("transcript_id", "")
            #print(gene_name, transcript_id)
            gtf_dict_f[transcript_id] = [egtf[0], int(egtf[3]), int(egtf[4]), gene_name, egtf[6]]
            exon_dict_f[transcript_id] = []
        elif (egtf[2] == "exon"):
            temp = egtf[8].split(';')
            transcript_id = [s for s in temp if "transcript_id" in s][0].replace("\"", "").replace(" ", "").replace("transcript_id", "")    
            exon_dict_f[transcript_id].append([int(egtf[3]), int(egtf[4])])
    return([gtf_dict_f, exon_dict_f])

out_function = parseGTF(gtf)
gtf_dict = out_function[0]
exon_dict = out_function[1]

#print(exon_dict["ENST00000641515"])

# now perform the blast step and store the results in an output file
blast_out = fasta + ".out"
s = os.system('module load blast/2.3.0+')
blast_command = "blastn -db " + db + " -query " + fasta + " -out " + blast_out + " -task blastn -outfmt 6" 
print(os.system(blast_command))

# now parse the BLAST output file and store the results into BED format
fin = open(blast_out).readlines()
for line in fin:
    each = line.strip().split("\t")
    gene = each[0].split("_")[0]
    transcript = each[1]
    if (transcript in gtf_dict.keys() and gtf_dict[transcript][3] == gene):
        chr = gtf_dict[transcript][0]
        strand = gtf_dict[transcript][4]
        score = each[11]
        start, end = int(each[8]), int(each[9])
        match_length = int(each[3])
        
        # sort and remove duplicates for the transcript exon list
        exon_transcript_sorted = sorted(exon_dict[transcript], key=itemgetter(0))
        transcript_length_array = [i[1]-i[0]+1 for i in exon_transcript_sorted]
        exon_index_transcript = np.cumsum(transcript_length_array).tolist()
        
        # if the gene is on the negative strand, swap the start and end coordinates by subtracting these from the transcript length
        if (strand == "-"):
            start = sum(transcript_length_array) - start
            end = sum(transcript_length_array) - end

        exon_index = bisect(exon_index_transcript, start)
        box1_length = exon_index_transcript[exon_index] - start
        box2_length = match_length - box1_length

        # convert to genome coordinates
        genome_start = exon_transcript_sorted[exon_index][1] - box1_length
        genome_end = exon_transcript_sorted[exon_index+1][0] + box2_length
        distance_exons = exon_transcript_sorted[exon_index+1][0] - exon_transcript_sorted[exon_index][1]
        start_box2 = box1_length + distance_exons - 1

        #print("###############")
        #print(gene, transcript, gtf_dict[transcript], exon_transcript_sorted, exon_transcript_sorted[0][0], gtf_dict[transcript][1], transcript_length_array, sum(transcript_length_array), len(fasta_sequence_db[transcript]))
        #print(exon_index_transcript, start, end, exon_index, box1_length, bisect(exon_index_transcript, end), box2_length)
        #print(genome_start, genome_end)

        bed_entry = [chr, str(genome_start), str(genome_end), each[0], "0", strand, str(genome_start), str(genome_end), "166,206,227", "2", str(box1_length)+","+str(box2_length), "0,"+str(start_box2)]
        print(bed_entry)
        fout.write("\t".join(bed_entry))
        fout.write("\n")
