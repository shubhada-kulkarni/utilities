# python script to perform BLAST on target sequences of Xenium and convert those results to a BED file
# for graphical representation in genome browser
#

import sys, os

fasta = sys.argv[1]
db = sys.argv[2]
gtf = sys.argv[3]

# parse the GTF file to store the name of the transcript and gene name (to later calculate gene length and convert co-ordinates)
print("Creating GTF dictionary")
gtf_dict = {}
for lgtf in open(gtf).readlines():
    if lgtf.startswith("#"):    continue
    egtf = lgtf.strip().split("\t")
    if (egtf[2] == "transcript"):
        temp = egtf[8].split(';')
        gene_name = [s for s in temp if "gene_name" in s][0].replace("\"", "").replace(" ", "").replace("gene_name", "")
        transcript_id = [s for s in temp if "transcript_id" in s][0].replace("\"", "").replace(" ", "").replace("transcript_id", "")
        #print(gene_name, transcript_id)
        gtf_dict[transcript_id] = [egtf[0], egtf[3], egtf[4], gene_name, egtf[6]]
#print(gtf_dict)

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
        print(gene, transcript, gtf_dict[transcript])
        start, end = each[8], each[9]
        start = int(gtf_dict[transcript][1]) + int(each[8])
        chr = gtf_dict[transcript][0]
        strand = gtf_dict[transcript][4]
        score = each[11]


        print(chr,  each[8]+","+each[9])




