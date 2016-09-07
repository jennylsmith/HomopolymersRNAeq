# HomopolymersRNAeq

#PS3 RNA-Seq Data Exploration
August 25, 2016


1.Files used:
	species: Gulf pipefish
	/research/bi610/data/rnaseq/PE_RNAseq_R1.fastq.gz
    /research/bi610/data/rnaseq/PE_RNAseq_R2.fastq.gz

2.  using Unix commands, determine the number of pairs of reads in the files 

determine the total number of read pairs by counting the number of header lines in the fastq
for file in *.fastq; do grep -c -H "^@HWI" $file; done
PE_RNAseq_R1.fastq:10000000
PE_RNAseq_R2.fastq:10000000

same number of reads because they are paired-end. So a total of 10 million pairs of reads

find the number of reads that failed Illumina's chastity filter.
Y = yes the read is bad and N = no, the read is not bad. 

cat PE_RNAseq_R1.fastq | grep "^@HWI" | grep -c "1:Y:" 
809609
cat PE_RNAseq_R2.fastq | grep "^@HWI" | grep -c "2:Y:" 
809609

find the proportion of reads that failed the chastity filter
809609/10000000 = 0.081 or 8% of reads 


3. Using Unix commands, extract and print a list of barcodes in R1 file. Count the
abundances of the barcodes in order from least to most abundant. Turn in the top 20 entries. 

unix command to pullout the barcodes from the header lines of fastq
 cat PE_RNAseq_R1.fastq | grep "^@HWI" | \
  sed -E 's/^@HWI.{28}:[0-9]+:[0-9]+ 1:[A-Z]{1}:0:([A-Z]+)/\1/' | sort \
  | uniq -c | sort -r -nk 1 > barcodes.txt

head -20 barcodes.txt 
4018067 CGATGT
1071176 TGACCT
 980546 CAGATT
 888640 GTCCGT
 862571 AGTTCT
 759964 CCGTCT
  64436 CCGCCT #notice that this barcode is very similar to the one above, indicating a seq error 
  48509 CAGCTT 
  42142 CCGACT
  39385 AGCTCT
  38547 AGATCT
  35690 CAGACT
  34640 TAACCT
  33312 CGCTGT
  32614 GCCCGT
  30782 TTACCT
  29319 GTCCCT
  27504 CGATCT
  27113 CGACCT
  27035 TCACCT


a. What type of barcodes were used? and what length are they?

These barcodes are index barcodes because they are found in the header of the fastq file. 
The length of the barcodes is 6 nucleotides. 


b. How many multiplexed libraries were likely submitted to sequencing?

There were most likely 6 libraries submitted to sequencing. The other barcodes are likely 
due to sequencing errors when the barcode was being read. This is because the top 6 most 
abundant barcodes are found 700,000 times or more, while the other barcodes are found much 
less often, and have similar sequences to the most abundant barcodes. This indicates that
there could have been a mismatch during the sequencing and caused single nucleotide 
difference in the index sequences. There were a total of 2,559 "indexes" that were pulled
out of the index header in the fastq files. 

wc of the barcodes. Has barcodes with N? 
wc -l barcodes.txt 
2559 barcodes.txt

c. Did the person that prepped these libraries get even coverage across them? 

No, there was not even coverage, because many barcodes are more abundant in the 
dataset than others. For example, the most abundant barcode is present in 4 million
reads, but the next most abundant is present in only 1 million reads. 


4. finding possible adaptor contamination in the sequence reads. 

Forward adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
Reverse adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

5. use grep to determine if there is forward or reverse adapter contamination. 

for file in *.fastq; do grep -c -H "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" $file ; done
PE_RNAseq_R1.fastq:9621
PE_RNAseq_R2.fastq:0

For file in *.fastq; do grep -c -H "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" $file; done 
PE_RNAseq_R1.fastq:0
PE_RNAseq_R2.fastq:9399

There are approximately 9,000 sequences that have adaptor contamination in this
dataset. 

a. execute grep using different subsets of the adaptor sequence. Is there variation
in the number of matches? 

make a variable to hold the forward and reverse adaptor sequences
forward='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
forwardHalf=$(echo $forward | cut -c 1-20) #first 20 nucleotides of F adaptor seq
forwardHalf2=$(echo $forward | cut -c 1-10)  #first 10 nucleotides of F adaptor seq

reverse='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
reverseHalf=$(echo $reverse | cut -c 1-20) #first 20 nucleotides of R adaptor seq
reverseHalf2=$(echo $reverse | cut -c 1-10) #first 10 nucleotides of R adaptor seq

For loop to grep for the partial adaptor sequences. 
ls -1 | for file in *.fastq ; do forwardContam1=$(grep $forwardHalf $file | wc -l); \
forwardContam2=$(grep $forwardHalf2 $file | wc -l) ; \ 
reverseContam1=$(grep $reverseHalf $file | wc -l); \
reverseContam2=$(grep $reverseHalf2 $file | wc -l); \
echo "the forward contamination is" $forwardContam1 "," $forwardContam2 "." \
"The reverse contamination is" $reverseContam1 "," $reverseContam2; done 

the forward contamination is 24572 , 66834 . The reverse contamination is 0 , 66834
the forward contamination is 1 , 70637 . The reverse contamination is 27238 , 70637

This result indicates that in PE_RNASeq_R1.fastq (which would be listed first due to the
ls -1 command), has 24,572 reads which may contain the first 20 nucleotides of the
forward adaptor sequence and an additional 66,834 reads may contain 10 nucleotides of the
forward adaptor sequence. There are 0 R1 reads with reverse adaptor contamination,
which  should not be present in the first read, which is also confirmed above when grepping
for the entire reverse sequence. However, 66,834 R1 reads with the first 10 nucleotides of
the reverse adaptor, which indicates the difficulty of resolving issues when there is only
a small amount of adaptor contamination in the reads. 

The reverse adaptor contamination shows 27,238 R2 reads which may contain partial adaptor 
sequences. But again, 70,637 sequences in R2 had potential contamination of the 10 
nucleotides of the reverse adaptor, and highlights difficulty of partial contamination. 
Also the the first 10 nucleotides of the forward and reverse adaptors are very similar,
as seen below. 

echo $forwardHalf2 $reverseHalf2
AGATCGGAAG AGATCGGAAG


b. what could happen to the assembly process if adaptor sequences are present in
the raw data? 

The assembly process could result in in accuracies when building contigs, because
the adaptor sequences would overlap with each other, forming an edge, and the sequences
from two different mRNA transcripts could be made into a chimeric contig. 


6. Use Gmap-Gsnap to align the R1 reads to the reference ribosomal RNA sequences.  

a. File used:
	/research/bi610/data/rnaseq/Gacu_rRNA.fasta

b. Command line used to build the ribosomal RNA database
module load gmap-gsnap 
gmap_build -d pipefish -D /home12/jsmith16/bi623/160825_PS3 -k 15 Gacu_rRNA.fasta &


c. Job submission for gsnap assembly:

cd /home12/jsmith16/bi623/160825_PS3/

module load gmap-gsnap/2013-11-23

gsnap --nthreads 12 -m 20 -d pipefish -D /home12/jsmith16/bi623/160825_PS3/pipefish \
-B 4 -O --split-output gsnap_out -A sam \
 /home12/jsmith16/bi623/160825_PS3/PE_RNAseq_R1.fastq 

The --split-output option is to split the output types into separate SAM files. The output
types are A) Nomapping, which are the reads that failed to align to the database, 
B) circular, which finds alignments that go around the origin of a circular chromosome, such 
as plasmids, C) unpaired_mult are reads with multiple alignments to the reference, 
D) the translocation file contains reads where part is mapped to a distant location in the
genome, and E) unpaired_uniq are reads that map to single location. 

d. Count the aligned reads by parsing the SAM output files, but remember not to count 
reads more than once if they align non-uniquely. Consider a read of rRNA origin if it 
aligns to the rRNA reference in any capacity (e.g. multiply, uniquely, translocally, etc.).

unix commands to remove the header from the sam files and then sort by the read_id field
to count each read only once. 

cat gsnap_out.nomapping | grep -v "^@" | cut -f 1 -d " " | sort | uniq | sort | wc -l
9932206

cat gsnap_out.unpaired_mult | grep -v "^@" | cut -f 1 -d "     " | sort | uniq | sort | wc -l 
51631

cat gsnap_out.unpaired_uniq | grep -v "^@" | cut -f 1 -d "   " | sort | uniq | wc -l  
16163

e. The nomapping file indicates all the unique RNA-seq reads that were not originating from
ribosomal RNA, since they did not align to the rRNA database. There were 9,932,206 reads
that did not align to the rRNA reference. Thus, these are the reads which should be
included in downstream analyses, since they most likely originated from mRNA transcripts 
during library preparation. 

There were a total of 51,631 reads which originated from rRNA contamination during library
preparation and had multiple alignments. Also, 16,163 reads were also rRNA and they aligned
uniquely to the reference database. There were no alignments found in circular 
and translocation output files. This indicates that the mRNA isolation was very efficient
for this library preparation, since the identified rRNA contamination was in a total of
67,794 reads out of 10,000,000. This is a total 0.67% of read which originated from rRNA.


7. create a shell script or python script to determine proportions of reads that have 
a polyA tail. PolyA tails have more than 15 consecutive As. 

shell script to grep out the sequence reads from the fastq file and then grep and count
the number of reads with 15 or more As. 

cat PE_RNAseq_R1.fastq | grep -A 1 "@HWI" | grep -v "\-" | grep -c  -E "[A]{15,}"  
39871

There are a total of 39,871 reads containing a polyA tail, which is approximately
0.4% of all reads. Reads with polyA tails can be useful for annotation of genes, since
it provides information about possible splice variants. 

8. Identify reads with 15 or more consecutive Cs, Gs and Ts. Explain why their observed
frequencies do or do not make sense. 

shell script to identify sequence reads with homopolymers composed of Cs, Ts, and Gs. 

cat PE_RNAseq_R1.fastq | grep -A 1 "@HWI" | grep -v "\-" | grep -c  -E "[C]{15,}"  
2063

cat PE_RNAseq_R1.fastq | grep -A 1 "@HWI" | grep -v "\-" | grep -c  -E "[T]{15,}"  
32851

cat PE_RNAseq_R1.fastq | grep -A 1 "@HWI" | grep -v "\-" | grep -c  -E "[G]{15,}"  
2523


Using a python script for this same application of PolyA tails and homopolymer identification

python3.5 PS3_homopolymers.py 
The number of polyA tails is 39871
The number of polyC is 2063
The number of polyT is 32851
The number of polyG is 2523

The observed frequencies of polyA sequences (39,871 reads) and polyT sequences (32,851)
do make sense. Library preparation used oligodT to enrich for mRNA transcripts, which 
bind only to the polyA tails. All mRNA sequences have a polyA tail and this was
reverse transcribed into cDNA, resulting in polyT regions in the cDNA. Thus, during 
sequencing, there would be similar numbers of reads with polyA and poly T regions. Also, 
it makes sense for homopolymers of 15 or more Cs and Gs (2,063 and 2,523 sequences respectively)
to be relatively rare in mRNA sequences. 

