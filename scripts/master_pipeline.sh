##################################################
## Chapter 3; DNA Divers; Bioiformatic pipeline ##
##################################################

# raw sequences (fastq files) to taxonomic assignment and read counts

#####
## 1. Quality control, trim, merge, quality control merged reads
#####

### 1.1 Quality control and trim

# generate qc plots

# sequence run 1
mkdir fasqc
fastqc -o /home/genome2/beseneav/divers_methods1/fastqc/ --extract -f fastq *.fastq

# sequence run 2
mkdir fasqc
gzip -d *fastq.gz
fastqc -o /home/genome2/beseneav/divers_methods2/fastq/ --extract -f fastq *.fastq

# sequence run 3 
mkdir fasqc
gzip -d *fastq.gz
fastqc -o /home/genome2/beseneav/elas_other/fasqc/ --extract -f fastq *.fastq

# activate OBITools environment
source /home/genome2/beseneav/applications/OBITools/OBITools-env/bin/activate

### 1.2 Trim sequences based on quality

# sequence run 1
obicut -e 150 lib-diversmethods1_S1_L001_R1_001.fastq > divmeth1_R1_trim150.fastq
obicut -e 150 lib-diversmethods1_S1_L001_R2_001.fastq > divmeth1_R2_trim150.fastq

# sequence run 2
obicut -e 148 elas02-diversmethods2_S2_L001_R1_001.fastq > e_divmeth2_R1_trim148.fastq
obicut -e 148 elas02-diversmethods2_S2_L001_R2_001.fastq > e_divmeth2_R2_trim148.fastq
obicut -e 148 tele02-diversmethods2_S1_L001_R1_001.fastq > t_divmeth2_R1_trim148.fastq
obicut -e 148 tele02-diversmethods2_S1_L001_R2_001.fastq > t_divmeth2_R2_trim148.fastq

# sequence run 3
obicut -e 150 other-elas02_S2_L001_R1_001.fastq > e_other_R1_trim150.fastq
obicut -e 150 other-elas02_S2_L001_R2_001.fastq > e_other_R2_trim150.fastq

### 1.3 Merge paired-end reads

# sequence run 1
illuminapairedend -r divmeth1_R1_trim150.fastq divmeth1_R2_trim150.fastq > paired_divmeth1.fastq

# sequence run 2
illuminapairedend -r e_divmeth2_R1_trim148.fastq e_divmeth2_R2_trim148.fastq > paired_e_divmeth2.fastq
illuminapairedend -r t_divmeth2_R1_trim148.fastq t_divmeth2_R2_trim148.fastq > paired_t_divmeth2.fastq

# sequence run 3
illuminapairedend -r e_other_R1_trim150.fastq e_other_R2_trim150.fastq > paired_e_other.fastq

### 1.4 Quality control the aligned reads

# sequence run 1
obiannotate -S goodali:'"Good_divmeth1" if score>30.00 else "Bad_divmeth1"' paired_divmeth1.fastq| obisplit -t goodali

export PATH=$PATH:~/applications/seqtk/
seqtk seq -a Good_divmeth1.fastq > Good_divmeth1.fasta

# sequence run 2
obiannotate -S goodali:'"Good_e_divmeth2" if score>30.00 else "Bad_e_divmeth2"' paired_e_divmeth2.fastq| obisplit -t goodali
obiannotate -S goodali:'"Good_t_divmeth2" if score>30.00 else "Bad_t_divmeth2"' paired_t_divmeth2.fastq| obisplit -t goodali

export PATH=$PATH:~/applications/seqtk/
seqtk seq -a Good_e_divmeth2.fastq > Good_e_divmeth2.fasta
seqtk seq -a Good_t_divmeth2.fastq > Good_t_divmeth2.fasta

# sequence run 3
obiannotate -S goodali:'"Good_e_other" if score>30.00 else "Bad_e_other"' paired_e_other.fastq| obisplit -t goodali

export PATH=$PATH:~/applications/seqtk/
seqtk seq -a Good_e_other.fastq > Good_e_other.fasta


#####
## 2. Demultiplex, filter aligned reads by length, dereplicate reads
#####

### 2.1 Demultiplex sample barcodes

# sequence run 1
ngsfilter -t /home/genome2/beseneav/divers_methods1/demulti/ngs_filter_divmeth1_v2.txt --fasta-output -u unidentified_divmeth.fasta /home/genome2/beseneav/divers_methods1/Good_divmeth1.fasta --DEBUG > divmeth1.filtered.fasta

# sequence run 2
ngsfilter -t /home/genome2/beseneav/divers_methods2/demulti/ngsfilter_e_dnadivers2.txt --fasta-output -u unidentified_e_divmeth2.fasta /home/genome2/beseneav/divers_methods2/Good_e_divmeth2.fasta --DEBUG > e_divmeth2.filtered.fasta
ngsfilter -t /home/genome2/beseneav/divers_methods2/demulti/ngsfilter_t_dnadivers2.txt --fasta-output -u unidentified_t_divmeth2.fasta /home/genome2/beseneav/divers_methods2/Good_t_divmeth2.fasta --DEBUG > t_divmeth2.filtered.fasta

# sequence run 3
ngsfilter -t /home/genome2/beseneav/elas_other/demulti/ngsfilter_elas_other.txt --fasta-output -u unidentified_e_other.fasta /home/genome2/beseneav/elas_other/Good_e_other.fasta --DEBUG > e_other.filtered.fasta

### 2.2 Filter reads by length & calculate per sample statistics

# sequence run 1
obigrep -p 'seq_length>130' -p 'seq_length<190' -s '^[ACGT]+$' divmeth1.filtered.fasta > divmeth1.filtered_length_noN.fasta

obistat -c sample -a seq_length divmeth1.filtered_length_noN.fasta > sample_stats_divmeth1.filtered_length.txt

# sequence run 2
obigrep -p 'seq_length>130' -p 'seq_length<210' -s '^[ACGT]+$' e_divmeth2.filtered.fasta > e_divmeth2.filtered_length_noN.fasta
obigrep -p 'seq_length>130' -p 'seq_length<190' -s '^[ACGT]+$' t_divmeth2.filtered.fasta > t_divmeth2.filtered_length_noN.fasta

obistat -c sample -a seq_length e_divmeth2.filtered_length_noN.fasta > sample_stats_e_divmeth2.filtered_length.txt
obistat -c sample -a seq_length t_divmeth2.filtered_length_noN.fasta > sample_stats_t_divmeth2.filtered_length.txt

# sequence run 3
obigrep -p 'seq_length>130' -p 'seq_length<210' -s '^[ACGT]+$' e_other.filtered.fasta > e_other.filtered_length_noN.fasta

obistat -c sample -a seq_length e_other.filtered_length_noN.fasta > sample_stats_e_other.filtered_length.txt

### 2.3 Dereplicate reads and add id

# sequence run 1
obiuniq -m sample divmeth1.filtered_length_noN.fasta > divemeth1.unique.fasta

obiannotate --seq-rank divemeth1.unique.fasta | obiannotate --set-identifier '"'dme1'_%0*d" %(9,seq_rank)' > divemeth1.new9.fasta

# sequence run 2
obiuniq -m sample e_divmeth2.filtered_length_noN.fasta > e_divmeth2.unique.fasta
obiuniq -m sample t_divmeth2.filtered_length_noN.fasta > t_divmeth2.unique.fasta

obiannotate --seq-rank e_divmeth2.unique.fasta | obiannotate --set-identifier '"'dme2'_%0*d" %(9,seq_rank)' > e_divmeth2.new9.fasta
obiannotate --seq-rank t_divmeth2.unique.fasta | obiannotate --set-identifier '"'dmt2'_%0*d" %(9,seq_rank)' > t_divmeth2.new9.fasta

# sequence run 3
obiuniq -m sample e_other.filtered_length_noN.fasta > e_other.unique.fasta

obiannotate --seq-rank e_other.unique.fasta | obiannotate --set-identifier '"'othe'_%0*d" %(9,seq_rank)' > e_other.new9.fasta


#####
## 3. Remove chimeras and cluster amplicons
#####

# R scripts from https://github.com/metabarpark/

### 3.1 Format fasta files for vsearch software

# sequence run 1
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_obifasta2vsearch -i divemeth1.new9.fasta -o divmeth1.vsearch.fasta
sed 's/ ;/;/g' divmeth1.vsearch.fasta > divmeth1.vsearch.mod.fasta

# sequence run 2
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_obifasta2vsearch -i e_divmeth2.new9.fasta -o e_divmeth2.vsearch.fasta
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_obifasta2vsearch -i t_divmeth2.new9.fasta -o t_divmeth2.vsearch.fasta
sed 's/ ;/;/g'  e_divmeth2.vsearch.fasta >  e_divmeth2.vsearch.mod.fasta
sed 's/ ;/;/g'  t_divmeth2.vsearch.fasta >  t_divmeth2.vsearch.mod.fasta

# sequence run 3
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_obifasta2vsearch -i e_other.new9.fasta -o e_other.vsearch.fasta
sed 's/ ;/;/g'  e_other.vsearch.fasta >  e_other.vsearch.mod.fasta

### 3.2 Chimera detection and removal

# sequence run 1
vsearch=~/applications/vsearch/bin/vsearch
$vsearch --uchime_denovo divmeth1.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/divmeth1.nonchimeras.fasta --chimeras vsearch_output/divmeth1.chimeras.fasta --threads 28 --uchimeout vsearch_output/divmeth1.uchimeout.txt &> vsearch_output/log.divmeth1_chimeras

# sequence run 2
vsearch=~/applications/vsearch/bin/vsearch
$vsearch --uchime_denovo e_divmeth2.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/e_divmeth2.nonchimeras.fasta --chimeras vsearch_output/e_divmeth2.chimeras.fasta --threads 28 --uchimeout vsearch_output/e_divmeth2.uchimeout.txt &> vsearch_output/log.e_divmeth2_chimeras
$vsearch --uchime_denovo t_divmeth2.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/t_divmeth2.nonchimeras.fasta --chimeras vsearch_output/t_divmeth2.chimeras.fasta --threads 28 --uchimeout vsearch_output/t_divmeth2.uchimeout.txt &> vsearch_output/log.t_divmeth2_chimeras

# sequence run 3
vsearch=~/applications/vsearch/bin/vsearch
$vsearch --uchime_denovo e_other.vsearch.mod.fasta --sizeout --nonchimeras vsearch_output/e_other.nonchimeras.fasta --chimeras vsearch_output/e_other.chimeras.fasta --threads 28 --uchimeout vsearch_output/e_other.uchimeout.txt &> vsearch_output/log.e_other_chimeras

#### 3.2 Clustering amplicons & recount abundance by samples

# sequence run 1
swarm=~/applications/swarm/src/swarm
$swarm -d 1 -f -z -t 20 -o swarm_output/divmeth1_SWARM1_output -s swarm_output/divmeth1_SWARM1_stats -w swarm_output/divmeth1_SWARM1_seeds.fasta vsearch_output/divmeth1.nonchimeras.fasta

obitab -o divemeth1.new9.fasta > divmeth1.new.tab
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_recount_swarm ./swarm_output/divmeth1_SWARM1_output divmeth1.new.tab
awk -F";" '{if(NR>1){print ">"$1"\n"$NF}}' divmeth1_SWARM1_output.counts.csv > divmeth1_SWARM1_output.counts.fasta

# sequence run 2
swarm=~/applications/swarm/src/swarm
$swarm -d 1 -f -z -t 20 -o swarm_output/e_divmeth2_SWARM1_output -s swarm_output/e_divmeth2_SWARM1_stats -w swarm_output/e_divmeth2_SWARM1_seeds.fasta vsearch_output/e_divmeth2.nonchimeras.fasta
$swarm -d 1 -f -z -t 20 -o swarm_output/t_divmeth2_SWARM1_output -s swarm_output/t_divmeth2_SWARM1_stats -w swarm_output/t_divmeth2_SWARM1_seeds.fasta vsearch_output/t_divmeth2.nonchimeras.fasta

obitab -o e_divmeth2.new9.fasta > e_divmeth2.new.tab
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_recount_swarm ./swarm_output/e_divmeth2_SWARM1_output e_divmeth2.new.tab
awk -F";" '{if(NR>1){print ">"$1"\n"$NF}}' e_divmeth2_SWARM1_output.counts.csv > e_divmeth2_SWARM1_output.counts.fasta
obitab -o t_divmeth2.new9.fasta > t_divmeth2.new.tab
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_recount_swarm ./swarm_output/t_divmeth2_SWARM1_output t_divmeth2.new.tab
awk -F";" '{if(NR>1){print ">"$1"\n"$NF}}' t_divmeth2_SWARM1_output.counts.csv > t_divmeth2_SWARM1_output.counts.fasta

# sequence run 3
swarm=~/applications/swarm/src/swarm
$swarm -d 1 -f -z -t 20 -o swarm_output/e_other_SWARM1_output -s swarm_output/e_other_SWARM1_stats -w swarm_output/e_other_SWARM1_seeds.fasta vsearch_output/e_other.nonchimeras.fasta

obitab -o e_other.new9.fasta > e_other.new.tab
Rscript /home/genome2/beseneav/applications/R_scripts_metabarpark/owi_recount_swarm ./swarm_output/e_other_SWARM1_output e_other.new.tab
awk -F";" '{if(NR>1){print ">"$1"\n"$NF}}' e_other_SWARM1_output.counts.csv > e_other_SWARM1_output.counts.fasta


#####
## 4. Taxonomic assignment
#####

### 4.1 Ecotag

# generate reference databases for ecotag using ecoPCR
# all vertebrates, no human, elas02 primer
ecoPCR -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 -e 3 -l 130 -L 210 -r 7742 -i 9606 GTTGGTHAATCTCGTGCCAGC CATAGTAGGGTATCTAATCCTAGTTTG > elas02_fish_nohuman.ecopcr

obigrep -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 --require-rank=species --require-rank=genus --require-rank=family elas02_fish_nohuman.ecopcr > elas02_fish_nohuman_clean.fasta
obiuniq -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 elas02_fish_nohuman_clean.fasta > elas02_fish_nohuman_clean_uniq.fasta
obiannotate --uniq-id elas02_fish_nohuman_clean_uniq.fasta > db_elas02_fish_nohuman_clean_uniq.fasta

# all vertebrates, no human, tele02 primer
ecoPCR -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 -e 3 -l 130 -L 210 -r 7742 -i 9606 AAACTCGTGCCAGCCACC GGGTATCTAATCCCAGTTTG > tele02_fish_nohuman.ecopcr

obigrep -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 --require-rank=species --require-rank=genus --require-rank=family tele02_fish_nohuman.ecopcr > tele02_fish_nohuman_clean.fasta
obiuniq -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 tele02_fish_nohuman_clean.fasta > tele02_fish_nohuman_clean_uniq.fasta
obiannotate --uniq-id tele02_fish_nohuman_clean_uniq.fasta > db_tele02_fish_nohuman_clean_uniq.fasta



# sequence run 1

# sequence run 2
ecotag -d /home/genome/marianilab/db_obitools/EMBL_Oct2021/EMBL_r143 -R ~/applications/db_obitools/EMBL_human/db_elas02_fish_nohuman_clean_uniq.fasta --sort=count e_divmeth2_SWARM1_output.counts.fasta > e_divemeth2_SWARM1_all_try3.ecotag.fasta


# sequence run 3



### 4.2 Sintax


# generate reference databases for sintax using meta-fish-lib
# https://github.com/genner-lab/meta-fish-lib


# sequence run 1

# sequence run 2
$vsearch --threads "$THREADS" --sintax /home/genome2/beseneav/divers_methods2/swarm_output/e_divmeth2_SWARM1_output.counts.fasta --db /home/genome2/beseneav/divers_methods2/blast_output/db/references.12s.miya.divers2.ALL.fasta --sintax_cutoff 0.7 --tabbedout /home/genome2/beseneav/divers_methods2/taxa/e_SWARM1_sintax_output_ALL.tsv


# sequence run 3



### 4.3 blastn

# sequence run 1

# sequence run 2
makeblastdb -in ./db/references.12s.miya.divers2.ALL.fasta -dbtype nucl -blastdb_version 5
blastn -task blastn -num_threads 1 -evalue 0.05 -word_size 7 -max_target_seqs 10 -db db/references.12s.miya.divers2.ALL.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out result/elas02-divers2-blast-v3.out -query db/e_divmeth2_SWARM1_output.counts.fasta
echo -e 'id\tblastDbid\tblastEvalue\tblastLength\tblastPident\tblastNident\tblastScore\tblastBitscore' | sed -e "s/-e //g" > result/headers
cat result/headers result/elas02-divers2-blast-v3.out > result/elas02_divers2_blast_v3.tsv


# sequence run 3




