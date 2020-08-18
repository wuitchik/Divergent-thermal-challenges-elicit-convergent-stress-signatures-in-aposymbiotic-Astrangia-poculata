## Astrangia MPCC 2017 transcriptome assembly ##

# working directory

/projectnb/coral/MPCC_2017/transcriptome_assembly

# cat all fastqs

cat ../*.trim > astrangia.trim.fastq


# running Trinity:
module load bowtie2
module load bowtie2
module load samtools
module load jellyfish
module load salmon
module load python3/3.6.5
module load trinity/2.8.4
module load java

## Trinity file "trin" ##

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trinity # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be
#$ -pe omp 24

Trinity --seqType fq --max_memory 750G --single astrangia.trim.fastq --CPU 24

# Qsub it, and the resulting file is Trinity.fasta

# getting uniprot_swissprot KB database
echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" >getuni2
nano getuni2
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N getuni # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be
qsub getuni2

# getting annotations (this file is over 3G, will take a while)
echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz" >getgo
nano getgo
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N getgo # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be
qsub getgo

# indexing the fasta database
module load blast+
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
nano mdb
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N mbd # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M wuitchik@bu.edu #your email
#$ -m be
qsub mbd

# splitting the transcriptome into 100 chunks
splitFasta.pl ast_MPCC2017_transcriptome.fasta 100

# blasting all 100 chunks to uniprot in parallel, 3 cores per chunk
module load blast+/
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 3 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
module load python3
/projectnb/coral/MPCC_2018/scc6_qsub_launcher.py -N blast -P coral -M wuitchik -j y -h_rt 24:00:00 -jobsfile bl
qsub blast_array.qsub

# combining all blast results
cat subset*br > ast_MPCC2017_transcriptome.br
rm -f subset*

# for trinity-assembled transcriptomes: annotating with "isogroup" (=component)
grep ">" ast_MPCC2017_transcriptome.fasta | perl -pe 's/>((\S+)_i\d+).+/$1\t$2/' >ast_MPCC2017_transcriptome_seq2iso.tab 
cat ast_MPCC2017_transcriptome.fasta | perl -pe 's/>((\S+)_i\d+)/>$1 gene=$2/' >ast_MPCC2017_transcriptome_iso.fasta

# extracting gene names (per isogroup) these are Misha perl scritps that can be found: https://github.com/z0on/annotatingTranscriptomes/blob/master/getGeneNameFromUniProtKB.pl
echo "getGeneNameFromUniProtKB.pl blast=ast_MPCC2017_transcriptome.br prefix=astrangia fastaQuery=ast_MPCC2017_transcriptome_iso.fasta" >getgn
module load python3
/projectnb/coral/MPCC_2018/scc6_qsub_launcher.py -N getgn -P coral -M wuitchik -j y -h_rt 24:00:00 -jobsfile getgn
qsub getgn_array.qsub

# extracting GO annotations (per isogroup): https://github.com/z0on/annotatingTranscriptomes/blob/master/getGOfromUniProtKB.pl
echo "getGOfromUniProtKB.pl blast=ast_MPCC2017_transcriptome.br prefix=astrangia fastaQuery=ast_MPCC2017_transcriptome_iso.fasta" >getgo
module load python3
/projectnb/coral/MPCC_2018/scc6_qsub_launcher.py -N getgo -P coral -M wuitchik -j y -h_rt 24:00:00 -jobsfile getgo
qsub getgo_array.qsub

