
#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=16G
#$ -l h_rt=24:0:0
#$ -cwd
#$ -j y
#$ -t 1-6

export OMP_NUM_THREADS=$NSLOTS

module load fastqc
module load python
source /data/home/hfx381/presto/bin/activate

### Need to edit the -t flag at the top based on number of samples to process

### EDITABLE VARIABLE - directory between "ngs_data" and ${SGE_TASK_ID} - CHOOSE WHICH DATABATCH TO PROCESS

data="/data/BCI-ESCS/richard/ngs_data/***MASTER-VARIABLE***/${SGE_TASK_ID}/"

#PATH to genome folder for bismark alignment

genome_folder="/data/BCI-ESCS/richard/processing_files/genome/clock_genome/"

#PATH to asm_primer sequences fasta file

primer_fasta="/data/BCI-ESCS/richard/processing_files/fasta_files/asm_primers_f_r.fasta"

#PATH to the python functions

pythonfiles="/data/home/hfx381/presto/bin/"

### EDITABLE VARIABLE - ASSIGN THE APPROPRIATE DATABATCH

databatch="***MASTER-VARIABLE***"

#Make new directory to store the CpG Text file outputs:

mkdir /data/BCI-ESCS/richard/ngs_data/output_meth_reports/methseq_$databatch

#PATH to deliver the CpG Text File outputs

output_path="/data/BCI-ESCS/richard/ngs_data/output_meth_reports/methseq_$databatch/"


#Processing Code Begins Here:

cd $data

sample_name=$(for folder in *; do echo $folder; done)

fastqc $sample_name/*.fastq.gz

gunzip $sample_name/*.fastq.gz

module unload fastqc

$pythonfiles/FilterSeq.py quality -s $sample_name/*R1_001.fastq -q 20 --outdir $sample_name --outname "$sample_name"_r1 --log $sample_name/"$sample_name"_r1_qual.log

$pythonfiles/FilterSeq.py quality -s $sample_name/*R2_001.fastq -q 20 --outdir $sample_name --outname "$sample_name"_r2 --log $sample_name/"$sample_name"_r2_qual.log

$pythonfiles/ParseLog.py -l $sample_name/*r1_qual.log $sample_name/*r2_qual.log -f ID QUALITY

rm $sample_name/*r1_qual.log $sample_name/*r2_qual.log

$pythonfiles/FilterSeq.py trimqual -s $sample_name/*r1_quality-pass.fastq -q 20 --outdir $sample_name --outname "$sample_name"_r1 --log $sample_name/"$sample_name"_r1_trimqual.log

$pythonfiles/FilterSeq.py trimqual -s $sample_name/*r2_quality-pass.fastq -q 20 --outdir $sample_name --outname "$sample_name"_r2 --log $sample_name/"$sample_name"_r2_trimqual.log

$pythonfiles/ParseLog.py -l $sample_name/*r1_trimqual.log $sample_name/*r2_trimqual.log -f ID QUALITY

rm $sample_name/*r1_trimqual.log $sample_name/*r2_trimqual.log

$pythonfiles/FilterSeq.py length -s $sample_name/*r1_trimqual-pass.fastq -n 100 --outdir $sample_name --outname "$sample_name"_r1 --log $sample_name/"$sample_name"-r1_length.log

$pythonfiles/FilterSeq.py length -s $sample_name/*r2_trimqual-pass.fastq -n 100 --outdir $sample_name --outname "$sample_name"_r2 --log $sample_name/"$sample_name"-r2_length.log

$pythonfiles/ParseLog.py -l $sample_name/*r1_length.log $sample_name/*r2_length.log -f ID LENGTH

rm $sample_name/*r1_length.log $sample_name/*r2_length.log

$pythonfiles/PairSeq.py -1 $sample_name/*r1_length-pass.fastq -2 $sample_name/*r2_length-pass.fastq --outdir $sample_name --coord illumina

module load seqtk

seqtk trimfq -b 2 $sample_name/*r2_length-pass_pair-pass.fastq > $sample_name/"$sample_name"_r2_barcode_cc_trimmed.fastq

$pythonfiles/MaskPrimers.py extract -s $sample_name/*r2_barcode_cc_trimmed.fastq --mode cut --len 21 --pf BARCODE --outdir $sample_name --outname "$sample_name"_r2 --log $sample_name/"$sample_name"_r2_primers-pass.log

$pythonfiles/ParseLog.py -l $sample_name/*r2_primers-pass.log --outdir $sample_name -f ID PRIMER BARCODE ERROR

rm $sample_name/*r2_primers-pass.log

$pythonfiles/PairSeq.py -1 $sample_name/*r1_length-pass_pair-pass.fastq -2 $sample_name/*r2_primers-pass.fastq --outdir $sample_name --2f BARCODE --coord illumina

seqtk trimfq -b 20 $sample_name/*r1_length-pass_pair-pass_pair-pass.fastq > $sample_name/"$sample_name"_r1_primers_trimmed.fastq

seqtk trimfq -b 20 $sample_name/*r2_primers-pass_pair-pass.fastq > $sample_name/"$sample_name"_r2_primers_trimmed.fastq

$pythonfiles/BuildConsensus.py -s $sample_name/*r1_primers_trimmed.fastq --outdir $sample_name --bf BARCODE --maxerror 0.1 --outname "$sample_name"_r1 --log $sample_name/"$sample_name"_r1_consensus.log

$pythonfiles/BuildConsensus.py -s $sample_name/*r2_primers_trimmed.fastq --outdir $sample_name --bf BARCODE --maxerror 0.1 --outname "$sample_name"_r2 --log $sample_name/"$sample_name"_r2_consensus.log

$pythonfiles/ParseLog.py -l $sample_name/*r1_consensus.log $sample_name/*r2_consensus.log --outdir $sample_name -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT

rm $sample_name/*r1_consensus.log $sample_name/*r2_consensus.log

$pythonfiles/PairSeq.py -1 $sample_name/*r1_consensus-pass.fastq -2 $sample_name/*r2_consensus-pass.fastq --outdir $sample_name --coord presto

$pythonfiles/SplitSeq.py group -s $sample_name/*r1_consensus-pass_pair-pass.fastq --outdir $sample_name -f CONSCOUNT --num 2 --outname "$sample_name"_r1

$pythonfiles/ParseHeaders.py table -s $sample_name/*r1_atleast-2.fastq --outdir $sample_name -f ID PRCONS CONSCOUNT DUPCOUNT

$pythonfiles/SplitSeq.py group -s $sample_name/*r2_consensus-pass_pair-pass.fastq --outdir $sample_name -f CONSCOUNT --num 2 --outname "$sample_name"_r2

$pythonfiles/ParseHeaders.py table -s $sample_name/*r2_atleast-2.fastq --outdir $sample_name -f ID PRCONS CONSCOUNT DUPCOUNT

$pythonfiles/PairSeq.py -1 $sample_name/*r1_atleast-2.fastq -2 $sample_name/*r2_atleast-2.fastq --outdir $sample_name --coord presto

deactivate

module unload python
module load bismark
module load bowtie2
module load samtools

bismark --bowtie2 --minins 300 --maxins 1000 --score_min L,0,-1 --genome $genome_folder -1 $sample_name/*r1_atleast-2_pair-pass.fastq -2 $sample_name/*r2_atleast-2_pair-pass.fastq -o $sample_name

bismark_methylation_extractor -p --gzip --bedgraph --buffer_size 10G --cytosine_report --genome $genome_folder $sample_name/*bismark_bt2_pe.bam -o $sample_name

samtools sort $sample_name/*bismark_bt2_pe.bam > $sample_name/"$sample_name"_bismark_bt2_pe_sorted.bam

samtools index $sample_name/*sorted.bam

gunzip $sample_name/CpG_OT*.txt.gz

mv $sample_name/CpG_OT*.txt $sample_name/"$sample_name"_meth_calls.txt

cp $sample_name/*meth_calls.txt $output_path

rm $sample_name/*primers-pass.fastq $sample_name/*primers-pass_pair-pass.fastq $sample_name/*length-pass.fastq $sample_name/*quality-pass.fastq $sample_name/*trimqual-pass.fastq $sample_name/*trimqual-pass_pair-pass.fastq $sample_name/*primers_trimmed.fastq $sample_name/*cc_trimmed.fastq $sample_name/*trimqual-pass_pair-pass_pair-pass.fastq $sample_name/*consensus-pass.fastq $sample_name/*consensus-pass_pair-pass.fastq $sample_name/*atleast-2.fastq $sample_name/*under-2.fastq $sample_name/*trimmed.fq

gzip $sample_name/*.fastq

exit