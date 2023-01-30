#!/bin/sh

#Set the directories for output, input, and various programs
PYTHON="/usr/bin/env python2"
TRANSIT=/home/will/baderlab/Projects/transit-mod
BWA=/usr/bin/bwa
REFS=/home/will/baderlab/Projects/Tnseq_Mtb/my_Tnseq/Results/PBS_Mtb_Tnseq/Analysis/GCF_000195955.2_ASM19595v2_genomic.fna
FASTQ_DIR=/run/media/will/phd_data1/2017-11-26_Mtb_TnSeq/raw
OUT_DIR=/run/media/will/phd_data1/2017-11-26_Mtb_TnSeq/2018-09-15_TPP
PREFIXES_OUTFILE=$OUT_DIR/`basename $FASTQ_DIR`_prefixes.txt

#For Genomes with multiple replicons you will need to set REPLICON_ID
REPLICON_ID="NC_000962.3"

# Settings for creating a CSV file from wig files
GENBANK_FILE=/home/will/baderlab/Projects/Tnseq_Mtb/my_Tnseq/Results/PBS_Mtb_Tnseq/Analysis/GCF_000195955.2_ASM19595v2_genomic.gbff
CSV_OUTFILE=$OUT_DIR/`basename $FASTQ_DIR`.csv
UNIQUE_FIELDS="locus_tag"
FIELDS="product regulatory_class bound_moiety"

#TPP Settings
PRIMER=AACCTGTTA
MISMATCHES=2
WINDOW_SIZE=6
#BWA_ALG="mem"
BWA_ALG="aln"

##############Beginning of Operations
mkdir -p $OUT_DIR

COUNTER=0
INITIAL_START_TIME=$SECONDS
for FASTQ in $FASTQ_DIR/*R1_001.fastq; do
  (( COUNTER += 1 ))
  echo "******** Run $COUNTER: $FASTQ ********"
  #Note: The next four lines assume a particular way of choosing names for forward and reverse reads. If your specific naming convention differs then you will need to change these lines.
  READS1=$FASTQ
  READS2=${FASTQ/R1_001.fastq/R2_001.fastq}
  OUTNAME=$(basename $FASTQ)
  OUTNAME=${OUTNAME/_R1_001.fastq/}

  ITERATION_START_TIME=$SECONDS
#  $PYTHON $TRANSIT/src/tpp.py -himar1 -bwa $BWA -bwa-alg $BWA_ALG -ref $REFS -replicon-id $REPLICON_ID -reads1 $READS1 -reads2 $READS2 -window-size $WINDOW_SIZE -primer $PRIMER -mismatches $MISMATCHES -output $OUT_DIR/$OUTNAME
  ITERATION_END_TIME=$SECONDS
  (( ITERATION_TIME = ITERATION_END_TIME - ITERATION_START_TIME ))
 
  (( TOTAL_RUN_TIME = SECONDS - INITIAL_START_TIME )) 
  (( CURRENT_AVG = TOTAL_RUN_TIME / COUNTER ))
  echo "******** TPP finished in $ITERATION_TIME seconds! Average iteration time over $COUNTER iterations:  $CURRENT_AVG seconds. ********"
done

echo "TPP iterations done! Performing post-processing operations:"
echo "Checking TA sites hit for all samples..."
NUM_REPLICONS=`echo "$REPLICON_ID" | wc -w`
TA_HITS_OUTFILE=$OUT_DIR/`basename $FASTQ_DIR`_TAs_hit.txt
$PYTHON $TRANSIT/mod_misc_scripts/get_runstats.py $OUT_DIR $NUM_REPLICONS TAs_hit > $TA_HITS_OUTFILE
$PYTHON $TRANSIT/mod_misc_scripts/get_runstats.py $OUT_DIR $NUM_REPLICONS TAs_hit
echo ""
echo "Creating prefixes file with all prefixes from all runs..."
basename -a $OUT_DIR/*.wig | cut -c-13 | uniq > $PREFIXES_OUTFILE
echo "Created '$PREFIXES_OUTFILE'."
echo ""
echo "Creating CSV file with all samples processed by TPP..."
$PYTHON $TRANSIT/mod_misc_scripts/wig_gb_to_csv.py -l $PREFIXES_OUTFILE -g $GENBANK_FILE -u $UNIQUE_FIELDS -f $FIELDS -o $CSV_OUTFILE
echo "Created '$CSV_OUTFILE'."
echo ""
echo "********** TPP driver script finished in a total of $TOTAL_RUN_TIME seconds **********"
