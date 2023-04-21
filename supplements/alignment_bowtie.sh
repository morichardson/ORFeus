#!/bin/bash


usage()
{
  echo "usage: $0 [-u] [-a adapter] [-q phred_cutoff] [-m min_length] [-M max_length]
  [-r riboseq_accession_1] [-r riboseq_accession_2] ... [-r riboseq_accession_n]
  ladder_fasta_file ncrna_fasta_file genome_fasta_file annotations_file output_dir"
  echo "  -u unique            require unique-mapping of a read"
  echo "  -a adapter           adapter sequence for trimming"
  echo "  -q phred_cutoff      quality score cutoff (phred score)"
  echo "  -m min_length        minimum read length to align"
  echo "  -M max_length        maximum read length to align"
  echo "  -M max_length        maximum read length to align"
  echo "  -r riboseq_accs      ribo-seq SRA accession"
  exit 1
}

# Set the argument defaults
unique=false
adapter=""
phred_cutoff=20
min_length=25
max_length=100

# Read the parameter arguments
while getopts "ua:q:m:M:r:" opt; do
  case $opt in
    u) unique=true ;;
    a) adapter=$OPTARG ;;
    q) phred_cutoff=$OPTARG ;;
    m) min_length=$OPTARG ;;
    M) max_length=$OPTARG ;;
    r) riboseq_accs+=("$OPTARG") ;;
  esac
done

shift $((OPTIND -1)) # Shift to the unread arguments

# Read the input/output arguments
ladder_seqs=$1
ncrna_seqs=$2
dna_seqs=$3
annotations_file=$4
output_dir=$5

# Download FASTQ raw reads file(s)
echo Starting fastq prep
touch ${output_dir}riboseq.fastq.gz

for riboseq_acc in "${riboseq_accs[@]}"; do

  # Check if this argument is an existing FASTQ file
  if [[ -f $riboseq_acc ]]; then
    echo "Copying file ${riboseq_acc} to ${output_dir}riboseq.fastq.gz"
    cat $riboseq_acc >> ${output_dir}riboseq.fastq.gz

  # Otherwise, download the FASTQ file for this accession
  else
    echo "Downloading ${riboseq_acc} data to ${output_dir}riboseq.fastq.gz"
    fastq-dump --gzip -O $output_dir $riboseq_acc
    cat ${output_dir}${riboseq_acc}.fastq.gz >> ${output_dir}riboseq.fastq.gz
    rm ${output_dir}${riboseq_acc}.fastq.gz

  fi

done

echo Finished fastq prep


# Create binary index files for the ladder, tRNA, rRNA, and genome fasta files
if [ -f $ladder_seqs ]; then
  echo Started building ladder binary
  ladder_idx=$(echo "$ladder_seqs" | sed -e 's/\.[^.]*$//')
  bowtie-build $ladder_seqs $ladder_idx
  echo Finished building ladder binary
fi

if [ -f $ncrna_seqs ]; then
  echo Started building ncRNA binary
  ncrna_idx=$(echo "$ncrna_seqs" | sed -e 's/\.[^.]*$//')
  bowtie-build $ncrna_seqs $ncrna_idx
  echo Finished building ncRNA binary
fi

echo Started building genome binary
dna_idx=$(echo "$dna_seqs" | sed -e 's/\.[^.]*$//')
bowtie-build $dna_seqs $dna_idx
echo Finished building genome binary

# Trim adapter sequences and filter reads by quality and length
echo Starting cutadapt filtering
cutadapt -a $adapter -q $phred_cutoff -m $min_length -M $max_length -o ${output_dir}filter_trimmed.fastq.gz ${output_dir}riboseq.fastq.gz
rm ${output_dir}riboseq.fastq.gz
gunzip ${output_dir}filter_trimmed.fastq.gz
echo Finished cutadapt filtering

# First, align to ladder index to subtract
if [ -f $ladder_seqs ]; then
  echo Started aligning reads to ladder
  bowtie -v 2 -y -m 1 -a --best --strata -S -p 2 --un ${output_dir}ladder_nomatch.fastq --max ${output_dir}ladder_multi.fastq --al ${output_dir}ladder_match.fastq ${ladder_idx} ${output_dir}filter_trimmed.fastq ${output_dir}ladder_match.sam
  rm ${output_dir}ladder_multi.fastq ${output_dir}ladder_match.fastq ${output_dir}filter_trimmed.fastq ${output_dir}ladder_match.sam
  echo Finished aligning reads to ladder
else
  mv ${output_dir}filter_trimmed.fastq ${output_dir}ladder_nomatch.fastq
fi

# Second, align to ncRNA to subtract
if [ -f $ncrna_seqs ]; then
  echo Started aligning reads to ncRNA
  bowtie -v 2 -y -m 1 -a --best --strata -S -p 2 --un ${output_dir}ncrna_nomatch.fastq --max ${output_dir}ncrna_multi.fastq --al ${output_dir}ncrna_match.fastq ${ncrna_idx} ${output_dir}ladder_nomatch.fastq ${output_dir}ncrna_match.sam
  rm ${output_dir}ncrna_multi.fastq ${output_dir}ncrna_match.fastq ${output_dir}ladder_nomatch.fastq ${output_dir}ncrna_match.sam
  echo Finished aligning reads to ncRNA
else
  mv ${output_dir}ladder_nomatch.fastq ${output_dir}ncrna_nomatch.fastq
fi

# Finally, align to the genome index
# Report all valid alignments (-a) with max 2 mismatches allowed (-v 2)
echo Started aligning reads to genome
if [ $unique ]; then
  bowtie -v 2 -y -m 1 -a --best --strata -S -p 2 --un ${output_dir}genome_nomatch.fastq --max ${output_dir}genome_multi.fastq --al ${output_dir}genome_match.fastq ${dna_idx} ${output_dir}ncrna_nomatch.fastq ${output_dir}genome_match.sam
  rm ${output_dir}genome_multi.fastq ${output_dir}genome_match.fastq ${output_dir}ncrna_nomatch.fastq
else
  bowtie -v 2 -y -a --best --strata -S -p 2 --un ${output_dir}genome_nomatch.fastq --al ${output_dir}genome_match.fastq ${dna_idx} ${output_dir}ncrna_nomatch.fastq ${output_dir}genome_match.sam
  rm ${output_dir}genome_match.fastq ${output_dir}ncrna_nomatch.fastq
fi
echo Finished aligning reads to genome

echo Started converting sam to bam
samtools view -S -b ${output_dir}genome_match.sam | samtools sort - -o ${output_dir}genome_match_sorted.bam
mv ${output_dir}genome_match_sorted.bam ${output_dir}all_hits.bam
rm ${output_dir}*.sam
samtools index ${output_dir}all_hits.bam
echo Finished converting sam to bam


# Calculate the read density by counting the 3’ ends mapped to each genomic position
#echo Started calculating read density
#bedtools genomecov -ibam ${output_dir}all_hits.bam -bg -5 -strand + > ${output_dir}plus.bg
#bedtools genomecov -ibam ${output_dir}all_hits.bam -bg -5 -strand - > ${output_dir}minus.bg
#echo Finished calculating read density
