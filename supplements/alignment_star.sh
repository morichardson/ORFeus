#!/bin/bash


usage()
{
  echo "usage: $0 [-a adapter] [-q phred_cutoff] [-m min_length] [-M max_length]
  [-r riboseq_accession_1] [-r riboseq_accession_2] ... [-r riboseq_accession_n]
  ladder_fasta_file ncrna_fasta_file genome_fasta_file annotations_file output_dir"
  echo "  -a adapter           adapter sequence for trimming"
  echo "  -q phred_cutoff      quality score cutoff (phred score)"
  echo "  -m min_length        minimum read length to align"
  echo "  -M max_length        maximum read length to align"
  echo "  -M max_length        maximum read length to align"
  echo "  -r riboseq_accs      ribo-seq SRA accession"
  exit 1
}

# Set the argument defaults
adapter=""
phred_cutoff=20
min_length=25
max_length=100

# Read the parameter arguments
while getopts "a:q:m:M:r:" opt; do
  case $opt in
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
echo Started aligning reads to genome [unique]
dna_idx=$(echo $(dirname "$dna_seqs"))
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $dna_idx --genomeFastaFiles $dna_seqs --sjdbGTFfile $annotations_file --outFileNamePrefix ${output_dir}
STAR --runThreadN 4 --genomeDir $dna_idx --readFilesIn ${output_dir}ncrna_nomatch.fastq --outFileNamePrefix ${output_dir} --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD --outSAMattributes NH HI AS nM NM MD --outSAMmultNmax 1 --outFilterMultimapNmax -1 --outFilterMismatchNmax 2
samtools view -h -q 255 ${output_dir}Aligned.sortedByCoord.out.bam > ${output_dir}Unique.sortedByCoord.out.bam
samtools sort ${output_dir}Unique.sortedByCoord.out.bam > ${output_dir}all_hits.bam
samtools index ${output_dir}all_hits.bam
rm ${output_dir}ncrna_nomatch.fastq ${output_dir}Aligned.sortedByCoord.out.bam ${output_dir}Unique.sortedByCoord.out.bam
rm ${output_dir}Log.progress.out ${output_dir}Log.final.out ${output_dir}Log.out ${output_dir}SJ.out.tab
rm -r ${output_dir}_STARtmp
echo Finished aligning reads to genome

# Calculate the read density by counting the 5’ ends mapped to each genomic position
#echo Started calculating read density
#bedtools genomecov -ibam ${output_dir}all_hits.bam -bg -5 -strand + > ${output_dir}plus.bg
#bedtools genomecov -ibam ${output_dir}all_hits.bam -bg -5 -strand - > ${output_dir}minus.bg
#echo Finished calculating read density
