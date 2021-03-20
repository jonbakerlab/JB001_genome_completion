# ------ Unicycler Assembly -------
# Use minimap2 to map reads from the metagenomic long-read set to the draft genome
minimap2 \
    -ax map-ont \
    -t 48 \
    --sam-hit-only \
    draft_mag_assembly.fasta nanopore_reads.fastq > mapped_reads.sam

# Convert SAM file to BAM file
samtools view \
    -b \
    -@ 48 \
    mapped_reads.sam > mapped_reads.bam

# Use anvio to sort and index BAM file
anvi-init-bam \
    assembly.bam \
    -o assembly.sorted.bam \
    -T 48

# Extract reads that mapped for use in assembly
samtools fastq \
    -F 4 \
    -@ 48 \
    assembly.sorted.bam > mapped_nanopore_reads.fastq

# Perform hybrid assembly using Unicycler
unicycler \
    -1 short_reads_mapped_to_draft_1.fastq \
    -2 short_reads_mapped_to_draft_2.fastq \
    -l mapped_nanopore_reads.fastq \
    --threads 48 \
    -o unicycler

# The resulting genome had 4 contigs, one ~650kb circular contig, and 3 short (<6000bp) contigs.  
# The contigs were manually examined using Anvi'o with the scripts below. 
# The full pipeline to get to the Anvi'o profile stage can be found at http://merenlab.org/software/anvio/.
#   -The 3 short contigs were discarded on the basis of different coverage and GC content compared to the main contig, as well as blast hits 
#   to other organisms within the 3 short contigs

anvi-profile \
    -i assembly2.sorted.bam \
    -c assembly.contigs.db \
    -o assembly.profile \
    -S JB001_final \
    --min-contig-length 500 \ # no contigs in the Unicycler assembly were below this threshold
    --cluster-contigs \
    -T 48

anvi-refine \
-p assembly.profile/PROFILE.db \
-c assembly.contigs.db \
-C DEFAULT

# ------ Flye Assembly -------
# Remove reads mapping to the human genome
minimap2 \
    -ax map-ont \
    -t 48 \
    refseq_human.fna nanopore_reads.fastq > mapped-to-human.sam

# Convert SAM file to BAM file
samtools view \
    -b \
    -@ 48 \
    mapped-to-human.sam > mapped-to-human.bam

# Use anvio to sort and index BAM file
anvi-init-bam \
    mapped-to-human.bam \
    -o mapped-to-human.sorted.bam \
    -T 48

# Extract reads that did not map to the human genome
samtools fastq \
    -f 4 \
    -@ 48 \
    mapped-to-human.sorted.bam > nanopore_reads_nohuman.fastq

# Assemble the remaining long reads with meta-flye
flye \
    --nano-raw nanopore_reads_nohuman.fastq \
    --genome-size 100m \
    --threads 48 \
    --meta \
    --out-dir flye-assembly

# Identify the JB001 contig(s) to extract from within the metagenome
megablast \
    -d draft_assembly.fasta \
    -i flye_metagenome_assembly.fasta \
    -D 3 \
    -a 48 \
    -o flye-blast-for-jb001.txt


# ------ Use Trycycler to create a consensus assembly using the flye and Unicycler draft assemblies ------
trycycler cluster \
    --threads 56 \
    --assemblies assemblies/*.fasta \
    --reads mapped_nanopore_reads.fastq \
    --out_dir trycycler

trycycler reconcile \
     -t 56 \
    --reads mapped_nanopore_reads.fastq \
    --min_identity 94 \
    --max_indel_size 410 \
    --cluster_dir cluster_001/

trycycler msa \
    --cluster_dir cluster_001/ \
    --threads 56

trycycler partition \
    --reads mapped_nanopore_reads.fastq \
    --cluster_dirs cluster_001/

trycycler consensus \
    --cluster_dir cluster_001/ \
    --verbose \
    --threads 48

# ------ Polish the consensus assembly using the long reads ------
medaka_consensus \
    -i mapped_nanopore_reads.fastq \
    -d cluster_001/7_final_consensus.fasta \
    -o medaka \
    -m r941_min_high_g360 \
    -t 48


# ------ Polish the consensus assembly using the Illumina short reads ------
# You should repeat Pilon  multiple times until Pilon stops making changes to the assembly
# For further instructions for this section, consult https://github.com/rrwick/Trycycler/wiki/Polishing-after-Trycycler 

bowtie2-build after_medaka.fasta after_medaka.fasta

bowtie2 \
    -1 short_reads_1.fastq \
    -2 short_reads_2.fastq \
    -x after_medaka.fasta \
    --threads 48 \
    -I "$insert_min" \      # You will need to determine the insert min and max for your dataset, instructions at https://github.com/rrwick/Trycycler/wiki/Polishing-after-Trycycler 
    -X "$insert_max" \
    --local \
    --very-sensitive-local | samtools sort > illumina_alignments.bam

samtools index illumina_alignments.bam

pilon \
    --genome before.fasta \
    --frags illumina_alignments.bam \
    --output after_pilon \
    --changes \
    --threads 48

# ------ Rotate the final assembly to start with the JB001 dnaA gene ------
circlator fixstart \
    --verbose \
    --genes_fa tm7-g6-dnaa.fasta \
    after-pilon.fasta \
    fixstart_prefix