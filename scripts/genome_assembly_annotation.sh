#!/bin/bash

#define sample name
id=APE-7

#assembly and binning
#for APE-1 to APE-6
spades.py -o $id'_spades_assembly_result' \
-1 $id'_R1.fq' \
-2 $id'_R2.fq' \
-k 23,43,63,83,99,113 \
--careful --cov-cutoff auto -m 800 -t 60

cd $id'_spades_assembly_result'
mv contigs.fasta $id'.contigs.fasta'
metawrap binning -o INITIAL_BINNING -t 30 -a $id'.contigs.fasta' \
--metabat2 --maxbin2 --concoct ../$id'_R1.fq' ../$id'_R2.fq'

metawrap bin_refinement -o BIN_REFINEMENT -t 30 -A INITIAL_BINNING/metabat2_bins/ \
-B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 95 -x 5

#for other strains
flye --nano-raw $id'_clean.fq' \
-o $id'_flye_meta_assembly_result' \
--meta -t 40 -i 3

cd $id'_flye_meta_assembly_result'
bwa index assembly.fasta
bwa mem -t $t assembly.fasta ../$id'_R1.fq' ../$id'_R2.fq' |\
samtools view - -Sb | samtools sort - -@ $t -o illumina.sorted.bam

samtools index illumina.sorted.bam
java -Xmx320G -jar ~/data2/software/pilon-1.23.jar --genome assembly.fasta \
--fix all --changes --frags illumina.sorted.bam --output assembly_pillon \
--outdir assembly_pillon --threads 20 --vcf 2

bwa index assembly_pillon/assembly_pillon.fasta
bwa mem -t $t assembly_pillon/assembly_pillon.fasta ../$id'_R1.fq' ../$id'_R2.fq' |\
samtools view - -Sb | samtools sort - -@ $t -o illumina.sorted.bam

samtools index illumina.sorted.bam
java -Xmx320G -jar ~/data2/software/pilon-1.23.jar --genome assembly_pillon/assembly_pillon.fasta \
--fix all --changes --frags illumina.sorted.bam --output assembly_pillon2 \
--outdir assembly_pillon --threads 20 --vcf 2

metawrap binning -o INITIAL_BINNING -t 50 -a assembly_pillon/assembly_pillon2.fasta \
--metabat2 --maxbin2 --concoct ../$id'_R1.fq' ../$id'_R2.fq'

metawrap bin_refinement -o BIN_REFINEMENT -t 50 -A INITIAL_BINNING/metabat2_bins/ \
-B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 95 -x 5

#select the bin classified as Cyanobacteria 
mkdir -p ~/genome/$id
mv bin1.fa ~/genome/$id/$id'_genomic.fna'

#annotation
cd ~/genome/$id

#16S rRNA
~/data2/software/rnammer-1.2/rnammer -S bac -gff $id'_16srrna.gff' \
-m ssu -f $id'.16srRNA.fasta' $id'_genomic.fna'

#gene prediction
prodigal -a $id'.protein.fa' -d $id'.gene.fa' -f gff -g 11 \
-i $id'_genomic.fna' -o $id'.gff' -s $id'.stat' -p single -q -m
sed 's/\s#\s.*//' $id'.gene.fa' > $id'_gene.fna'
sed 's/\s#\s.*//' $id'.protein.fa' > $id'_protein.faa'
sed -i 's/*//' $id'_protein.faa'

#KEGG database
diamond blastp -d ~/data2/database/kegg/diamond/Metagenomes \
-q $id'_protein.faa' -o $id'_kegg' \
-k 1 -e 1e-5 --id 40 --query-cover 50 --subject-cover 50 -f 100 -p 40

diamond view -a $id'_kegg.daa' -o $id'_kegg.table' -f 6 -p 40
perl ~/Scripts/blast_out_annotation.pl $id'_kegg.table' ~/data2/database/kegg/Metagenomes.txt $id'_kegg_anno.table'

#eggNOG5 database
python3 ~/data2/software/eggnog-mapper-2.1.8/emapper.py -i $id'_protein.faa' \
-m diamond --pident 40 --query_cover 50 --subject_cover 50 \
--output $id'_nog' --cpu 30

#codon usage
~/data2/software/codonW/codonw $id'_gene.fna' -all_indices $id'_gene.out' $id'_gene.blk' -nomenu -nowarn -silent -coa_cu -rscu