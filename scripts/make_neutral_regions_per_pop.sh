#!/usr/bin/env bash
set -euo pipefail
#
# make_neutral_regions_per_pop.sh
# input：
#   1.gff
#   SingleCopyGene.txt  (orthogroup <tab> gene_id)
#   <pop>_BSM_sig_omega.txt  (one orthogroup id per line)
# output（each pop）:
#   neutral_regions.<POP>.bed   (0-based, half-open)
#   neutral_regions.<POP>.gff   (gff3-like: CDS lines for neutral single-copy genes + artificial intergenic lines)
#
# usage: ./make_neutral_regions_per_pop.sh POPLIST "G2 G4 G5 G7 G8"  
#

# ------ config ------
GFF="1.gff"
SINGLE_COPY="SingleCopyGene.txt"
POPS=(G2 G4 G5 G7 G8)   # pop list
OUTDIR="neutral_regions"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${OUTDIR}" "${TMPDIR}"
# -----------------------------

for cmd in awk sed bedtools sort bgzip tabix bcftools; do
  if ! command -v ${cmd} >/dev/null 2>&1; then
    echo "ERROR: command '${cmd}' not found in PATH. Please install. (bedtools, bgzip/tabix from htslib, bcftools, awk, sed, sort)" >&2
    exit 1
  fi
done

# extract CDS/gene row from GFF 
GFF_CDS="${TMPDIR}/cds_lines.tsv"
awk 'BEGIN{FS="\t";OFS="\t"} $1!~/^#/ && ($3=="CDS" || $3=="gene" || $3=="CDS_region"){ 
    gid="."; 
    # try to extract ID or Name from column 9
    if($9 ~ /ID=/){ match($9, /ID=([^;]+)/, m); gid=m[1]; }
    else if($9 ~ /Name=/){ match($9, /Name=([^;]+)/, m); gid=m[1]; }
    else if($9 ~ /gene=/){ match($9, /gene=([^;]+)/, m); gid=m[1]; }
    print $1,$4-1,$5,$7,$9,gid;
}' "${GFF}" > "${GFF_CDS}"

# If no gid found, try to use column 9 token as gene id fallback
awk 'BEGIN{FS="\t";OFS="\t"} { if($6=="."){ # parse first token of attr
      split($5,a,"[;= ]"); $6=a[1];
    } print }' "${GFF_CDS}" > "${GFF_CDS}.tmp" && mv "${GFF_CDS}.tmp" "${GFF_CDS}"

# Create BED of all CDS (0-based half-open: chrom, start, end, geneid, ., strand)
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$6,".",$4}' "${GFF_CDS}" > "${TMPDIR}/all_cds.bed"

# Merge CDS per contig to get occupied coding regions
sort -k1,1 -k2,2n "${TMPDIR}/all_cds.bed" > "${TMPDIR}/all_cds.sorted.bed"
bedtools merge -i "${TMPDIR}/all_cds.sorted.bed" -c 4,6 -o distinct,distinct > "${TMPDIR}/all_cds.merged.bed"

# Compute intergenic regions = complement of merged CDS within each contig.
# Need contig lengths -> derive from GFF (max end per seqid)
awk 'BEGIN{FS="\t";OFS="\t"} $1!~/^#/ {if($1!=""){if($1 in max){ if($5+0>max[$1]) max[$1]=$5+0 } else max[$1]=$5+0 }} END{for(i in max) print i, max[i]}' "${GFF}" > "${TMPDIR}/contig_maxlen.tsv"
# convert to genome file for bedtools (chrom length)
awk 'BEGIN{FS=OFS="\t"} {print $1,0,$2}' "${TMPDIR}/contig_maxlen.tsv" > "${TMPDIR}/genome.size.bed"

# bedtools complement (gives intergenic regions)
bedtools complement -i "${TMPDIR}/all_cds.merged.bed" -g "${TMPDIR}/genome.size.bed" > "${TMPDIR}/intergenic.bed"

# Standardize intergenic BED columns (chrom,start,end,.,.,.)
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,"intergenic",".","+"}' "${TMPDIR}/intergenic.bed" > "${TMPDIR}/intergenic.formatted.bed"

# Preprocess SingleCopyGene.txt -> mapping orthogroup -> gene_id(s)
# Expected format: OG0001<TAB>geneX (one per line per gene). We'll group gene IDs per orthogroup.
SINGLE=${SINGLE_COPY}
if [[ ! -f "${SINGLE}" ]]; then
  echo "ERROR: SingleCopyGene.txt not found: ${SINGLE}" >&2
  exit 1
fi
# create mapping file orthogroup<TAB>geneid
awk 'BEGIN{FS=OFS="\t"} NF>=2{print $1,$2}' "${SINGLE}" > "${TMPDIR}/orth2gene.tsv"

# Index gene -> BED line for quick lookup
# all_cds.bed has gene_id in col4
awk 'BEGIN{FS=OFS="\t"} {print $4,$0}' "${TMPDIR}/all_cds.bed" > "${TMPDIR}/gene2bed.tsv"  # geneid TAB bedcols...

# For each population: produce neutral bed/gff
for pop in "${POPS[@]}"; do
  selfile="${pop}_BSM_sig_omega.txt"
  echo "Processing ${pop} ..."
  selected_ogs="/dev/null"
  if [[ -f "${selfile}" ]]; then
    # get list of selected orthogroup ids (first column)
    awk 'NF>=1{print $1}' "${selfile}" > "${TMPDIR}/${pop}.selected_ogs.txt"
    selected_ogs="${TMPDIR}/${pop}.selected_ogs.txt"
  else
    echo "  Warning: selection file ${selfile} not found -> assume no selected orthogroups for ${pop}."
    > "${TMPDIR}/${pop}.selected_ogs.txt"
    selected_ogs="${TMPDIR}/${pop}.selected_ogs.txt"
  fi

  # get list of single-copy orthogroups present (from orth2gene) but excluding selected ones
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{sel[$1]=1; next} { if(!( $1 in sel )) print $1,$2 }' "${selected_ogs}" "${TMPDIR}/orth2gene.tsv" > "${TMPDIR}/${pop}.orth2gene.neutral.tsv"
  # This file contains orthogroup<TAB>geneid for single-copy OGs NOT under positive selection for this pop

  # From geneids get corresponding CDS bed lines
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{g[$1]=1; next} { if($1 in g) { print $2,$3,$4,$5,$6,$7 } }' "${TMPDIR}/${pop}.orth2gene.neutral.tsv" "${TMPDIR}/gene2bed.tsv" > "${TMPDIR}/${pop}.neutral_coding.bed"

  # merge coding neutral bed (may be overlapping due to paralog coordinates etc)
  sort -k1,1 -k2,2n "${TMPDIR}/${pop}.neutral_coding.bed" > "${TMPDIR}/${pop}.neutral_coding.sorted.bed"
  if [[ -s "${TMPDIR}/${pop}.neutral_coding.sorted.bed" ]]; then
    bedtools merge -i "${TMPDIR}/${pop}.neutral_coding.sorted.bed" -c 4,6 -o distinct,distinct > "${TMPDIR}/${pop}.neutral_coding.merged.bed"
    awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,"coding_neutral",".","+"}' "${TMPDIR}/${pop}.neutral_coding.merged.bed" > "${TMPDIR}/${pop}.neutral_coding.formatted.bed"
  else
    # empty neutral coding (no single-copy neutral genes)
    > "${TMPDIR}/${pop}.neutral_coding.formatted.bed"
  fi

  # Combine neutral coding + intergenic -> full neutral regions
  cat "${TMPDIR}/${pop}.neutral_coding.formatted.bed" "${TMPDIR}/intergenic.formatted.bed" \
    | sort -k1,1 -k2,2n > "${TMPDIR}/${pop}.neutral_combined.sorted.bed"

  # Merge to remove overlaps and produce final neutral bed
  bedtools merge -i "${TMPDIR}/${pop}.neutral_combined.sorted.bed" -c 4,5,6 -o distinct,distinct,distinct > "${OUTDIR}/neutral_regions.${pop}.bed"

  # Create neutral GFF: include original CDS GFF lines for neutral coding genes + synthetic intergenic GFF lines
  # Extract original GFF lines for coding genes in neutral set (we use gene IDs from orth2gene.neutral)
  # gene IDs list:
  awk 'BEGIN{FS=OFS="\t"} {print $2}' "${TMPDIR}/${pop}.orth2gene.neutral.tsv" | sort -u > "${TMPDIR}/${pop}.neutral_geneids.txt" || true
  # find original GFF lines containing these IDs (simple grep -F)
  if [[ -s "${TMPDIR}/${pop}.neutral_geneids.txt" ]]; then
    # build grep pattern file
    awk '{print $0}' "${TMPDIR}/${pop}.neutral_geneids.txt" > "${TMPDIR}/${pop}.neutral_geneids.grep"
    # grep -F -f will match attributes
    grep -F -f "${TMPDIR}/${pop}.neutral_geneids.grep" "${GFF}" > "${OUTDIR}/neutral_regions.${pop}.coding.gff" || true
  else
    > "${OUTDIR}/neutral_regions.${pop}.coding.gff"
  fi

  # create intergenic gff entries from neutral bed: fields: seqid source feature start end . strand phase attributes
  awk 'BEGIN{FS=OFS="\t"} {print $1,"make_neutral","intergenic",$2+1,$3,".",$6,".", "ID=intergenic_"NR";Note=intergenic_region"}' "${OUTDIR}/neutral_regions.${pop}.bed" > "${OUTDIR}/neutral_regions.${pop}.intergenic.gff"

  # concatenate coding + intergenic -> final gff (simple)
  (echo "##gff-version 3"; cat "${OUTDIR}/neutral_regions.${pop}.coding.gff"; cat "${OUTDIR}/neutral_regions.${pop}.intergenic.gff") > "${OUTDIR}/neutral_regions.${pop}.gff"

  echo "  Wrote: ${OUTDIR}/neutral_regions.${pop}.bed  and  ${OUTDIR}/neutral_regions.${pop}.gff"
done

echo "All populations processed. Neutral files in ${OUTDIR}/"
