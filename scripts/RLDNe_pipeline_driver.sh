#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  ./RLDNe_pipeline_driver.sh --vcf FILE --s2p FILE --outdir DIR [--neut-bed FILE] [--thin-list "0,1000,5000,10000"] [--single THIN DATASET]
Options:
  --vcf         Original or compressed VCF (.vcf.gz)
  --s2p         sample->population mapping file (sample<TAB>pop)
  --outdir      Output directory
  --neut-bed    (Optional) neutral regions BED file (if wanting to generate neutral data)
  --thin-list   Comma-separated list of thinning values (bp) (default "0,1000,5000,10000"; 0 means no thinning)
  --single THIN DATASET
                Single run mode: process only this combination (suitable for parallel calls). DATASET is "all" or "neutral"
Examples:
  # All combinations (serial)
  ./RLDNe_pipeline_driver.sh --vcf in.vcf.gz --s2p sample2pop.tsv --outdir out

  # Specify lists
  ./RLDNe_pipeline_driver.sh --vcf in.vcf.gz --s2p sample2pop.tsv --outdir out --thin-list "0,5000"

  # Single mode (suitable for parallel submission)
  ./RLDNe_pipeline_driver.sh --vcf in.vcf.gz --s2p sample2pop.tsv --outdir out --single 1000 all

  # GNU parallel example (need to generate combos.txt first)
  # Each line in combos.txt: THIN DATASET
  parallel -a combos.txt -j 8 --colsep ' ' ./RLDNe_pipeline_driver.sh --vcf in.vcf.gz --s2p sample2pop.tsv --outdir out --neut-bed neutral.bed --single {1} {2}
EOF
  exit 1
}

# defaults
THIN_LIST_DEFAULT="0,1000,5000,10000"

VCF=""
S2P=""
OUTDIR=""
NEUT_BED=""
THIN_LIST="$THIN_LIST_DEFAULT"
SINGLE_MODE=0
SINGLE_THIN=""
SINGLE_DS=""

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --s2p) S2P="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --neut-bed) NEUT_BED="$2"; shift 2;;
    --thin-list) THIN_LIST="$2"; shift 2;;
    --single) SINGLE_MODE=1; SINGLE_THIN="$2"; SINGLE_DS="$3"; shift 3;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

if [[ -z "$VCF" || -z "$S2P" || -z "$OUTDIR" ]]; then
  usage
fi

command -v vcftools >/dev/null || { echo "vcftools is required"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools is required"; exit 1; }
command -v bgzip >/dev/null || { echo "bgzip is required"; exit 1; }
if [[ ! -x ./convert_vcf_to_RLDNe.sh ]]; then
  echo "convert_vcf_to_RLDNe.sh needs to be in current directory and executable"
  exit 1
fi

mkdir -p "$OUTDIR"

# helper: create neutral vcf if not exists
NEUT_VCF="$OUTDIR/filtered_neutral.vcf.gz"
prepare_neutral() {
  if [[ -z "$NEUT_BED" ]]; then
    echo "ERROR: neutral dataset requested but --neut-bed not provided"
    exit 1
  fi
  if [[ ! -f "$NEUT_VCF" ]]; then
    echo "Generating neutral VCF to $NEUT_VCF"
    bcftools view -R "$NEUT_BED" -m2 -M2 -v snps "$VCF" -Oz -o "$NEUT_VCF"
    bcftools index -f "$NEUT_VCF"
  else
    echo "Reusing existing neutral VCF: $NEUT_VCF"
  fi
}

# worker: process one combination (dataset in {all,neutral})
process_one() {
  local thin="$1"
  local ds="$2"

  local base_vcf="$VCF"
  if [[ "$ds" == "neutral" ]]; then
    prepare_neutral
    base_vcf="$NEUT_VCF"
    dsname="neutral"
  else
    dsname="all"
  fi

  tag="ds-${dsname}_thin-${thin}"
  work_v="$OUTDIR/${tag}.vcf"
  outprefix="$OUTDIR/${tag}"

  if [[ "$thin" == "0" ]]; then
    vcftools --gzvcf "$base_vcf" --recode --recode-INFO-all --out "$OUTDIR/${tag}" >/dev/null
    mv "$OUTDIR/${tag}.recode.vcf" "$work_v"
  else
    vcftools --gzvcf "$base_vcf" --thin "$thin" --recode --recode-INFO-all --out "$OUTDIR/${tag}" >/dev/null
    mv "$OUTDIR/${tag}.recode.vcf" "$work_v"
  fi

  bgzip -f "$work_v"
  work_vgz="${work_v}.gz"
  bcftools index -f "$work_vgz"

  ./convert_vcf_to_RLDNe.sh "$work_vgz" "$S2P" "$outprefix"
  echo "Generated: $outprefix.RLDNe.tsv"
}

# main
if [[ $SINGLE_MODE -eq 1 ]]; then
  # single combo mode
  if [[ -z "$SINGLE_THIN" || -z "$SINGLE_DS" ]]; then
    echo "Single mode requires 2 parameters: THIN DATASET"
    exit 1
  fi
  process_one "$SINGLE_THIN" "$SINGLE_DS"
  exit 0
fi

# otherwise: parse lists and run all combos (serial)
IFS=',' read -r -a THINLIST_ARR <<< "$THIN_LIST"

DATASETS=("all")
if [[ -n "$NEUT_BED" ]]; then DATASETS+=("neutral"); fi

for ds in "${DATASETS[@]}"; do
    for thin in "${THINLIST_ARR[@]}"; do
      echo "Processing ds=$ds thin=$thin"
      process_one "$thin" "$ds"
    done
done


echo "All completed, output in $OUTDIR"
