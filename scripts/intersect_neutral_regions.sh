#!/usr/bin/env bash
set -euo pipefail

# intersect_neutral_regions.sh
# Usage:
#   ./intersect_neutral_regions.sh out_intersect.gff neutral_regions.*.gff[.gz]
# Produces:
#   - GFF3 with regions present (overlapping) in ALL input neutral_regions.*.gff files.
#   - BED (0-based, half-open) file of the same intersected regions.
# Requires: bedtools, awk, sort, zcat/gunzip

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 out_intersect.gff neutral_regions.*.gff[.gz]"
  exit 1
fi

OUT_GFF="$1"; shift
INPUTS=( "$@" )

# derive OUT_BED from OUT_GFF (replace .gff suffix, else append .bed)
if [[ "$OUT_GFF" == *.gff ]]; then
  OUT_BED="${OUT_GFF%.gff}.bed"
else
  OUT_BED="${OUT_GFF}.bed"
fi

# expand globs (if any were quoted)
expanded=()
for p in "${INPUTS[@]}"; do
  if compgen -G "$p" >/dev/null; then
    for f in $p; do expanded+=("$f"); done
  else
    expanded+=("$p")
  fi
done
INPUTS=( "${expanded[@]}" )

if [[ ${#INPUTS[@]} -lt 1 ]]; then
  echo "No input files found."
  exit 1
fi

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found in PATH"; exit 1; }

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

BED_FILES=()
i=0
for f in "${INPUTS[@]}"; do
  i=$((i+1))
  if [[ ! -f "$f" ]]; then
    echo "File not found: $f"
    exit 1
  fi
  base="$(basename "$f")"
  bed="$TMPDIR/${base%.gz}.bed"
  # convert GFF -> BED (0-based half-open)
  if [[ "$f" == *.gz ]]; then
    zcat "$f" | \
      awk 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {
        start = $4 - 1; end = $5;
        if (end>start) print $1, start, end;
      }' | sort -k1,1 -k2,2n > "$bed"
  else
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {
      start = $4 - 1; end = $5;
      if (end>start) print $1, start, end;
    }' "$f" | sort -k1,1 -k2,2n > "$bed"
  fi
  BED_FILES+=( "$bed" )
  # record mapping idx -> basename (for later)
  echo -e "$i\t$base" >> "$TMPDIR/file_index_map.tsv"
done

NFILES=${#BED_FILES[@]}

# run bedtools multiinter
MULTI="$TMPDIR/multiinter.bed"
bedtools multiinter -i "${BED_FILES[@]}" > "$MULTI"

# Keep only regions present in all input files (column 4 == NFILES)
FILTERED="$TMPDIR/intersect_filtered.bed"
awk -v N="$NFILES" 'BEGIN{OFS="\t"} $4==N{print $1,$2,$3,$4,$5}' "$MULTI" > "$FILTERED" || true

# If no intersection, write minimal outputs and exit
if [[ ! -s "$FILTERED" ]]; then
  echo "##gff-version 3" > "$OUT_GFF"
  echo "# No regions present in all $NFILES files" >> "$OUT_GFF"
  # create an empty/annotated BED
  {
    echo "# No regions present in all $NFILES files"
  } > "$OUT_BED"
  echo "Wrote empty GFF to $OUT_GFF and placeholder BED to $OUT_BED"
  exit 0
fi

# Translate index list (col5 like "1,2,3") to comma-separated basenames using awk:
# build the GFF body and the BED file in one awk pass
awk -v OFS="\t" '
  BEGIN{ FS="\t" }
  NR==FNR {
    idx=$1; name=$2;
    map[idx]=name;
    next
  }
  {
    chrom=$1; start0=$2; end=$3; support=$4; idxlist=$5;
    start = start0 + 1;             # 1-based inclusive for GFF start
    # build SupportFiles by mapping indices
    n=split(idxlist, a, ",");
    files="";
    for(i=1;i<=n;i++){
      if(i>1) files=files","map[a[i]];
      else files=map[a[i]];
    }
    id = "intersect_"NR;
    # print GFF line
    printf "%s\tintersection\tregion\t%d\t%d\t.\t.\t.\tID=%s;SupportCount=%s;SupportFiles=%s\n", chrom, start, end, id, support, files;
    # print BED line (0-based half-open): chrom, start0, end, name, score
    # use name=id and score=support
    printf "%s\t%d\t%d\t%s\t%s\n", chrom, start0, end, id, support;
  }
' "$TMPDIR/file_index_map.tsv" "$FILTERED" > "$TMPDIR/intersect_outputs.tmp"

# Split the combined file into GFF body and BED
# The awk output has GFF line then BED line alternating; separate them
awk 'NR%2==1{print > "'"$TMPDIR"'/intersect.gff.body"} NR%2==0{print > "'"$TMPDIR"'/intersect.bed.body"}' "$TMPDIR/intersect_outputs.tmp"

# write final GFF3
{
  echo "##gff-version 3"
  cat "$TMPDIR/intersect.gff.body"
} > "$OUT_GFF"

# write final BED (ensure sorted)
sort -k1,1 -k2,2n "$TMPDIR/intersect.bed.body" > "$OUT_BED"

echo "Wrote intersection GFF3 to $OUT_GFF and BED to $OUT_BED (regions present in all $NFILES input files)."
