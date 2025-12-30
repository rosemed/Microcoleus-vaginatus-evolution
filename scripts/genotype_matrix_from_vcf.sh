#!/usr/bin/env bash
set -euo pipefail

# usage: ./genotype_matrix_from_vcf.sh input.vcf.gz out.tsv "<pop_pattern>" [neutral.bed]

VCF="$1"
OUT="$2"
POP_PATTERN="$3"
NEUT="${4:-}"

command -v bcftools >/dev/null || { echo "bcftools required"; exit 1; }
command -v python3 >/dev/null || { echo "python3 required"; exit 1; }

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

if [ -n "$NEUT" ]; then
  bcftools view -R "$NEUT" -m2 -M2 -v snps -f PASS "$VCF" -Ou \
    | bcftools view -i 'F_MISSING<=0.1 && QUAL>=30' -Oz -o "$TMPDIR/filtered.vcf.gz"
else
  bcftools view -m2 -M2 -v snps -f PASS "$VCF" -Ou \
    | bcftools view -i 'F_MISSING<=0.1 && QUAL>=30' -Oz -o "$TMPDIR/filtered.vcf.gz"
fi
bcftools index -f "$TMPDIR/filtered.vcf.gz"
bcftools query -f '%CHROM_%POS[\t%GT]\n' "$TMPDIR/filtered.vcf.gz" > "$TMPDIR/genotypes.tsv"
bcftools query -l "$TMPDIR/filtered.vcf.gz" > "$TMPDIR/samples.txt"

python3 - "$TMPDIR/genotypes.tsv" "$TMPDIR/samples.txt" "$OUT" "$POP_PATTERN" <<'PY'
import sys,glob
geno_tsv = sys.argv[1]
samples_file = sys.argv[2]
out_tsv = sys.argv[3]
pop_pattern = sys.argv[4]

with open(samples_file) as f:
    samples = [x.strip() for x in f if x.strip()]
locus_ids = []
matrix = []
with open(geno_tsv) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 2: continue
        locus = parts[0]
        gts = parts[1:]
        locus_ids.append(locus)
        matrix.append(gts)

s2i = {s:i for i,s in enumerate(samples)}
# build sample rows
rows = {s:[] for s in samples}
for col_idx in range(len(locus_ids)):
    col = [matrix[col_idx][i] for i in range(len(samples))]
    for i,s in enumerate(samples):
        gt = col[i]
        if gt is None or gt=='.' or gt=='./.' or gt=='./.':
            val = "NA"
        else:
            gt_token = gt.split(':')[0].replace('|','/')
            if '/' in gt_token:
                parts = gt_token.split('/')
                if '1' in parts:
                    val = "1"
                elif '0' in parts:
                    val = "0"
                else:
                    val = "NA"
            else:
                if gt_token == '.':
                    val = "NA"
                else:
                    try:
                        if int(gt_token) >= 1:
                            val = "1"
                        else:
                            val = "0"
                    except:
                        val = "NA"
        rows[samples[i]].append(val)

with open(out_tsv,'w') as O:
    O.write("sample\t" + "\t".join(locus_ids) + "\n")
    for s in samples:
        O.write(s + "\t" + "\t".join(rows[s]) + "\n")
print("Writing matrix: %s" % out_tsv)
PY

echo "Done: $OUT"
