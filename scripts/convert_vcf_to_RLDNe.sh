#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./convert_vcf_to_RLDNe.sh filtered.vcf.gz sample2pop.tsv out_prefix
# sample2pop.tsv format: sample_name<TAB>pop_name   (header optional, will be ignored if text)
VCF="$1"
SAMPLE2POP="$2"
OUTPREFIX="$3"

command -v bcftools >/dev/null || { echo "bcftools required"; exit 1; }
command -v python3   >/dev/null || { echo "python3 required"; exit 1; }

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# index if needed
if [ ! -f "${VCF}.csi" ] && [ ! -f "${VCF}.tbi" ]; then
  bcftools index -f "$VCF"
fi

# export CHR_POS and GTs (biallelic assumed)
bcftools query -f '%CHROM\_%POS[\t%GT]\n' "$VCF" > "$TMPDIR/genotypes.tsv"
bcftools query -l "$VCF" > "$TMPDIR/samples.txt"

# call python converter
python3 - "$TMPDIR/genotypes.tsv" "$TMPDIR/samples.txt" "$SAMPLE2POP" "$OUTPREFIX" <<'PY'
import sys
geno_tsv = sys.argv[1]
samples_file = sys.argv[2]
sample2pop_file = sys.argv[3]
outprefix = sys.argv[4]
out_file = outprefix + ".RLDNe.tsv"

# read samples (VCF order)
with open(samples_file) as f:
    samples = [x.strip() for x in f if x.strip()]

# read sample->pop mapping
s2pop = {}
with open(sample2pop_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            s2pop[parts[0]] = parts[1]

# read loci and matrix (rows loci)
loci = []
matrix = []
with open(geno_tsv) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        loci.append(parts[0])
        matrix.append(parts[1:])

n_samples = len(samples)
n_loci = len(loci)

# mapping: REF -> 1, ALT -> 2, missing -> 0
# For haploid: duplicate allele into .A1 and .A2
def gt_to_alleles(gt):
    if gt is None or gt == '.' or gt == './.' or gt == '.|.':
        return ("0","0")
    gt_tok = gt.split(':')[0].replace('|','/')
    if '/' in gt_tok:
        parts = gt_tok.split('/')
    else:
        parts = [gt_tok]
    # try to find an allele: if any allele == 1 -> ALT (code 2), else if 0 -> REF (code 1)
    if any(p == '1' for p in parts):
        a = "2"
    elif any(p == '0' for p in parts):
        a = "1"
    else:
        # other values (shouldn't happen after biallelic filter) => mark missing
        a = "0"
    # duplicate for haploid
    return (a,a)

# Build per-sample rows
rows = {s: [] for s in samples}
for li in range(n_loci):
    col = matrix[li]
    for i,s in enumerate(samples):
        gt = col[i] if i < len(col) else '.'
        a1,a2 = gt_to_alleles(gt)
        rows[s].append((a1,a2))

# Write output table: header then rows
with open(out_file,'w') as O:
    # header
    hdr = ["pop","ind_id"]
    for idx, locus in enumerate(loci, start=1):
        hdr.append(f"Locus_{idx}.A1")
        hdr.append(f"Locus_{idx}.A2")
    O.write("\t".join(hdr) + "\n")
    # rows
    for s in samples:
        pop = s2pop.get(s,"NA")
        alleles = []
        for (a1,a2) in rows[s]:
            alleles.append(a1); alleles.append(a2)
        O.write(pop + "\t" + s + "\t" + "\t".join(alleles) + "\n")

print("Writing:", out_file)
PY
