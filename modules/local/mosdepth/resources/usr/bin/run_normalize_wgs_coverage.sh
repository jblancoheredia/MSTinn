#!/usr/bin/env bash

# By blancoj@mskcc.org on 15FEB26 for CMO Technology Innovation Lab

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <mosdepth.summary.txt> [out.tsv]" >&2
  exit 1
fi

IN="$1"
OUT="${2:-${IN%.txt}.autosomes.norm_wholegenome.tsv}"

awk -F'\t' -v OFS='\t' '
  function is_auto_chrom(c) { return (c ~ /^([1-9]|1[0-9]|2[0-2])$/) }

  NR==1 { header=$0; next }

  {
    chrom=$1; len=$2; bases=$3; mean=$4

    rows[++n]=$0
    row_chrom[n]=chrom
    row_len[n]=len
    row_bases[n]=bases
    row_mean[n]=mean

    if (is_auto_chrom(chrom)) {
      L += len
      B += bases
    }
  }

  END {
    if (L == 0) {
      print "ERROR: no autosome whole-genome rows (1-22) found for baseline." > "/dev/stderr"
      exit 2
    }
    baseline = B / L

    print header, "norm_mean_to_autosomes_1_22"

    for (i=1; i<=n; i++) {
      chrom=row_chrom[i]
      if (is_auto_chrom(chrom)) {
        norm = row_mean[i] / baseline
        print rows[i], norm
      }
      # skip non-autosome whole-genome rows entirely
      # (so output is chr1..chr22 only)
    }
  }
' "$IN" > "$OUT"

echo "Wrote: $OUT" >&2
