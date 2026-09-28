#!/bin/bash

# Usage: intron_counter.sh <input.bam> <output.txt>
#
# Tallies cellranger RE:A region tags (E = exonic, N = intronic, I = intergenic)
# per cell barcode. Only records cellranger counted in its UMI matrix are kept:
#   - xf tag has bit 8 set (one representative read per UMI; drops PCR
#     duplicates, multimappers, antisense and non-gene-assigned reads)
#   - unpaired, or first-in-pair (both mates of a PE pair carry xf bit 8)
#   - CB present

INPUT_BAM="$1"
OUTPUT="$2"

echo "Input BAM: $INPUT_BAM"
echo "Output: $OUTPUT"

# --- Run extraction ---
samtools view -@ 4 -F 0x904 "$INPUT_BAM" | awk '
BEGIN { FS = "\t" }
{
    flag = $2
    # paired (0x1) but not first-in-pair (0x40)
    if (flag % 2 == 1 && int(flag / 64) % 2 == 0) next
    cb = ""; re = "NA"; xf = 0
    for (i = 12; i <= NF; i++) {
        if (substr($i, 1, 5) == "CB:Z:") cb = substr($i, 6)
        else if (substr($i, 1, 5) == "RE:A:") re = substr($i, 6)
        else if (substr($i, 1, 5) == "xf:i:") xf = substr($i, 6) + 0
    }
    if (cb == "" || int(xf / 8) % 2 == 0) next
    print cb "\t" re
}' | sort | uniq -c > "$OUTPUT"
