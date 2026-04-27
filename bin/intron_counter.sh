#!/bin/bash

# Usage: intron_counter.sh <input.bam> <output.txt>

INPUT_BAM="$1"
OUTPUT="$2"

echo "Input BAM: $INPUT_BAM"
echo "Output: $OUTPUT"

# --- Run extraction ---
samtools view -@ 4 "$INPUT_BAM" | awk '
{
    cb = "NA"; re = "NA"
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) { split($i, a, ":"); cb = a[3] }
        if ($i ~ /^RE:A:/) { split($i, a, ":"); re = a[3] }
    }
    print cb "\t" re
}' | sort | uniq -c > "$OUTPUT"
